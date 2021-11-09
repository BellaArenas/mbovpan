// Author: Noah Austin Legall
// Note to self: do NOT run this program within the github folder
// test outside of the github folder. we do not want to make the github too full with run information 

/*
Finish the pipeline
Extend Figure 2 with runtime reports
Testing in sapelo2!!!
*/

log.info """ 
    M B O V P A N (v0.1)    
=============================
A pangenomic pipeline for the analysis of
Mycobacterium bovis isolates 

Project : $workflow.projectDir
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
Manifest's pipeline version: $workflow.manifest.version
=============================
"""

// Where are we grabbing our sequence data from?     
input = ""

// Where should results be published to?
output = ""

// Are these paired-end short reads or longreads?
mode = "short"

// Are we computing SNPs or computing the pangenome? 
run_mode = "snp"

// How many threads will be available to run the pipeline. 
// Automatically uses all the cpus that are available 
threads = Runtime.getRuntime().availableProcessors()

reads = ""

// record the path for the M. bovis reference genome
ref = "$workflow.projectDir/ref/mbovAF212297_reference.fasta"
range = "$workflow.projectDir/auxilary/chrom_ranges.txt" 

// are default parameters included?
if(params.input == null || params.output == null){
    println "necessary paramaters aren't supplied - supply values for input and/or output"
    exit 0
}

else {
    println "necessary paramaters supplied"
    input = params.input
    output = params.output
}

// record the mode the program should be ran in
if(params.mode == "short" || params.mode == "long"){
    mode = params.mode
    println "mbovpan will run in ${mode} read mode"
}
else {
    println "mbovpan will run in default ${mode} read mode"
}

// what part of the pipeline should be ran?
if(params.run == "all" ){
    println "mbovpan will infer snps and pangenome"
    run_mode = "all"
}

else if(params.run == "snp"){
    println "mbovpan will only infer snps"
    run_mode = "snp"
}

else if(params.run == "pan"){
    println "mbovpan will only infer the pangenome"
    run_mode = "pan"
}

else {
    println "mbovpan will infer both snps and pangenome by default"
    run_mode = "all"
}

// how many threads will be utilized
if(params.threads != null){
    println "mbovpan will run using ${threads} threads"
    threads = params.threads
}
else{
    println "mbovpan will run using ${threads} threads by default"
}

// how we download the reads should be based on what mode we run mbovpan with
if(mode == "longread"){
    reads = Channel.fromPath("input*.fastq.gz").ifEmpty { error "Cannot find the read files" }
}
else{
    reads = Channel.fromFilePairs("$input*_{1,2}.fastq.gz").ifEmpty { error "Cannot find the read files" }
}

reads.into {
    reads_process
    reads_trim
}

 // PART 1: Processing 

 log.info """ 
Summary of pipeline run

mode: $mode
run type: $run_mode
reference location: $ref
input: $input
output: $output
no. of threads: $threads
=====================================
"""

/* AUTOMATIC QC of read data */

if(mode == "short"){
    process pre_fastqc {

    publishDir = output

    input:
    tuple sample_id, file(reads_file) from reads_process

    output:
    file("pre_fastqc_${sample_id}_logs") into fastqc_ch1

    script:
    """
    mkdir  pre_fastqc_${sample_id}_logs
    fastqc -o  pre_fastqc_${sample_id}_logs -f fastq -q ${reads_file}
    """
    }


//change this around to get rid of the weird naming
    process fastp {

    publishDir = output

    cpus Math.floor(threads/2)

    input:
    tuple sample_id, file(reads_file) from reads_trim

    output:
    file("${sample_id}_trimmed_R*.fastq") into fastp_ch
 
    script:
    """
    fastp -w ${task.cpus} -q 30 --detect_adapter_for_pe -i  ${reads_file[0]} -I  ${reads_file[1]} -o  ${sample_id}_trimmed_R1.fastq -O  ${sample_id}_trimmed_R2.fastq
    """

    }

    fastp_ch.into {
        fastp_reads1
        fastp_reads2
        fastp_reads3
    }


    process post_fastqc {

    publishDir = output

    input:
    tuple file(trim1), file(trim2) from fastp_reads1

    output:
    file("post_fastqc_${trim1.baseName - ~/_trimmed_R*/}_logs") into fastqc_ch2

    script:
    """
    mkdir  post_fastqc_${trim1.baseName - ~/_trimmed_R*/}_logs
    fastqc -o  post_fastqc_${trim1.baseName - ~/_trimmed_R*/}_logs -f fastq -q ${trim1} ${trim2}
    """
    }
}
else {
    process pre_nanostat {
    publishDir = output

    conda "$workflow.projectDir/envs/longread.yaml"

    cpus Math.floor(threads/2)

    input:
    tuple sample_id, file(read_file) from reads_process

    output:
    file("pre_nanoplot_${sample_id}_logs") into nanoplot_ch

    script:
    """
    mkdir pre_nanoplot_${sample_id}_logs
    NanoPlot -t ${tasks.cpu} --verbose -o ./pre_nanoplot_${sample_id}_logs -p ${sample_id} --N50 --fastq ${read_file}
    """
    }

    process filtlong {
    publishDir = output

    conda "$workflow.projectDir/envs/longread.yaml"

    input:
    tuple sample_id, file(read) from reads_trim

    output:
    file("${sample_id}_trimmed.fastq") into filtlong_ch

    script:
    """
    filtlong --keep_percent 90 --min_length 1000 --min_window_q 9 ${read} > ${sample_id}_trimmed.fastq
    """
    }

    filtlong_ch.into {
        filtlong_reads1
        filtlong_reads2
        filtlong_reads3
    }

    process post_nanostat {
    publishDir = output

    conda "$workflow.projectDir/envs/longread.yaml"

    cpus Math.floor(threads/2)

    input:
    file(trim) from filtlong_reads1

    output:
    file("post_nanoplot_${trim.baseName}_logs") into nanoplot_ch

    script:
    """
    mkdir post_nanoplot_${trim.baseName}_logs
    NanoPlot -t ${tasks.cpu} --verbose -o ./post_nanoplot_${trim.baseName}_logs -p ${trim.baseName} --N50 --fastq ${trim}
    """
    }
}

/*
 process multiqc {
    publishDir = output

    input:
    file(pre) from fastqc_ch1.collect()
    file(post) from fastqc_ch2.collect()

    output:
    file("mbovpan_seq_qual*")

    script:
    """
    multiqc -n mbovpan_seq_qual ${pre} ${post}
    """
}
*/


// MODE 1: Variant Calling 
bam = Channel.create()

if(run_mode == "snp" || run_mode == "all"){
    if(mode == "short"){

    process read_map {
    publishDir = output 

    cpus Math.floor(threads/2)
    
    conda "$workflow.projectDir/envs/samtools.yaml"
   
    input:
    tuple file(trim1), file (trim2) from fastp_reads2 
    //path(reference) from ref

    output:
    file("${trim1.baseName - ~/_trimmed_R*/}.bam")  into map_ch 

    script:
    """
    bwa mem -t ${task.cpus}-M -R "@RG\\tID:${trim1.baseName - ~/_trimmed_R*/}\\tSM:${trim1.baseName - ~/_trimmed_R*/}\\tPL:ILLUMINA\\tPI:250" ${ref} ${trim1} ${trim2} | samtools view -Sb | samtools sort -o ${trim1.baseName - ~/_trimmed_R*/}.bam
    """
    }

    bam = map_ch
}

else{
    process longread_map {
    publishDir = output 

    cpus Math.floor(threads/2)
    
    conda "$workflow.projectDir/envs/samtools.yaml"
   
    input:
    tuple file(trim1), file (trim2) from filtlong_reads2 
    //path(reference) from ref

    output:
    file("${trim1.baseName - ~/_trimmed_R*/}.bam")  into map_ch 

    script:
    """
    minimap2 -ax map-pb --MD ${ref} ${trim} | samtools view -Sb | samtools sort -o ${trim1.baseName - ~/_trimmed_R*/}.bam
    """
    }

    bam = map_ch

}



    

// picard MarkDuplicates INPUT={sorted_bamfile} OUTPUT={rmdup_bamfile} ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE=dup_metrics.csv

// Important to have 'USE_JDK_DEFLATER=true USE_JDK_INFLATER=true'
    process mark_dups {
    publishDir = output 

    conda "$workflow.projectDir/envs/picard.yaml"    

    input:
    file(bam) from bam

    output:
    file("${bam.baseName}.nodup.bam") into nodup_ch

    script:
    """
    picard MarkDuplicates INPUT=${bam} OUTPUT=${bam.baseName}.nodup.bam ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE=dup_metrics.csv USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
    samtools index ${bam.baseName}.nodup.bam
    cp ${bam.baseName}.nodup.bam.bai $workflow.launchDir
    """
    }

    process freebayes {
    publishDir = output 

    cpus Math.floor(threads/2)

    conda "$workflow.projectDir/envs/freebayes.yaml"

    input:
    file(bam) from nodup_ch
    //file(range) from chrom_range
    path(reference) from ref

    output:
    file("${bam.baseName - ~/.nodup/}.vcf") into freebayes_ch 

    script:
    """
    cp $workflow.launchDir/${bam.baseName - ~/.nodup/}.nodup.bam.bai ./
    freebayes-parallel ${range} ${task.cpus} -f ${ref} ${bam} > ${bam.baseName - ~/.nodup/}.vcf
    """
    }

//Start with a basic QUAL > 150, later add a parameter for changing this
    process vcf_filter {
    publishDir = output 

    conda "$workflow.projectDir/envs/vcflib.yaml"

    input:
    file(vcf) from freebayes_ch

    output:
    file("${vcf.baseName}.filtered.vcf") into filter_ch

    script:
    """
    vcffilter -f "QUAL > 150" ${vcf} > ${vcf.baseName}.filtered.vcf
    """

    } 
}

// vcf filtering + generate alignment? 
// pangenome steps (might need to separately create environment for pangenome software. inclusion in main environment might lead to conflicts)
// roary, quast OR busco

// PART 3: Pangenome 
assembly = Channel.create()
if(run_mode == "pan" || run_mode == "all"){
    if(mode == "short"){
    
    process assembly {
    publishDir = output 

    cpus Math.floor(threads/2)

    input:
    tuple file(trim1), file(trim2) from fastp_reads3

    output:
    file("${trim1.baseName}.scaffold.fasta") into shortassembly_ch

    script:
    """
    spades.py --threads ${task.cpus} --only-assembler --careful -1 ${trim1} -2 ${trim2} -o ${trim1.baseName}
    cd ${trim1.baseName}
    mv scaffolds.fasta  ../${trim1.baseName}.scaffold.fasta
    """
}
assembly_ch = shortassembly_ch
}
else{
    
    process longread_assembly {
    publishDir = output 

    conda "$workflow.projectDir/envs/raven.yaml"

    cpus Math.floor(threads/2)

    input:
    file(trim) from filtlong_reads3

    output:
    file("${trim1.baseName}.scaffold.fasta") into longassembly_ch

    script:
    """
    raven --threads ${task.cpus} ${trim} > ${trim.baseName}.scaffold.fasta"
    """
}  

assembly_ch = longassembly_ch
}

process annotate {
    publishDir = output

    conda "$workflow.projectDir/envs/prokka.yaml"

    input:
    file(assembly) from assembly_ch

    output:
    file("${assembly.baseName}.annot.gff") into annotate_ch

    script:
    """
    prokka  --outdir ./${assembly.baseName} --prefix ${assembly.baseName}.annot ${assembly}
    cp ./${assembly.baseName}/${assembly.baseName}.annot.gff ./
    """
}

process roary {
    publishDir = output

    conda "$workflow.projectDir/envs/roary.yaml"

    cpus Math.floor(threads/2)

    input:
    file(gff) from annotate_ch.collect()

    output:
    file("*") into roary_ch

    script:
    """
    roary -p ${task.cpus} $gff
    """
}
}


