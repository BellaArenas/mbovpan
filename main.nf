// Author: Noah Austin Legall
nextflow.enable.dsl=2
     
input = ""
output = ""
run_mode = "snp"

params.version = false 
if(params.version){
    println("v0.1")
    exit(0)
}

params.help = false
if(params.help){
    println(
"""
    M B O V P A N (v1)    
=============================
A pangenomic pipeline for the analysis of
Mycobacterium bovis isolates 


usage: nextflow run mbovpan/mbovpan.nf [options] --input ./path/to/input --output ./path/to/output
  options:
    --run [all|snp|pan]: 
        Specifies in what mode to run mbovpan in [DEFAULT:all]
    --qual [INT]:
        The minimum QUAL score for a SNP to be considered [DEFAULT:20]
    --depth [INT]:
        The minimum DP score for a SNP to be considered [DEFAULT:25]
    --mapq [INT]:
        The minimum MQ score for a SNP to be considered [DEFAULT:40]
    --threads [INT]:
        How many threads to use for the programs [DEFAULT:(number of avail. threads)/2]
    --help
        Prints this help message
    --version
        Prints the current version 
=============================
"""

    )
    exit(0)
}

// Automatically uses all the cpus that are available 
// If not specified, use 50% of available resources 
params.threads = Math.floor(Runtime.getRuntime().availableProcessors()/2)

reads = ""

// inital values based on previous benchmark work from dissertation
// can be changed by user
params.qual = 20

params.depth = 25 

params.mapq = 40


if(params.qual){
    qual = params.qual as Integer
    }
if(params.depth){
    depth = params.depth as Integer
    }
if(params.mapq){
    mapq = params.mapq as Integer
    }


// record the path for key files in the repo
ref = "$workflow.projectDir/ref/mbovAF212297_reference.fasta" // reference genome
range = "$workflow.projectDir/chrom_ranges.txt" // for freebayes-parallel
spotyping = "$workflow.projectDir/scripts/SpoTyping/SpoTyping.py" //spoligotyping
check = "$workflow.projectDir/scripts/lineage_check.py" 
lineage_table = "$workflow.projectDir/scripts/lineage_table.py"

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

if(params.threads){
    println "mbovpan will run using ${params.threads} threads"
    threads = params.threads
}
else{
    println "mbovpan will run using ${threads} threads by default"
}


// finding the reads from the input file location 
reads = Channel.fromFilePairs("$input*{1,2}*.f*q*").ifEmpty { error "Cannot find the read files" }


println(""" 
    M B O V P A N (v1.0.0)    
=============================
A pangenomic pipeline for the analysis of
Mycobacterium bovis isolates 


Project : $workflow.projectDir
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
Manifest's pipeline version: $workflow.manifest.version
=============================

Summary of pipeline run

run type: $run_mode
reference location: $ref
input: $input
output: $output
no. of threads: $threads
QUAL: $qual
MAPQ: $mapq
DEPTH: $depth
=====================================
""")

// Start of the nextflow processes

workflow {

    if (params.version) {
        println("v0.1")
        System.exit(0)
    }

    if (params.help) {
        println("""
            M B O V P A N (v1.0.0)

            =============================

            A pangenomic pipeline for the analysis of Mycobacterium bovis isolates

            Usage: nextflow run mbovpan --input ./path/to/input --output ./path/to/output --run [all|snp|pan]

            =============================

        """)
        System.exit(0)
    }

    if (params.input == null || params.output == null) {
        error "necessary parameters aren't supplied - supply values for input and/or output"
    }

    println("necessary parameters supplied")

    println("mbovpan will run using ${params.threads} threads by default")

    // Define the channel for read pairs

    reads = Channel.fromFilePairs("${params.input}*{1,2}*.f*q*", size: 2, flat: true).ifEmpty { error "Cannot find the read files" }

    spotyping(reads)
    pre_fastqc(spotyping.out.spoligo_ch)
    fastp(spotyping.out.spoligo_ch)
    post_fastqc(fastp.out.trimmed_reads)
    lineage(fastp.out.trimmed_reads)
    if(run_mode == "snp" || run_mode == "all"){
        read_map(fastp.out.trimmed_reads)
        mark_dups(read_map.out.bam)
    }

}

process spotyping {

    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/spotyping", mode: 'copy'

    conda "${workflow.projectDir}/envs/spotyping.yaml"

    errorStrategy 'ignore'

    debug true

    input:

    tuple val(sample_id), path(read1), path(read2)

    output:

    tuple val(sample_id), path("${read1.baseName}.fastq"), path("${read2.baseName}.fastq"), emit: spoligo_ch

    path("${read1.baseName - ~/_1/}.out"), emit: spoligo_multi

    script:

    """

    python3 ${workflow.projectDir}/scripts/SpoTyping/SpoTyping.py ${read1} ${read2} -o ${read1.baseName - ~/_1/}.out > stdout.txt 

    """

}   

process pre_fastqc {
    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/fastqc", mode: 'copy'

    conda 'bioconda::fastqc'

    input:

    tuple val(sample_id), path(read_one), path(read_two)

    output:

    path("pre_fastqc_${read_one.baseName - ~/_1*/}_logs/*"), emit: fastqc_ch1

    script:
    """
    mkdir pre_fastqc_${read_one.baseName - ~/_1*/}_logs

    fastqc -o pre_fastqc_${read_one.baseName - ~/_1*/}_logs -f fastq -q ${read_one} ${read_two}
    """

}

process fastp {
    tag "${sample_id}"
    conda 'bioconda::fastp'

    publishDir "${params.output}/mbovpan_results/read_trimming", mode: 'copy'

    input:
    tuple val(sample_id), path(read_one), path(read_two)

    output:
    tuple val(sample_id), path("${read_one.baseName.replaceAll('_1$', '')}_trimmed_R1.fastq"), path("${read_one.baseName.replaceAll('_1$', '')}_trimmed_R2.fastq"), emit: trimmed_reads

    script:
    """
    fastp -w ${params.threads} -q 30 --detect_adapter_for_pe -i ${read_one} -I ${read_two} -o ${read_one.baseName.replaceAll("_1\$", "")}_trimmed_R1.fastq -O ${read_one.baseName.replaceAll("_1\$", "")}_trimmed_R2.fastq    
    """
}

process post_fastqc {
    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/fastqc", mode: 'copy'

    conda 'bioconda::fastqc'

    input:
    tuple val(sample_id), path(trim_one), path(trim_two)

    output:
    path("post_fastqc_${trim_one.baseName - ~/_1*/}_logs/*"), emit: fastqc_ch2

    script:
    """
    mkdir post_fastqc_${trim_one.baseName - ~/_1*/}_logs

    fastqc -o post_fastqc_${trim_one.baseName - ~/_1*/}_logs -f fastq -q ${trim_one} ${trim_two}
    """
}

process lineage {
    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/lineage", mode: 'copy'

    conda "${workflow.projectDir}/envs/tbprofile.yaml"

    input:
    tuple val(sample_id), path(trim_one), path(trim_two)

    output:
    path("${trim_one.baseName.replaceAll('_trimmed_R.*', '')}.results.json"), emit: tbprofile_ch

    script:
    """
    tb-profiler profile --read1 ${trim_one} --read2 ${trim_two} --no_delly --prefix ${trim_one.baseName.replaceAll('_trimmed_R.*', '')}
    cp ./results/${trim_one.baseName.replaceAll('_trimmed_R.*', '')}.results.json ./
    """
}

process read_map {
    tag "${sample_id}"

    publishDir "${params.output}", mode: 'copy'

    conda "samtools bioconda::bowtie2"

    input:
    tuple val (sample_id), path(trim_one), path(trim_two)

    output:
    tuple val(sample_id), path("${trim_one.baseName.replaceAll('_trimmed_R.*', '')}.bam"), emit: bam

    script:
    """
    bowtie2 --threads ${task.cpus} -x ${workflow.projectDir}/ref/mbov_bowtie_index -1 ${trim_one} -2 ${trim_two} | samtools view -Sb | samtools sort -o ${trim_one.baseName.replaceAll('_trimmed_R.*', '')}.bam
    """
}

process mark_dups {
    tag "${bam.baseName}"

    publishDir "${params.output}/mbovpan_results/readmapping", mode: 'copy'

    conda "${workflow.projectDir}/envs/picard.yaml"

    input:
    tuple val(sample_id), file(bam)

    output:
    tuple val(sample_id), file("${bam.baseName}.nodup.bam"), emit: nodup_bam
    path "${bam.baseName}.nodup.bam.bai", emit: nodup_bai

    script:
    """
    picard MarkDuplicates INPUT=${bam} OUTPUT=${bam.baseName}.nodup.bam ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE=dup_metrics.csv USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
    samtools index ${bam.baseName}.nodup.bam
    """
}

//     process freebayes {
//     publishDir = "$output/mbovpan_results/variant_calling" 

//     cpus threads

//     conda "$workflow.projectDir/envs/freebayes.yaml"

//     input:
//     file(bam) from nodup1_ch
//     path(reference) from ref

//     output:
//     file("${bam.baseName - ~/.nodup/}.vcf") into freebayes_ch 

//     script:
//     """
//     cp $workflow.launchDir/${bam.baseName - ~/.nodup/}.nodup.bam.bai ./
//     freebayes-parallel ${range} ${task.cpus} -f ${ref} ${bam} > ${bam.baseName - ~/.nodup/}.vcf
//     """
//     }

//     process vcf_filter {
//     publishDir = "$output/mbovpan_results/variant_calling" 

//     conda "$workflow.projectDir/envs/vcflib.yaml"

//     debug true

//     input:
//     file(vcf) from freebayes_ch

//     output:
//     file("${vcf.baseName}.filtered.vcf") into filter_ch

//     script:
//     """
//     vcffilter -f "QUAL > ${qual}" ${vcf} | vcffilter -f "DP > ${depth}" | vcffilter -f "MQM > ${mapq}" |  vcffilter -f "TYPE = snp" | bedtools intersect -header -a - -b $workflow.projectDir/ref/pe_ppe_regions.gff3 -v > ${vcf.baseName}.filtered.vcf
//     """

//     } 

//     filter_ch.into {
//         filter1_ch
//         filter2_ch
//     }

//     stats_ch = fastp_reads4.merge(nodup2_ch).merge(filter2_ch)


//     process psuedo_assembly {
//         publishDir = "$output/mbovpan_results/assemblies"

//         conda "$workflow.projectDir/envs/consensus.yaml"

//         input:
//         file(vcf) from filter1_ch 

//         output:
//         file("${vcf.baseName - ~/.filtered/}.consensus.fasta") into fasta_ch

//         script:
//         """
//         bgzip ${vcf} 
//         bcftools index ${vcf}.gz
//         bcftools norm --check-ref s --fasta-ref $ref -Ov ${vcf}.gz > ${vcf.baseName}.norm.vcf
//         bgzip ${vcf.baseName}.norm.vcf
//         bcftools index ${vcf.baseName}.norm.vcf.gz
//         cat ${ref} | vcf-consensus ${vcf.baseName}.norm.vcf.gz > ${vcf.baseName - ~/.filtered/}.dummy.fasta
//         sed 's|LT708304.1|${vcf.baseName - ~/.filtered/}}|g' ${vcf.baseName - ~/.filtered/}.dummy.fasta > ${vcf.baseName - ~/.filtered/}.consensus.fasta
//         """
//     }
    
//     process iqtree_phylo {
//         publishDir = "$output/mbovpan_results/phylogeny"
        
//         conda "$workflow.projectDir/envs/iqtree.yaml"

//         errorStrategy "ignore"
        
//         cpus threads 

//         input:
//         file(aln) from fasta_ch.collect()
        
//         output:
//         file("*") into phylo_ch 
//         file("phylo_out.txt") into statistics_ch
        
//         script:
//         """
//         cat *.fasta > mbovpan_align.fasta
//         cat $workflow.projectDir/ref/mbov_reference.fasta >> mbovpan_align.fasta
//         snp-sites -o mbovpan_align.snp_only.fasta mbovpan_align.fasta
//         iqtree -s mbovpan_align.snp_only.fasta -m MFP -nt ${task.cpus} -bb 1000 -pre mbovpan_align -o "LT708304.1"
//         echo "file created to make statistics file build conda env last" > phylo_out.txt
//         """
//     }

//     process stats {
//         publishDir = "$output/mbovpan_results/statistics"

//         conda "$workflow.projectDir/envs/statistics.yaml"

//         input:
//         file(nec_files) from stats_ch
//         file(dummy_variable) from statistics_ch

//         output:
//         file("${nec_files[0].baseName}.stats") into output_stat_ch


//         script:
//         """
//         python $workflow.projectDir/scripts/statistics.py ${nec_files[0]} ${nec_files[1]} ${nec_files[2]} ${nec_files[3]} > ${nec_files[0].baseName}.stats
//         """

//     }
// }

// assembly = Channel.create()
// if(run_mode == "pan" || run_mode == "all"){
    
//     process assembly {
//     publishDir = output 
    
//     conda "bioconda::spades"
    
//     errorStrategy "ignore"

//     cpus threads

//     input:
//     tuple file(trim1), file(trim2) from fastp_reads3

//     output:
//     file("${trim1.baseName - ~/_trimmed_R*/}.scaffold.fasta") into assembly_ch

//     script:
//     """
//     mkdir ${trim1.baseName}
//     spades.py -1 ${trim1} -2 ${trim2} --careful -o ${trim1.baseName} -t ${task.cpus} --only-assembler
//     cd ${trim1.baseName}
//     mv scaffolds.fasta  ../${trim1.baseName - ~/_trimmed_R*/}.scaffold.fasta
//     """
// }

// assembly_ch = assembly_ch

// assembly_ch.into {
//     assembly_ch1
//     assembly_ch2
// }

// process quast {
//     publishDir = "$output/mbovpan_results/statistics"

//     conda "$workflow.projectDir/envs/quast.yaml"
    
//     cpus threads

//     input:
//     file(assemblies) from assembly_ch1.collect()
    
//     output:
//     file("*") into quast_ch
    
//     script:
//     """
//     quast -o assembly_stats ${assemblies} -r ${ref} -t ${task.cpus}
//     """
// }


// process annotate {
//     publishDir = "$output/mbovpan_results/annotations"
    
//     cpus threads

//     conda "$workflow.projectDir/envs/prokka.yaml"

//     errorStrategy "ignore"

//     input:
//     file(assembly) from assembly_ch2

//     output:
//     file("${assembly.baseName - ~/.scaffold/}.annot.gff") into annotate_ch

//     script:
//     """
//     prokka  --outdir ./${assembly.baseName - ~/.scaffold/} --cpus ${task.cpus} --prefix ${assembly.baseName - ~/.scaffold/}.annot ${assembly} 
//     cp ./${assembly.baseName - ~/.scaffold/}/${assembly.baseName - ~/.scaffold/}.annot.gff ./
//     """
// }

// process panaroo {
//     publishDir = "$output/mbovpan_results/pangenome"

//     cpus threads

//     input:
//     file(gff) from annotate_ch.collect()

//     output:
//     file("*") into roary_ch

//     script:
//     """
//     panaroo -i *.gff -o ./ -t ${task.cpus} -a core --core_threshold 0.98 --clean-mode strict
//     """
// }

// roary_ch.into{
//     roary_ch2
//     roary_ch3
//     roary_ch4
//     roary_ch5
//     roary_ch_filter
// }


// process iqtree_core {
//         publishDir = "$output/mbovpan_results/phylogeny"
        
//         conda "$workflow.projectDir/envs/iqtree.yaml"
        
//         cpus threads 
        
//         errorStrategy 'ignore'
        
        
//         input:
//         file(input) from roary_ch.collect()
    
//         output:
//         file("*") into iqtreecore_ch
        
//         script:
//         """
//         snp-sites -o mbovpan_align.snp_only.fasta core_gene_alignment.aln
//         iqtree -s mbovpan_align.snp_only.fasta -m MFP -nt ${task.cpus} -bb 1000 -pre mbovis_core 
//         """
//     }


// process filter_pan {
//     publishDir = "./"
    
//     conda "bioconda::blast=2.9.0 pandas"

//     debug 'true'
    
//     input:
//     file("*") from roary_ch_filter.collect()

//     output:
//     file('mbovis_filtered_cogs.csv')
    
    
//     script:
//     """
//     # first determine the genes that are present and use blast to find their true annotation
//     blastn -query pan_genome_reference.fa -max_target_seqs 1 -db $workflow.projectDir/ref/mbovis_reference -out mb.out -outfmt "6 delim=, qseqid sseqid length qstart qend sstart send qlen slen"
//     echo "qseqid,sseqid,length,qstart,qend,sstart,send,qlen,slen" > heading.txt 
//     cat heading.txt mb.out >> mb.out.csv

//     # once mb.out is created, we can use python to 1) create length %, 2) filter by 75% or higher, and 3) reduce repetitive blast results
//     python $workflow.projectDir/scripts/blast_result_filter.py gene_presence_absence.csv
//     """
// }

// }

// process multiqc {
//     publishDir = "$output/mbovpan_results/statistics"
    
//     conda "$workflow.projectDir/envs/multiqc.yaml"

//     input:
//     file(pre) from fastqc_ch1.collect().ifEmpty([])
//     file(post) from fastqc_ch2.collect().ifEmpty([])
    
//     output:
//     file("mbovpan_report*")

//     script:
//     """
//     multiqc -n mbovpan_report .
//     """
// }

// process mbovis_verification {
//     publishDir = "$output/mbovpan_results/lineage_info"

//     input:
//     file(spoligotype_info) from spoligo_multi.collect().ifEmpty([])
//     file(lineage_info) from tbprofile_ch.collect().ifEmpty([])

//     output:
//     file("mbovpan_lineage_info.csv")

//     script:
//     """
//     python ${lineage_table} > mbovpan_lineage_info.csv
//     """
// }
