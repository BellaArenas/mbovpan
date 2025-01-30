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

// reads = Channel.fromFilePairs("$input/*_{1,2}.f*q*", size: 2, flat: true).ifEmpty { error "Cannot find the read files" }

// Define channels for BAM and BAI files
// mark_dups_ch = Channel.fromPath("/xdisk/lilianasalvador/isabellaarenas/mbovpan/outputs/snp_testing_0119/mbovpan_results/readmapping/*.nodup.bam")
//                      .map { path -> tuple(path.baseName.replaceAll('.nodup', ''), path) }

// mark_dups_bai_ch = Channel.fromPath("/xdisk/lilianasalvador/isabellaarenas/mbovpan/outputs/snp_testing_0119/mbovpan_results/readmapping/*.nodup.bam.bai")
//                      .map { path -> tuple(path.baseName.replaceAll('.nodup.bam', ''), path) }

freebayes_ch = Channel.fromPath('/xdisk/lilianasalvador/isabellaarenas/mbovpan/outputs/snp_testing_0121/mbovpan_results/variant_calling/*.vcf')
                     .filter { !it.name.endsWith('.filtered.vcf') }
                     .map { path -> tuple(path.baseName.replaceFirst(/\.vcf$/, ''), path) }

// trimmed_reads = Channel.fromFilePairs("/xdisk/lilianasalvador/isabellaarenas/mbovpan/outputs/snp_testing_0119/mbovpan_results/read_trimming/*_trimmed_R{1,2}.fastq", flat: true)

// nodup_bam = Channel.fromPath("/xdisk/lilianasalvador/isabellaarenas/mbovpan/outputs/snp_testing_0119/mbovpan_results/readmapping/*.nodup.bam")

// assembly_ch = Channel.fromPath("/xdisk/lilianasalvador/isabellaarenas/mbovpan/outputs/pilon_01/mbovpan_results/assembly/*/*.scaffold.fasta")

// contigs_ch = Channel.fromPath("/xdisk/lilianasalvador/isabellaarenas/mbovpan/outputs/pilon_01/mbovpan_results/assembly/*/*.contigs.fasta")


// assembly_ch = Channel.fromPath("/xdisk/lilianasalvador/isabellaarenas/mbovpan/outputs/pilon_01/mbovpan_results/assembly/*/*.scaffold.fasta")
//     .map { path -> tuple(path.baseName.replaceFirst(/\.scaffold$/, ''), path) }

// contigs_ch = Channel.fromPath("/xdisk/lilianasalvador/isabellaarenas/mbovpan/outputs/pilon_01/mbovpan_results/assembly/*/*.contigs.fasta")
//     .map { path -> tuple(path.baseName.replaceFirst(/\.contigs$/, ''), path) }

// sorted_bam = Channel.fromPath("/xdisk/lilianasalvador/isabellaarenas/mbovpan/outputs/pilon/mbovpan_results/pilon/*.sorted.bam")



// bai = Channel.fromPath("/xdisk/lilianasalvador/isabellaarenas/mbovpan/outputs/pilon/mbovpan_results/pilon/*.sorted.bam.bai")

// panaroo_ch = Channel.fromPath("/xdisk/lilianasalvador/isabellaarenas/mbovpan/outputs/galway_girl/mbovpan_results/pangenome/*")

// filter_ch = Channel.fromPath("/xdisk/lilianasalvador/isabellaarenas/mbovpan/outputs/snp_testing_0119/mbovpan_results/variant_calling/*.filtered.vcf.gz")

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

    spotyping(reads)
    pre_fastqc(spotyping.out.spoligo_ch)
    fastp(spotyping.out.spoligo_ch)
    post_fastqc(fastp.out.trimmed_reads)
    lineage(fastp.out.trimmed_reads)
    if(run_mode == "snp" || run_mode == "all"){
        // read_map(fastp.out.trimmed_reads)
        // mark_dups(read_map.out.bam)
        // freebayes(mark_dups.out.nodup_bam, mark_dups.out.nodup_bai)
        // bcftools_norm(freebayes.out.freebayes_ch)
        bcftools_norm(freebayes_ch)
        vcf_filter(bcftools_norm.out.norm_vcf_ch)
        // vcf_filter(freebayes_ch)
        // stats_ch = fastp.out.trimmed_reads.merge(mark_dups.out.nodup_bam).merge(vcf_filter.out.filter_ch)
        // stats_ch = trimmed_reads.merge(nodup_bam).merge(freebayes_ch)
        psuedo_assembly(vcf_filter.out.filter_ch)
        // psuedo_assembly(filter_ch)
        // iqtree_phylo(psuedo_assembly.out.fasta_ch.collect())
        // stats(stats_ch, iqtree_phylo.out.statistics_ch)
        // iqtree_phylo(psuedo_assembly.out.fasta_ch.collect())
        // stats(stats_ch, iqtree_phylo.out.statistics_ch)
    }

    if(run_mode == "pan" || run_mode == "all"){
        assembly(fastp.out.trimmed_reads)
        quast_assembly(assembly.out.assembly_ch)
        bwa_index(assembly.out.contigs_ch)
        bwa_mem_inputs = bwa_index.out.bwa_index_out
        .join(assembly.out.contigs_ch, by: 0)
        .join(fastp.out.trimmed_reads, by: 0)
        bwa_mem(bwa_mem_inputs)
        samtools_index(bwa_mem.out.sorted_bam)
        pilon_inputs = assembly.out.contigs_ch
            .join(bwa_mem.out.sorted_bam, by: 0)
            .join(samtools_index.out.bai, by: 0)
        pilon(pilon_inputs)
        quast_pilon(pilon.out.pilon_fasta)
        annotate(pilon.out.pilon_fasta)
        panaroo(annotate.out.annotate_ch.collect())
        iqtree_core(panaroo.out.panaroo_ch.collect())
        filter_pan(panaroo.out.panaroo_ch.collect())
    }
    multiqc(pre_fastqc.out.fastqc_ch1.collect(), post_fastqc.out.fastqc_ch2.collect())
    mbovis_verification(spotyping.out.spoligo_multi.collect(), lineage.out.tbprofile_ch.collect())
}
    
process spotyping {
    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/spotyping", mode: 'copy'

    conda "${workflow.projectDir}/envs/spotyping.yaml"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_1.fastq"), path("${sample_id}_2.fastq"), emit: spoligo_ch
    path("${sample_id}.out"), emit: spoligo_multi
    path("SITVIT_ONLINE.*.xls"), emit: sitvit_xls

    script:
    """
    echo "Starting spotyping process for sample: ${sample_id}"
    echo "Read1: ${read1}"
    echo "Read2: ${read2}"

    if [ -z "${read1}" ] || [ -z "${read2}" ]; then
        echo "Error: One of the input files is null"
        exit 1
    fi

    # Decompress the input files if they are gzipped
    if [[ "${read1}" == *.gz && "${read2}" == *.gz ]]; then
        echo "Decompressing ${read1} and ${read2}"
        gunzip -c ${read1} > ${sample_id}_1.fastq
        gunzip -c ${read2} > ${sample_id}_2.fastq
        python3 ${workflow.projectDir}/scripts/SpoTyping/SpoTyping.py -o ${sample_id}.out --debug -O . ${sample_id}_1.fastq ${sample_id}_2.fastq
    else
        python3 ${workflow.projectDir}/scripts/SpoTyping/SpoTyping.py -o ${sample_id}.out --debug -O . ${read1} ${read2}
    fi

    if [ \$? -ne 0 ]; then
        echo "Error in SpoTyping.py script"
        exit 1
    fi

    echo "spotyping process completed successfully for sample: ${sample_id}"
    echo "Checking output files:"
    ls -lh .
    """
}

process pre_fastqc {
    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/fastqc", mode: 'copy'

    conda "${workflow.projectDir}/envs/fastqc.yaml"

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

    publishDir "${params.output}/mbovpan_results/read_trimming", mode: 'copy'
    
    conda 'bioconda::fastp'

    input:
    tuple val(sample_id), path(read_one), path(read_two)

    output:
    tuple val(sample_id), path("${sample_id}_trimmed_R1.fastq"), path("${sample_id}_trimmed_R2.fastq"), emit: trimmed_reads

    script:
    """
    fastp -w ${task.cpus} -q 30 --detect_adapter_for_pe -i ${read_one} -I ${read_two} -o ${sample_id}_trimmed_R1.fastq -O ${sample_id}_trimmed_R2.fastq    
    """
}

process post_fastqc {
    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/fastqc", mode: 'copy'

    conda "${workflow.projectDir}/envs/fastqc.yaml"

    input:
    tuple val(sample_id), path(trim_one), path(trim_two)

    output:
    path("post_fastqc_${trim_one.baseName.replaceAll('_R[12]$', '')}_logs/*"), emit: fastqc_ch2

    script:
    """
    mkdir post_fastqc_${trim_one.baseName.replaceAll('_R[12]$', '')}_logs

    fastqc -o post_fastqc_${trim_one.baseName.replaceAll('_R[12]$', '')}_logs -f fastq -q ${trim_one} ${trim_two}
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
    // path("${trim_one.baseName.replaceAll('_trimmed_R.*', '')}.docx"), emit: tbprofile_docx

    script:
    """
    tb-profiler profile --read1 ${trim_one} --read2 ${trim_two} --no_delly --prefix ${trim_one.baseName.replaceAll('_trimmed_R.*', '')} --docx --docx_template /groups/lilianasalvador/barenas/scripts/mbovpan/ref/lineage_template.docx
    cp ./results/${trim_one.baseName.replaceAll('_trimmed_R.*', '')}.results.json ./
    """
}

process read_map {
    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/alignment", mode: 'copy'

    conda "samtools bioconda::bowtie2"

    input:
    tuple val (sample_id), path(trim_one), path(trim_two)

    output:
    tuple val(sample_id), path("${trim_one.baseName.replaceAll('_trimmed_R.*', '')}.bam"), emit: bam

    script:
    """
    bowtie2 --threads ${params.threads} -x ${workflow.projectDir}/ref/mbov_bowtie_index -1 ${trim_one} -2 ${trim_two} | samtools view -Sb | samtools sort -o ${trim_one.baseName.replaceAll('_trimmed_R.*', '')}.bam
    """
}

process mark_dups {
    tag "${bam.baseName}"

    publishDir "${params.output}/mbovpan_results/readmapping", mode: 'copy'

    conda "${workflow.projectDir}/envs/picard.yaml"

    input:
    tuple val(sample_id), file(bam)

    output:
    tuple val(sample_id), path("${bam.baseName}.nodup.bam"), emit: nodup_bam
    tuple val(sample_id), path("${bam.baseName}.nodup.bam.bai"), emit: nodup_bai

    script:
    """
    picard MarkDuplicates INPUT=${bam} OUTPUT=${bam.baseName}.nodup.bam ASSUME_SORTED=true REMOVE_DUPLICATES=true METRICS_FILE=dup_metrics.csv USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
    samtools index ${bam.baseName}.nodup.bam
    """
}

process freebayes {
    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/variant_calling", mode: 'copy'

    conda "${workflow.projectDir}/envs/freebayes.yaml"
    // conda "${workflow.projectDir}/envs/consensus.yaml"

    input:
    tuple val(sample_id), path(nodup_bam)
    tuple val(sample_id), path(nodup_bai)

    output:
    tuple val(sample_id), path("${sample_id}.vcf"), emit: freebayes_ch
 
    script:
    """
    freebayes-parallel ${range} ${params.threads} -f ${ref} ${nodup_bam} > ${sample_id}.vcf
    """
}

// bcftools mpileup -f ${ref} ${nodup_bam} | bcftools call -mv -Ov -o ${nodup_bam.baseName.replaceAll('.nodup', '')}.vcf

process bcftools_norm {
    tag "${vcf}"

    publishDir "${params.output}/mbovpan_results/variant_calling", mode: 'copy'

    conda "${workflow.projectDir}/envs/consensus.yaml"

    input:
    tuple val(sample_id), path(vcf)

    output:
    tuple val(sample_id), path("${sample_id}.norm.vcf.gz"), emit: norm_vcf_ch
    tuple val(sample_id), path("${sample_id}.vcf.gz.csi"), emit: vcf_index_ch

    script:
    """
    set -euxo pipefail

    echo "Compressing VCF file: ${vcf}"
    bgzip ${vcf}

    echo "Indexing VCF file"
    bcftools index ${vcf}.gz

    echo "Normalizing VCF file"

    bcftools norm --threads ${params.threads} --check-ref s --fasta-ref ${ref} -Oz ${vcf}.gz -o ${sample_id}.norm.vcf.gz
    """
}

process vcf_filter {
    tag "${normalized_vcf}"

    publishDir "${params.output}/mbovpan_results/variant_calling", mode: 'copy'

    conda "${workflow.projectDir}/envs/vcflib.yaml"

    input:
    tuple val(sample_id), path(normalized_vcf)

    output:
    tuple val(sample_id), path("${sample_id}.norm.filtered.vcf"), emit: filter_ch

    // script:
    // """    
    // vcffilter -f "QUAL > ${qual}" ${vcf} | vcffilter -f "DP > ${depth}" | vcffilter -f "MQM > ${mapq}" | vcffilter -f "TYPE = snp" | bedtools intersect -header -a - -b ${workflow.projectDir}/ref/pe_ppe_regions.gff3 -v | bgzip > ${vcf.baseName.replaceAll(".vcf","")}.filtered.vcf.gz
    // """

    script:
    """
    echo "Filtering VCF file using bcftools filter: ${normalized_vcf}"
    bcftools filter -i 'QUAL>${qual} & INFO/DP>${depth} & INFO/MQM>${mapq} & TYPE="snp"' --IndelGap 2 ${normalized_vcf} | bedtools intersect -header -a - -b ${workflow.projectDir}/ref/pe_ppe_regions.gff3 -v > ${sample_id}.norm.filtered.vcf
    """
}

process psuedo_assembly {
    tag "${filtered_vcf}"

    publishDir "${params.output}/mbovpan_results/assemblies", mode: 'copy'

    conda "${workflow.projectDir}/envs/consensus.yaml"

    input:
    tuple val(sample_id), path(filtered_vcf)

    output:
    path("${sample_id}.consensus.fasta"), emit: fasta_ch

    script:
    """
    set -euxo pipefail

    echo "Starting psuedo_assembly process for ${filtered_vcf}"
    bgzip ${filtered_vcf}

    echo "Indexing"
    bcftools index ${filtered_vcf}.gz

    echo "Generating consensus sequence"
    bcftools consensus -f ${ref} ${filtered_vcf}.gz > ${sample_id}.dummy.fasta

    echo "Replacing sequence name in consensus sequence"
    sed 's|LT708304.1|${sample_id}|g' ${sample_id}.dummy.fasta > ${sample_id}.consensus.fasta
    """
}



process iqtree_phylo {
    publishDir "${params.output}/mbovpan_results/phylogeny", mode: 'copy'

    conda "${workflow.projectDir}/envs/iqtree.yaml"

    input:
    path(aln)

    output:
    path("*"), emit: phylo_ch 
    path("phylo_out.txt"), emit: statistics_ch

    script:
    """
    cat *.fasta > mbovpan_align.fasta
    cat ${workflow.projectDir}/ref/mbov_reference.fasta >> mbovpan_align.fasta
    snp-sites -o mbovpan_align.snp_only.fasta mbovpan_align.fasta
    iqtree -s mbovpan_align.snp_only.fasta -m MFP -nt ${params.threads} -bb 1000 -pre mbovpan_align -o "LT708304.1"
    echo "file created to make statistics file build conda env last" > phylo_out.txt
    """
}

process stats {
    publishDir "${params.output}/mbovpan_results/statistics"

    conda "$workflow.projectDir/envs/statistics.yaml"

    input:
    tuple val(sample_id1), path(fastq1), path(fastq2), val(sample_id2), path(bam_file), path(vcf_file)
    path(dummy_variable)

    output:
    path("${sample_id1}.stats"), emit: output_stat_ch

    script:
    """
    python $workflow.projectDir/scripts/statistics.py ${fastq1} ${fastq2} ${bam_file} ${vcf_file} > ${sample_id1}.stats
    """
}

process assembly {
    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/assembly/${sample_id}", mode: 'copy'

    conda "${workflow.projectDir}/envs/spades.yaml"

    input:
    tuple val(sample_id), path(trim_one), path(trim_two)

    output:
    tuple val(sample_id), path("${sample_id}.scaffold.fasta"), emit: assembly_ch
    tuple val(sample_id), path("${sample_id}.contigs.fasta"), emit: contigs_ch

    script:
    """
    mkdir ${sample_id}
    spades.py -1 ${trim_one} -2 ${trim_two} --careful -o ${sample_id} -t ${params.threads}
    cp ${sample_id}/contigs.fasta ${sample_id}.contigs.fasta
    cp ${sample_id}/scaffolds.fasta ${sample_id}.scaffold.fasta
    """
}

process quast_assembly {
    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/quast/${sample_id}", mode: 'copy'

    conda "${workflow.projectDir}/envs/quast.yaml"

    input:
    tuple val(sample_id), path (assembly_fasta)

    output:
    path "*", emit: quast_assembly_output

    script:
    """
    quast -o raw ${assembly_fasta} -r ${ref} -t ${task.cpus}
    """
}

process bwa_index {
    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/assembly/${sample_id}", mode: 'copy'

    conda "${workflow.projectDir}/envs/samtools.yaml"

    input:
    tuple val(sample_id), path(contigs_ch)

    output:
    tuple val(sample_id), path("*"), emit: bwa_index_out

    script:
    """
    echo "Sample ID: ${sample_id}"
    echo "Contigs FASTA: ${contigs_ch}"
    
    # Index the assembly FASTA file
    bwa index ${contigs_ch}
    """
}

process bwa_mem {
    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/assembly/${sample_id}", mode: 'copy'

    conda "${workflow.projectDir}/envs/samtools.yaml"

    input:
    tuple val(sample_id), path(bwa_index_out), path(contigs_ch), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}.contigs.sorted.bam"), emit: sorted_bam

    script:
    """
    echo "Sample ID: ${sample_id}"
    echo "Assembly FASTA: ${contigs_ch}"
    echo "Read 1: ${read1}"
    echo "Read 2: ${read2}"
    
    # Align reads to the indexed assembly FASTA and generate sorted BAM file
    bwa mem -t ${task.cpus} ${contigs_ch} ${read1} ${read2} | samtools view -Sb | samtools sort -@${task.cpus} -o ${sample_id}.contigs.sorted.bam
    """
}

process samtools_index {
    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/assembly/${sample_id}", mode: 'copy'

    conda "${workflow.projectDir}/envs/samtools.yaml"

    input:
    tuple val(sample_id), path(sorted_bam)

    output:
    tuple val(sample_id), path("${sorted_bam}.bai"), emit: bai

    script:
    """
    samtools index ${sorted_bam}
    """
}

process pilon {
    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/pilon", mode: 'copy'

    conda "bioconda::pilon"

    input:
    tuple val(sample_id), path (contigs_ch), path (sorted_bam), path (bai)

    output:
    tuple val(sample_id), path("${sample_id}.contigs.pilon.fasta"), emit: pilon_fasta

    script:
    """
    echo "Assembly FASTA: ${contigs_ch}"
    echo "Sorted BAM: ${sorted_bam}"
    echo "BAM Index: ${bai}"
    
    pilon --genome ${contigs_ch} --bam ${sorted_bam} --output ${sample_id}.contigs.pilon
    """
}

process quast_pilon {
    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/quast/${sample_id}", mode: 'copy'

    conda "${workflow.projectDir}/envs/quast.yaml"

    input:
    tuple val(sample_id), path (pilon_fasta)

    output:
    path "*", emit: quast_pilon_results

    script:
    """
    quast -o polished ${pilon_fasta} -r ${ref} -t ${task.cpus}
    """
}

process annotate {
    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/annotations/${sample_id}", mode: 'copy'

    conda "${workflow.projectDir}/envs/prokka.yaml"

    input:
    tuple val(sample_id), path(pilon_fasta)

    output:
    path("${sample_id}.gff"), emit: annotate_ch
    path("prokka_output/*")

    script:
    """
    echo "Pilon FASTA: ${pilon_fasta}"
    echo "Sample ID: ${sample_id}"
    
    prokka --centre X --compliant ${pilon_fasta} --proteins ${workflow.projectDir}/ref/af2122.97.gb --prefix ${sample_id} --outdir prokka_output
    mv prokka_output/${sample_id}.gff ${sample_id}.gff
    """
}

process panaroo {
    publishDir "${params.output}/mbovpan_results/pangenome"

    conda "${workflow.projectDir}/envs/panaroo.yaml"

    input:
    path gff

    output:
    path("*"),  emit: panaroo_ch

    script:
    """
    panaroo -i *.gff -o ./ -t ${params.threads} -a core --core_threshold 0.98 --clean-mode strict --remove-invalid-genes
    """
}

process iqtree_core {
    publishDir "${params.output}/mbovpan_results/phylogeny", mode: 'copy'

    conda "${workflow.projectDir}/envs/iqtree.yaml"

    input:
    path(input)

    output:
    path("*"), emit: iqtreecore_ch

    script:
    """
    snp-sites -o mbovpan_align.snp_only.fasta core_gene_alignment.aln
    iqtree -s mbovpan_align.snp_only.fasta -m MFP -nt ${params.threads} -bb 1000 -pre mbovis_core 
    """
}

process filter_pan {
    publishDir "${params.output}/mbovpan_results/filtered_pan", mode: 'copy'

    conda "bioconda::blast=2.9.0 pandas"

    input:
    path(input)

    output:
    path('mbovis_filtered_cogs.csv'), emit: filtered_pan_ch

    script:
    """
    # first determine the genes that are present and use blast to find their true annotation
    blastn -query pan_genome_reference.fa -max_target_seqs 1 -db ${workflow.projectDir}/ref/mbovis_reference -out mb.out -outfmt "6 delim=, qseqid sseqid length qstart qend sstart send qlen slen"
    echo "qseqid,sseqid,length,qstart,qend,sstart,send,qlen,slen" > heading.txt 
    cat heading.txt mb.out >> mb.out.csv

    # once mb.out is created, we can use python to 1) create length %, 2) filter by 75% or higher, and 3) reduce repetitive blast results
    python ${workflow.projectDir}/scripts/blast_result_filter.py gene_presence_absence.csv
    """
}

process multiqc {
    publishDir "${params.output}/mbovpan_results/statistics", mode: 'copy'

    conda "${workflow.projectDir}/envs/multiqc.yaml"

    input:
    path pre
    path post

    output:
    path("mbovpan_report*"), emit: multiqc_ch

    script:
    """
    export LC_ALL=en_US.utf8
    export LANG=en_US.utf8
    multiqc -n mbovpan_report .
    """
}

process mbovis_verification {
    publishDir "${params.output}/mbovpan_results/lineage_info", mode: 'copy'

    conda "${workflow.projectDir}/envs/pandas.yaml"

    input:
    path spoligotype_info
    path lineage_info

    output:
    path("mbovpan_lineage_info.csv"), emit: mbovis_verification_ch

    script:
    """
    python ${lineage_table} > mbovpan_lineage_info.csv
    """
}