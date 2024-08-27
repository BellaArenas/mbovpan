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

//reads = Channel.fromFilePairs("${params.input}*{1,2}*.f*q*", size: 2, flat: true).ifEmpty { error "Cannot find the read files" }

// mark_dups_ch = Channel.fromPath("/groups/lilianasalvador/barenas/scripts/mbovpan/process_09_new/mbovpan_results/readmapping/*.bam")
// mark_dups_bai_ch = Channel.fromPath("/groups/lilianasalvador/barenas/scripts/mbovpan/process_09_new/mbovpan_results/readmapping/*.bam.bai")

// freebayes_ch = Channel.fromPath('/groups/lilianasalvador/barenas/scripts/mbovpan/process_09_new/mbovpan_results/variant_calling/*.vcf')

// assembly_ch = Channel.fromPath("/groups/lilianasalvador/barenas/scripts/mbovpan/process_13/*scaffold.fasta")

// panaroo_ch = Channel.fromPath("/groups/lilianasalvador/barenas/scripts/mbovpan/process_17/mbovpan_results/pangenome/*")

filter_ch = Channel.fromPath("/groups/lilianasalvador/barenas/scripts/mbovpan/process_09_new/mbovpan_results/variant_calling/*.vcf")
    .filter { it.name.endsWith('.filtered.vcf') }

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

    // spotyping(reads)
    // pre_fastqc(spotyping.out.spoligo_ch)
    // fastp(spotyping.out.spoligo_ch)
    // post_fastqc(fastp.out.trimmed_reads)
    // lineage(fastp.out.trimmed_reads)
    if(run_mode == "snp" || run_mode == "all"){
        // read_map(fastp.out.trimmed_reads)
        // mark_dups(read_map.out.bam)
        // freebayes(mark_dups.out.nodup_bam, mark_dups.out.nodup_bai)
        // vcf_filter(freebayes.out.freebayes_ch)
        // psuedo_assembly(vcf_filter.out.filter_ch)
        psuedo_assembly(filter_ch)
        iqtree_phylo(psuedo_assembly.out.fasta_ch.collect())
        // stats(psuedo_assembly.out.fasta_ch, iqtree_phylo.out.statistics_ch)
    }

    if(run_mode == "pan" || run_mode == "all"){
        assembly(fastp.out.trimmed_reads)
        quast(assembly.out.assembly_ch)
        annotate(assembly.out.assembly_ch)
        panaroo(annotate.out.annotate_ch.collect())
        iqtree_core(panaroo.out.panaroo_ch.collect())
        filter_pan(panaroo.out.panaroo_ch.collect())
        multiqc(pre_fastqc.out.fastqc_ch1.collect(),post_fastqc.out.fastqc_ch2.collect())
        mbovis_verification(spotyping.out.spoligo_multi.collect(), lineage.out.tbprofile_ch.collect())
    }

}

process spotyping {
    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/spotyping", mode: 'copy'

    conda "${workflow.projectDir}/envs/spotyping.yaml"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${read1.baseName}.fastq"), path("${read2.baseName}.fastq"), emit: spoligo_ch
    path("${read1.baseName - ~/_1/}.out"), emit: spoligo_multi

    script:
    """
    echo "Starting spotyping process for sample: ${sample_id}"
    echo "Read1: ${read1}"
    echo "Read2: ${read2}"

    if [ -z "${read1}" ] || [ -z "${read2}" ]; then
        echo "Error: One of the input files is null"
        exit 1
    fi

    python3 ${workflow.projectDir}/scripts/SpoTyping/SpoTyping.py ${read1} ${read2} -o ${read1.baseName - ~/_1/}.out > stdout.txt 2> stderr.txt
    if [ \$? -ne 0 ]; then
        echo "Error in SpoTyping.py script"
        cat stderr.txt
        exit 1
    fi

    echo "spotyping process completed successfully for sample: ${sample_id}"
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

    publishDir "${params.output}/mbovpan_results/read_trimming", mode: 'copy'
    
    conda 'bioconda::fastp'

    input:
    tuple val(sample_id), path(read_one), path(read_two)

    output:
    tuple val(sample_id), path("${read_one.baseName.replaceAll('_1$', '')}_trimmed_R1.fastq"), path("${read_two.baseName.replaceAll('_2$', '')}_trimmed_R2.fastq"), emit: trimmed_reads

    script:
    """
    fastp -w ${params.threads} -q 30 --detect_adapter_for_pe -i ${read_one} -I ${read_two} -o ${read_one.baseName.replaceAll("_1\$", "")}_trimmed_R1.fastq -O ${read_two.baseName.replaceAll("_2\$", "")}_trimmed_R2.fastq    
    """
}

process post_fastqc {
    tag "${sample_id}"

    publishDir "${params.output}/mbovpan_results/fastqc", mode: 'copy'

    conda 'bioconda::fastqc'

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

    input:
    tuple val(sample_id), path(nodup_bam)
    tuple val(sample_id), path(nodup_bai)

    output:
    path("${nodup_bam.baseName.replaceAll('.nodup', '')}.vcf"), emit: freebayes_ch
 
    script:
    """
    freebayes-parallel ${range} ${params.threads} -f ${ref} ${nodup_bam} > ${nodup_bam.baseName.replaceAll('.nodup', '')}.vcf
    """
}
    


process vcf_filter {
    tag "${vcf.baseName}"

    publishDir "${params.output}/mbovpan_results/variant_calling", mode: 'copy'

    conda "${workflow.projectDir}/envs/vcflib.yaml"

    input:
    path(vcf)

    output:
    path("${vcf.baseName}.filtered.vcf"), emit: filter_ch

    script:
    """
    vcffilter -f "QUAL > ${qual}" ${vcf} | vcffilter -f "DP > ${depth}" | vcffilter -f "MQM > ${mapq}" |  vcffilter -f "TYPE = snp" | bedtools intersect -header -a - -b ${workflow.projectDir}/ref/pe_ppe_regions.gff3 -v > ${vcf.baseName}.filtered.vcf
    """
}
process psuedo_assembly {
    tag "${vcf.baseName}"

    publishDir "${params.output}/mbovpan_results/assemblies", mode: 'copy'

    conda "${workflow.projectDir}/envs/consensus.yaml"

    input:
    path(vcf)

    output:
    path("${vcf.baseName.replaceAll('.filtered', '')}.consensus.fasta"), emit: fasta_ch

    script:
    """
    set -euxo pipefail

    echo 'Compressing VCF file'
    bgzip ${vcf}
    
    echo 'Indexing compressed VCF file'
    bcftools index ${vcf}.gz
    
    echo 'Creating output directory if it does not exist'
    mkdir -p ${params.output}/mbovpan_results/assemblies
    
    echo 'Generating initial mismatch report'
    bcftools norm --check-ref s --fasta-ref ${ref} -Ov ${vcf}.gz > ${vcf.baseName}.norm.vcf 2> ${params.output}/mbovpan_results/assemblies/${vcf.baseName.replaceAll('.filtered', '')}_initial_mismatch_report.txt
    
    echo 'Force fixing mismatches'
    bcftools norm --check-ref s --fasta-ref ${ref} -Ov ${vcf}.gz > ${vcf.baseName}.norm.vcf
    
    echo 'Generating post-normalization mismatch report'
    bcftools norm --check-ref s --fasta-ref ${ref} -Ov ${vcf.baseName}.norm.vcf > /dev/null 2> ${params.output}/mbovpan_results/assemblies/${vcf.baseName.replaceAll('.filtered', '')}_post_normalization_mismatch_report.txt
    
    echo 'Compressing normalized VCF file'
    bgzip ${vcf.baseName}.norm.vcf
    
    echo 'Indexing normalized VCF file'
    bcftools index ${vcf.baseName}.norm.vcf.gz
    
    echo 'Generating consensus sequence using bcftools consensus'
    bcftools consensus -f ${ref} ${vcf.baseName}.norm.vcf.gz > ${vcf.baseName.replaceAll('.filtered', '')}.dummy.fasta
    
    echo 'Replacing sequence name in consensus sequence'
    sed "s|LT708304.1|${vcf.baseName.replaceAll('.filtered', '')}|g" ${vcf.baseName.replaceAll('.filtered', '')}.dummy.fasta > ${vcf.baseName.replaceAll('.filtered', '')}.consensus.fasta
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
    iqtree -s mbovpan_align.snp_only.fasta -m MFP -nt ${task.cpus} -bb 1000 -pre mbovpan_align -o "LT708304.1"
    echo "file created to make statistics file build conda env last" > phylo_out.txt
    """
}

process stats {
    publishDir "${params.output}/mbovpan_results/statistics", mode: 'copy'

    conda "${workflow.projectDir}/envs/statistics.yaml"

    input:
    path(nec_files) 
    path(dummy_variable)

    output:
    file("${nec_files[0].baseName}.stats"), emit: output_stat_ch

    script:
    """
    python ${workflow.projectDir}/scripts/statistics.py ${nec_files[0]} ${nec_files[1]} ${nec_files[2]} ${nec_files[3]} > ${nec_files[0].baseName}.stats
    """
}

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

process assembly {
    tag "${sample_id}"

    publishDir "${params.output}", mode: 'copy'

    conda "bioconda::spades python=3.9"

    input:
    tuple val(sample_id), path(trim_one), path(trim_two)

    output:
    path("${trim_one.baseName.replaceAll('_trimmed_R.*', '')}.scaffold.fasta"), emit: assembly_ch

    script:
    """
    mkdir ${trim_one.baseName.replaceAll('_trimmed_R.*', '')}
    spades.py -1 ${trim_one} -2 ${trim_two} --careful -o ${trim_one.baseName.replaceAll('_trimmed_R.*', '')} -t ${params.threads} --only-assembler
    cd ${trim_one.baseName.replaceAll('_trimmed_R.*', '')}
    mv scaffolds.fasta ../${trim_one.baseName.replaceAll('_trimmed_R.*', '')}.scaffold.fasta
    """
}

process quast {
    tag "${assemblies.baseName.replaceAll('.scaffold', '')}"

    publishDir "${params.output}/mbovpan_results/statistics"

    conda "${workflow.projectDir}/envs/quast.yaml"
    
    input:
    path (assemblies)
    
    output:
    path ("*"), emit: quast_output
    
    script:
    """
    quast -o assembly_stats ${assemblies} -r ${ref} -t ${params.threads}
    """
}

process annotate {
    tag "${assembly.baseName.replaceAll('.scaffold', '')}"

    publishDir "${params.output}/mbovpan_results/annotations"
    
    conda "$workflow.projectDir/envs/prokka.yaml"

    input:
    path assembly
    
    output:
    path "${assembly.baseName.replaceAll('.scaffold', '')}.annot.gff", emit: annotate_ch
    
    script:
    """
    prokka --outdir ./${assembly.baseName.replaceAll('.scaffold', '')} --prefix ${assembly.baseName.replaceAll('.scaffold', '')}.annot ${assembly}
    ls ./${assembly.baseName.replaceAll('.scaffold', '')}  # Add this line
    cp ./${assembly.baseName.replaceAll('.scaffold', '')}/${assembly.baseName.replaceAll('.scaffold', '')}.annot.gff ./
    """
}

process panaroo {
    publishDir "${params.output}/mbovpan_results/pangenome"

    conda "${workflow.projectDir}/envs/panaroo.yaml"

    input:
    path gff

    output:
    path (input),  emit: panaroo_ch

    script:
    """
    panaroo -i *.gff -o ./ -t ${params.threads} -a core --core_threshold 0.98 --clean-mode strict
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
    iqtree -s mbovpan_align.snp_only.fasta -m MFP -nt ${task.cpus} -bb 1000 -pre mbovis_core 
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
