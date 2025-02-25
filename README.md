[![image search api](https://user-images.githubusercontent.com/40697188/193487621-a4b91a1c-19b6-42df-9e63-7fcff0658be0.png)](https://github.com/noahaus)

# mbovpan
Mbovpan is a nextflow bioinformatic pipeline for _Mycobacterium bovis_ pangenome analysis. The goal of Mbovpan is to make the insights from _M. bovis_ genomics easily accesible to reseachers.  

Mbovpan can be ran in three separate modes: SNP mode for only inferring Single Nucleotide Polymorphisms, PAN mode for assessing gene presence absences, or ALL to do a full analysis (default). Below you can witness the main workflow encapsulated by mbovpan. 
![image](https://github.com/salvadorlab/mbovpan/assets/40697188/315e9533-1567-48c1-aa0c-f1b5c12e2589)


### Installation  

The first step would be to download the pipeline using the following git command.
```
$git clone https://github.com/salvadorlab/mbovpan.git
```

We encourage the user to add the 'mbovpan' directory to the PATH variable. This makes running the pipeline much easier due to not needing to specify the absolute path of the 'mbovpan.nf' file.

```
#if NOT in the directory, give the absolute path
export PATH=$PATH:/path/to/mbovpan

#if already in the directory
export PATH=$PATH:$(pwd)
```  

After downloading, the user will need to create the mbovpan environment that will make it possible to run the pipeline. 


```
#create the environment with mamba as the first package, then activate
$conda create -p </path/to/>mbovpan 
$conda activate </path/to/>mbovpan
(</path/to/>mbovpan)$conda install -c conda-forge mamba
(</path/to/>mbovpan)$conda install -c bioconda -c conda-forge nextflow=22.10.6
(</path/to/>mbovpan)$conda install -c bioconda -c conda-forge pandas
(</path/to/>mbovpan)$conda install -c bioconda -c conda-forge panaroo

#now test that everything downloaded appropriately with a simple help command
(</path/to/>mbovpan)$ nextflow run mbovpan --help 
```

### Quickstart  
#### Downloading example sequence data 
```
#create a new directory that will house NGS files, navigate to the newly created file
(</path/to/>mbovpan)$mkdir mbovis_input
(</path/to/>mbovpan)$cd mbovis_input

#using sratoolkit download 5 M. bovis sequences extracted from United Kingdom badgers
(</path/to/>mbovpan)$fasterq-dump --verbose --split-3 SRR10482974
(</path/to/>mbovpan)$fasterq-dump --verbose --split-3 SRR10482944
(</path/to/>mbovpan)$fasterq-dump --verbose --split-3 ERR11893527
(</path/to/>mbovpan)$fasterq-dump --verbose --split-3 SRR23174187
(</path/to/>mbovpan)$fasterq-dump --verbose --split-3 SRR23174144

#once the downloads are complete, exit into the previous directory
(</path/to/>mbovpan)$cd ../
```
#### Example run
```
(</path/to/>mbovpan)$nextflow run mbovpan --input ./mbovis_input/ --run snp --output ./ 
```
In this command, **'nextflow run'** is the command used to look at and execute the pipeline instructions in **'mbovpan'**. This looks into the mbovpan directory and runs the workflow instructions that are  contained in 'main.nf'

With the test data already downloaded, the parameter **'--input ./mbovis_input/'** looks for if paired end sequences are present in the input directory to initiate the pipeline. Raw sequence data is analyzed to generate the spoligotyping pattern of the isolates, and if the pattern that is generated does not present as _M. bovis_, the file will be filtered out from analysis. 

**'--run snp'** signifies what analysis mode mbovpan utilizes. **'snp'** mode maps the paired end files to the reference genome while, **'pan'** mode creates de novo genomes from scratch. if no option is supplied, the pipeline will run both. 

**--output ./** stipulates where the output directory "mbovpan_results" will be created

### Quickstart with Docker
#### Building and running Dockerfile
```
# Download the Dockerfile in a directory where you would want to run mbovpan
docker build -t mbovpan:image .
docker run -it mbovpan:image /bin/bash
```
#### Example run
```
(base)# nextflow run ./main.nf --input ./mbovis_input/ --run snp --output ./ --threads 16
```
### Additional Usages

```
# boost the number of threads utilized
(</path/to/>mbovpan)$nextflow run mbovpan/mbovpan.nf --input ./mbovpan/seqs/ --run snp --output ./ --threads 16

# modulate the minimum quality and maximum SNP depth required
(</path/to/>mbovpan)$nextflow run mbovpan/mbovpan.nf --input ./mbovpan/seqs/ --run snp --output ./ --qual 20 --depth 20

# Using most of the parameters that mbovpan has to offer to decipher the pangenome
(</path/to/>mbovpan)$nextflow run mbovpan/mbovpan.nf --input ./mbovpan/seqs/ --run snp --output ./ --qual 100 --depth 5 --mapq 50 --threads 30 --run pan
```

### Inputs

mbovpan requires as input paired end FASTQ files originating from _Mycobacterium bovis_. mbovpan runs the tool **spotyping** that can use the reads to deduce what MTBC member the sequences originate from. 

### Outputs

_from the mbovpan manuscript_

 #### all mode - Maximum Likelihood Phylogenies  
For each mode in Mbovpan, a maximum likelihood phylogeny based on SNPs is produced in order to visualize the genomic variation amongst the isolates. If run in ‘pan’ mode, the phylogeny produced will be based on SNPs that were located in core genes inferred from Roary. Otherwise, if the pipeline is run in ‘snp’ mode, SNPs will be inferred directly from the sequences aligned to the reference genome. Given user supplied metadata, the phylogenies will be annotated with the metadata at the tips of the tree. Phylogenies will be generated using Iqtree with 1000 Ultrafast bootstraps for nodal support, and the nucleotide substitution model being inferred through use of the implemented ModelFinder approach (Kalyaanamoorthy et al., 2017; Hoang et al., 2018; Minh et al., 2020).  
 
#### all mode - Spoligotyping and Lineage Classification   
Computation of the spoligotype patterns and lineage of M. bovis isolates are oftentimes important aspects of outbreak investigations and analysis of global variation (Fuente et al., 2015; Zimpel et al., 2017). Mbovpan will output data tables as CSV files that link an isolate with their spoligotyping pattern and associated Lineage information through the use of SpoTyping and TB-Profiler (Xia et al., 2016; Phelan et al., 2019), respectively.

#### snp mode – Read mapping files and Variant Calling Format  
If the user is concerned with computing the use of SNPs only, mbovpan will produce files associated with read-mapping to the integrated reference genome (as well as a duplicate read removal version in BAM format). Freebayes is then used to infer SNP sites and then output the results in a VCF file. A subsequent filtering of SNPs occurs where the user specifies threshold values for mapping quality, SNP depth, and SNP quality.  
 
#### pan mode - Accessory Genome PCA and PanGWAS   
To further assess the similarity of the accessory genome based on user provided metadata, mbovpan produces a Principal Component Analysis score plot that will reduce the large dimensionality of the accessory genome. Given multiple fields of metadata in
CSV format as input, there will be as many PCA plots as metadata available. Using the same metadata, Mbovpan will also produce a genomic variation visualization, and a pangenome Genome Wide Association Study (panGWAS) implemented through the software Scoary (Brynildsrud et al., 2016) to check if accessory genes can discriminate the inputted binary traits. The program requires as input the pangenome inferred through mbovpan alongside a traits file that links a particular trait with an isolate. This analysis will output genes that have a significant association (based on a Chi-Squared test) to the trait of interest that is not based on lineage alone. These results are returned as a CSV formatted table.  

#### pan mode – M. bovis virulence gene presence/absence matrix  
If the ‘pan’ mode is chosen, after the pangenome is inferred from Roary, Mbovpan takes the gene presence absence data and only keeps the accessory genes that are present in a predetermined list of M. bovis virulent genes. This results in the creation of an output file with only these genes, and additionally a visualization to show how the patterns in virulence gene presence and absence content relates to the user supplied metadata. The matrix will be paired with a hierarchical clustering dendrogram made from applying the Ward clustering algorithm on a Jaccard Similarity matrix generated from the virulent gene presence/absence matrix (Ceres et al., 2022).   
 


### Help
```
    M B O V P A N (v0.1)    
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
```

### Authors
Noah Legall, Ph. D. (conceptualization)
Bella Arenas (maintaince) 
Liliana C. M. Salvador, Ph. D. (supervisor)
 


