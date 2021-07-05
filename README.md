## Table of contents  
1. [About this Repository](About-this-Repository)  
2. [readQC, alginment and read counting ](readQC,-alginment-and-read-counting)  
3. [Analysis of RNAseq count matrix ](Analysis-of-RNAseq-count-matrix )  
4. [Miscellaneous](Miscellaneous)  

## About this Repository  
This is a collection of scripts to analyse RNA sequencing data from start (FASTQ or unaligned BAM file) to end (differential gene expression, heatmap, enrichment analysis). Some chunks of code were adopted from this [limma tutorial](https://bioconductor.org/packages/release/workflows/html/RNAseq123.html)
<br/><br/>  
## readQC, alginment and read counting  
These scripts were developed on a HPC with a [SLURM scheduler](https://slurm.schedmd.com/quickstart.html). Those are array jobs, however the backbone of the scripts should work on any linux based operating system. There is a script to align single end reads and a seperate one for paired end reads. The data per sample can distributed over several input FASTQ files. The scripts will run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [qualimap](http://qualimap.conesalab.org) alongside with [STAR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) based alignment of the reads. QC metrics will be collected in a seperate script. Finally, we count the reads using [Rsubread](https://pubmed.ncbi.nlm.nih.gov/30783653/).  
>* this [script](FASTQ_to_COUNTs/prep_GRCh38_refGenome_for_STAR.sh) helps to index the reference genome for STAR 
>* this [script](FASTQ_to_COUNTs/STAR_alignment_PE.sh) performs QC steps and STAR alignment of reads for paired end reads. For single end reads use [this one](FASTQ_to_COUNTs/STAR_alignment_SE.sh). The scripts create folders for fastQC and qualimap analysis also log files are processed to be easy accessible for R.  
>* this [script](FASTQ_to_COUNTs/2_collect_QC_input.sh) will create a text file with an overview of colleted QC metrics and clean up the analysis directory. It will put all bam files into one folder for the next step of the analsis.  
>* this [script](FASTQ_to_COUNTs/3_PE_count_reads.sh) counts the reads and creates the input matrix for subsequent analysis with DESeq2 or voom limma. 

## Analysis of RNAseq count matrix  
The [DESeq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8) analysis pipeline does not require a lot of coding effort. It is couple of R-commands from data ingest to differential gene expression results. Thus, this [script](analyseCOUNTs/DESeq2/deseq_differential_expression.R) also includes initial data wrangling to get the data in shape for the DESeq2 data object and conversion of EntrezIDs to NCBI Symbol. You can use the vst matrix in case you want to do analyse the data outside of DESeq2.  
<br/>
The [limma voom](https://genomebiology.biomedcentral.com/articles/10.1186/gb-2014-15-2-r29) analysis pipeline as implemented in this [scirpt](analyseCOUNTs/limma_voom/voom_incl_batch_removal.R) includes several QC plots to find outlying patterns in your count data and a batch removal procedure of a known batches in the dataset. It includes setting up a design and a contrast matrix for more multifacetted datasets alongside with a wrapper to convert from EntrezID to Gene Symbol and export differential gene expression results.

## Miscellaneous
