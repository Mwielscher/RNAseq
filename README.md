## Table of contents  
1. [About this Repository](About-this-Repository)  
2. [readQC, alginment and read counting ](readQC,-alginment-and-read-counting)  
3. [Analysis of RNAseq count matrix ](Analysis-of-RNAseq-count-matrix )  
4. [Miscellaneous](Miscellaneous)  

## About this Repository  
This is collection of scripts to analyse RNA sequencing data from start (FASTQ or unaligned BAM file) to end (differential gene expression, heatmap, enrichment analysis). Some chunks of code are adopted from this [limma tutorial](https://bioconductor.org/packages/release/workflows/html/RNAseq123.html)
<br/><br/>  
## readQC, alginment and read counting  
These scripts were developed on a HPC with a [SLURM scheduler](https://slurm.schedmd.com/quickstart.html). Those are array jobs, however the backbone of the scripts should work on any linux based operating system. There is a script to align single end reads and a seperate one for paired end reads. The data per sample can distributed over several input FASTQ files. The scripts will run [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [qualimap](http://qualimap.conesalab.org) alongside with [STAR](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3530905/) based alignment of the reads. QC metrics will be collected in a seperate script. Finally, we count the reads using [Rsubread](https://pubmed.ncbi.nlm.nih.gov/30783653/).  
>* this [script](FASTQ_to_COUNTs/prep_GRCh38_refGenome_for_STAR.sh) helps to index the reference genome for STAR 
>* this [script](FASTQ_to_COUNTs/STAR_alignment_PE.sh) performs QC steps and STAR alignment of in reads for paired end reads and [this one](FASTQ_to_COUNTs/STAR_alignment_SE.sh) for sinlge end reads 

## Analysis of RNAseq count matrix  

## Miscellaneous
