#!/bin/bash

module load intel/16 htslib/1.3.2 samtools/1.3.1


dat='/home/lv71395/mwielsch/binfl/BLUEPRINT/analysis/'
out='/home/lv71395/mwielsch/binfl/BLUEPRINT/analysis/result/'
#id=$(awk -v var="$SLURM_ARRAY_TASK_ID" 'FNR==var{print $1}' ${dat}PE_samples)
# ------------   data wrangling for mulitple PE fastq files
id=$(awk -v var="$SLURM_ARRAY_TASK_ID" 'FNR==var{print $1}' ${dat}PE_samples_feb)

grep "${id}" ${dat}input_align_2 | awk '{ print $2}' > ${dat}test_${id}
sort ${dat}test_${id} > ${dat}test1_${id}
check_file=$(awk 'NR==1{print $1}' ${dat}test1_${id})
awk '0 == (NR + 1) % 2' <${dat}test1_${id} | tr "\n" "," | sed 's/,$//'>${dat}odd_line_${id}
awk '0 == NR % 2' <${dat}test1_${id} | tr "\n" "," | sed 's/,$//'>${dat}even_line_${id}
files=$(paste -d" " ${dat}odd_line_${id} ${dat}even_line_${id})

rm ${dat}test_${id}
rm ${dat}test1_${id}
rm ${dat}odd_line_${id}
rm ${dat}even_line_${id}


## QC fast Q files
/home/lv71395/mwielsch/FastQ-Screen-0.14.1/./fastq_screen ${check_file} --subset 1000000 --outdir ${out}
/home/lv71395/mwielsch/FastQC/./fastqc ${check_file} --extract --outdir=${out}
sed -i '1,2d' ${out}${id}_screen.txt 

###  align reads to GRCH38 - sort and index

/home/lv71395/mwielsch/./STAR --runThreadN 12 --genomeDir /home/lv71395/mwielsch/binfl/STAR_genome\
 --readFilesCommand zcat --readFilesIn ${files} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 31000000000\
 --outFileNamePrefix ${out}${id}_

samtools index -b ${out}${id}_Aligned.sortedByCoord.out.bam 


module purge
module load intel/18 intel-mkl/2018 pcre2/10.35 R/4.0.2

/home/lv71395/mwielsch/qualimap_v2.2.1/./qualimap rnaseq -bam ${out}${id}_Aligned.sortedByCoord.out.bam\
 --java-mem-size=15g -gtf /home/lv71395/mwielsch/binfl/STAR_genome/input/Homo_sapiens.GRCh38.101.gtf\
 -s -outdir ${out} -outformat PDF -p non-strand-specific -outfile ${id}
 



# make STAR LOG readable for R
sed '/:/d' ${out}${id}_Log.final.out > ${out}${id}_summary
sed -i '/^$/d' ${out}${id}_summary

###  ------------    header QC file
#Genome Reads_processed Unmapped Unmapped_inPERC One_hit_one_genome One_hit_one_genome_inPERC Multiple_hits_one_genome Multiple_hits_one_genome_inPERC One_hit_multiple_genomes One_hit_multiple_genomes_in_PERC Multiple_hits_multiple_genomes Multiple_hits_multiple_genomes_in_PERC


