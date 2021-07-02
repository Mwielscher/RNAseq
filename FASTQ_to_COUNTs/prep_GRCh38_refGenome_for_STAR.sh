#!/bin/bash
#BATCH -J STAR_ref_build
#SBATCH -N 1

/home/lv71395/mwielsch/./STAR --runThreadN 12 --runMode genomeGenerate --genomeDir /home/lv71395/mwielsch/work/STAR_genome\
 --genomeFastaFiles /home/lv71395/mwielsch/work/STAR_genome/input/Homo_sapiens.GRCh38.dna.primary_assembly.fa\
 --sjdbGTFfile /home/lv71395/mwielsch/work/STAR_genome/input/Homo_sapiens.GRCh38.101.gtf
