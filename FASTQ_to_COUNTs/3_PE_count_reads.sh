#!/bin/bash
#BATCH -J COUNT
#SBATCH -N 1

module load intel/18 intel-mkl/2018 pcre2/10.35 R/4.0.2

in='/home/lv71395/mwielsch/binfl/BLUEPRINT/analysis/result/bam_files/'
out='/home/lv71395/mwielsch/binfl/BLUEPRINT/analysis/result/count_matrix/'

R --vanilla << EOF

library(Rsubread)

pe=read.table("/home/lv71395/mwielsch/binfl/BLUEPRINT/analysis/PE_for_counting")
gtffile = file.path("/home/lv71395/mwielsch/work/STAR_genome/input/Homo_sapiens.GRCh38.101.gtf")
#filenames = file.path(paste0("${in}",list.files(path = "${in}",pattern=".bam$")))

filenames = file.path(paste0("${in}",pe\$V1,"_Aligned.sortedByCoord.out.bam"))

fc = featureCounts(files=filenames,annot.inbuilt = "hg38",
                    isPairedEnd=TRUE,nthreads=12)

save(fc,file="${out}BLUEPRINT_PE_feature_count_object.RData")

write.table(as.data.frame(fc\$counts),file="${out}BLUEPRINT_PE_count_matrix.txt",col.names=T,row.names=T,quote=F)

write.table(as.data.frame(fc\$annotation),file="${out}BLUEPRINT_PE_feature_annot.txt",col.names=T,row.names=T,quote=F)


EOF
