#!/bin/bash



module load intel/18 intel-mkl/2018 pcre2/10.35 R/4.0.2

out='/home/lv71395/mwielsch/binfl/BLUEPRINT/analysis/result/'
dat='/home/lv71395/mwielsch/binfl/BLUEPRINT/analysis/'

R --vanilla <<EOF

infiles = system(paste0("awk '{print \$1}' ${dat}all_samples_ok"),intern=T)
metric=c("ID", "Number_of_inputREADS","average_mapped_READ_length","unmapped_READs","Uniquely_mapped_reads","reads_mapped_to_multiple_loci_perc","rRNA_content","mtRNA_content",
	"mouse_content","Ecoli_content","Adapter_content","FASTQC_basic_STATs","FASTQC_perBase_Quality","GC_content_perRead",
	"Mismatch_rate_per_base_perc","chimeric reads_perc","reads_mapped_to_too_many_loci")
res=as.data.frame(matrix(NA,nrow=length(infiles),ncol=length(metric)))
colnames(res)=c(metric)
rownames(res)=as.character(infiles)

for (k in as.character(infiles)) {
	star.dat=read.table(paste0("${out}",k,"_summary"),sep="\t")
	
	system(paste("grep ",k," /home/lv71395/mwielsch/binfl/BLUEPRINT/analysis/input_align | awk '{ print $2}' > /home/lv71395/mwielsch/binfl/BLUEPRINT/analysis/test_",k,sep=""))
	system(paste0("sort /home/lv71395/mwielsch/binfl/BLUEPRINT/analysis/test_",k, " > /home/lv71395/mwielsch/binfl/BLUEPRINT/analysis/test1_",k))
	tt1=read.table(paste0("/home/lv71395/mwielsch/binfl/BLUEPRINT/analysis/test1_",k))
	t=tt1[1,] 
	system(paste0("rm /home/lv71395/mwielsch/binfl/BLUEPRINT/analysis/test_",k))
	system(paste0("rm /home/lv71395/mwielsch/binfl/BLUEPRINT/analysis/test1_",k))
	t1=gsub(".gz","",as.character(t))
	t=gsub(".fastq","_fastqc",as.character(t1))
	print(c(t))
	ts=gsub("_fastqc","",as.character(t))

	src.dat=read.table(paste0("/home/lv71395/mwielsch/binfl/BLUEPRINT/analysis/result/",ts,"_screen.txt"),sep="\t",fill=T)
	colnames(src.dat)=c("Genome", "Reads_processed", "Unmapped", "Unmapped_inPERC", "One_hit_one_genome", "One_hit_one_genome_inPERC", "Multiple_hits_one_genome", "Multiple_hits_one_genome_inPERC", "One_hit_multiple_genomes", "One_hit_multiple_genomes_in_PERC", "Multiple_hits_multiple_genomes", "Multiple_hits_multiple_genomes_in_PERC")
	print(c("SRC_DAT_works"))
	fast.sum=read.table(paste0("${out}",t,"/summary.txt"),sep="\t")
	print(c("fast_sum_works"))
	GC=system(paste0("sed -n '10p' ${out}",t, "/fastqc_data.txt | awk '{print \$2}'"),intern=T)

	res[k,1:6]=c(as.character(k),star.dat\$V2[2], star.dat\$V2[6],src.dat\$Unmapped_inPERC[1],star.dat\$V2[5],star.dat\$V2[13] )
	res[k,7:11]=c(100-src.dat\$Unmapped_inPERC[4], 100-src.dat\$Unmapped_inPERC[5],src.dat\$One_hit_one_genome_inPERC[2],src.dat\$One_hit_one_genome_inPERC[3],src.dat\$One_hit_one_genome_inPERC[6])
	res[k,12:17]=c(fast.sum\$V1[1], fast.sum\$V1[2],GC,star.dat\$V2[7], star.dat\$V2[17], star.dat\$V2[15])
}

write.table(res,file='${out}QC_alignment_overview.txt',sep="\t",col.names=T,row.names=T,quote=F)

EOF

sed -i 's/%//g' ${out}QC_alignment_overview.txt

###   clean  up out directory !!!

#mkdir ${out}bam_files
#mv ${out}*Aligned.sortedByCoord.out.bam* ${out}bam_files
#mkdir ${out}fastQC_reports
#mv ${out}*.html ${out}fastQC_reports
#mv ${out}*.pdf ${out}fastQC_reports
#mkdir ${out}QC_and_LOGs
#mv ${out}*.zip ${out}QC_and_LOGs
#ämv ${out}*.out ${out}QC_and_LOGs
#mv ${out}*screen.txt ${out}QC_and_LOGs
#mv ${out}*.tab ${out}QC_and_LOGs
#mv ${out}*summary ${out}QC_and_LOGs

