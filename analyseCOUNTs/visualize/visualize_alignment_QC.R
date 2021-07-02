qc16=read.table("QC_alignment_2016_overview.txt",header=T)
qc20=read.table("QC_alignment_2020_overview.txt",header=T,row.names=1)
qc20$ID=as.character(rownames(qc20))
qc20=qc20[,c("ID","Number_of_inputREADS")]
qc16=qc16[,c("ID","Number_of_inputREADS")]
qc=as.data.frame(rbind(qc20,qc16))
##  ---------------------------------
all.n=read.table("final_count_matrix_head.txt")
colnames(all.n)=c("ID","eff_counts")
plot.dat=merge(all.n,qc,by.x="ID",by.y="ID",all.x=T)

library(ggplot2)
library(tidyr)
##  QC20 had mtRNA and rRNA swapped
in.dat=qc16[,c("ID","reads_mapped_to_multiple_loci_perc","rRNA_content","mtRNA_content")]

data_long <- gather(in.dat, QC_metric, percent_of_total_reads, reads_mapped_to_multiple_loci_perc:mtRNA_content, factor_key=TRUE)


ggplot(data_long, aes(fill=QC_metric, y=percent_of_total_reads, x=ID)) + 
     geom_bar(position="dodge", stat="identity") +
     coord_cartesian(ylim=c(0,50))+
    theme(axis.text.x = element_text(size=11,angle = 90, vjust = 0.5, hjust=1),
                     axis.text.y = element_text(face="bold"),
                    axis.title.x = element_blank())
#################
######    QC plot 2
in.dat2=qc16[,c("ID","Number_of_inputREADS")]
colnames(in.dat2)[2]=c("total_reads")

ggplot(in.dat2, aes( y=total_reads/1000000, x=ID)) + 
  geom_bar(stat="identity",fill="steelblue") +
  theme(axis.text.x = element_text(size=11,angle = 90, vjust = 0.5, hjust=1),
        axis.title.x = element_blank())
