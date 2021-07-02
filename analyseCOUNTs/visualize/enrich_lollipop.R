### ------------------------------------
##    http://yulab-smu.top/clusterProfiler-book/chapter12.html 
ora.dat=read.table("FINAL_RESULT_RNA/ORA_PRE_TREATMENT_MSigDB_Hallmark_2020_noFILTER.txt",sep="\t",header=T)

ora.dat$rich=(ora.dat$in_overlap/ora.dat$sigGenes)/(ora.dat$set_size/backgroundN)
ora.dat1=ora.dat[ora.dat$Adjusted.P.value <0.05,]
colnames(ora.dat1)[2]=c("count")

library(ggplot2)
library(forcats)
library(enrichplot)

ggplot(ora.dat1, showCategory = 20, 
       aes(rich, fct_reorder(Term, rich))) + 
  geom_segment(aes(xend=0, yend = Term)) +
  geom_point(aes(color=Adjusted.P.value, size = count)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 8)) +
  theme_minimal() + 
  xlab("enrichment") +
  ylab(NULL) + 
  ggtitle("GO Biological Process")
