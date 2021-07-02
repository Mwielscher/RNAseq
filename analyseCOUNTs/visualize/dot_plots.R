###   --------------------------    make nice dot plots 
genes=as.character(res$SYMBOL[1:50])
library(ggsci)
library(ggplot2)
for (k in genes) {
pdf(paste0("FINAL_RESULT_RNA/PRE_TREATMENT_dot_plots/",k,"_dotplot.pdf"))
  inter=as.data.frame(cbind(lcpm3[k,],as.character(samp$Treatment)))
  colnames(inter)=c("normalized_expression","samp")
  inter$samp <- factor(inter$samp , levels=c("LPS", "ApoA1_LPS","contr"))
  inter$normalized_expression=as.numeric(as.character(inter$normalized_expression))
  p<-ggplot(inter, aes(x=samp, y=normalized_expression, fill=samp)) +
      geom_dotplot(binaxis='y', stackdir='center',dotsize = 2) +
      stat_summary(fun.y = median, fun.ymin = median, fun.ymax = median,
                 geom = "crossbar", width = 0.5,lwd=0.2) +
      scale_fill_aaas() +
    theme_bw() +
    theme(legend.position = "none") +
    theme(axis.title.x = element_blank()) +
    theme(axis.title.y = element_text(size = 12)) +
    theme(axis.text.y = element_text(size=12)) +
    theme(axis.text.x = element_text(angle = 45, hjust=1,size=14)) +
    ggtitle(paste(k, "gene expression")) +
    theme(plot.title = element_text(hjust = 0.5,size = 18, face = "bold"))
  print(p)
dev.off()
}
