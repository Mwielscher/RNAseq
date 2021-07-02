plot.d=as.data.frame(dat.voom1$E[,colnames(dat.voom1)%in% phe$comb_ID[phe$outcome %in%c("contr","contr_LPS_24",
                                                                                        "LPS_ApoA1_10_24","LPS_ApoA1_100_24",
                                                                                       "LPS_ApoA1_1000_24")]])
plot.d1=entrez_to_symbol(plot.d)
phe1=phe[phe$comb_ID %in% colnames(plot.d1),]


traits=c("contr_vs_contr_LPS_24","LPS_ApoA1_10_24_vs_contr_LPS_24","LPS_ApoA1_100_24_vs_contr_LPS_24","LPS_ApoA1_1000_24_vs_contr_LPS_24")
for (j in traits){
  print (paste ("this is ",j,"now"))
  tt=read.table(paste0("ApoA1_new_controls_DEG_lists_",j,".txt"),sep="\t",header=T)
  tt1=tt[tt$adj.P.Val < 0.05,]
  genes=as.character(tt1$SYMBOL)

  library(ggsci)

  pdf(paste0("summary_plots",j,"_dotplot.pdf"))
  for (k in genes) {
      inter=as.data.frame(cbind(as.numeric(plot.d1[k,]),as.character(phe1$outcome)))
      colnames(inter)=c("normalized_expression","samp")
       inter$normalized_expression=as.numeric(as.character(inter$normalized_expression))
       inter$samp <- factor(inter$samp , levels=c("contr","contr_LPS_24","LPS_ApoA1_10_24","LPS_ApoA1_100_24","LPS_ApoA1_1000_24"))
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
  }
  dev.off()
}
