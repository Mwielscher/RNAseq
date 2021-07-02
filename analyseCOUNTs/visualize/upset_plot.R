traits=c("contr_vs_contr_LPS_24","LPS_ApoA1_10_24_vs_contr_LPS_24","LPS_ApoA1_100_24_vs_contr_LPS_24","LPS_ApoA1_1000_24_vs_contr_LPS_24")

for (k in traits){
  print (paste ("this is ",k,"now"))
  tt=read.table(paste0("ApoA1_new_controls_DEG_lists_",k,".txt"),sep="\t",header=T)
  tt1=tt[tt$adj.P.Val < 0.05,]
  assign(k,as.character(tt1$SYMBOL))
}

library(UpSetR)
listInput2 <- list(contr_vs_contr_LPS_24 = as.character(contr_vs_contr_LPS_24),
                   LPS_ApoA1_10_24_vs_contr_LPS_24 = as.character(LPS_ApoA1_10_24_vs_contr_LPS_24),
                   LPS_ApoA1_100_24_vs_contr_LPS_24 = as.character(LPS_ApoA1_100_24_vs_contr_LPS_24),
                   LPS_ApoA1_1000_24_vs_contr_LPS_24 = as.character(LPS_ApoA1_1000_24_vs_contr_LPS_24))
                  

upset(fromList(listInput2), nsets=7,mainbar.y.label = "Intersections of \ndifferential expressed Genes", 
      keep.order=T,order.by = "freq",
      sets.x.label = "total DEG",text.scale = c(1.3, 1.3, 1, 1, 1.3, 1),cutoff = 2)                  
