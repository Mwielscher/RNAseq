####   enrichments incl ballon plots
###    ----------------------------------------------------------------------------------------------------------
################--------------------   --------------- overlaps
library(enrichR)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)

enrich.list=c('MSigDB_Hallmark_2020','WikiPathways_2019_Human','KEGG_2019_Human' ,'Reactome_2016',
              "HumanCyc_2016",'BioPlanet_2019',"MGI_Mammalian_Phenotype_Level_4_2019","GO_Biological_Process_2018",
              "Human_Phenotype_Ontology","ARCHS4_Tissues")

traits=c("contr_vs_contr_LPS_24","LPS_ApoA1_10_24_vs_contr_LPS_24","LPS_ApoA1_100_24_vs_contr_LPS_24","LPS_ApoA1_1000_24_vs_contr_LPS_24")

for(k in enrich.list){
  rm(ora.final)
  for (j in traits){
    print(paste("this is ",k,"for ",j,"now!!",sep=" "))
    ###        -------------------   prep input gene list
    tt=read.table(paste0("ApoA1_new_controls_DEG_lists_",j,".txt"),sep="\t",header=T)
    tt1=tt[tt$adj.P.Val < 0.05,]
    genes1=as.character(tt1$SYMBOL)
    ###    -----------------   calulate overlap
    enr.result = enrichr(as.character(genes1),k)[[k]]
    ##    ---------------------------   flter and prep output datasets
    #inter2=enr.result[enr.result$Adjusted.P.value <0.05,]  
    inter2=enr.result
    inter3=inter2 %>% separate(Overlap, c("in_overlap", "set_size"))
    inter3$overlap_PERC=as.numeric(inter3$in_overlap)/length(genes1)
    inter3$overlap_SET=as.numeric(inter3$in_overlap)/as.numeric(inter3$set_size)
    inter3$sigGenes=length(genes1)
    inter3$trait=rep(j,nrow(inter3))
    inter3=inter3[,c("Term","in_overlap", "set_size","Adjusted.P.value","overlap_PERC","trait", "sigGenes","Genes")]
    if(!exists("ora.final")){
      ora.final=inter3
    }else{
      ora.final=rbind(ora.final,inter3)
    }
  }
  write.table(ora.final,paste0("ORA_list_",k,"_all_taits_noFILTER.txt"),sep="\t",col.names=T,row.names=F,quote=F)
}

####   --------------------------------------------------------------------------------------------------
#######    ------------------------------------------------------------------------------------
####     ----------------------     make balloon plots

enrich.list=c('MSigDB_Hallmark_2020','WikiPathways_2019_Human','KEGG_2019_Human' ,'Reactome_2016',
              "HumanCyc_2016",'BioPlanet_2019',"MGI_Mammalian_Phenotype_Level_4_2019","GO_Biological_Process_2018",
              "Human_Phenotype_Ontology","ARCHS4_Tissues")


k=c("GO_Biological_Process_2018")


in.dat=read.table(paste0("ORA_list_",k,"_all_taits_noFILTER.txt"),sep="\t",header=T)
#traits=unique(in.dat$trait)
rm(fin.ord)
for (i in traits){
  inter=in.dat[in.dat$trait %in% i,]
  inter=inter[order(inter$Adjusted.P.value),]
  inter=inter[1:5,]
  if(!exists("fin.ord")){
    fin.ord=inter
  } else { fin.ord=rbind(fin.ord,inter)}
}
fin.ord=fin.ord[!is.na(fin.ord$Term),]
fin.ord=fin.ord[!duplicated(fin.ord$Term),]
fin.ord=fin.ord[order(fin.ord$Adjusted.P.value), ]
fin.ord$ord=c(1:nrow(fin.ord))
## -----------------------------------------
fin.dat=fin.ord[,c("Term","ord")]
colnames(fin.dat)=c("ID","ORD")
rm(plot.dat)
for (t in traits) {
  
  enr=in.dat[in.dat$trait %in% t,]
  enr$p.adjust=-log10(enr$Adjusted.P.value)
  enr$genelist=rep(as.character(t),nrow(enr))
  clus.dat1=enr[,c("Term","overlap_PERC","p.adjust","trait")]
  clus.dat=merge(fin.dat,clus.dat1,all.x=T,by.x="ID",by.y="Term")
  clus.dat=clus.dat[order(as.numeric(as.character(clus.dat$ORD))),]
  clus.dat$p.adjust[is.na(clus.dat$p.adjust)]=0.9
  clus.dat$trait[is.na(clus.dat$trait)]=as.character(t)
  clus.dat$overlap_PERC[is.na(clus.dat$overlap_PERC)]=0 
  if (!exists("plot.dat")){
    plot.dat=clus.dat
  }else {
    plot.dat=rbind(plot.dat,clus.dat)
  }
} # cluster loop
plot.dat=plot.dat[!is.na(plot.dat$trait),]
colnames(plot.dat)[3]=c("Genes.in.Overlap")
## -- order X axis
plot.dat$trait = factor(plot.dat$trait,levels = c("contr_vs_contr_LPS_24","LPS_ApoA1_10_24_vs_contr_LPS_24",
                                                  "LPS_ApoA1_100_24_vs_contr_LPS_24","LPS_ApoA1_1000_24_vs_contr_LPS_24"))

plot.dat$p.adjust=ifelse(plot.dat$p.adjust>10,10,plot.dat$p.adjust)   ## winsorize to 10
plot.dat$ID=strtrim(as.character(plot.dat$ID),55)      ### cut after 35 characters
p1 <- ggballoonplot(plot.dat,x="trait",y="ID", size="Genes.in.Overlap",fill="p.adjust") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold",angle = 45, hjust = 1)) +
  ggtitle(paste(k,sep=" ")) +
  theme(plot.title = element_text(face="bold",color="darkblue",hjust = 0.5))
p1 = p1 + scale_fill_gradient(low = "#4DBBD5FF", high = "#DC0000FF")
p1
