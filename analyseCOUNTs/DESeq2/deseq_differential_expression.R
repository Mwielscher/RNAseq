rm(list=ls())
samp=read.csv("Donor_raw_sel_pooled_med.csv",header=T)
master=read.table("MASTER_SAMPLE_LIST_V2.txt",sep="\t",header=T)
master=master[,colnames(master) %in% c("BARCODE_NAME","Janas_donors","Contaminated_cells")]
samp1=merge(master,samp,all.x=T,by.x="BARCODE_NAME",by.y="BARCODE_NAME")
## exclude and harmonize 4 samples are contaminated colname == Contaminated_cells
## a few samples need to be forwarded to Jana colname == "Janas_donors"
samp1=samp1[!samp1$Contaminated_cells %in% c("Yes (to be removed)"),]
rownames(samp1)=as.character(samp1$BARCODE_NAME)
dat=read.table("count_matrix_anton_ok.txt",header=T)
### check samples
dat=dat[,colnames(dat) %in% samp$BARCODE_NAME]
samp$BARCODE_NAME[!samp$BARCODE_NAME %in% colnames(dat)]
rownames(samp1)=as.character(samp1$BARCODE_NAME)
samp1=samp1[as.character(colnames(dat)),]
###----  check RNA seq data and TRIM
samp_read=colSums(dat)
table(samp_read < 2000000)
dat=dat[,!samp_read < 2000000 ]
## ------------  also kick out PCA outlier and possible SWOPS
#swops=c("GM_Donor_1_S71937","IL2_Donor_12_S71899","IL2_Donor_4_S71910")
outlier=c("IL2_Donor_16_S71943")
dat=dat[,!colnames(dat) %in% c(outlier) ]
#dat=dat[,!colnames(dat) %in% c("GM_Donor_1_S71937","IL2_Donor_12_S71899","IL2_Donor_4_S71910","IL2_Donor_16_S71943") ]
## -----------------  for now kick out JANA samples and outlier
lana=as.character(samp1$BARCODE_NAME)[samp1$Janas_donors %in% c("Yes")]
dat=dat[,!colnames(dat) %in% c(lana) ]
##
samp1=samp1[samp1$BARCODE_NAME %in% colnames(dat),]
samp1=samp1[as.character(colnames(dat)),]
gene_read=rowSums(dat)
table(gene_read>50)
dat=dat[gene_read>50, ]


## ---------------------------------------------------------------------------------------------------------------------------------------------
##############   ----------------------------------------------------   compare IL2 to GM now
##  --------------------------------------------------------------------------------------------------------

##---------------- inspect pheno data
##  ------------------   for now -----------
samp1=samp1[,1:19]
table(samp1$comb_pheno)
#table(samp1$Genotype,samp1$comb_pheno)
###   -------------------   now exclude all GM cells
#table(grepl("IL2_Donor",as.character(samp1$BARCODE_NAME)))
#phe=samp1[grepl("IL2_Donor",as.character(samp1$BARCODE_NAME)),]
#dat=dat[,colnames(dat) %in% as.character(phe$BARCODE_NAME)]
#phe=phe[as.character(colnames(dat)),]
#############     ------------ do selection
phe=samp1
phe=phe[as.character(colnames(dat)),]

library(DESeq2)
dds = DESeqDataSetFromMatrix(countData = dat,
                             colData = phe,
                             design = ~ Genotype)

# design = ~ IL2_group + Genotype)
dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)
mat.vsd=assay(vsd)
#write.table(mat.vsd,file="DESEQ_normalized_data_full_IL2_donors_22_04_2021.txt",sep="\t",col.names=T,row.names=T,quote=F)
#####  ------------   quick look
pcs=plotPCA(vsd, intgroup=c("comb_pheno"))   ### that 
outlier_find=as.data.frame(cbind(pcs$data$PC1,pcs$data$PC2,pcs$data$name))
outlier_find=outlier_find[order(-as.numeric(outlier_find$V1)),]

head(vsd@colData$Genotype)
vsd1=vsd[,vsd@colData$Genotype %in% c("ARPC1B", "HEM1", "ND", "PAK2", "WAS", "WDR1")]
plotPCA(vsd1, intgroup=c("Genotype"))   ### that 



### restrict  
# 1. go for ND
phe1=phe[phe$Genotype %in% c("WDR1"),]
dat1=dat[,colnames(dat) %in% phe1$BARCODE_NAME]
phe1=phe1[colnames(dat1),]
###  ---------------------  funny outliyeng sample
#phe1=phe1[!phe1$BARCODE_NAME %in% c("IL2_Donor_3_S71917"),]
#dat1=dat[,colnames(dat) %in% phe1$BARCODE_NAME]
#phe1=phe1[colnames(dat1),]

table(phe1$comb_pheno)

##  -----------------------------
library(DESeq2)
dds = DESeqDataSetFromMatrix(countData = dat1,
                             colData = phe1,
                             design = ~ comb_pheno)

dds$comb_pheno <- relevel(dds$comb_pheno, ref = "WDR1_GM")

dds <- DESeq(dds)
vsd <- vst(dds)
mat.vsd=assay(vsd)
#plotPCA(vsd, intgroup=c("IL2_group", "Genotype"))
plotPCA(vsd, intgroup=c("comb_pheno"))
#write.table(mat.vsd,file="DESEQ_normalized_data_Lana_FINAL.txt",sep="\t",col.names=T,row.names=T,quote=F)
##### Regressions 
res =results(dds)
resultsNames(dds)
###     ----------------- get result  
library(apeglm)
library(annotate)
require('org.Hs.eg.db')
traits =c("comb_pheno_WDR1_IL2_vs_WDR1_GM" ) 
for (k in traits){
  print(paste("this is ",k, "now !!!!"))
  resLFC <- lfcShrink(dds, coef=as.character(k), type="apeglm")
  resLFC=resLFC[!is.na(resLFC$baseMean), ]
  resLFC=resLFC[!is.na(resLFC$padj), ]
  res2=as.data.frame(resLFC)
  genes=as.character(rownames(res2))
  
  test=getSYMBOL(genes,'org.Hs.eg.db')
  test1=as.data.frame(test)
  test1$ID=as.character(rownames(test1))
  colnames(test1)=c("SYMBOL","ID")
  res2=as.data.frame(res2)
  res2$ID=as.character(rownames(res2))
  res.fin1=merge(res2,test1,all.x=T)
  colnames(res.fin1)[1]=c("ENTREZ_ID")
  res.fin1=res.fin1[order(res.fin1$padj),]
  write.table(res.fin1, file=paste("IL2_vs_GM_",k,".txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
  
  res.fin2=res.fin1[res.fin1$padj <0.05 & abs(res.fin1$log2FoldChange) > 1 ,]   # change this 1
  write.table(res.fin2, file=paste("IL2_vs_GM_",k,"_FILTERED.txt",sep=""),sep="\t",col.names=T,row.names=F,quote=F)
  
}
####  ------------------------------------------------------------------------------------
##----------------------------------------------------------------------

overview=read.table("RNA_vs_morphology_2/OVERVIEW_all_regressions.txt",header=T,sep="\t")
colnames(overview)=c("feature","contr_only","incl_PID")
overview$incl_PID[is.na(overview$incl_PID)]=c(0)
overview$contr_only[is.na(overview$contr_only)]=c(0)
overview$target=overview$contr_only- overview$incl_PID
overview=overview[order(-overview$target),]

###    ----------------------------------------------------------------------------------------------------------
################--------------------   --------------- overlaps
library(enrichR)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)




traits =c("ARPC1B_IL2_vs_ARPC1B_GM","HEM1_IL2_vs_HEM1_GM","ND_IL2_vs_ND_GM","PAK2_IL2_vs_PAK2_GM","WAS_IL2_vs_WAS_GM","WDR1_IL2_vs_WDR1_GM")

enrich.list=c("Aging_Perturbations_from_GEO_down","Aging_Perturbations_from_GEO_up","GO_Biological_Process_2018",
              "GO_Molecular_Function_2018","GWAS_Catalog_2019","KEGG_2019_Human",'MSigDB_Hallmark_2020')



k=c("Aging_Perturbations_from_GEO_up")

for(k in enrich.list){
  rm(ora.final)
  for (j in traits){
    print(paste("this is ",k,"for ",j,"now!!",sep=" "))
    ###        -------------------   prep input gene list
    tt_nd=read.table(paste("IL2_vs_GM/IL2_vs_GM_comb_pheno_",j,"_FILTERED.txt",sep=""),header=T,row.names = 1,sep="\t")
    genes1=as.character(tt_nd$SYMBOL)
    genes1=genes1[!is.na(genes1)]
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
    inter3=inter3[inter3$Adjusted.P.value < 0.05, ]
    if(!exists("ora.final")){
      ora.final=inter3
    }else{
      ora.final=rbind(ora.final,inter3)
    }
  }
  write.table(ora.final,paste0("IL2_vs_GM/ORA_list_IL2_vs_GM_",k,"_all_taits.txt"),sep="\t",col.names=T,row.names=F,quote=F)
}

#####    ----------------------------------------------------------- ballon plot


enrich.list=c("Aging_Perturbations_from_GEO_down","Aging_Perturbations_from_GEO_up","GO_Biological_Process_2018",
              "GO_Molecular_Function_2018","GWAS_Catalog_2019","KEGG_2019_Human",'MSigDB_Hallmark_2020')

k=c("KEGG_2019_Human")



in.dat=read.table(paste0("IL2_vs_GM/ORA_list_IL2_vs_GM_",k,"_all_taits.txt"),sep="\t",header=T)
#traits=unique(in.dat$trait)
traits =c("ARPC1B_IL2_vs_ARPC1B_GM","HEM1_IL2_vs_HEM1_GM","ND_IL2_vs_ND_GM","PAK2_IL2_vs_PAK2_GM","WAS_IL2_vs_WAS_GM","WDR1_IL2_vs_WDR1_GM")

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
plot.dat$trait = factor(plot.dat$trait,levels = c("ND_IL2_vs_ND_GM","ARPC1B_IL2_vs_ARPC1B_GM","HEM1_IL2_vs_HEM1_GM"
                                                  ,"PAK2_IL2_vs_PAK2_GM","WAS_IL2_vs_WAS_GM","WDR1_IL2_vs_WDR1_GM"))

plot.dat$p.adjust=ifelse(plot.dat$p.adjust>10,10,plot.dat$p.adjust)   ## winsorize to 10
plot.dat$ID=strtrim(as.character(plot.dat$ID),55)      ### cut after 35 characters
p1 = ggballoonplot(plot.dat,x="trait",y="ID", size="Genes.in.Overlap",fill="p.adjust") +
  theme_bw() +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(face = "bold",angle = 45, hjust = 1)) +
  ggtitle(paste(k,sep=" ")) +
  theme(plot.title = element_text(face="bold",color="darkblue",hjust = 0.5))
p1 = p1 + scale_fill_gradient(low = "#4DBBD5FF", high = "#DC0000FF")
p1































