rm(list=ls())
### --------  data wrangling select project related data
samp=read.csv("Donor_raw_sel_pooled_med.csv",header=T)
master=read.table("MASTER_SAMPLE_LIST_V2.txt",sep="\t",header=T)
master=master[,colnames(master) %in% c("BARCODE_NAME","Janas_donors","Contaminated_cells")]
samp1=merge(master,samp,all.x=T,by.x="BARCODE_NAME",by.y="BARCODE_NAME")
##  4 samples are contaminated colname == Contaminated_cells
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
## ------------  also kick out PCA outlier 
outlier=c("IL2_Donor_16_S71943")
dat=dat[,!colnames(dat) %in% c(outlier) ]
## -----------------  for now kick out JANA samples and outlier
lana=as.character(samp1$BARCODE_NAME)[samp1$Janas_donors %in% c("Yes")]
dat=dat[,!colnames(dat) %in% c(lana) ]
##
samp1=samp1[samp1$BARCODE_NAME %in% colnames(dat),]
samp1=samp1[as.character(colnames(dat)),]
gene_read=rowSums(dat)
table(gene_read>50)
dat=dat[gene_read>50, ]
##---------------- inspect pheno data
samp1=samp1[,1:19]
table(samp1$comb_pheno)
phe=samp1
## -------------   align countmatrix an phenotype file
phe=phe[as.character(colnames(dat)),]

#### -------------------------------  initial round to get an overview
library(DESeq2) 
dds = DESeqDataSetFromMatrix(countData = dat,
                             colData = phe,
                             design = ~ Genotype)

dds <- DESeq(dds)
vsd <- vst(dds, blind=FALSE)
mat.vsd=assay(vsd)
#write.table(mat.vsd,file="DESEQ_normalized_data_full_IL2_donors_22_04_2021.txt",sep="\t",col.names=T,row.names=T,quote=F)
#####  ------------   quick look
pcs=plotPCA(vsd, intgroup=c("comb_pheno"))  
head(vsd@colData$Genotype)
vsd1=vsd[,vsd@colData$Genotype %in% c("ARPC1B", "HEM1", "ND", "PAK2", "WAS", "WDR1")]
plotPCA(vsd1, intgroup=c("Genotype"))   

### restrict for one GT
phe1=phe[phe$Genotype %in% c("WDR1"),]
dat1=dat[,colnames(dat) %in% phe1$BARCODE_NAME]
phe1=phe1[colnames(dat1),]

##  -----------------------------   run regressions
library(DESeq2)
dds = DESeqDataSetFromMatrix(countData = dat1,
                             colData = phe1,
                             design = ~ comb_pheno)

dds$comb_pheno = relevel(dds$comb_pheno, ref = "WDR1_GM")
dds = DESeq(dds)
vsd = vst(dds)
mat.vsd=assay(vsd)
plotPCA(vsd, intgroup=c("comb_pheno"))
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
