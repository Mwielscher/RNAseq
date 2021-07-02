###   -----------   ApoA1 RNAseq fresh batch of controls
rm(list=ls())
library(limma)
library(Glimma)
library(edgeR)
###    data in 
phe=read.table("merged_phe_new_controls.txt",header=T,sep="\t")
load("edgeR_DEG_input_data_ready_to_use.RData")
old=as.data.frame(dat1$counts)
old=old[,colnames(old) %in% phe$comb_ID]
new=read.table("ApoA1_count_matrix_REP_CONTR_SE.txt",header=T)
new=new[,colnames(new)%in%phe$comb_ID]
###   -------------------  merge datasets
new$ID=as.character(rownames(new))
old$ID=as.character(rownames(old))
dat=merge(new,old,by.x="ID",by.y="ID")
rownames(dat)=as.character(dat$ID)
dat=as.matrix(dat[,-1])
##  ------------------
coverage=apply(dat,2,sum)
table(coverage < 5000000)
##  ----------  make DGE object
dat=dat[,as.character(phe$comb_ID)]
outcome=as.factor(phe$outcome)
dat1 = DGEList(counts=dat, samples =phe)
#  ----------  check genes
table(rowSums(dat1$counts==0)==ncol(dat1))   ### genes 0 in all samples
keep.exprs = filterByExpr(dat1, group=outcome)
table(keep.exprs)
dat1=dat1[keep.exprs,]
dat1 <- calcNormFactors(dat1, method = "TMM")
dat1$samples$norm.factors
####  check samples
cpm = cpm(dat1,length=dat1$length)
lcpm = cpm(dat1, log=TRUE)
L = mean(dat1$samples$lib.size) * 1e-6
M = median(dat1$samples$lib.size) * 1e-6
c(L, M)
QC_metric=summary(lcpm)
QC_metric
###     ------------------------------------   plotting
lcpm.cutoff <- log2(10/M + 2/L)
## ------------------------- this will use the lcpm from before removal of low expressed genes
library(RColorBrewer)
nsamples = ncol(dat1)
col = brewer.pal(nsamples, "Paired")
par(mfrow=c(1,2))
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="A. Raw data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}
##   ------------------  here we rerun lcpm on cleaned up dataset 
lcpm <- cpm(dat1, log=TRUE)
plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
title(main="B. Filtered data", xlab="Log-cpm")
abline(v=lcpm.cutoff, lty=3)
for (i in 2:nsamples){
  den <- density(lcpm[,i])
  lines(den$x, den$y, col=col[i], lwd=2)
}

###   -------------------------------------

####   ------------------------------------------------------  NORMALISATION
dat1 <- calcNormFactors(dat1, method = "TMM")
dat1$samples$norm.factors
#####   ======================================================
### ---------------------------------------------------------- boxplot
x2 <- dat1
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5
## ------------------------------
par(mfrow=c(1,2))
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="RNAseq ApoA1: Unnormalised data",ylab="Log-cpm")
#x2 <- calcNormFactors(x2)  
#x2$samples$norm.factors
lcpm <- cpm(dat1, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="RNAseq ApoA1: Normalised data",ylab="Log-cpm")
#####  there are mean outliers in the normalized data
find_median_outliers=apply(lcpm,2,median)
as.data.frame(find_median_outliers)


####-----------------------------------
###  --------------------------    now exclude median outliers
dat1=dat1[,!colnames(dat1) %in% c("S76706","S76668","S76665")]
lcpm <- cpm(dat1, log=TRUE)
dat1 <- calcNormFactors(dat1, method = "TMM")

####  ------------------------
dev.off()
#SIGMA
#LPS_treatment
cols=as.character(dat1$samples$LPS_treatment)
cols=gsub("1","red",cols)
cols=gsub("0","blue",cols)
plotMDS(dat1,  pch=16,col=cols,dim.plot=c(1,2))
##  ---------
cols=as.character(dat1$samples$batch)
cols=gsub("1","red",cols)
cols=gsub("2","blue",cols)
cols=gsub("3","black",cols)
cols=gsub("4","green",cols)
#cols=gsub("3","green",cols)
plotMDS(dat1,  pch=16,col=cols,dim.plot=c(1,2))
## -------------------------------------------------------
#####   ggplot
mds = plotMDS(dat1)
pc12= as.data.frame(cbind(names(mds$x),mds$x,mds$y))
colnames(pc12)=c("comb_ID","PC1","PC2")
pc12a=merge(pc12,phe,by.x="comb_ID",by.y="comb_ID")
library(dplyr)
library(tidyr)
library(ggplot2)
ggplot(pc12a, aes(x=as.numeric(PC1), y=as.numeric(PC2), shape=as.factor(LPS_treatment),color=as.factor(Batch.nr))) +
  geom_point(size=3.5) +
  theme_bw()
###  ---------------------------------------------------    batch removal
outcome=as.factor(dat1$samples$outcome)
Batch.nr=as.factor(dat1$samples$Batch.nr)
batch=as.factor(dat1$samples$batch)
design = model.matrix(~0+outcome + Batch.nr + batch )
colnames(design) <- gsub("outcome", "", colnames(design))
colnames(design)

####
dat.voom1= voomWithQualityWeights(dat1, design, var.design= design,plot = T)
#dat.voom.batch= removeBatchEffect(dat.voom1,batch=Batch.nr,batch2=batch2)

dat.voom.batch= removeBatchEffect(dat.voom1,batch=Batch.nr)
dat.voom1$E=dat.voom.batch
mds=plotMDS(dat.voom1)
#####   ggplot
mds = plotMDS(dat.voom1)
pc12= as.data.frame(cbind(names(mds$x),mds$x,mds$y))
colnames(pc12)=c("comb_ID","PC1","PC2")
pc12a=merge(pc12,phe,by.x="comb_ID",by.y="comb_ID")
library(dplyr)
library(tidyr)
library(ggplot2)
ggplot(pc12a, aes(x=as.numeric(PC1), y=as.numeric(PC2), shape=as.factor(LPS_treatment),color=as.factor(Batch.nr))) +
  geom_point(size=3.5) +
  theme_bw()
###   ----------------------------------------------------------------
##  ---------------   explore
##  ------- remove SIGMA ApoA1
dat.voom1=dat.voom1[,!colnames(dat.voom1)%in% phe$comb_ID[phe$SIGMA %in% c(1)]]
mds = plotMDS(dat.voom1)
pc12= as.data.frame(cbind(names(mds$x),mds$x,mds$y))
colnames(pc12)=c("comb_ID","PC1","PC2")
pc12a=merge(pc12,phe,by.x="comb_ID",by.y="comb_ID")
library(dplyr)
library(tidyr)
library(ggplot2)
ggplot(pc12a, aes(x=as.numeric(PC1), y=as.numeric(PC2), shape=as.factor(LPS_treatment),color=as.factor(Batch.nr))) +
  geom_point(size=3.5) +
  theme_bw()
#####  AopA1 time
ggplot(pc12a, aes(x=as.numeric(PC1), y=as.numeric(PC2), shape=as.factor(ApoA1_time),color=as.factor(Batch.nr))) +
  geom_point(size=3.5) +
  theme_bw()


ggplot(pc12a, aes(x=as.numeric(PC1), y=as.numeric(PC2), shape=as.factor(ApoA1_conc),color=as.factor(Batch.nr))) +
  geom_point(size=3.5) +
  theme_bw()

###  -------------------------------------------------------------------------------------------------------------
### -----------------------------------------------------------------------------
#####   -----------    LPS ARM 
pc12b=pc12a[pc12a$outcome %in% c("contr_LPS_4","LPS_ApoA1_10_4","LPS_ApoA1_100_4","LPS_ApoA1_1000_4"),]

ggplot(pc12b, aes(x=as.numeric(PC1), y=as.numeric(PC2), shape=as.factor(ApoA1_conc),color=as.factor(pc12b$outcome))) +
  geom_point(size=3.5) +
  theme_bw()

## --------------------------
pc12b=pc12a[pc12a$outcome %in% c("contr_LPS_24","LPS_ApoA1_10_24","LPS_ApoA1_100_24","LPS_ApoA1_1000_24"),]

ggplot(pc12b, aes(x=as.numeric(PC1), y=as.numeric(PC2), shape=as.factor(ApoA1_conc),color=as.factor(pc12b$outcome))) +
  geom_point(size=3.5) +
  theme_bw()

## ----------------------  preTreatment
pc12b=pc12a[pc12a$outcome %in% c("contr_LPS_preT","ApoA1_LPS_1000_24"),]

ggplot(pc12b, aes(x=as.numeric(PC1), y=as.numeric(PC2), shape=as.factor(ApoA1_conc),color=as.factor(pc12b$outcome))) +
  geom_point(size=3.5) +
  theme_bw()

###  --------------------------------------------------------------------
###  ------------------ ApoA1   ARM
pc12a=pc12a[!pc12a$inference %in% c("contr"),]
pc12b=pc12a[pc12a$outcome %in% c("contr","ApoA1 _10_4","ApoA1 _100_4","ApoA1_1000_4"),]
ggplot(pc12b, aes(x=as.numeric(PC1), y=as.numeric(PC2), shape=as.factor(ApoA1_conc),color=as.factor(pc12b$outcome))) +
  geom_point(size=3.5) +
  theme_bw()

pc12b=pc12a[pc12a$outcome %in% c("contr","ApoA1 _10_24","ApoA1 _100_24","ApoA1_1000_24"),]
ggplot(pc12b, aes(x=as.numeric(PC1), y=as.numeric(PC2), shape=as.factor(ApoA1_conc),color=as.factor(pc12b$outcome))) +
  geom_point(size=3.5) +
  theme_bw()
################## ----------------------------------------------------------------
### ========================================             REGRESSIONS
## --------------------------------------------------------------------------------
## exclude SIGMA and controls
dat.voom1=dat.voom1[,!colnames(dat.voom1)%in% phe$comb_ID[phe$inference %in% c("contr")]]
##   -----------------------------------------------------   create data object
save(phe,dat.voom1,file="APOA1_data_new_controls_batch_removed_for_REGRESSION.RData")
##  ------------ start here
rm(list=ls())
load("APOA1_data_new_controls_batch_removed_for_REGRESSION.RData")
###-----
colnames(dat.voom1$design)
contr.matrix <- makeContrasts(
  LPS_ApoA1_10_24_vs_contr_LPS_24 = LPS_ApoA1_10_24-contr_LPS_24, 
  LPS_ApoA1_100_24_vs_contr_LPS_24 = LPS_ApoA1_100_24-contr_LPS_24,
  LPS_ApoA1_1000_24_vs_contr_LPS_24 = LPS_ApoA1_1000_24-contr_LPS_24, 
  LPS_ApoA1_10_4_vs_contr_LPS_4 = LPS_ApoA1_10_4-contr_LPS_4,
  LPS_ApoA1_100_4_vs_contr_LPS_4 = LPS_ApoA1_100_4-contr_LPS_4,
  LPS_ApoA1_1000_4_vs_contr_LPS_4 = LPS_ApoA1_1000_4-contr_LPS_4,
  contr_vs_contr_LPS_4 = contr-contr_LPS_4,
  contr_vs_contr_LPS_24 = contr-contr_LPS_24,
  ApoA1_10_24_vs_contr = ApoA1_10_24-contr,
  ApoA1_100_24_vs_contr = ApoA1_100_24-contr,
  ApoA1_1000_24_vs_contr = ApoA1_1000_24-contr,
  ApoA1_10_4_vs_contr = ApoA1_10_4-contr,
  ApoA1_100_4_vs_contr = ApoA1_100_4-contr,
  ApoA1_1000_4_vs_contr = ApoA1_1000_4-contr,
  ApoA1_LPS_1000_24_vs_contr_LPS_preT = ApoA1_LPS_1000_24-contr_LPS_preT,
 levels = colnames(dat.voom1$design))

vfit = limma::lmFit(dat.voom1, dat.voom1$design)
vfit = contrasts.fit(vfit, contrasts=contr.matrix)
efit = eBayes(vfit)
plotSA(efit, main="Final model: Mean-variance trend")
summary(decideTests(efit))
write.table(as.data.frame(summary(decideTests(efit))),file="initial_overview.txt",sep="\t",col.names = T,row.names = T,quote=F)
#######   --------------------
traits=colnames(contr.matrix)
library(annotate)
require('org.Hs.eg.db')

for (k in 1:length(traits)){
  inter = topTreat(efit, coef=k, n=Inf)
  genes=as.character(rownames(inter))
  test=getSYMBOL(genes,'org.Hs.eg.db')
  test1=as.data.frame(test)
  test1$ID=as.character(rownames(test1))
  colnames(test1)=c("SYMBOL","ID")
  inter=as.data.frame(inter)
  inter$ID=as.character(rownames(inter))
  res.fin1=merge(inter,test1,all.x=T)
  colnames(res.fin1)[1]=c("ENTREZ_ID")
  res.fin1=res.fin1[order(res.fin1$adj.P.Val),]
  
  k.trait=as.character(traits[k])
  write.table(res.fin1,file=paste0("ApoA1_new_controls_DEG_lists_",k.trait,".txt"),sep="\t",col.names = T,row.names = T,quote=F)
}

###    -----------------------------------------------------------------------------------
#######    -----------------------------      enrichments and boxplots now
## contr -- contr_LPS --- LPS_ApoA1_10 ---  LPS_ApoA1_100 ---LPS_ApoA1_1000

###   ------------------------------------   organise datasets
rownames(phe)=as.character(phe$comb_ID)
phe=phe[as.character(colnames(dat.voom1)),]
###  -----------------------------  entrez to SYMBOL function --only use for visualizations
entrez_to_symbol=function(dat){
  require(org.Hs.eg.db)
  require(annotate)
  genes=as.character(rownames(dat))
  test=getSYMBOL(genes,'org.Hs.eg.db')
  test1=as.data.frame(test)
  test1$ID=as.character(rownames(test1))
  colnames(test1)=c("SYMBOL","ID")
  dat=as.data.frame(dat)
  dat$ID=as.character(rownames(dat))
  res.fin1=merge(dat,test1,all.x=T)
  res.fin1=res.fin1[!duplicated(res.fin1$SYMBOL),]
  res.fin1=res.fin1[!is.na(res.fin1$SYMBOL),]
  rownames(res.fin1)=as.character(res.fin1$SYMBOL)
  res.fin1=res.fin1[,!colnames(res.fin1)%in% c("SYMBOL","ID")]
  return(res.fin1)
}
#### -----------------------------------------------------------------------------
## --------------------------    boxplots 
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
###     ----------------------------------------------------------------------------------------------------
###    ------------------------------------------------------------------------------------------------
#####   ---------------    upset plot
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



sets = c("contr_vs_contr_LPS_24","LPS_ApoA1_10_24_vs_contr_LPS_24","LPS_ApoA1_100_24_vs_contr_LPS_24","LPS_ApoA1_1000_24_vs_contr_LPS_24")
#####  ----------------------------------------------------------------------------------------------------------
######    -------------------------------------------------                                              heatmap
##############################################
## use plot.d1 and phe1
library(gplots)
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x,method="euclidean")
mycol <- colorpanel(1000,"blue","white","red")
# LPS_ApoA1_100_24_vs_contr_LPS_24 is a vector for creating the upset plot

heatmap.2(as.matrix(plot.d1[LPS_ApoA1_1000_24_vs_contr_LPS_24,]), Colv=FALSE, labCol=as.character(phe1$outcome),key=F,col=mycol,scale="row",
          dendrogram="row",trace="none",  margins=c(6,9), lwid = c(1,10), lhei = c(3,17), 
          hclust=hclustfunc,distfun=distfunc)

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
#######3    ------------------------------------------------------------------------------------
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







•	Chief Editor: Robert O. Bonow, MD, MS 

Deputy Editors: Patrick T. O’Gara, MD, Marc S. Sabatine, MD,, Clyde W. Yancy, MD,


Marta Koch, Editor-in-Chief:  MKoch@lancet.com | David Holmes, Deputy Editor  David.Holmes@lancet.com


General editorial enquiries:editorial@lancet.com 
Lancet Journal Office: +44 (0) 207 424 4950 
