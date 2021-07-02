rm(list=ls())
library(limma)
library(Glimma)
library(edgeR)
###    data in 
phe=read.table("merged_phe_new_controls.txt",header=T,sep="\t")
load("edgeR_DEG_input_data_ready_to_use.RData")   ### countmatrix
old=as.data.frame(dat1$counts)
old=old[,colnames(old) %in% phe$comb_ID]
new=read.table("ApoA1_count_matrix_REP_CONTR_SE.txt",header=T). ## new countmatrix
new=new[,colnames(new)%in%phe$comb_ID]
###   -------------------  merge datasets   --- in a "regular" RNAseq we would start of with one countmatrix and the merge would not be necessary
new$ID=as.character(rownames(new))
old$ID=as.character(rownames(old))
dat=merge(new,old,by.x="ID",by.y="ID")
rownames(dat)=as.character(dat$ID)
dat=as.matrix(dat[,-1])
##  ------------------
coverage=apply(dat,2,sum)
table(coverage < 5000000)
##  ----------  make DGE object  --- edgeR
## ---- large chunks of the code are adopted from https://bioconductor.org/packages/release/workflows/html/RNAseq123.html 
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


####-----------------------------------    in a "regular" RNA seq experiments we would not spot irregularities at this point
###  --------------------------    now exclude median outliers
dat1=dat1[,!colnames(dat1) %in% c("S76706","S76668","S76665")]
lcpm <- cpm(dat1, log=TRUE)
dat1 <- calcNormFactors(dat1, method = "TMM")
####  ------------------------
dev.off()
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
###.  -----------------------------------------------------------------------
###  --------------------------------------------------- !!!!!   batch removal  !!!!!!!!!!!
outcome=as.factor(dat1$samples$outcome)
Batch.nr=as.factor(dat1$samples$Batch.nr)
batch=as.factor(dat1$samples$batch)
design = model.matrix(~0+outcome + Batch.nr + batch )
colnames(design) <- gsub("outcome", "", colnames(design))
colnames(design)
####
dat.voom1= voomWithQualityWeights(dat1, design, var.design= design,plot = T)
dat.voom.batch= removeBatchEffect(dat.voom1,batch=Batch.nr)
dat.voom1$E=dat.voom.batch
mds=plotMDS(dat.voom1)
#### ---- check process  
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
##  ---------------   explore and subset
##  ------- subset
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
################## ----------------------------------------------------------------
### ========================================             REGRESSIONS
## --------------------------------------------------------------------------------
## exclude wet-lab-experiment-controls, replicates, and samples unrelated to the research question
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
#######   --------------------  export results 
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

