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

######       ---------------------------------------------------------------------
###   ----------------  collapse replicastes
rownames(samp)=as.character(samp$ID)

to.colapse=unique(samp$ID.1[duplicated(samp$ID.1)])

for (k in to.colapse){
  
  i=as.character(samp$ID[samp$ID.1 %in% k])
  dat[,k]=rowSums(dat[,i])
  dat=dat[,!colnames(dat) %in% i]
  
}
write.table(dat,file="EMPROR_count_matrix_repl_collapsed.txt",col.names=T,row.names=T,quote=F)
samp=samp[!duplicated(samp$ID.1),]
write.table(samp,file="phenotype_repl_collapsed.txt",col.names=T,row.names=T,quote=F)


###  ----------------------------------------  start edgeR limma
####     --------------------------
dat=read.table("EMPROR_count_matrix_repl_collapsed.txt",header=T)
genes=read.table("EMPROR_feature_annot.txt",header=T)
samp=read.table("phenotype_repl_collapsed.txt",header=T)
samp$ID.1=gsub("-",".",as.character(samp$ID.1))
rownames(samp)=as.character(samp$ID.1)
samp=samp[as.character(colnames(dat)),]



### minimum reads accross all samples
dat[is.na(dat)] = 0
sumUP2=apply(dat,2,sum)

sumUP=apply(dat,1,sum)
table(sumUP>50)  ## initial with 50              
dat=dat[sumUP>50,]  ## that leaves 16739  behind

##  -- 2nd iteration 
dat=dat[rownames(dat) %in% keep2,]
genes=genes[genes$GeneID%in% as.character(row.names(dat)),]
rownames(genes)=as.character(genes$GeneID)
genes=genes[as.character(rownames(dat)), ]
### ------------------- align for DEG object
library(limma)
library(Glimma)
library(edgeR)

outcome=samp$celltype
dat1 = DGEList(counts=dat, group=outcome)
dat1$batch = samp$batch
dat1$celltype = samp$celltype
dat1$age_group =samp$age_group
dat1$genes = genes$GeneID
dat1$length = genes$Length
######
###   ---------------------  object set up complete
cpm = cpm(dat1,gene.length=length)
lcpm <- cpm(dat1,gene.length=length, log=TRUE)
L <- mean(dat1$samples$lib.size) * 1e-6
M <- median(dat1$samples$lib.size) * 1e-6
c(L, M)
QC_metric=summary(lcpm)
QC_metric
table(rowSums(dat1$counts==0)==9)
keep.exprs = filterByExpr(dat1, group=outcome)
table(keep.exprs)
##  record and go back to start !! ?
keep2=rownames(dat1)[keep.exprs]
###  -------------------------
####   ------------------------------------------------------  NORMALISATION
dat1 <- calcNormFactors(dat1, method = "TMM")
dat1$samples$norm.factors

### ---------------------------------------------------------- boxplot
x2 <- dat1
x2$samples$norm.factors <- 1
x2$counts[,1] <- ceiling(x2$counts[,1]*0.05)
x2$counts[,2] <- x2$counts[,2]*5
## ------------------------------
par(mfrow=c(1,2))
nsamples <- ncol(x2)
col <- brewer.pal(nsamples, "Paired")
lcpm <- cpm(x2, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="A. Example: Unnormalised data",ylab="Log-cpm")
#x2 <- calcNormFactors(x2)  
#x2$samples$norm.factors
lcpm <- cpm(dat1, log=TRUE)
boxplot(lcpm, las=2, col=col, main="")
title(main="B. Example: Normalised data",ylab="Log-cpm")

##   ----------------------------------------------------------
pr2=plotMDS(lcpm, labels=dat1$samples$group)
p.lcpm=lcpm
#colnames(p.lcpm)=as.character(dat1$samples$group)

pr <- prcomp(t(p.lcpm), center=TRUE, scale=TRUE )
#plot(pr)
p.dat=pr$x
assoc.dat=as.data.frame(cbind(p.dat, samp))

pc12=as.data.frame(cbind(p.dat[,2:20], samp))
pc12$batch=as.factor(pc12$batch)
pc12$p.batch=as.factor(pc12$p.batch)
pc12$age_group=as.factor(pc12$age_group)

#pc12$comb_ID=as.character(rownames(pc12))
library(dplyr)
library(tidyr)

#pc12a=pc12 %>% separate(comb_ID, c("exposure", "time"))

###   
ggplot(pc12, aes(x=PC6, y=PC7, shape=celltype,color=age_group)) +
  geom_point(size=4) +
  scale_color_manual(values=c("#999999", "#E69F00","blue")) +
  theme_bw()


###############   check PC assoc to batch

sig_PCs=c(15)
mt_vs_PCs=data.frame(ID= rep(NA,sig_PCs),Estimate=rep(NA,sig_PCs),std_err=rep(NA,sig_PCs),t_val=rep(NA,sig_PCs),P_val=rep(NA,sig_PCs))

for (k in 1:sig_PCs){
  fit=summary(glm(paste("PC",k,"~batch",sep=""),data=assoc.dat))$coefficients[2,]
  mt_vs_PCs[k,1]=paste("PC_",k,sep="")
  mt_vs_PCs[k,2:5]=fit
}
print(paste("multiple testing threshold is ",0.05/sig_PCs,sep=""))
print("++++++++++ associations between mitochnodrial RNA content and PCs in analysis +++++++++++++")
mt_vs_PCs

write.table(mt_vs_PCs,file="PC_assoc_to_BATCH_repl_collapsed.txt",sep="\t",col.names=T)

####

library("factoextra")

pr2=prcomp(t(p.lcpm), center=TRUE, scale=TRUE )
fviz_screeplot(pr2, addlabels = TRUE, ylim = c(0, 30))
###  -------   UMAP

library(umap)

init.umap = umap(p.dat[,2:6])
plot.it=as.data.frame(init.umap$layout)
colnames(plot.it) =c("UMAP_1","UMAP_2")
u.plot=as.data.frame(cbind(pc12,plot.it))
u.plot$age_group[u.plot$age_group %in% c(2016)]=NA
#p.batch

ggplot(u.plot, aes(x=UMAP_1, y=UMAP_2, shape=celltype,color=p.batch)) +
  geom_point(size=4) +
  scale_color_manual(values=c("#999999", "#E69F00")) +
  theme_bw()


ggplot(u.plot, aes(x=UMAP_1, y=UMAP_2, shape=celltype,color=age_group)) +
  geom_point(size=4) +
  theme_bw()

###########    -------------------------   dendrogramm
##we used values adjusted by variance stabilizing transformation from DESeq2 
# in which batch-effects had been corrected for with ComBat (47).

#http://www.sthda.com/english/wiki/beautiful-dendrogram-visualizations-in-r-5-must-known-methods-unsupervised-machine-learning

init.dendro = p.dat[,2:6]
batch=as.factor(samp$batch)
library(sva)
combat_edata = ComBat(dat=t(init.dendro), batch=batch)
dd <- dist(t(combat_edata), method = "euclidean")
hc <- hclust(dd, method = "ward.D2")

# vanilla
plot(hc, hang = -1, cex = 0.6)
##  
hcd <- as.dendrogram(hc)
par(mar=c(7,4,2,0))
plot(hcd, type = "rectangle", ylab = "Height",tip.color = coloring)

library(dendextend)
library(ggsci)
dend=as.dendrogram(hc)

all.l=labels(dend)
samp1=samp[as.character(all.l),]
mypal = pal_npg("nrc", alpha=1)(length(unique(samp$celltype)))
coloring=rep("#E64B35FF",nrow(samp1))
coloring[samp1$celltype %in% c("EMPROR")] = c("#4DBBD5FF")
coloring[samp1$celltype %in% c("MONO")] =c("#00A087FF")
coloring[samp1$celltype %in% c("PDC")] =c("#3C5488FF")


dend %>% set("labels_col", coloring) %>% # cset("labels_cex", 2) %>% # Change size
  plot() # plot


##  fan needs an outer margin
#omar=
par(mar=c(0,0,0,0))
plot(as.phylo(hc), type = "fan")
########


#dd <- dist(scale(init.dendro), method = "euclidean")
dd <- dist(init.dendro, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
###  ----------------------------- data 
dd <- dist(t(lcpm), method = "euclidean")
#dd <- dist(init.dendro, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")


#plot(hc)
library(ape)
plot(hc, hang = -1, cex = 0.6)

plot(as.phylo(hc))

plot(as.phylo(hc), type = "fan")

###
dd <- dist(u.plot[,9:10], method = "euclidean")
hc <- hclust(dd, method = "ward.D2")
plot(hc, hang = -1, cex = 0.6)
####

