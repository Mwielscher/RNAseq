library(gplots)
hclustfunc <- function(x) hclust(x, method="complete")
distfunc <- function(x) dist(x,method="euclidean")
mycol <- colorpanel(1000,"blue","white","red")
# LPS_ApoA1_100_24_vs_contr_LPS_24 is a vector for creating the upset plot

heatmap.2(as.matrix(plot.d1[LPS_ApoA1_1000_24_vs_contr_LPS_24,]), Colv=FALSE, labCol=as.character(phe1$outcome),key=F,col=mycol,scale="row",
          dendrogram="row",trace="none",  margins=c(6,9), lwid = c(1,10), lhei = c(3,17), 
          hclust=hclustfunc,distfun=distfunc)
