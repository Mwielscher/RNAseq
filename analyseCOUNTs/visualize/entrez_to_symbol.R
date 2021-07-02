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
