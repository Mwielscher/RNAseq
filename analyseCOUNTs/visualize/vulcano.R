library(EnhancedVolcano)
res=read.table("FINAL_RESULT_RNA/FINAL_LIST_PRETREATMENT_ApoA1.txt",header=T,sep="\t")

colnames(res)[3]=c('log2FoldChange')
colnames(res)[6]=c('pvalue')
rownames(res)=as.character(res$SYMBOL)
v.labels=as.character(res$SYMBOL[1:50])

p=EnhancedVolcano(res,
                lab = rownames(res),
                selectLab=v.labels,
                drawConnectors = TRUE,
                arrowheads = FALSE,
                labSize = 2,
                x = 'log2FoldChange',
                y = 'pvalue')


###       ------------------------------ VULCOANO used
options(ggrepel.max.overlaps = Inf)
EnhancedVolcano(res,
                lab = rownames(res),
                drawConnectors = TRUE,
                arrowheads = FALSE,
                ylim = c(0, max(-log10(res$pvalue), na.rm = TRUE) + 2),
                labSize = 2,
                pCutoff = 0.002020671,
                FCcutoff = 0.9,
                cutoffLineType = 'blank',
                widthConnectors = 0.3,
                hline = c(0.0000598),
                hlineCol = c("black"),
                hlineType = c("longdash"),
                vline =c(-1,1),
                vlineCol = c("black"),
                vlineType = c("longdash"),
                x = 'log2FoldChange',
                y = 'pvalue')
