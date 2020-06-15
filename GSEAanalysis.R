#.libPaths("")
## try http:// if https:// URLs are not supported
#source("https://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade") ## you may need this
#biocLite("DOSE")
#biocLite("clusterProfiler")
#biocLite("ggjoy")
#cnetplot,dotplot,joyplot,enrichMap,barplot
#biocLite("GSEABase")
#biocLite("graph")

library(DOSE)
library(clusterProfiler)
library(ggplot2)
library(ggjoy)

setwd("./GSEA")

#Monocytes
#ranked_genelist by log2Foldchange
gene <- read.table("Mono_ATACfakerexpr_pvalue_FD.txt",header = T)
gene <- gene[!duplicated(gene$genesymbol), ]
gene2 <- gene[,2]
names(gene2) <- gene$genesymbol
#genelist <- read.table("tumor_versus_normal.txt",sep ="\t") #
library("GSEABase")
library("annotate")
kegg <- read.gmt("RA.gmt") #kegg
res <- GSEA(gene2,TERM2GENE = kegg,pvalueCutoff = 1)
write.table(res@result,"Mono_pvalue_RA.txt",quote = FALSE,row.names = T)
#gseaplot
pdf(file="Mono_pvalue_RAgsea.pdf")
gseaplot(res,geneSetID = "Rheumatoid_Arthritis",title = "Rheumatoid_Arthritis")
dev.off()
#dotplot
pdf(file="Mono_pvalue_RAdotplot.pdf",width = 6)
dotplot(res,showCategory = 8,colorBy="pvalue",split=".sign",font.size = 6) + facet_grid(.~.sign) #
dev.off()
