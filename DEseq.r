encoding="utf-8"

library(DESeq)
file="CRPstim.gene.txt"
CountTable = read.table (file,header =T ,row.names = 1 )
condition = factor (c("con","con","con","stim","stim","stim"))
cds = newCountDataSet(CountTable,condition)
cds=estimateSizeFactors(cds)
count_norm = counts(cds,normalized=TRUE)
write.table(count_norm,file="CRPstim_Norm.txt",row.names=TRUE,col.names=TRUE,quote=FALSE,sep="\t")
