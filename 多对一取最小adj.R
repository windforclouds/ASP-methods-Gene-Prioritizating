###按照adj.P最小来进行取值
options(stringsAsFactors = F)
GSE30784 <- read.table("GSE30784new.txt", header = T, row.names = 1, sep = "\t", na.strings = "")

TCGARna <- read.csv("TCGAresultsRNA.csv", header = T, row.names = 1)
#按照padj和logFC提取差异基因
TCGARnaDE <- subset(TCGARna, TCGARna$padj < 0.01 & abs(TCGARna$log2FoldChange) > 1 )
write.csv(TCGARnaDE, "TcgaDeRna.csv")

GSE30784 <- GSE30784[!is.na(GSE30784$Gene.symbol) & !grepl('//', GSE30784$Gene.symbol),]
GSE30784 <- subset(GSE30784, GSE30784$adj.P.Val < 0.01 & abs(GSE30784$logFC) > 1)

tail(GSE38010$logFC)
#这里写了一个函数来提取对应同一个的多个探针里adj.P.Val最小的那个基因
probeMerge <- function(x){
  #使用了aggregate函数根据Gene.symbol进行分组，对每一个分组取adj.P.Val的最小值
  probe2symbol= aggregate(x[, c("adj.P.Val", "symbol")], by=list(x[,"symbol"]),min)
  #使用merge函数，根据c("adj.P.Val", "Gene.symbol")进行匹配
  return(merge(x, probe2symbol, by = c("adj.P.Val", "Gene.symbol")))
}
GSE38010uniq <- probeMerge(GSE30784)[, -c(7, 8)]#调用函数，同时去除基因注释


#查看了一下TCGA RNASeq分析的脚本
#> levels(dds$Group)
#[1] "01" "11"
#所以原始分析的比较是正常配对-肿瘤样本，因为DESeq2分析运行起来非常慢（1-2h），不再改代码重新运算
#直接对这里logFC取相反数
upGeneTCGA <- rownames(TCGARnaDE[TCGARnaDE$log2FoldChange < 0,])
dnGeneTCGA <- rownames(TCGARnaDE[TCGARnaDE$log2FoldChange > 0,])

upGene30784 <- GSE30784uniq[GSE30784uniq$logFC > 0,]$Gene.symbol
dnGene30784 <- GSE30784uniq[GSE30784uniq$logFC < 0,]$Gene.symbol

upGene <- intersect(upGene30784, upGeneTCGA)
dnGene <- intersect(dnGene30784, dnGeneTCGA)

write.csv(upGene, "upGene.csv")
write.csv(dnGene, "dnGene.csv")


iqrs <- apply(exprs(eset), 1, IQR)
library(hgu133plus2.db)
library(genefilter)
prbs <- findLargest(featureNames(eset), testStat = iqrs, data = "hgu133plus2.db")
eset2 <- eset[prbs, ]
dim(eset2)
exprSet <- exprs(eset2)