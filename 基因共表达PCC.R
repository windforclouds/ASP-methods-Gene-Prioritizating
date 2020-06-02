#2019/6/13 yj
#基于芯片数据找到27genes的表达值做皮尔森相关系数
setwd("E:/BaiduNetdiskDownload/RStudio/Ryj/GSE26484_RAW") #设置工作路径
gse26484<-read.csv("GSE26484.txt",header = T,sep = "\t") #载入标准化去重注释好gene symbol的表达矩阵
gene27<-read.csv("gene27.txt",header = T,sep = "\t")  #载入27genes
gse26484a<-merge(gse26484,gene27,by="X")
rownames(gse26484a)<-gse26484a[,1]
gse26484a<-t(gse26484a[,-1])
corgse26484<-as.data.frame(cor(gse26484a,method = "pearson"))

library(Hmisc)
pgse26484 <- rcorr(as.matrix(gse26484a), type = "pearson")
#library(psych)  ##psych包中提供的corr.test()函数可以一次算完
#pgse26484<-corr.test(gse26484a,method = "pearson",use = "complete")
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat) 
  data.frame(row = rownames(cormat)[row(cormat)[ut]],column = rownames(cormat)[col(cormat)[ut]], cor =(cormat)[ut], p = pmat[ut] )
}##用来将相关系数表和p值表整合在一起
pgse26484<-as.data.frame(pgse26484$P)
PCC26484<-flattenCorrMatrix(corgse26484,pgse26484)
write.csv(PCC26484,"PCC26484.csv")


#a<-(gse26484a[,1]-mean(gse26484a[,1]))/sd(gse26484a[,1])
#b<-(gse26484a[,2]-mean(gse26484a[,2]))/sd(gse26484a[,2])
#pcc26484<-sum(a*b)/19