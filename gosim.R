#计算gene GO语义相似性
#date:2019/6/20
#author:yj
setwd("D:/Science/MS/data/MS_PPI") #设置工作路径
library(igraph)
edge61_162<-read.csv("string_interactions 61+162.tsv",header = T,sep = "\t")
ms<- make_graph(t(edge61_162[,1:2]),directed = FALSE)
gene26<-read.csv("gene26.txt",header = F)
gene26a<-as.vector(t(gene26))
gene27<-read.csv("gene27.txt",header = F)
#install.packages("GOSemSim")
library(GOSemSim)
hsGO <- godata('org.Hs.eg.db', ont="BP")
hsGO2 <- godata('org.Hs.eg.db', keytype = "SYMBOL", ont="BP", computeIC=FALSE)

gene27_gosim<-gene27
assign(paste("MS","TNFAIP3",sep = "_"),ms-vertices(gene26a))
neib_TNFAIP3<-names(neighbors(MS_TNFAIP3,"TNFAIP3"))  #获取name即为邻居
sim<-mgeneSim(c("TNFAIP3",neib_TNFAIP3),semData=hsGO2, measure="Wang")
gene27_gosim[1,2]<-sum(sim[1,])-1#特殊节点
  
for (i in 1:26) {
  ms_minus<-ms-vertices(gene26a[-i])
  graph<-assign(paste("MS",gene26a[i],sep = "_"),ms_final)
  if (degree(ms_minus,gene26a[i])==0) {
    gene27_gosim[i+1,2]<-NA
    }else{
    neib<-names(neighbors(ms_minus,c(gene26a[i])))
    sim<-mgeneSim(c(gene26a[i],neib),semData=hsGO2, measure="Wang")
    gene27_gosim[i+1,2]<-sum(sim[1,])-1
   }
}
write.csv(gene27_gosim,"gene27_gosim.csv")
