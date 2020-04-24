#date:2019/5/27
#author:yj
library(igraph)
edge61_162<-read.csv("string_interactions 61+162.tsv",header = T,sep = "\t") ##读取上一步得到的PPIN文件
ms<- make_graph(t(edge61_162[,1:2]),directed = FALSE) ##使用igraph导入PPIN
#ave1<-all_shortest_paths(ms, 1, 1:110); 计算ms所有最短路径
#spms<-as.data.frame(shortest.paths(ms))##dataframe of all shortest paths of the network(spms)
gene26<-read.csv("gene26.txt",header = F) ##读取网络中已知MS相关基因(去除与数据库得到的基因重复部分)
gene26a<-as.vector(t(gene26))  ##转为向量格式
gene27<-read.csv("gene27.txt",header = F) ##读取网络中已知全部MS相关基因

##循环得到小网络中每个节点与大网络的平均最短路径
gene27_avesp<-gene27
for (i in 1:26) {
  assign(paste("MS","TNFAIP3",sep = "_"),ms-vertices(gene26a))  
  SPMS_TNFAIP3<-as.data.frame(shortest.paths(MS_TNFAIP3))
  #gene27_inf<-gene27
  #gene27_0<-gene27
  #gene27_0inf<-gene27
  #gene27_inf[1,2]<-sum(SPMS_TNFAIP3$TNFAIP3==Inf)
  #gene27_0inf[1,2]<-sum(SPMS_TNFAIP3$TNFAIP3==Inf)
  gene27_avesp[1,2]<-sum(SPMS_TNFAIP3[SPMS_TNFAIP3$TNFAIP3!=Inf,1])/(nrow(SPMS_TNFAIP3)-1) ##特殊节点
  ms_minus<-ms-vertices(gene26a[-i])
  ms_final<-ms_minus-V(ms_minus)[degree(ms_minus)==0]
  assign(paste("MS",gene26a[i],sep = "_"),ms_final)
  SPMS<-assign(paste("SPMS",gene26a[i],sep = "_"),as.data.frame(shortest.paths(ms_final)))##we get dataframes of shortestpaths
  if (gene26a[i] %in% rownames(SPMS)) {
    spfinal<-SPMS[,gene26a[i]]
    spsum<-sum(spfinal[spfinal!=Inf])
    if (sum(spfinal!=Inf)<=40) {          ##Inf可以自定义，我这里为40
      gene27_avesp[i+1,2]<-Inf}else{
      gene27_avesp[i+1,2]<-spsum/(nrow(SPMS)-1)
      rm(spfinal,spsum)
    }
  }else{
    gene27_avesp[i+1,2]<-Inf
  }
}
write.csv(gene27_avesp,"gene27_avesp.csv")
