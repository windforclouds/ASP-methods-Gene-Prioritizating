#date:2020/12/21
#author:yj
'''1.ready for one txt/csv file(such as MSPPIN_vertex.csv), first 
column is gene name, and the second colname is the type 
identification code 0, 1, and 2
(0 means groupB, 1 maens groupA, 2 means share genes);
2.Get your edges-network file ready(such as MSPPIN_edge.csv) '''

#set as working directory
#setwd("~/PS-V2N")

#if you have not installed igraph, please print install.package("igraph")
library(igraph)

node_all<-read.csv("MSPPIN_vertex.csv",header = T,sep = ",")
groupB<-as.vector(subset(node_all,aorb == 1|aorb == 2)[,1])
groupA<-as.vector(subset(node_all,aorb == 0|aorb == 2)[,1]) 
group_share<-as.vector(subset(node_all,aorb == 2)[,1])
groupA_noshare<-as.vector(subset(node_all,aorb == 0)[,1])
#build up graph by your edge informations of network
edge_all<-read.csv("MSPPIN_edge.csv",header = T,sep = ",")
net_all<- make_graph(t(edge_all),directed = FALSE)

#calulate share gene's self value, namely use ave_closenessB to measure
net_groupB<-net_all-vertices(groupA_noshare)
net_groupB<-net_groupB-V(net_groupB)[degree(net_groupB)==0]
SPgroupB<-as.data.frame(shortest.paths(net_groupB))
closenessB<-data.frame(symbol=c(NA),closeness=c(NA))
i=0  
for (a in rownames(SPgroupB)) {
  i=i+1
  closenessB[i,1]<-a
  CBi<-SPgroupB[a,]
  for (m in colnames(CBi)) {
    CBi[2,m]<-CBi[1,m]
  }
  closenessB[i,2]<-1/sum(CBi[2,][CBi[2,]!=Inf])
}
ave_closenessB<-mean(closenessB$closeness)


#calculate by PS-V2N methods
groupA_avesp<-data.frame(symbol=c(NA),num_V_newnet=c(NA),sum_SPL=c(NA),sum_Inf=c(NA),num_neighbor_in_B=c(NA),PSV2N=c(NA))
i=0
for (node in groupA) {
  i=i+1
  groupA_avesp[i,1]<-node
  if (node %in% group_share) {  
    net_minus<-net_all-vertices(groupA_noshare) ##rebuild net_minus
  }else{
    net_minus<-net_all-vertices(setdiff(groupA_noshare,node))
  }
  net_final<-net_minus-V(net_minus)[degree(net_minus)==0] ##remove outliers of net_minus
  #assign(paste("MS",node,sep = "_"),net_minus) ##生成net_minus网络
  SPMS<-assign(paste("SPMS",node,sep = "_"),as.data.frame(shortest.paths(net_minus)))##we get net_minus's dataframes of shortestpaths
  if (node %in% rownames(SPMS)) {
    spfinal<-SPMS[node,]
    ave_degree<-sum(degree(net_minus))/length(spfinal)  ##calculate net_minus's average degree
    if (node %in% group_share) { 
      for (m in colnames(spfinal)) {
        if (m == node) {
          spfinal[2,node]<-(1/ave_closenessB)
        }else{
          spfinal[2,m]<-(1/spfinal[1,m])
        }
      }
    }else{
      for (m in colnames(spfinal)) {
        spfinal[2,m]<-(1/spfinal[1,m])
      }
    }
    spsum<-sum(spfinal[2,][spfinal[2,]!=Inf])
    groupA_avesp[i,2]<-length(spfinal)
    groupA_avesp[i,3]<-spsum
    groupA_avesp[i,4]<-sum(spfinal[1,]==Inf)
    groupA_avesp[i,5]<-length(intersect(neighbors(net_minus,node)$name,groupB))
    if (node %in% group_share) {
      groupA_avesp[i,6]<-spsum/(nrow(SPMS))+(1/(sum(spfinal[1,]==Inf)+1))*(1/diameter(net_minus))/nrow(SPMS) ##nonshare genes
    }else{
      groupA_avesp[i,6]<-spsum/(nrow(SPMS)-1)+(1/(sum(spfinal[1,]==Inf)+1))*(1/diameter(net_minus))/(nrow(SPMS)-1) ##share genes
    }
    rm(spsum,SPMS)
  }else{
    groupA_avesp[i,6]<-Inf
  }
}

#the results is ordered by PS-V2N
groupA_avesp<-groupA_avesp[order(-groupA_avesp$PSV2N),]
order_by_PSV2N<-c(1:i)
groupA_avesp<-cbind(order_by_PSV2N,groupA_avesp)

write.csv(groupA_avesp,"groupA_avesp.csv",row.names = F)
