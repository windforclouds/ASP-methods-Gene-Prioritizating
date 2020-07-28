#date:2019/5/27
#author:yj
library(igraph)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RCurl)
library(stringr)
library(XML)
node_all<-read.csv("61+156_01.csv",header = T,sep = ",")
groupA<-as.vector(subset(node_all,aorb == 1|aorb == 2)[,1])
groupB<-as.vector(subset(node_all,aorb == 0|aorb == 2)[,1]) 
group_share<-as.vector(subset(node_all,aorb == 2)[,1])
groupB_noshare<-as.vector(subset(node_all,aorb == 0)[,1])

edge_all<-read.csv("61+162_edge.csv",header = T,sep = ",")
net_all<- make_graph(t(edge_all),directed = FALSE)

##计算share节点的ASPvalue
groupB_avesp<-data.frame(symbol=c(NA),sum_SPL=c(NA),sum_Inf=c(NA),ASPLvalue=c(NA))
i=0
for (node in groupB) {
  i=i+1
  groupB_avesp[i,1]<-node
  if (node %in% group_share) {  ##计算share节点的ASPvalue
    net_minus<-net_all-vertices(groupB_noshare)}else{
      net_minus<-net_all-vertices(setdiff(groupB_noshare,node))
    }
  net_final<-net_minus-V(net_minus)[degree(net_minus)==0]
  assign(paste("MS",node,sep = "_"),net_minus)
  SPMS<-assign(paste("SPMS",node,sep = "_"),as.data.frame(shortest.paths(net_minus)))##we get dataframes of shortestpaths
  if (node %in% rownames(SPMS)) {
    spfinal<-SPMS[,node]
    spsum<-sum(spfinal[spfinal!=Inf])
    groupB_avesp[i,2]<-spsum
    groupB_avesp[i,3]<-sum(spfinal==Inf)
    groupB_avesp[i,4]<-(spsum+sum(spfinal==Inf)*diameter(net_minus))/(nrow(SPMS)-1) ##给与未能达到的点乘以该网络直径
    rm(spfinal,spsum,SPMS)
  }else{
    groupB_avesp[i,4]<-Inf
  }
}

groupgene<- bitr(groupB_avesp$symbol, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
groupB_avesp$ENTREZID<-groupgene$ENTREZID
groupB_avesp$NCBI_url <- paste("https://www.ncbi.nlm.nih.gov/gene/",groupB_avesp$ENTREZ)
groupB_avesp$NCBI_url <- gsub(" ","",groupB_avesp$NCBI_url)

#构建getNodesTxt，dealNodeTxt函数爬取NCBI
# 根据xpath获取节点内容：
getNodesTxt <- function(html_txt1,xpath_p){
  els1 = getNodeSet(html_txt1, xpath_p)
  # 获得Node的内容，并且去除空字符：
  els1_txt <- sapply(els1,xmlValue)[!(sapply(els1,xmlValue)=="")]
  # 去除\n：
  str_replace_all(els1_txt,"(\\n )+","")
}

# 处理节点格式，为character且长度为0的赋值为NA：
dealNodeTxt <- function(NodeTxt){
  ifelse(is.character(NodeTxt)==T && length(NodeTxt)!=0 , NodeTxt , NA)
}

# xpath精确定位：
for(i in 1:nrow(groupB_avesp)){
  # 获得网址：
  doc <- getURL(groupB_avesp[i,"NCBI_url"])
  cat("成功获得网页！\t")
  # 获得网页内容
  html_txt1 = htmlParse(doc, asText = TRUE)
  # 获得Full Name:
  groupB_avesp[i,"FullName"] <- dealNodeTxt(getNodesTxt(html_txt1,'//*[@id="summaryDl"]/dd[2]/text()'))
  # 获得HGNC ID:
  groupB_avesp[i,"HGNC_ID"] <- str_replace_all(getNodesTxt(html_txt1,'//*[@id="summaryDl"]/dd[preceding-sibling::dt[text()="Primary source" and position()=1 ] ]')," |HGNC|:","")
  cat("写入HGNC_ID\t")
  # 获得Gene type:
  groupB_avesp[i,"GeneType"] <- dealNodeTxt(getNodesTxt(html_txt1,'//*[@id="summaryDl"]/dd[preceding-sibling::dt[text()="Gene type" and position()=1 ] ]'))
  cat("写入GeneType\t")
  # 获得summary：
  groupB_avesp[i,"Summary"] <- dealNodeTxt(getNodesTxt(html_txt1,'//*[@id="summaryDl"]/dd[preceding-sibling::dt[text()="Summary" and position()=1 ] ]'))
  cat("写入Summary\n")
  print(paste("完成第",i,"个了！"))
}

write.csv(groupB_avesp,"groupB_avesp.csv",row.names = F)
