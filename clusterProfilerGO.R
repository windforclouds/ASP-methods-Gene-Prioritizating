##运用clusterProfiler来对获取到的gene list做GO,KEGG分析并绘图
library("clusterProfiler")
library("org.Hs.eg.db")
gene<-read.csv("string 61+10neighbor.txt",header = F)
gene<-bitr(t(gene),fromType = "SYMBOL",toType = "ENTREZID",OrgDb="org.Hs.eg.db")

##得到MF,CC,BP以及汇总的GO信息
ego_MF <- enrichGO(OrgDb="org.Hs.eg.db",gene=gene$ENTREZID,keyType = "ENTREZID",pvalueCutoff = 0.05,ont = "MF")
ego_CC <- enrichGO(OrgDb="org.Hs.eg.db",gene=gene$ENTREZID,keyType = "ENTREZID",pvalueCutoff = 0.05,ont = "CC")
ego_BP <- enrichGO(OrgDb="org.Hs.eg.db",gene=gene$ENTREZID,keyType = "ENTREZID",pvalueCutoff = 0.05,ont = "BP")
ego_ALL <- enrichGO(OrgDb="org.Hs.eg.db",gene=gene$ENTREZID,keyType = "ENTREZID",pvalueCutoff = 0.05,ont = "ALL")

ego_MF<-as.data.frame(ego_MF@result)
ego_CC<-as.data.frame(ego_CC@result)
ego_BP<-as.data.frame(ego_BP@result)
ego_ALL<-as.data.frame(ego_ALL@result)

##将GO,KEGG信息写出
write.table(ego_BP,"ego_BP61.txt",quote = F,row.names = F, sep = "\t" )
write.table(ego_CC,"ego_CC61.txt",quote = F,row.names = F, sep = "\t" )
write.table(ego_MF,"ego_MF61.txt",quote = F,row.names = F, sep = "\t" )
write.table(ego_ALL,"ego_ALL61.txt",quote = F,row.names = F, sep = "\t" )

MSKEGG <- enrichKEGG(gene = gene$ENTREZID,organism = "hsa", pvalueCutoff = 0.05)
SUMKEGG <-MSKEGG@result
symbolKEGG<- setReadable(MSKEGG, OrgDb = org.Hs.eg.db, keyType="ENTREZID")  #将gene ID转为symbol
symbolKEGG <-symbolKEGG@result 
write.table(enrichKEGG,"enrichKEGG61.txt",quote = F,row.names = F ,sep = "\t" )
  
ego_MF <- as.data.frame(ego_MF)[1:display_number[1], ]
ego_CC <- as.data.frame(ego_CC)[1:display_number[2], ]
ego_BP <- as.data.frame(ego_BP)[1:display_number[3], ]

ego_MF <- ego_MF[order(ego_MF$Count),]
ego_CC <- ego_CC[order(ego_CC$Count),]
ego_BP <- ego_BP[order(ego_BP$Count),]
display_number = c(15, 10, 15)

##简历一个分组用数据框,后续绘图可用
go_enrich_df <- data.frame(ID=c(ego_BP$ID, ego_CC$ID, ego_MF$ID),
                           Description=c(ego_BP$Description, ego_CC$Description, ego_MF$Description),
                           GeneNumber=c(ego_BP$Count, ego_CC$Count, ego_MF$Count),
                           type=factor(c(rep("biological process", display_number[1]), rep("cellular component", display_number[2]),
                           rep("molecular function", display_number[3])), levels=c("molecular function", "cellular component", "biological process")))

go_enrich_df$number <- factor(rev(1:nrow(go_enrich_df)))

##处理GO_term方便后续的绘图
shorten_names <- function(x, n_word=4, n_char=40){
  if (length(strsplit(x, " ")[[1]]) > n_word || (nchar(x) > 40))
  {
    if (nchar(x) > 40) x <- substr(x, 1, 40)
    x <- paste(paste(strsplit(x, " ")[[1]][1:min(length(strsplit(x," ")[[1]]), n_word)],
                     collapse=" "), "...", sep="")
    return(x)
  } 
  else
  {
    return(x)
  }
}

labels=(sapply(
  levels(go_enrich_df$Description)[as.numeric(go_enrich_df$Description)],
  shorten_names))
names(labels) = rev(1:nrow(go_enrich_df))

##绘图
CPCOLS <- c("#8DA1CB", "#FD8D62", "#66C3A5")
p <- ggplot(data=go_enrich_df, aes(x=number, y=GeneNumber, fill=type)) +
  geom_bar(stat="identity", width=0.8) + coord_flip() + 
  scale_fill_manual(values = CPCOLS) + theme_bw() + 
  scale_x_discrete(labels=labels) +
  xlab("GO term") + 
  theme(axis.text=element_text(face = "bold", color="gray50")) +
  labs(title = "The Most Enriched GO Terms")

p
