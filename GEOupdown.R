setwd('E:/BaiduNetdiskDownload/RStudio/Ryj/GSE38010_RAW')
GSE38010_output<- read.table("GSE38010_output.txt",header = T)
GSE38010up_p0.05_1.5 <- as.data.frame(GSE38010_output[GSE38010_output$logFC > log(1.5)/log(2) & GSE38010_output$P.Value<0.05,c(1,4)])
GSE38010dn_p0.05_1.5 <- as.data.frame(GSE38010_output[GSE38010_output$logFC < -log(1.5)/log(2) & GSE38010_output$P.Value<0.05,c(1,4)])

GSE38010up_p0.01 <- as.data.frame(GSE38010_output[GSE38010_output$logFC > 1 & GSE38010_output$P.Value<0.01,c(1,4)])
GSE38010dn_p0.01 <- as.data.frame(GSE38010_output[GSE38010_output$logFC < -1 & GSE38010_output$P.Value<0.01,c(1,4)])

GSE38010up_p0.05 <- as.data.frame(GSE38010_output[GSE38010_output$logFC > 1 & GSE38010_output$P.Value<0.05,c(1,4)])
GSE38010dn_p0.05 <- as.data.frame(GSE38010_output[GSE38010_output$logFC < -1 & GSE38010_output$P.Value<0.05,c(1,4)])

GSE38010up_p0.01_1.5 <- as.data.frame(GSE38010_output[GSE38010_output$logFC > log(1.5)/log(2) & GSE38010_output$P.Value<0.01,c(1,4)])
GSE38010dn_p0.01_1.5 <- as.data.frame(GSE38010_output[GSE38010_output$logFC < -log(1.5)/log(2) & GSE38010_output$P.Value<0.01,c(1,4)])

write.table(GSE38010up_p0.05_1.5, "GSE38010up_p0.05_1.5.txt",quote = F,sep = "\t")
write.table(GSE38010dn_p0.05_1.5, "GSE38010dn_p0.05_1.5.txt",quote = F,sep = "\t")


GSE38010up_p0.25 <- as.data.frame(GSE38010_output[GSE38010_output$logFC > log(1.5)/log(2) & GSE38010_output$P.Val<0.25,c(1,5)])
GSE38010dn_p0.25 <- as.data.frame(GSE38010_output[GSE38010_output$logFC < -log(1.5)/log(2) & GSE38010_output$P.Val<0.25,c(1,5)])

GSE38010up_0.05 <- as.data.frame(GSE38010_output[GSE38010_output$logFC > 1 & GSE38010_output$adj.P.Val<0.05,c(1,5)])
GSE38010dn_0.05 <- as.data.frame(GSE38010_output[GSE38010_output$logFC < -1 & GSE38010_output$adj.P.Val<0.05,c(1,5)])
write.table(GSE38010up_0.05, "GSE38010up_0.05.txt",quote = F,sep = "\t")
write.table(GSE38010dn_0.05, "GSE38010dn_0.05.txt",quote = F,sep = "\t")