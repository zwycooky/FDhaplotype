options(stringAsFactor=F)

svdat <- read.table("syri.sniffiles.merged.sv.final.txt",header=F)
dupdat <- read.table("syri.DUP.sv.txt",header=F)
invdat <- read.table("syri.INV.sv.txt",header=F)

DELsv <- svdat[svdat[,8] == "DEL",]
INSsv <- svdat[svdat[,8] == "INS",]

DEL_len <- sum(DELsv[,7])
INS_len <- sum(INSsv[,7])
DUP_len <- sum(dupdat[,7])
INV_len <- sum(invdat[,7])

final_merged_sv <- rbind(as.matrix(svdat),as.matrix(dupdat),as.matrix(invdat))
final_merged_sv <- final_merged_sv[order(final_merged_sv[,1],as.numeric(final_merged_sv[,2])),]

final_merged_DID_sv <- rbind(as.matrix(svdat),as.matrix(dupdat))
final_merged_DID_sv <- final_merged_DID_sv[order(final_merged_DID_sv[,1],as.numeric(final_merged_DID_sv[,2])),]

write.table(final_merged_sv,"FDDB_AB_all4type_SV.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(final_merged_DID_sv,"FDDB_AB_DEL.INS.DUP_SV.txt",row.names=F,col.names=F,sep="\t",quote=F)


## start plot ##
mycol <- colorRampPalette(c("#80A984","#125560"))
mycol_vec <- mycol(4)

#pie(c(DEL_len,INS_len,DUP_len,INV_len),border=NA,col=mycol_vec,labels=c("DEL","INS","DUP","INV"))

SVname <- c("DEL","INS","DUP","INV")
SVlen <- c(DEL_len,INS_len,DUP_len,INV_len)
SVnum <- c(nrow(DELsv), nrow(INSsv), nrow(dupdat), nrow(invdat))

A <- data.frame(SVlen,SVname)
A$SVname <- factor(A$SVname,levels=A$SVname)

library(ggplot2)
library(ggforce)
 
ggplot()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(), 
        axis.text.y = element_blank(),
        axis.text.x = element_blank(),
        legend.title=element_blank(), 
        panel.border = element_blank(),
        panel.background = element_blank())+#去除没用的ggplot背景，坐标轴
  xlab("")+ylab('')+#添加颜色
  scale_fill_manual(values = mycol_vec) +
  geom_arc_bar(data=A,
               stat = "pie",
               aes(x0=0,y0=0,r0=1,r=2,
                   amount=SVlen,fill=SVname)
  )

