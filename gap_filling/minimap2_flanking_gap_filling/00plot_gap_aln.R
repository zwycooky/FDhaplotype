options(stringsAsFactors=F)

Args <- commandArgs(T)
input <- Args[1]
output <- Args[2]

dat <- read.table(input,header=F)

flanking_len <- 100000

pdf(output,width=5,height=2.5)
par(mar=c(3,1,3,1))
plot(0,type="n",axes=F,xlab="",ylab="",xlim=c(0,flanking_len*2 + 2000),ylim=c(0,1))


min_aln <- min(dat[,8:9])
max_aln <- max(dat[,8:9])

aln_len <- max_aln - min_aln
format_index <- (flanking_len*2+2000) / aln_len

for (i in 1:nrow(dat)) {
	if (grepl("L",dat[i,1])) {
		if (dat[i,5] == "+") {
			x1 <- dat[i,3]
			x2 <- dat[i,4]
			x3 <- (dat[i,8] - min_aln) * format_index
			x4 <- (dat[i,9] - min_aln) * format_index
			polygon(x=c(x1,x2,x4,x3),y=c(1,1,0,0),border=NA,col=grey(0.7,alpha=0.5),xpd=NA)
		}else{
			x1 <- dat[i,3]
			x2 <- dat[i,4]
			x3 <- (dat[i,8] - min_aln) * format_index
			x4 <- (dat[i,9] - min_aln) * format_index
			polygon(x=c(x1,x2,x3,x4),y=c(1,1,0,0),border=NA,col=grey(0.7,alpha=0.5),xpd=NA)
		}
	}else{
		if (dat[i,5] == "+") {
			x1 <- dat[i,3] + flanking_len + 2000
			x2 <- dat[i,4] + flanking_len + 2000
			x3 <- (dat[i,8] - min_aln) * format_index
			x4 <- (dat[i,9] - min_aln) * format_index
			polygon(x=c(x1,x2,x4,x3),y=c(1,1,0,0),border=NA,col=grey(0.7,alpha=0.5),xpd=NA)
		}else{
			x1 <- dat[i,3] + flanking_len + 2000
			x2 <- dat[i,4] + flanking_len + 2000
			x3 <- (dat[i,8] - min_aln) * format_index
			x4 <- (dat[i,9] - min_aln) * format_index
			polygon(x=c(x1,x2,x3,x4),y=c(1,1,0,0),border=NA,col=grey(0.7,alpha=0.5),xpd=NA)
		}
	}
}
lines(x=c(0,flanking_len),y=c(1,1),xpd=NA,lwd=3)
lines(x=c(flanking_len + 2000,flanking_len*2 + 2000),y=c(1,1),xpd=NA,lwd=3)
lines(x=c((min_aln-min_aln)*format_index,(max_aln-min_aln)*format_index),y=c(0,0),xpd=NA,lwd=3)

text(x=flanking_len/2,y=1.1,labels=paste("gap left ",flanking_len," Kb",sep=""),xpd=NA)
text(x=flanking_len/2+flanking_len+2000,y=1.1,labels=paste("gap right ",flanking_len," Kb",sep=""),xpd=NA)

text(x=(flanking_len*2+2000)/2,y=-0.1,labels=paste(dat[1,6]," ",round(aln_len/1000)," Kb",sep=""), xpd=NA)
dev.off()
