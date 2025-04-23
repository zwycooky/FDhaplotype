options(stringsAsFactors=F)

library("progress")

met_dat <- read.table("met_Trait.forR.txt")
fuding_family <- read.table("fuding_family_pred_haps.merged.txt",header=T)
bins <- read.table("bins.txt",header=F)
bin_geno <- read.table("bin_geno.txt")
rownames(bins) <- rownames(bin_geno)
bin_geno <- t(bin_geno)

fuding_off <- unique(sub("\\.","-",sub("X","",fuding_family$accession)))
rownames(bin_geno) <- sub("\\.","-",sub("X","",rownames(bin_geno)))
met_dat <- met_dat[rownames(met_dat) %in% fuding_off,]

## order ##
bin_geno <- bin_geno[order(rownames(bin_geno)),]
met_dat <- met_dat[order(rownames(met_dat)),]

bin_met_res <- matrix(0,nrow=ncol(bin_geno)*ncol(met_dat),ncol=3)
line_count <- 0
pb <- progress_bar$new(total=ncol(bin_geno))
for (i in 1:ncol(bin_geno)) {
	geno <- bin_geno[,i]
	for (j in 1:ncol(met_dat)) {
		metTrait <- met_dat[,j]
		if (is.na(table(geno[!is.na(metTrait)])[1]) | is.na(table(geno[!is.na(metTrait)])[2])) {
			pvalue <- NA
		}else if (table(geno[!is.na(metTrait)])[1] < 3 | table(geno[!is.na(metTrait)])[2] < 3) {
			pvalue <- NA
		}else{
			pvalue <- summary(aov(metTrait~geno))[[1]][,5][1]
		}
		line_count = line_count + 1
		bin_met_res[line_count,] <- c(colnames(bin_geno)[i],colnames(met_dat)[j],pvalue)
	}
	pb$tick()
}
bin_met_res <- as.data.frame(bin_met_res)
bin_met_res[,3] <- as.numeric(bin_met_res[,3])
bin_met_sig_res <- bin_met_res[which(bin_met_res[,3] < 0.01),]
bin_met_sig2_res <- bin_met_res[which(bin_met_res[,3] < 0.001),]
write.table(bin_met_sig_res,file="bin_met_sig0.01_res.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(bin_met_sig2_res,file="bin_met_sig0.001_res.txt",row.names=F,col.names=F,sep="\t",quote=F)
write.table(bin_met_res,file="bin_met_res.txt",row.names=F,col.names=F,sep="\t",quote=F)

## -------------- plot -------------- ##
## -------------- plot -------------- ##
chr_uni <- unique(bins[,1])
chr_len_vec <- rep(0,length(chr_uni))
for (i in 1:length(chr_uni)) {
	tmp <- bins[bins[,1] == chr_uni[i],]
	chr_len_vec[i] <- as.numeric(tail(tmp,n=1)[3])
}
chr_len_vec <- c(0,chr_len_vec)

## format pos ##
format_bins <- bins
for (i in 1:nrow(format_bins)) {
	chr <- as.numeric(format_bins[i,1])
	format_bins[i,2] <- as.numeric(format_bins[i,2]) + sum(chr_len_vec[1:chr])
	format_bins[i,3] <- as.numeric(format_bins[i,3]) + sum(chr_len_vec[1:chr])
}

## get chr text pos ##
chr_text_pos_vec <- rep(0,length(chr_len_vec)-1)
for (i in 2:length(chr_len_vec)) {
	chr_text_pos_vec[i-1] <- chr_len_vec[i]/2 + sum(chr_len_vec[1:i-1])
}

## start plot ##
plot_dat <- format_bins[bin_met_res[,1],]
chr_class_col <- ifelse(as.logical(plot_dat[,1] %% 2),"#80A984","#125560")
ymax <- max(-log10(bin_met_res[,3]), na.rm=T)

library(Cairo)
Cairo::CairoPNG( 
  filename = "Met_binMarker.png",
  width = 10,
  height = 3.3,
  units = "in",
  dpi = 400)

x <- (plot_dat[,2] + plot_dat[,3]) / 2
y <- -log10(bin_met_res[,3])

plot(x, y, bty="l", pch=20, xlab="", ylab="-log10(pvalue)", xaxt="n", cex=0.8, col=chr_class_col)
#segments(x0=plot_dat[,2], x1=plot_dat[,3], y0=y, y1=y, col=chr_class_col)
#text(x=chr_text_pos_vec, y=-0.08*ymax, labels= paste("chr",1:15,sep=""),xpd=NA)
abline(h=c(2,3),col=c(grey(0.6),"black"),lty=2)

# plot mQTL number #
qtl_n <- rep(0,nrow(format_bins))
names(qtl_n) <- rownames(format_bins)
qtl_sig <- table(bin_met_sig2_res[,1])
qtl_n[names(qtl_sig)] <- qtl_sig

mycol_fun <- colorRampPalette(c("#80A984","#125560"))
mycol <- mycol_fun(length(qtl_n))[as.numeric(cut(as.numeric(qtl_n),breaks = length(qtl_n)))]

ymax_d <- ymax * 0.3
ra <- ymax_d / max(qtl_n)
for (i in 1:length(qtl_n)) {
	bin_id <- names(qtl_n)[i]
	b_s <- format_bins[bin_id,2]
	b_e <- format_bins[bin_id,3]
	qn <- qtl_n[i]
	
	rect(b_s, -qn*ra -0.1*ymax, b_e, -0.1 * ymax, col = mycol[i], xpd=NA, border= NA)
}

lines(x=c(0,0), y=c(-0.1 * ymax, -ymax_d),xpd=NA)
tick_vec <- seq(-0.1*ymax,-ymax_d,length=3)
tick_text <- round(seq(0,max(qtl_n),length=3))
for (i in 1:length(tick_vec)) {
	lines(x=c(0,-10000000), y=c(tick_vec[i], tick_vec[i]), xpd=NA)
	text (x=-15000000, y=tick_vec[i], labels=tick_text[i], xpd=NA, adj=1, cex=0.8)
}

text(x=chr_text_pos_vec,y=-0.5*ymax, labels=paste("chr",c(1:15),sep=""), xpd=NA)

dev.off()

x <- read.table("sig_bin.txt")
#met_ann <- read.table("met_ann.txt")
sig_qtl_16 <- bin_met_sig2_res[bin_met_sig2_res[,1] %in% names(qtl_n)[qtl_n >= 15],]
write.table(sig_qtl_16,"sig_bin_met.txt",row.names=F,col.names=F,sep="\t",quote=F)

