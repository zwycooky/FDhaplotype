options(stringsAsFactors=F)

bins <- read.table("bins.txt",header=F)
pred_hap <- read.table("fuding_family_pred_haps.merged.txt",header=T)

accession <- unique(pred_hap$accession)
bin_geno <- matrix(0,nrow = nrow(bins),ncol=length(accession))
for (i in 1:nrow(bins)) {
	chr <- bins[i,1]
	b_s <- bins[i,2]
	b_e <- bins[i,3]
	
	tmp_geno <- rep(NA,length(accession))
	names(tmp_geno) <- accession
	
	tmp_pred_hap <- pred_hap[pred_hap$chr == chr & pred_hap$start <= b_s & pred_hap$end >= b_e,]
	tmp_geno[tmp_pred_hap$accession] <- tmp_pred_hap$haps
	tmp_geno <- ifelse(tmp_geno==1,'A','B')
	bin_geno[i,] <- tmp_geno
}

rownames(bin_geno) <- paste('bin',1:nrow(bin_geno),sep="")
colnames(bin_geno) <- accession
write.table(bin_geno,"bin_geno.txt",row.names=T,col.names=T,sep="\t",quote=F)

