options(stringAsFactor=F)

svdat <- read.table("sniffiles.filtered.SV.txt",header=F)
svdat <- svdat[order(svdat[,1],svdat[,2]),]

chr_vec <- unique(svdat[,1])

filter_overlap <- function(dat) {
	remove <- NULL
	for (i in 1:nrow(dat)) {
		s <- dat[i,1]
		e <- dat[i,2]
		if (sum((s >= dat[,1] & s <= dat[,2]) | (e >= dat[,1] & e <= dat[,2]) | (s <= dat[,1] & e >= dat[,2])) > 1) {
			remove <- c(remove,i)
		}
	}
	return(remove)
}

total_iterations <- length(chr_vec)
pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)

filtered_svdat <- NULL
for (i in 1:length(chr_vec)) {
	tmp_dat <- svdat[svdat[,1] == chr_vec[i],]
	overlap_remove_vec1 <- filter_overlap(tmp_dat[,2:3])
	overlap_remove_vec2 <- filter_overlap(tmp_dat[,5:6])
	overlap_remove_vec <- unique(overlap_remove_vec1,overlap_remove_vec2)
	
	filtered_svdat <- rbind(filtered_svdat,tmp_dat[-overlap_remove_vec,])
	setTxtProgressBar(pb, i)
}

close(pb)

write.table(filtered_svdat,"sniffiles.filtered.SV.1.txt", row.names=F,col.names=F,sep="\t",quote=F)
