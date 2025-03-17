options(stringAsFactor=F)

syri_dat <- read.table("syri.minimap2.filtered.sv.1.txt",header=F)
sniffiles_dat <- read.table("sniffiles.filtered.SV.1.txt",header=F)

chr_vec <- unique(syri_dat[,1])

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

find_overlap_between_syri_sniffiles <- function(syri,sniffiles) {
	
	remove <- NULL
	new <- NULL
	overlaped <- NULL
	
	for (i in 1:nrow(sniffiles)) {
		s <- sniffiles[i,1]
		e <- sniffiles[i,2]
		overlap_nums <- which((s >= syri[,1] & s <= syri[,2]) | (e >= syri[,1] & e <= syri[,2]) | (s <= syri[,1] & e >= syri[,2]))
		if (length(overlap_nums) > 1) {
			remove <- c(remove, i)
		}else if(length(overlap_nums) == 0) {
			new <- c(new, i)
		}else if (length(overlap_nums) == 1) {
			overlaped <- c(overlaped, overlap_nums)
		}
	}
	
	res <- list(remove = remove, new = new, overlaped = overlaped)
	return(res)
	
}


total_iterations <- length(chr_vec)
pb <- txtProgressBar(min = 0, max = total_iterations, style = 3)

merged_svdat <- NULL
for (i in 1:length(chr_vec)) {
	tmp_syri <- syri_dat[syri_dat[,1] == chr_vec[i],]
	tmp_sniffiles <- sniffiles_dat[sniffiles_dat[,1] == chr_vec[i],]
	
	tmp_syri_ins <- tmp_syri[tmp_syri[,8] == "INS",]
	tmp_syri_del <- tmp_syri[tmp_syri[,8] == "DEL",]
	
	tmp_sniffiles_ins <- tmp_sniffiles[tmp_sniffiles[,8] == "INS",]
	tmp_sniffiles_del <- tmp_sniffiles[tmp_sniffiles[,8] == "DEL",]
	
	overlaped_syri_sniffiles_ins <- find_overlap_between_syri_sniffiles(tmp_syri_ins[,5:6],tmp_sniffiles_ins[,5:6])
	overlaped_syri_sniffiles_del <- find_overlap_between_syri_sniffiles(tmp_syri_del[,2:3],tmp_sniffiles_del[,2:3])
	
	# merge ins #
	tmp_syri_ins <- cbind(tmp_syri_ins,1,0)
	tmp_syri_ins[unique(overlaped_syri_sniffiles_ins$overlaped),10] <- 1
	sniffiles_ins_new <- tmp_sniffiles_ins[overlaped_syri_sniffiles_ins$new,]
	sniffiles_ins_new <- cbind(sniffiles_ins_new,0,1)
	
	# merge del #
	tmp_syri_del <- cbind(tmp_syri_del,1,0)
	tmp_syri_del[unique(overlaped_syri_sniffiles_del$overlaped),10] <- 1
	sniffiles_del_new <- tmp_sniffiles_del[overlaped_syri_sniffiles_del$new,]
	sniffiles_del_new <- cbind(sniffiles_del_new,0,1)
	
	merged_svdat <- rbind(merged_svdat, as.matrix(tmp_syri_ins), as.matrix(sniffiles_ins_new), as.matrix(tmp_syri_del), as.matrix(sniffiles_del_new))
	setTxtProgressBar(pb, i)
}

close(pb)

merged_svdat <- merged_svdat[order(merged_svdat[,1],merged_svdat[,2]),]

write.table(merged_svdat,"syri.sniffiles.merged.sv.txt", row.names=F,col.names=F,sep="\t",quote=F)
