options(stringsAsFactors=F)

library(plotrix)

Args <- commandArgs(T)

input <- Args[1]
output <- Args[2]

dat <- read.table(input, header=F)
dat <- dat[order(dat[,1]),]

contigs <- unique(dat[,1])

tar_length_vec <- NULL
for (i in 1:length(contigs)) {
	tmp_dat <- dat[dat[,1] == contigs[i],]
	tmp_tar_len <- tmp_dat[1,3]
	tar_start <- min(tmp_dat[,5:6])
	tar_end <- max(tmp_dat[,5:6])
	
	tar_len <- tar_end - tar_start
	tar_length_vec <- c(tar_length_vec,tar_len)
}
names(tar_length_vec) <- contigs

## set plot len ##
plot_len <- 6

if (tar_length_vec[1] / tar_length_vec[2] > 5) {
	left_ratio <- 5
	right_ratio <- 1
}else if (tar_length_vec[1] / tar_length_vec[2] <= 5){
	left_ratio <- 1
	right_ratio <- 5
}else{
	left_ratio <- plot_len * (tar_length_vec[1] / (tar_length_vec[1] + tar_length_vec[2]))
	right_ratio <- plot_len - left_ratio
}

# query plot len #
q_plot_len <- 6.6
q_start <- min(dat[,7:8])
q_end <- max(dat[,7:8])
q_len <- q_end - q_start

# function format #
format_pos <- function(start,pos,len,ratio) {
	return(ratio * ((pos - start) / len))
}

## start plot ##
pdf(output,width=5.5,height=2.8)

plot (0,type="n",axes=F,xlab="",ylab="",xlim=c(0,6.6),ylim=c(1,2))

for (i in 1:2) {
	tmp_dat <- dat[dat[,1] == contigs[i],]
	tar_start <- min(tmp_dat[,5:6])
	tar_end <- max(tmp_dat[,5:6])
	
	contig_len <- tmp_dat[1,3]
	
	contigs_id <- names(tar_length_vec)[i]
	
	for (j in 1:nrow(tmp_dat)) {
		tmp_tar_start <- tmp_dat[j,5]
		tmp_tar_end <- tmp_dat[j,6]
		tmp_q_start <- tmp_dat[j,7]
		tmp_q_end <- tmp_dat[j,8]
		
		# format tar pos #
		if (i == 1) {
			tmp_tar_start <- format_pos(tar_start, tmp_tar_start, tar_length_vec[i], left_ratio)
			tmp_tar_end <- format_pos(tar_start, tmp_tar_end, tar_length_vec[i], left_ratio)
		}else{
			tmp_tar_start <- left_ratio + 0.6 + format_pos(tar_start, tmp_tar_start, tar_length_vec[i], right_ratio)
			tmp_tar_end <- left_ratio + 0.6 + format_pos(tar_start, tmp_tar_end, tar_length_vec[i], right_ratio)
		}
		
		# format q pos #
		tmp_q_start <- format_pos(q_start, tmp_q_start, q_len, q_plot_len)
		tmp_q_end <- format_pos(q_start, tmp_q_end, q_len, q_plot_len)
		
		polygon(x=c(tmp_tar_start,tmp_q_start,tmp_q_end,tmp_tar_end), y=c(2,1,1,2), col=grey(0.8,alpha=0.5), xpd=NA, border=NA)
		
	}
	
	if (i == 1) {
		tar_start_plot <- format_pos(tar_start, tar_start, tar_length_vec[i], left_ratio)
		tar_end_plot <- format_pos(tar_start, tar_end, tar_length_vec[i], left_ratio)
	}else{
		tar_start_plot <- left_ratio + 0.6 + format_pos(tar_start, tar_start, tar_length_vec[i], right_ratio)
		tar_end_plot <- left_ratio + 0.6 + format_pos(tar_start, tar_end, tar_length_vec[i], right_ratio)
	}
	
	lines(x=c(tar_start_plot,tar_end_plot), y=c(2,2), lwd=3)
	ymult <- getYmult()
	if (i == 1) {
		lines(x=c(tar_end_plot-0.07,tar_end_plot+0.07), y=c(2+0.07*ymult,2-0.07*ymult), lwd=3, xpd=NA)
		lines(x=c(tar_end_plot+0.07,tar_end_plot+0.21), y=c(2+0.07*ymult,2-0.07*ymult), lwd=3, xpd=NA)
		
		left_remain <- contig_len - tar_end
		if (left_remain > 1000000) {
			left_remain <- paste(round(left_remain / 1000000,3), "Mb", sep=" ")
		}else{
			left_remain <- paste(round(left_remain / 1000,3), "Kb", sep=" ")
		}
		text(x=tar_end_plot, y=1.8, labels=left_remain, adj=1, xpd=NA)
	}else{
		lines(x=c(tar_start_plot-0.07,tar_start_plot+0.07), y=c(2+0.07*ymult,2-0.07*ymult), lwd=3, xpd=NA)
		lines(x=c(tar_start_plot-0.21,tar_start_plot-0.07), y=c(2+0.07*ymult,2-0.07*ymult), lwd=3, xpd=NA)
		
		right_remain <- tar_start
		if (right_remain > 1000000) {
			right_remain <- paste(round(right_remain / 1000000,3), "Mb", sep=" ")
		}else{
			right_remain <- paste(round(right_remain / 1000,3), "Kb", sep=" ")
		}
		text(x=tar_start_plot, y=1.8, labels=right_remain, adj=0, xpd=NA)
	}
	
	if (tar_length_vec[i] > 1000000) {
		tar_len_format <- paste(paste(contigs_id, ":",sep=""), round(tar_length_vec[i] / 1000000,3), "Mb", sep=" ")
	}else{
		tar_len_format <- paste(paste(contigs_id, ":",sep=""), round(tar_length_vec[i] / 1000,3), "Kb", sep=" ")
	}
	
	text(x=(tar_start_plot + tar_end_plot)/2, y=2.2, labels=tar_len_format, xpd=NA)
	
}

lines(x=c(0,6.6), y=c(1,1), lwd=3, xpd=NA)
query_id <- dat[1,2]
if (q_len > 1000000) {
	q_len_format <- paste(paste(query_id, ":",sep=""), round(q_len / 1000000,3), "Mb", sep=" ")
}else{
	q_len_format <- paste(paste(query_id, ":",sep=""), round(q_len / 1000,3), "Kb", sep=" ")
}
text(x=3.3, y=0.8, labels=q_len_format, xpd=NA)

dev.off()
