#!/usr/bin/perl

use strict;
use Cwd 'abs_path';

my ($sam,$phasing_data,$outprefix) = @ARGV[0,1,2];
my $Usage = "\n\t$0 <bam> <phased hap file> <output prefix>
\n";
die $Usage unless (@ARGV == 3);

#my $time1 = time();
## read hap data ##
my ($hap,$hap_pos);
open HAP,'<',"$phasing_data" or die "Error: Cannot open hap file:$!";
while (<HAP>) {
	chomp;
	my ($contig,$pos,$snp1,$snp2) = (split)[0,1,2,3];
	push @{$hap->{$contig}},"$pos\t$snp1\t$snp2";
	push @{$hap_pos->{$contig}}, $pos;
}
close HAP;

## start phasing pacbio/NGS reads ##
my ($pre_reads_id,$fir,$reads_tmp);
open my $sam_file,"samtools view -q 45 $sam|" or die "Error: Cannot open sam/bam file:$!";
while (<$sam_file>) {
	my ($reads_id,$contig) = (split /\t/,$_)[0,2];
	my $key = "$reads_id\t$contig";
	next if (!exists $hap_pos->{$contig});
	push @{$reads_tmp->{$key}},$_;
}
close $sam_file;

#open NOPHASING,'>',"$outprefix.nophase.txt";
open(my $phase1,"|gzip >$outprefix.hap1.fq.gz");
open(my $phase2,"|gzip >$outprefix.hap2.fq.gz");
#open PHASE1,'>',"$outprefix.hap1.txt";
#open PHASE2,'>',"$outprefix.hap2.txt";
#open LOWACC,'>',"$outprefix.low.accuracy.reads.txt";

foreach (keys %{$reads_tmp}) {
	chomp;

	my ($reads_id,$contig) = (split /\t/,$_)[0,1];
	my @tmp = @{$reads_tmp->{$_}};
	

	my $phasing_res = &phasing_reads(\@tmp, \@{$hap_pos->{$contig}}, \@{$hap->{$contig}});
	my $ifphased = (split /\t/,$phasing_res)[0];
	#my ($reads_seq,$reads_q) = (split /\t/,$tmp[0])[9,10];
		
	if ($ifphased > 0) {
		my ($phased_hap,$acc) = (split /\t/,$phasing_res)[0,1];
		if ($acc >= 90 && $phased_hap == 1) {
			#print PHASE1 "$reads_chr\t$reads_start\t$reads_end\t$phased_hap\t$hap_block\n";
			#print "$reads_chr\t$reads_start\t$reads_end\t$phased_hap\n";
			foreach (@tmp) {
				my ($rid,$rseq,$rq) = (split)[0,9,10];
				print $phase1 "\@$rid\n$rseq\n+\n$rq\n";
			}
		}elsif ($acc >= 90 && $phased_hap == 2) {
			foreach (@tmp) {
                                my ($rid,$rseq,$rq) = (split)[0,9,10];
                                print $phase2 "\@$rid\n$rseq\n+\n$rq\n";
                        }
			#print PHASE2 "$reads_chr\t$reads_start\t$reads_end\t$phased_hap\t$hap_block\n";
			#print "$reads_chr\t$reads_start\t$reads_end\t$phased_hap\n";
		}else {
			#print LOWACC "$reads_chr\t$reads_start\t$reads_end\t$phased_hap\t$hap_block\n";
		}
	}

}

close $phase1;
close $phase2;
#close NOPHASING;
#close PHASE1;
#close PHASE2;
#close LOWACC;

#my $time2 = time();
#my $time = ($time2 - $time1);
#print "$time s\n";

## sub progrem ##

sub phasing_reads {
	
	my ($tmp,$hap_pos,$hap) = @_;
	my $full_reads = 0;
	my ($reads_name,$not_phase);
	my ($accuracy,$snp1_count,$snp2_count,$total_count,$reads_len,$reads_chr,$reads_start,$reads_end);
	
	foreach (@{$tmp}) {
		my ($reads_id,$contig,$start,$MQ,$cigar,$seq) = (split /\t/,$_)[0,2,3,4,5,9];
		#print "$reads_id\n";
		#print length($seq) . "\n";
		# get full reads #
		if ($full_reads == 0) {
			$full_reads = $seq;
			$reads_name = $reads_id;
			$reads_len = length($full_reads);
		}
		$reads_start = $start;
		$reads_chr = $contig;

		# find overlap with phasing region #
		# get start and end position in contig with an alignment #
		my @cigar = &get_s_e_pos($start,$cigar);
		my $end = shift @cigar;
		$reads_end = $end;

		#print "start: $start\tend: $end\n";
			
		my $search_s = &binarySearch($start,$hap_pos);
		my $search_e = &binarySearch($end,$hap_pos);
		my @tmp_hap = @{$hap}[($search_s-1)..($search_e+1)];
			
		my $overlap = 0;
		foreach (@tmp_hap) {
			my ($pos,$snp1,$snp2) = (split);
			if ($pos >= $start && $pos <= $end) {
				$overlap = 1;
				## get base in reads ##
				my $relative_pos = $pos - $start + 1;
				my $target_base_num = &get_base_in_target($relative_pos,@cigar);
				if ($target_base_num eq 'no') { next };
				my $target_base = substr($seq,$target_base_num-1,1);
				#print "$contig\t$pos\t$target_base\t$snp1\t$snp2\n";
				
				if ($target_base eq $snp1) {
					$snp1_count ++;
					$total_count ++;
				}elsif ($target_base eq $snp2){
					$snp2_count ++;
					$total_count ++;
				}else {
					$total_count ++;
				}
				
			}elsif ($pos > $end) {
				last;
			}
		}
		
		if ($overlap == 0) {
			$not_phase += 1;
		}
		
		#print "$contig\t$start\t$end\n"
	}
	
	if ($total_count > 0) {
		if ($snp1_count > $snp2_count) {
			$accuracy = sprintf("%.2f",$snp1_count / $total_count * 100);
			return("1\t$accuracy\t$reads_chr\t$reads_start\t$reads_end");
		}elsif ($snp1_count < $snp2_count) {
			$accuracy = sprintf("%.2f",$snp2_count / $total_count * 100);
			return("2\t$accuracy\t$reads_chr\t$reads_start\t$reads_end");
		}else{
			return("0\t0\t$reads_chr\t$reads_start\t$reads_end");
		}
	}else{
		return("0\t0\t$reads_chr\t$reads_start\t$reads_end");
	}
}

sub get_base_in_target {
	my $pos = shift @_;
	my @cigar = @_;
	
	#print "$pos\n";
	
	## position on reads ##
	my $base_num = 0;
	## position on contig ##
	my $contig_num = 0;
	
	foreach (@cigar) {
		if (/S/) {
			my $num = (split /S/,$_)[0];
			$base_num += $num;
		}elsif (/H/) {
			next;
		}elsif (/I/) {
			my $num = (split /I/,$_)[0];
			
			$base_num += $num;
		}elsif (/D/) {
			my $num = (split /D/,$_)[0];
			if ($pos > $contig_num && $pos < $contig_num + $num) {
				## pos in the Deletion region: is an ERROR ##
				## warn "$pos $contig_num $num in the Deletion region: is an ERROR:$!";
				return("no")
			}
			$contig_num += $num;
		}elsif (/N/) {
			my $num = (split /N/,$_)[0];
                        if ($pos > $contig_num && $pos < $contig_num + $num) {
                                ## pos in the Deletion region: is an ERROR ##
                                ## warn "$pos $contig_num $num in the Deletion region: is an ERROR:$!";
                                return("no")
                        }
                        $contig_num += $num;
		}elsif (/M/) {
			my $num = (split /M/,$_)[0];
			## find the pos ##
			if ($pos > $contig_num && $pos < $contig_num + $num) {
				my $add_num = $pos - $contig_num;
				my $reads_pos = $add_num + $base_num;
				#print "$reads_pos\n";
				return($reads_pos);
			}else{
				$base_num += $num;
				$contig_num += $num;
			}
		}else{
			die "Cannot recongnize $_:$!";
		}
	}
}

sub get_s_e_pos {
	my ($start_pos,$cigar) = @_;
	## split cigar ##
	my @cigar = (split /(\d+\w)/,$cigar);
	@cigar = grep (!/\A\s*\z/,@cigar);
	
	my $total_len = 0;
	foreach (@cigar) {
		if (/M/) {
			my $num = (split /M/,$_)[0];
			$total_len += $num;
		}elsif (/D/) {
			my $num = (split /D/,$_)[0];
			$total_len += $num;
		}elsif (/N/) {
			my $num = (split /N/,$_)[0];
                        $total_len += $num;
		}
	}
	
	my $end_pos = $start_pos + $total_len - 1;
	
	unshift @cigar,$end_pos;
	return(@cigar);
}

sub binarySearch {
	my ($pos,$arr) = @_;
	my ($low, $high) = (0, $#$arr);
	
	while ($low <= $high) {
        	my $mid = int(($low + $high) / 2);
        	my $mid_pos = $arr->[$mid];
        	if ($mid_pos == $pos) {
            		return $mid;
        	} elsif ($mid_pos < $pos) {
            		$low = $mid + 1;
        	} else {
            		$high = $mid - 1;
        	}
    	}
    	return $low;
}



