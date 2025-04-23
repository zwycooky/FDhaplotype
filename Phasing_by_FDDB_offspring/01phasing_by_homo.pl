#!/usr/bin/perl

use strict;

my ($pred_hap,$formated_snps) = @ARGV[0,1];
my $Usage = "\n\t$0 <pred_haps> <formated snps>
\n";
 die $Usage unless (@ARGV == 2);

 my ($hapdat);
open IN,'<',"$pred_hap" or die;
while (<IN>) {
	chomp;
	my ($start,$end,$hap,$chr,$len,$accession) = split;
	$accession =~ s/X//;
	$accession =~ s/\./-/;
	push @{$hapdat->{$accession}->{$chr}}, "$start\t$end\t$len\t$hap";
}
close IN;

## open header ##
my @sample;
open HEA,'<',"header.txt" or die;
while (<HEA>) {
	chomp;
	@sample = split;
	@sample = @sample[4..(@sample-1)];	
	last;
}
close HEA;

## start phasing ##
open OUT,'>',"fd_all_phasing_res.txt";
open FIL,'>',"fd_filter_phasing_res.txt";
open SNP,'<',"$formated_snps" or die;
while (<SNP>) {
	chomp;
	my ($chr,$pos,$snp_id,$ref,$alt) = (split)[0,1,2,3,4];
	my @snp = split;
	@snp = @snp[5..(@snp-1)];
	
	my %tmp_hap;
	foreach (0..(@snp-1)) {
		my $accession = $sample[$_];
		my $geno = $snp[$_];
		if ($geno eq '0/1' || $geno eq './.') { next };
		
		## pred hap for target accession in tartget chr ##
		my @tmp_pred_hap = @{$hapdat->{$accession}->{$chr}};
		foreach (@tmp_pred_hap) {
			my ($start,$end,$hap) = (split)[0,1,3];
			if ($pos >= $start && $pos <= $end) {
				if ($hap == 1) {
					my ($hapA,$hapB);
					if ($geno eq '0/0') {
                        			$hapA = $ref;
						$hapB = $alt;
                			}elsif ($geno eq '1/1'){
                        			$hapA = $alt;
						$hapB = $ref;
                			}else{
                        			die "ERROR with $geno:$!";
                			}
					my $tmp_hap = "$hapA\t$hapB";
					$tmp_hap{$tmp_hap} ++;
				}elsif ($hap == 2) {
					my ($hapA,$hapB);
					if ($geno eq '0/0') {
                                                $hapA = $alt;
                                                $hapB = $ref;
                                        }elsif ($geno eq '1/1') {
                                                $hapA = $ref;
                                                $hapB = $alt;
                                        }else{
                                                die "ERROR with $geno:$!";
                                        }
					my $tmp_hap = "$hapA\t$hapB";
					$tmp_hap{$tmp_hap} ++;
				}
				last;
			}
		}
	}

	my ($final_hap,$hap_count,$total_hap_count,$ratio);
        foreach (sort keys %tmp_hap) {
        	if ($tmp_hap{$_} > $hap_count) {
                	$final_hap = $_;
                        $hap_count = $tmp_hap{$_};
                        $total_hap_count += $hap_count;
                }else{
                        $total_hap_count += $tmp_hap{$_};
                }
        }
	if ($total_hap_count > 0) {
		$ratio = sprintf("%.2f", $hap_count / $total_hap_count);
	}else{
		next;
	}
	
	print OUT "$chr\t$pos\t$snp_id\t$final_hap\t$hap_count\t$total_hap_count\t$ratio\n";
	if (length($ref) == 1 && length($alt) == 1) {
		if ( ($hap_count == $total_hap_count && $hap_count > 1) || ($total_hap_count - $hap_count == 1 && $hap_count > 3) ) {
			print FIL "$chr\t$pos\t$snp_id\t$final_hap\t$hap_count\t$total_hap_count\t$ratio\n";
		}
	}else{ ## for indel ##
		if ( ($hap_count == $total_hap_count && $hap_count > 2) ) {
			#print FIL "$chr\t$pos\t$snp_id\t$final_hap\t$hap_count\t$total_hap_count\t$ratio\n";
                }
	}
}
close SNP;

