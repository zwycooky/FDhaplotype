#!/usr/bin/perl

use strict;

my $pred_haps = "fuding_family_pred_haps.merged.txt";

my (%chr_s,%chr_e);
foreach (1..15) {
	$chr_s{$_} = 500000000;
	$chr_e{$_} = 0;
}

my $predHap;
open IN,'<',"$pred_haps" or die;
while (<IN>) {
	chomp;
	if (/\Astart/) {
		next;
	}
	my ($start,$end,$hap,$chr,$accession) = (split)[0,1,2,3,5];
	push @{$predHap->{$accession}->{$chr}},"$start\t$end\t$hap";
}
close IN;

## get co ##
my (@co,$co);
foreach (sort keys %{$predHap}) {
	my $accession = $_;
	foreach my $chr (1..15) {
		my @tmp = @{$predHap->{$accession}->{$chr}};
		my ($fir,$pre_hap,$pre_end,$s,$e,$hap);
		foreach (@tmp) {
			($s,$e,$hap) = split;
			if ($fir == 0) {
				$fir = 1;
				$pre_end = $e;
				$pre_hap = $hap;
				
				if ($s < $chr_s{$chr}) {
					$chr_s{$chr} = $s;
				}
				
			}elsif ($hap ne $pre_hap) {
				my $co_lft = $pre_end;
				my $co_rht = $s;
				push @co,"$accession\t$chr\t$co_lft\t$co_rht";
				push @{$co->{$chr}},"$co_lft\t$co_rht";
				
				$pre_end = $e;
				$pre_hap = $hap;
			}else{
				$pre_end = $e;
			}
		}
		
		if ($e > $chr_e{$chr}) {
			$chr_e{$chr} = $e;
		}
		
	}
}
open COOUT,'>',"co.txt";
foreach (@co) {
	print COOUT "$_\n";
}
close COOUT;

## get bins ##
my @bins;
foreach my $chr (sort {$a <=> $b} keys %{$co}) {
	my @tmp = @{$co->{$chr}};
	my ($fir,$pre_co_rht);
	foreach (sort {(split /\t/,$a)[0] <=> (split /\t/,$b)[0]} @tmp) {
		my ($co_lft,$co_rht) = split;
		if ($fir == 0) {
			my $bins = $chr_s{$chr};
			my $bine = $co_lft;
			push @bins,"$chr\t$bins\t$bine";
			$fir = 1;
			$pre_co_rht = $co_rht;
		}elsif ($pre_co_rht < $co_lft){
			my $bins = $pre_co_rht;
			my $bine = $co_lft;
			push @bins,"$chr\t$bins\t$bine";
			$pre_co_rht = $co_rht;
		}else{
			$pre_co_rht = $co_rht;
		}
	}
	my $bins = $pre_co_rht;
	my $bine = $chr_e{$chr};
	push @bins,"$chr\t$bins\t$bine";
}

foreach (@bins) {
	print "$_\n";
}
