#!/usr/bin/perl

use strict;
my $hap_pred = "fuding_family_pred_haps.txt";

my $pred;
open PRED,'<',"$hap_pred" or die;
while (<PRED>) {
        chomp;
        if (/\Astart/) { next };
        my ($start,$end,$hap,$chr,$sample) = (split)[0,1,2,3,5];

        if ($hap != 1 && $hap != 2) { next };
        push @{$pred->{$sample}->{$chr}}, "$start\t$end\t$hap";

}
close PRED;

print "start\tend\thaps\tchr\tlength\taccession\n";
foreach (sort keys %{$pred}) {
	my $sam_id = $_;
	foreach (1..15) {
		my $chr = $_;
		my @tmp = @{$pred->{$sam_id}->{$chr}};
		my @merged_hap = &merge_hap(@tmp);
		foreach (@merged_hap) {
			my ($start,$end,$hap) = split;
			my $len = $end - $start + 1;
			print "$start\t$end\t$hap\t$chr\t$len\t$sam_id\n";
		}
	}
}

sub merge_hap {
	my ($pre_hap,$pre_start,$pre_end,$fir,@merged_hap);
	foreach (@_) {
		my ($start,$end,$hap) = split;
		if ($fir == 0) {
			$pre_start = $start;
			$pre_end = $end;
			$pre_hap = $hap;
			$fir = 1;
			next;
		}
		if ($pre_hap == $hap) {
			$pre_end = $end;
		}elsif ($pre_hap != $hap && $end - $start < 1000000) {
			$pre_end = $end;
		}elsif ($pre_hap != $hap && $end - $start >= 1000000) {
			push @merged_hap,"$pre_start\t$pre_end\t$pre_hap";
			$pre_start = $start;
			$pre_end = $end;
			$pre_hap = $hap;
		}
	}
	
	push @merged_hap,"$pre_start\t$pre_end\t$pre_hap";
	return(@merged_hap);
}
