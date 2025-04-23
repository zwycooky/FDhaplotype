#!/usr/bin/perl

use strict;

my ($correct_hap,$homo_phasing) = @ARGV[0,1];
my $Usage = "\n\t$0 <correct_hap> <homo_phasing>
\n";
die $Usage unless (@ARGV == 2);


my (%hapdat,@merged_haps);
open IN,'<',"$correct_hap" or die;
while (<IN>) {
	chomp;
	my ($chr,$pos,$hap1,$hap2) = (split)[0,1,2,3];
	my $snp_id = "$chr\_$pos";
	$hapdat{$snp_id} = "$hap1\t$hap2";
	push @merged_haps,"$chr\t$pos\t$hap1\t$hap2";
}
close IN;

open NOTMATCH,'>',"SC_notmatch_FDoff.hap.txt";
open MATCH,'>',"SC_match_FDoff.hap.txt";
open IN,'<',"$homo_phasing" or die;
while (<IN>) {
	chomp;
	my ($snp_id,$hap1,$hap2) = (split)[2,3,4];
	my ($chr,$pos) = (split /_/,$snp_id)[0,1];
	if (exists $hapdat{$snp_id}) {
		my ($SChap1,$SChap2) = (split /\t/,$hapdat{$snp_id})[0,1];
		if ($SChap1 ne $hap1) {
			print NOTMATCH "$snp_id\t$SChap1\t$SChap2\t$hap1\t$hap2\n";
		}else{
			print MATCH "$snp_id\t$SChap1\t$SChap2\t$hap1\t$hap2\n";
		}
	}else{
		push @merged_haps,"$chr\t$pos\t$hap1\t$hap2";
	}
}
close IN;
close NOTMATCH;
close MATCH;

open OUT,'>',"FD_SC_7off_merged.hap.txt";
#open OUT2,'>',"FD_7off_INDEL.hap.txt";
foreach (sort {(split /lg/,(split /\t/,$a)[0])[1] <=> (split /lg/,(split /\t/,$b)[0])[1] || (split /\t/,$a)[1] <=> (split /\t/,$b)[1]} @merged_haps) {
	my ($hap1,$hap2) = (split)[2,3];
	if (length($hap1) == 1 && length($hap2) == 1) {
		print OUT "$_\n";
	}else{
	#	print OUT2 "$_\n";
	}
}
close OUT;
close OUT2;

