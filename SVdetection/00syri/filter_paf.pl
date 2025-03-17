#!/usr/bin/perl

use strict;
my ($paf) = @ARGV[0];
my $Usage = "\n\t$0 <paf>
\n";
die $Usage unless (@ARGV == 1);

open IN,'<',"$paf" or die;
while (<IN>) {
	chomp;
	my ($qal_len,$ral_len,$mapq) = (split)[9,10,11];
	if ($mapq > 20) {
		print "$_\n";
	}
}
close IN;
