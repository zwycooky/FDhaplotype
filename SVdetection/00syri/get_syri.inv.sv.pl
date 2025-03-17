#!/usr/bin/perl

use strict;
use Math::CDF qw(:all);
use Parallel::ForkManager;

my ($input_vcf) = @ARGV[0];
my $Usage = "\n\t$0 <syri.vcf>
\n";
die $Usage unless (@ARGV == 1);

## ForkManager setup ##

open IN,'<',"$input_vcf" or die;
while (<IN>) {
	chomp;
	
	if (/\A#/) { next };
	my ($Achr,$Astart,$type,$info) = (split)[0,1,4,7];
	my ($Aend,$Bchr,$Bstart,$Bend);
	
	my $A = (split /_/,$Achr)[0];

	if ($type eq '<INV>') {
		if ($info =~ /END=(\d+?);ChrB=(\w+?);StartB=(\d+?);EndB=(\d+?);/) {
			($Aend,$Bchr,$Bstart,$Bend) = ($1,$2,$3,$4);
		}

		my $B = (split /_/,$Bchr)[0];
		#if ($A ne $B) { next };

		#print "$Achr\t$Astart\t$Aend\t$Bchr\t$Bstart\t$Bend\tINS\n";
		my $sv_len = $Aend - $Astart;

		if ($sv_len < 50) { next };

		print "$Achr\t$Astart\t$Aend\t$Bchr\t$Bstart\t$Bend\t$sv_len\tINV\t-\t-\n";

	}

}
close IN;

