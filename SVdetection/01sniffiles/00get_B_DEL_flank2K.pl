#!/usr/bin/perl

my ($vcf,$genome) = @ARGV[0,1];
my $Usage = "\n\t$0 <vcf> <genome>
\n";
die $Usage unless (@ARGV == 2);

my (%seq,$id);
open IN,'<',"$genome" or die;
while (<IN>) {
	chomp;
	if (/>/) {
		$id = (split />/,(split)[0])[1];
		$seq{$id} = '';
	}else{
		$seq{$id} .= $_;
	}
}
close IN;

open IN,'<',"$vcf" or die;
while (<IN>) {
	chomp;
	if (/\A#/) { next };
	my ($chr,$start,$svid,$svtype,$pass,$info,$gt) = (split)[0,1,2,4,6,7,9];
	
	next unless ($svtype eq '<DEL>' && $gt =~ /0\/1/ && $info =~ /\APRECISE/);

	my $end;
	if ($info =~ /END=(\d+?;)/) {
		$end = $1;
	}

	my $flank1 = substr($seq{$chr},$start-2001,2000);
	my $flank2 = substr($seq{$chr},$end-1,2000);
	
	my $flankseq = $flank1 . $flank2;

	print ">$chr\_$svid\n$flankseq\n";

}
close IN;

