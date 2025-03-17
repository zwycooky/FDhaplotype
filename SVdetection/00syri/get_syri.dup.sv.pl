#!/usr/bin/perl

use strict;
use Math::CDF qw(:all);
use Parallel::ForkManager;

my ($input_vcf) = @ARGV[0];
my $Usage = "\n\t$0 <syri.vcf>
\n";
die $Usage unless (@ARGV == 1);

#samtools depth -a -q 50 -l 2000 -r lg10_FDB:3257-3529 pacbio_bam/FD.hifi.hapB.sorted.bam
#lg10_FDA        1       SYNAL1  N       <SYNAL> .       PASS    END=1414;ChrB=lg10_FDB;StartB=1;EndB=1414;Parent=SYN1;VarType=.;DupType=.       GT      1

my $fda = "/home/zwy/workspace/Fudingdabai/14SVdetection/2025SV/pacbio_bam/FD.hifi.hapA.sorted.bam";
my $fdb = "/home/zwy/workspace/Fudingdabai/14SVdetection/2025SV/pacbio_bam/FD.hifi.hapB.sorted.bam";

## ForkManager setup ##

open IN,'<',"$input_vcf" or die;
while (<IN>) {
	chomp;
	
	if (/\A#/) { next };
	my ($Achr,$Astart,$type,$info) = (split)[0,1,2,7];
	my ($Aend,$Bchr,$Bstart,$Bend);
	
	my $A = (split /_/,$Achr)[0];

	if ($info =~ /DupType=copygain/) {
		if ($info =~ /END=(\d+?);ChrB=(\w+?);StartB=(\d+?);EndB=(\d+?);/) {
			($Aend,$Bchr,$Bstart,$Bend) = ($1,$2,$3,$4);
		}

		my $B = (split /_/,$Bchr)[0];
		#if ($A ne $B) { next };

		#print "$Achr\t$Astart\t$Aend\t$Bchr\t$Bstart\t$Bend\tINS\n";
		my $sv_len = $Aend - $Astart;

		if ($sv_len < 50) { next };

		print "$Achr\t$Astart\t$Aend\t$Bchr\t$Bstart\t$Bend\t$sv_len\tDUP\t-\t-\n";

	}elsif ($info =~ /DupType=copyloss/){
		if ($info =~ /END=(\d+?);ChrB=(\w+?);StartB=(\d+?);EndB=(\d+?);/) {
                        ($Aend,$Bchr,$Bstart,$Bend) = ($1,$2,$3,$4);
                }
                #print "$Achr\t$Astart\t$Aend\t$Bchr\t$Bstart\t$Bend\tINS\n";
                my $sv_len = $Bend - $Bstart;
                
		my $B = (split /_/,$Bchr)[0];
		#if ($A ne $B) { next };

                if ($sv_len < 50) { next };
                        
                print "$Achr\t$Astart\t$Aend\t$Bchr\t$Bstart\t$Bend\t$sv_len\tDUP\t-\t-\n";
	}

}
close IN;

