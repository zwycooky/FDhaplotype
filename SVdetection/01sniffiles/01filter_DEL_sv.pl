#!/usr/bin/perl

use strict;

my ($paf,$vcf,$synfile) = @ARGV[0,1,2];
my $Usage = "\n\t$0 <paf> <vcf> <syri.out>
\n";
die $Usage unless (@ARGV == 3);

# read syn region #
my ($syn_pos,$syn);
open SYN,'<',"$synfile" or die;
while (<SYN>) {
	chomp;
	my ($Achr,$As,$Ae,$Bchr,$Bs,$Be,$type) = (split)[0,1,2,5,6,7,10];
	next unless ($type eq 'SYN');
	$Achr = $Achr . '_FDA';
	$Bchr = $Bchr . '_FDB';
	push @{$syn->{$Achr}}, "$Achr\t$As\t$Ae\t$Bchr\t$Bs\t$Be";
	push @{$syn_pos->{$Achr}}, $As;
}
close SYN;

# read sv pos #
my %svinfo;
open IN,'<',"$vcf" or die;
while (<IN>) {
	chomp;
	my ($chr,$start,$sv_id,$info) = (split)[0,1,2,7];
	if ($sv_id =~ /DEL/) {
		my $end;
		if ($info =~ /END=(\d+?);/) {
			$end = $1;
		}
		my $new_id = "$chr\_$sv_id";
		$svinfo{$new_id} = "$chr\t$start\t$end";
	}
}
close IN;

## read paf ##
my $paf_dat;
open IN,'<',"$paf" or die;
while (<IN>) {
	chomp;
	my $line = $_;
	my ($schr,$al_len) = (split)[5,10];
	if ($al_len > 3000) {
		push @{$paf_dat->{$schr}}, $line;
	}else{
		#print "$al_len\n";
	}
}
close IN;

foreach my $chr (sort keys %{$paf_dat}) {
	my @paf_chr = @{$paf_dat->{$chr}};
	my @syn_pos_chr = @{$syn_pos->{$chr}};
	my @syn_chr = @{$syn->{$chr}};

	foreach (sort {(split /\t/,$a)[7] <=> (split /\t/,$b)[7]} @paf_chr) {
		my ($svid,$sqs,$sqe,$strand,$Schr,$Ss,$Se) = (split)[0,2,3,4,5,7,8];
		my ($Qchr,$Qs,$Qe) = (split /\t/,$svinfo{$svid})[0,1,2];

		my $search_s = binarySearch($Ss, \@syn_pos_chr);
		my $search_e = binarySearch($Se, \@syn_pos_chr);
		
		my @tmp_syn = @syn_chr[$search_s-1,$search_e+1];
		my $find = 0;
		foreach (@tmp_syn) {
			my ($Achr,$As,$Ae,$Bchr,$Bs,$Be) = split;
			next unless ($Achr eq $Schr && $Bchr eq $Qchr);
			my $Aoverlap = overlap($As,$Ae,$Ss,$Se);
			my $Boverlap = overlap($Bs,$Be,$Qs,$Qe);
			if ( ($Aoverlap > 0.9 && $Boverlap > 0.9) || abs($Ss - $Qs) < 5000000 ) {
				$find = 1;
				last;
			}
		}
		if ($find == 1) {
			my $Ains_pos;
			if ($strand eq '+') {
				$Ains_pos = $Ss + 2000 - $sqs;
			}else{
				$Ains_pos = $Se - (2000 - $sqs);
			}
			my $svlen = abs($Qe - $Qs);
			if ($svlen < 50000) {
				print "$Schr\t$Ains_pos\t$Ains_pos\t$Qchr\t$Qs\t$Qe\t$svlen\tINS\n";
			}
		}else{
			#print "$Qchr\t$Qs\t$Qe\aINS\n";
		}
	}
}




## sub from here ##
sub binarySearch {
    my ($pos, $dat) = @_;
    my $high = @$dat;
    my $low = 0;

    # 检查 $pos 是否大于数组的最后一个元素
    if ($pos > $dat->[$high - 1]) {
        return $high;
    }

    while ($low < $high) {
        my $mid = int(($low + $high) / 2);
        if ($pos == $dat->[$mid]) {
            return $mid;
        } elsif ($pos > $dat->[$mid]) {
            $low = $mid + 1;
        } else {
            $high = $mid;
        }
    }

    return $low;
}

sub overlap {
        my ($s1,$e1,$s2,$e2) = @_;
        my @sort = sort {$a <=> $b} @_;
        my $sub1 = $e1 - $s1;
        my $sub2 = $e2 - $s2;
        my $sub_all = $sort[-1] - $sort[0];
        my $overlap_len = $sub1 + $sub2 - $sub_all;
	my $overlap_ratio = $overlap_len / $sub2;
        return($overlap_ratio);
}

