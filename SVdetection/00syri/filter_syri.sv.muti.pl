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
my $cpu = 15;
my $pm = Parallel::ForkManager->new($cpu);
my $vcf_dat;
open IN,'<',"$input_vcf" or die;
while (<IN>) {
	chomp;
	if (/\A#/) { next };
	my $line = $_;
	my ($Achr,$type,$info) = (split)[0,2,7];
	
	if ($info =~ /END=(\d+?);ChrB=(\w+?);StartB=(\d+?);EndB=(\d+?);/) {
        	my $Bchr = $2;
		$Bchr = (split /_/,$Bchr)[0];
		my $Atmp_chr = (split /_/,$Achr)[0];
		if ($Atmp_chr ne $Bchr) {
			next;
		}
        }

	if ($type =~ /INS/ || $type =~ /DEL/) {
		push @{$vcf_dat->{$Achr}}, $line;
	}

}
close IN;

foreach (sort keys %{$vcf_dat}) {
	$pm->start and next;
	&filter_vcf(\@{$vcf_dat->{$_}},$fda,$fdb);
	$pm->finish;
}
$pm->wait_all_children();

## -- sub -- ##
sub filter_vcf {

my ($vcf,$fda,$fdb) = @_;
foreach (@{$vcf}) {
	
	my ($Achr,$Astart,$type,$info) = (split)[0,1,2,7];
	my ($Aend,$Bchr,$Bstart,$Bend);
	
	if ($type =~ /\AINS\d+\z/) {
		if ($info =~ /END=(\d+?);ChrB=(\w+?);StartB=(\d+?);EndB=(\d+?);/) {
			($Aend,$Bchr,$Bstart,$Bend) = ($1,$2,$3,$4);
		}
		#print "$Achr\t$Astart\t$Aend\t$Bchr\t$Bstart\t$Bend\tINS\n";
		my $sv_len = $Bend - $Bstart;

		if ($sv_len < 50) { next };
		my (@depth,@flank,@sv);
		if ($sv_len > 1000) {
			my $s = $Bstart - 500;
			my $e = $Bstart + 500;
			chomp (my @depth = `samtools depth -a -l 2000 -r $Bchr:$s-$e $fdb|cut -f 3`);
			@flank = @depth[0..499];
			@sv = @depth[500..999];
			
			my $s = $Bend - 500;
                        my $e = $Bend + 500;
                        chomp (my @depth = `samtools depth -a -l 2000 -r $Bchr:$s-$e $fdb|cut -f 3`);
			
			push @flank, @depth[500..999];
		        push @sv, @depth[0..499];
		}else{
			my $flank_len = int($sv_len / 2);
			my $s = $Bstart - $flank_len;
                        my $e = $Bstart;
                        chomp (my @depth = `samtools depth -a -l 2000 -r $Bchr:$s-$e $fdb|cut -f 3`);
			push @flank, @depth;

			my $s = $Bstart;
                        my $e = $Bend;
                        chomp (my @depth = `samtools depth -a -l 2000 -r $Bchr:$s-$e $fdb|cut -f 3`);
			push @sv, @depth;

			my $s = $Bend;
                        my $e = $Bend + $flank_len;
                        chomp (my @depth = `samtools depth -a -l 2000 -r $Bchr:$s-$e $fdb|cut -f 3`);
                        push @flank, @depth;
		}

		my ($t_value, $p_value, $mean_flank, $mean_sv) = independent_t_test(\@flank, \@sv);
		if ($p_value == -1 && abs($mean_flank - $mean_sv) >= 5) {
			print "$Achr\t$Astart\t$Aend\t$Bchr\t$Bstart\t$Bend\t$sv_len\tINS\n";
		}elsif ($p_value < 0.0001 && abs($mean_flank - $mean_sv) >= 5) {
			print "$Achr\t$Astart\t$Aend\t$Bchr\t$Bstart\t$Bend\t$sv_len\tINS\n";
		}

	}elsif ($type =~ /\ADEL\d+\z/) {
		if ($info =~ /END=(\d+?);ChrB=(\w+?);StartB=(\d+?);EndB=(\d+?);/) {
			($Aend,$Bchr,$Bstart,$Bend) = ($1,$2,$3,$4);
                }
		#print "$Achr\t$Astart\t$Aend\t$Bchr\t$Bstart\t$Bend\tDEL\n";
		
		$Achr = (split /_/,$Achr)[0];
		my $sv_len = $Aend - $Astart;

                if ($sv_len < 50) { next };
                my (@depth,@flank,@sv);
                if ($sv_len > 1000) {
                        my $s = $Astart - 500;
                        my $e = $Astart + 500;
                        chomp (my @depth = `samtools depth -a -l 2000 -r $Achr:$s-$e $fda|cut -f 3`);
                        @flank = @depth[0..499];
                        @sv = @depth[500..999];

                        my $s = $Aend - 500;
                        my $e = $Aend + 500;
                        chomp (my @depth = `samtools depth -a -l 2000 -r $Achr:$s-$e $fda|cut -f 3`);

                        push @flank, @depth[500..999];
                        push @sv, @depth[0..499];
                }else{
                        my $flank_len = int($sv_len / 2);
                        my $s = $Astart - $flank_len;
                        my $e = $Astart;
                        chomp (my @depth = `samtools depth -a -l 2000 -r $Achr:$s-$e $fda|cut -f 3`);
                        push @flank, @depth;

                        my $s = $Astart;
                        my $e = $Aend;
                        chomp (my @depth = `samtools depth -a -l 2000 -r $Achr:$s-$e $fda|cut -f 3`);
                        push @sv, @depth;

                        my $s = $Aend;
                        my $e = $Aend + $flank_len;
                        chomp (my @depth = `samtools depth -a -l 2000 -r $Achr:$s-$e $fda|cut -f 3`);
                        push @flank, @depth;
                }
		
                my ($t_value, $p_value, $mean_flank, $mean_sv) = independent_t_test(\@flank, \@sv);
		if ( $p_value == -1 && abs($mean_flank - $mean_sv) >= 5 ) {
                	print "$Achr\_FDA\t$Astart\t$Aend\t$Bchr\t$Bstart\t$Bend\t$sv_len\tDEL\n";
		} elsif ( $p_value < 0.0001 && abs($mean_flank - $mean_sv) >= 5 ) {
                        print "$Achr\_FDA\t$Astart\t$Aend\t$Bchr\t$Bstart\t$Bend\t$sv_len\tDEL\n";
                }

	}

}

}

## --- ##
sub independent_t_test {
    my ($group1, $group2) = @_;

    # 计算每组的样本大小
    my $n1 = scalar(@$group1);
    my $n2 = scalar(@$group2);
    
    # 计算每组的均值
    my $mean1 = mean($group1);
    my $mean2 = mean($group2);
    
    if (is_constant($group1) || is_constant($group2)) {
        return (-1,-1,$mean1,$mean2);
    }

    # 计算每组的方差
    my $var1 = variance($group1, $mean1);
    my $var2 = variance($group2, $mean2);

    # 计算合并方差
    my $pooled_var = (($n1 - 1) * $var1 + ($n2 - 1) * $var2) / ($n1 + $n2 - 2);

    # 计算 t 值
    my $t_value = ($mean1 - $mean2) / sqrt($pooled_var * (1 / $n1 + 1 / $n2));

    # 计算自由度
    my $df = $n1 + $n2 - 2;

    # 计算 p 值（双侧检验）
    my $p_value = 2 * (1 - pt(abs($t_value), $df));

    return ($t_value, $p_value, $mean1, $mean2);
}

sub is_constant {
    my ($data) = @_;
    my $first_value = $data->[0];
    foreach my $value (@$data) {
        if ($value != $first_value) {
            return 0;
        }
    }
    return 1;
}

# 计算均值的子程序
sub mean {
    my ($data) = @_;
    my $sum = 0;
    foreach my $value (@$data) {
        $sum += $value;
    }
    return $sum / scalar(@$data);
}

# 计算方差的子程序
sub variance {
    my ($data, $mean) = @_;
    my $sum_squared_diff = 0;
    foreach my $value (@$data) {
        $sum_squared_diff += ($value - $mean) ** 2;
    }
    return $sum_squared_diff / (scalar(@$data) - 1);
}

