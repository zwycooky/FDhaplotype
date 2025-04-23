#!/usr/bin/perl

use strict;

my ($vcf,$contig_pos) = @ARGV[0];
my $Usage = "\n\t$0 <vcf>
\n";
die $Usage unless (@ARGV == 1);


my (@new_haps,$header);
open VCF,'<',"$vcf" or die;
while (<VCF>) {
	chomp;
	if (/\A##/) {
                next;
        }elsif (/\A#CHROM/) {
                my @line = split;
                my @sample = @line[9..(@line-2)];
                @sample = map {(split /\./,(split /\//,$_)[-1])[0]} @sample;
                $header = join ("\t",@sample);
		open HEADER,'>',"header.txt";
		print HEADER "chr\tpos\tref\talt\t$header\n";
		close HEADER;
        }else{
                my ($chr,$pos,$ref,$alt) = (split)[0,1,3,4];
		if ($ref =~ /[^ACGT]/ || $alt =~ /[^ACGT]/) { next };
		my @line = split;
		my $fd_geno = $line[@line-1];
		$fd_geno = (split /:/,$fd_geno)[0];
		if ($fd_geno ne '0/1') { next };
		
		my $snp_id = "$chr\_$pos";

                my @snp = @line[9..(@line-2)];
                @line = ();

                my @snp2;
                foreach (@snp) {
                        $_ = (split /:/,$_)[0];
                        if ($_ eq '0/1') {
                                push @snp2,'0/1';
                        }elsif ($_ eq '1/1') {
                                push @snp2,'1/1';
                        }elsif ($_ eq '0/0') {
                                push @snp2,'0/0';
                        }else{
                                push @snp2,$_;
                        }
                }
                @snp = ();
                my $snp = join ("\t",@snp2);
		my $chr_num = (split /lg/,$chr)[1];
		print "$chr_num\t$pos\t$snp_id\t$ref\t$alt\t$snp\n";
        }
}
close VCF;





