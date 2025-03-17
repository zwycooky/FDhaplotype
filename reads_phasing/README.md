# phasing NGS/Pacbio/ONT reads of FDDB
 A script of reads phasing

# Usage
```
Usage:
  separating_reads_by_haplotype.binarySearch.pl <bam> <phased hap file> <output prefix>
```
Test hap file is named as "lg1.10Mb.hap.txt", test bam file can be downlowded from 10.6084/m9.figshare.28609913
Notice: hap file is obtained by PollenSeq:
Zhang, W., Tariq, A., Jia, X. et al. Plant sperm cell sequencing for genome phasing and determination of meiotic crossover points. Nat Protoc 20, 690–708 (2025). https://doi.org/10.1038/s41596-024-01063-2
## Format of input file
Example of hap file:
```
lg1	29741	C	A
lg1	126320	T	G
lg1	126404	C	T
lg1	126432	A	C
lg1	139327	T	C
lg1	139365	C	T
lg1	511733	T	G
lg1	534000	G	A
lg1	534050	C	T
lg1	534095	T	C
lg1	534101	A	G
lg1	534133	A	T
lg1	534364	T	C
lg1	534802	C	G
```
The 1st col is chromosome id, 2nd is position, 3rd is hap1, 4th is hap2
