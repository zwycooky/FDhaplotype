# Identification of SV between the two haplotype genomes

Before the pipeline start you need to map HiFi reads of FDDB to hapA and hapB genome, respectively, which will be used for filtering SVs identified by syri and SV calling of sniffles.

# detect SV by syri 
```
# Alignment between the two genome
minimap2 -x asm5 -t 20 -c --eqx -o FDA.FDB.paf $FDA $FDB
# filter paf file by mapq (> 20)
perl 00syri/filter_paf.pl FDA.FDB.paf > FDA.FDB.filter.paf
# call SV by syri 
syri -c FDA.FDB.filter.paf -r $FDA -q $FDB -F P -f --nc 15 --cigar
# filter SV identified by syri
perl 00syri/filter_syri.sv.muti.pl syri.vcf > syri.minimap2.filtered.sv.txt
# get duplication
perl 00syri/get_syri.dup.sv.pl syri.vcf > syri.DUP.sv.txt
# get invertion
perl 00syri/get_syri.inv.sv.pl syri.vcf > syri.INV.sv.txt
```
$FDA and $FDB are the two haplotype genome of FDDB  
You need to modify the path of bam file mapped to FDA and FDB genome in line 15 and 16 of 00syri/filter_syri.sv.muti.pl.
# detect SV by sniffles
```
# detect deletions in FDA genome
sniffles --input FD.hifi.hapA.sorted.bam --vcf FD.hapA.sniffile.vcf
# detect deletions in FDB genome
sniffles --input FD.hifi.hapB.sorted.bam --vcf FD.hapB.sniffile.vcf
# paste 2 Kb flanking sequence of deletion identified in FDB into a 4 Kb sequence
perl 01sniffles/00get_B_DEL_flank2K.pl FD.hapB.sniffile.vcf $FDB > hapB_DEL_flank2K.fa
# Align the 4 Kb sequence to FDA genome
minimap2 -x asm5 -t 10 -N 1 -o hapB_DEL_flank2K.paf $FDA hapB_DEL_flank2K.fa
# filter deletions identified by sniffile in FDB genome
perl 01sniffles/01filter_DEL_sv.pl hapB_DEL_flank2K.paf FD.hapB.sniffile.vcf syri.out > hapB.filtered.sniffiles.DEL.txt
# paste 2 Kb flanking sequence of deletion identified in FDA into a 4 Kb sequence
perl 01sniffles/02get_A_DEL_flank2K.pl FD.hapA.sniffile.vcf $FDA > hapA_DEL_flank2K.fa
# Align the 4 Kb sequence to FDB genome
minimap2 -x asm5 -t 10 -N 1 -o hapA_DEL_flank2K.paf $FDB hapA_DEL_flank2K.fa
# filter deletions identified by sniffile in FDB genome
perl 01sniffles/03filter_DEL_sv.pl hapA_DEL_flank2K.paf FD.hapA.sniffile.vcf syri.out > hapA.filtered.sniffiles.DEL.txt
# merge SVs identified by sniffile
cat hapA.filtered.sniffiles.DEL.txt hapB.filtered.sniffiles.DEL.txt > sniffiles.filtered.SV.txt
```
# merge and filter SV identified by syri and sniffles
```
# filter out overlaped SV identified by sniffile
Rscript 02merge_and_filter/00filter_sniffiles_overlapedSV.R
# filter out overlaped SV identified by syri
Rscript 02merge_and_filter/00filter_syri_overlapedSV.R
# merge sniffile and syri results
Rscript 02merge_and_filter/01merge_sniffiles_syri.R
# final filter overlaped SVs
Rscript 02merge_and_filter/02final_filtered_overlap.R
# final merge into 4 types of SVs (DEL, INS, DUP, INV)
Rscript 02merge_and_filter/03stat_SV.R
```
The final SV file will named as "FDDB_AB_all4type_SV.txt"
