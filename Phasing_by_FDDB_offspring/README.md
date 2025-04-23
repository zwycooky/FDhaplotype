# phasing SNPs by FDDB offspring
 scripts for phasing SNPs of FDDB by offspring

# Usage
```
#1. format vcf of offspring:
  perl 00format_vcf_for_phasing.pl offspring.vcf > FD.7off.formated.snps
#2. phasing SNPs by offspring
  perl 01phasing_by_homo.pl fuding_family_pred_haps.txt FD.7off.formated.snps
#3. merge phased SNP from sperm cell and offspring into one file
  perl merge_SC_homo_phasing.pl FD.withUl5.chr.fill10.hap.txt fd_filter_phasing_res.txt
```
'fuding_family_pred_haps.txt' can be obtained by IBDHap-detector, see https://github.com/zwycooky/IBDHap-detector for details.  
'FD.withUl5.chr.fill10.hap.txt' is phased SNPs identified by sperm cell sequencing using PollenSeq pipeline, see https://github.com/zwycooky/PollenSeq for details.  
The final out file will be FD_SC_7off_merged.hap.txt
