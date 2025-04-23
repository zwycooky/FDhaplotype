# QTL detection in FDDB half-sib population
 scripts for QTL detectin in half-sib population

# Usage
```
#1. merge and format IBD block detected in offspring
  perl merge_pred_haps.pl fuding_family_pred_haps.txt > fuding_family_pred_haps.merged.txt
#2. construct bin maps in halp-sib population
  perl 00find_bin_from_FDhap.pl > bins.txt
#3. Genotypes of each bin were assigned in each individual based on whether the bin matched haplotype A or B of FDDB.
  Rscript 01assign_bin_geno.R
#4. QTL detection
  Rscript 02QTLdetect.R
```
'fuding_family_pred_haps.txt' can be obtained by IBDHap-detector, see https://github.com/zwycooky/IBDHap-detector for details.  
The final QTL will be stored at 'bin_met_res.txt'
