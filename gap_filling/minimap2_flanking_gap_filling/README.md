# Usage
```
Usage:
  00plot_gap_aln.R lg1.1.txt lg1.1.pdf
```
To assist in gap filling, 10 Kb flanking sequences on both sides of each gap were extracted. For example, ‘lg1.1.L’ and ‘lg1.1.R’ represent the 10 Kb sequences upstream and downstream of a gap on chromosome 1, respectively. These flanking sequences were aligned to the contigs of the verkko and Necat assemblies using minimap2 with the parameter -x asm5, and the output was saved in PAF format. The first 12 columns of the PAF file were extracted and saved as a text file (e.g., lg1.1.txt), which served as input for generating alignment plots.