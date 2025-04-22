# Usage
```
Usage:
  00plot_mummer_gap_filling.R lg2.1.n.txt lg2.1.n.pdf
```
Contigs from the verkko and Necat assemblies were aligned to the FDDB genome using nucmer with default parameters. The resulting alignment files were filtered using delta-filter with the parameters -r -q -l 10000 to retain reliable matches. The filtered delta files were then converted into tabular alignment format using show-coords, generating output files such as lg2.1.n.txt.