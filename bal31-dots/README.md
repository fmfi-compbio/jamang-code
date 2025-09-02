This folder contains data and scripts needed to make chromosome plots for BAL31 experiments.

Files:
* `Snakefile` script for snakemake
* `fisher.py` a Python script for computing Fisher test p-values
* `chromfigPaired.R` a script in R for creating the final images
* `jamAngA5.sizes` chromosome sizes
* `orig-nanopore-*.tsv.gz` input data files, one for each nanopore sample sequenced in the two experiments

Other tools/languages/libraries needed:
* bedtools (v2.31.1 used)
* perl (v5.38.2 used)
* R (v4.3.3 used)
  * ggplot2 library for R (v3.4.4 used)
* Python (v3.12.3) and libraries:
  * numpy==1.26.4
  * pandas==2.1.4+dfsg
  * SciPy==1.11.4
* Snakemake (v.7.32.4 used)


How to run:
```bash
# unzip data files
gunzip *.gz
# run snakemake on 1 thread
snakemake -c 1
```

Each of the input tsv files contains 4 columns:
* name of contig (uses an older naming scheme with `contig` instead of `chr`; this is renamed later in the Snakemake)
* position along the chromosome
* the number of occurrences of 21-mer at this position in the genome
* the number of occurrences of 21-mer at this position in the nanopore reads