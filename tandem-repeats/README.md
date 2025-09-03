This folder contains files needed to produce a scatter plot of tandem repeats in five genomes.

* `Snakefile` snakemake pipeline for preprocessing tantan results for visualization
* `draw.py` Python script for producing the plot
* `jamAng`, `jamPal`, `jamPhy`, `symKan` folders with source data for individual genomes. Each folder contains two files:
  * `{genome}-telo-self.bed` Position of telomeric repeats in the assembly
  * `tandem-repeats-all.txt` Output of tantan tandem repeat fidnign tool
* `ustMay-tandem-repeats-all.txt` the results of tantan for M. maydis (telomeric repeats not considered for this genome as the motif length is much shorter)
* `ustMay-tandem-repeats-proc.tsv`, `jamAng-tandem-repeats-proc.tsv`, `jamPal-tandem-repeats-proc.tsv`, `jamPhy-tandem-repeats-proc.tsv`, `symKan-tandem-repeats-proc.tsv` Processed tables of tandem repeats after filtering (keeping only repeats of length at least 30 and with at least 2 full repeats) and marking repeats overlapping telomeres. These are produced by the snakemake.


The columns in tantan output (`*tandem-repeats-all.txt`) are 0: chromosome, 1: start, 2: end (BED-style), 3: motif length, 4: number of repetitions, 5: motif, 6: occurrences.

The columns in postproccesed output: 0: repeat id, 1: motif length, 2: number of repetition, 3: tandem array length, 4: is telomeric?

Commands to run:
```bash
snakemake -c 1 all
python3 draw.py
```

Tantan outputs were created by commands of the form:
```bash
tantan -f 4 -w 250 jamAngA5.fa  > jamAng/tandem-repeats-all.txt
```
where jamAngA5.fa is the J. angkorensis assembly.

Software needed:
* [tantan](https://gitlab.com/mcfrith/tantan) (v.22 used)
* [overlapSelect](https://github.com/ucscGenomeBrowser/kent) from Kent utilities
* Snakemake (v.7.32.4 used)
* perl (v5.38.2 used)
* Python 3.12.3 with libraries:
  * pandas==2.1.4+dfsg
  * matplotlib==3.6.3
  * seaborn==0.13.2


