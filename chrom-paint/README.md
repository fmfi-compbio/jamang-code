This folder contains files needed to create chromosome paint figure where chromosomes of J. angkorensis are painted by nucleotide alignments with the three target genomes.

Files:
* `chrompaint2.py` the script to create the figure
* `jamAngA5.sizes` the sizes of J. angkorensis chromosomes
* `jamPalA1-LASTSPLIT-jamAngA5.psl.view2` alignments of J. pallidilutea to J. angkorensis (reformatted from psl output, see below)
* `jamPhyA2-LASTSPLIT-jamAngA5.psl.view2` alignments of P. phylloscopi to J. angkorensis
* `symKanA1-LASTSPLIT-jamAngA5.psl.view2` alignments of S. kandeliae to J. angkorensis

Python 3.12.3 with libraries:
* colorcet=3.1.0
* numpy==1.26.4
* pandas==2.1.4+dfsg
* matplotlib==3.6.3
* seaborn==0.13.2

Command to create the figure:
```bash
python3 chrompaint2.py
```

Commands to produce alignment file `jamPalA1-LASTSPLIT-jamAngA5.psl.view2` from assemblies `jamPalA1.fa` and `jamAngA5.fa` (using tool [last](https://gitlab.com/mcfrith/last) and its subtools `lastdb`, `lastal`, `maf-convert` as well as Perl one-liners). Analogous commands were used also for the two other pairs of species.

```bash
# format one genome as a database for Last aligner
lastdb jamPalA1-LASTSPLIT-jamAngA5.psl-tmp jamPalA1.fa
# run Last aligner and last-split
lastal jamPalA1-LASTSPLIT-jamAngA5.psl-tmp jamAngA5.fa -E1e-10 | last-split > jamPalA1-LASTSPLIT-jamAngA5.psl.maf
# convert results to psl
maf-convert psl jamPalA1-LASTSPLIT-jamAngA5.psl.maf > jamPalA1-LASTSPLIT-jamAngA5.psl
# reformat psl alignment file to psl.view2 with selected columns
perl -lane 'next unless /^[0-9]/; $g=$F[0]+$F[2]; $b=$F[1]+$F[3]+$F[5]+$F[7]; $s=sprintf "%.1f", 100*$g/($g+$b); print join("\t", @F[0..16], $s)' jamPalA1-LASTSPLIT-jamAngA5.psl > jamPalA1-LASTSPLIT-jamAngA5.psl.view
perl -lane 'print join("\t", @F[0,8..17])' jamPalA1-LASTSPLIT-jamAngA5.psl.view > jamPalA1-LASTSPLIT-jamAngA5.psl.view2
# remove temporary files
rm jamPalA1-LASTSPLIT-jamAngA5.psl.maf jamPalA1-LASTSPLIT-jamAngA5.psl-tmp.* jamPalA1-LASTSPLIT-jamAngA5.psl.view
```