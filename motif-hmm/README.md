This folder contains files for running profile HMM of telomeric motifs on selected long reads spanning the unique portion near chromosome end. The results are used in the [violin plots](../telo-violin-plots/README.md) of telomeric length.

For each of the four studied species, there is a folder containing the following files:
* `hmm-motifs.tsv` input file listing all telomeric motifs. Distal motifs are marked as `telo`, proximal as `sub` in the first column. The second column contains id of the chromosomal end, where suffix `s` means left end. suffix `e` right end.
* `hmm.count` The counts of bases emitted by different types of states for individual reads (see desciption in [telo-violin-plots](../telo-violin-plots/README.md) folder. This file is created from `hmm.count` files in individual subdirectories.
* `hmm.count.summary` Quartiles of counts from `hmm.count` for each chromosomal end (plus the number of reads)
* For each chrosomosomal end we also have a folder contingin the following files:
  * `chr_end/reads.fa`  Read parts used for analysis of this chromosomal end.
  * `chr_end/hmm.states.gz` The results of running Viterbi algorithm - state emitting each base of the read.
  * `chr-end/hmm.stats` Statistics of state usage per read. For each each telomeric match state we count how many times it was used. For the background state we distinguish "before", "middle" and "end".
  * `chr_end/hmm.count` Even more summarize state usage per read, where all telomeric states from one motif are grouped together.

Scripts:
* `hmm/check_chrom_ends.mk` Makefile for running analysis on all chromosomal ends. It needs file `hmm-motifs.tsv` for the species and for each chromsomal end it needs `chr_end/reads.fa`. It will produce all other files listed here.
* `hmm/rep_hmm.pl` Perl program running the Viterbi algorithm.
* `hmm/rep_hmm.prob` Parameters of the profile HMM used by `hmm/rep_hmm.pl`.
* `hmm/hmm_summary.py` Scripts computing quartiles of lengths from `hmm.count`.

Running scripts:
```bash
export DIR=jamAng; make -f hmm/check_chrom_ends.mk jamAng/hmm.count.summary
# similarly for the other 3 species
```

Software needed:
* GNU make (v4.3 was used)
* Python 3.12.3 with libraries:
  * pandas==2.1.4+dfsg
* Perl (v5.38.2 used)
