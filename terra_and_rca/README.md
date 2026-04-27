Folder containing code and data needed for RCA and TERRA analysis.

Files and folders:
* main snakemake file `telo-ends.smk`
* folders for experiments 
  * `rca1`: treated sample, experiment 1
  * `rca2`: control sample, experiment 1
  * `rca3`: treated sample, experiment 2
  * `rca4`: control sample, experiment 2
  * `rna-seq`: The first RNA-seq illumina library
  * `rna-seq2`: The second RNA-seq illumina library, stranded
  * `rna-seqN2`: nanopore direct RNA library
* Each folder for experiment contains:
  * links to several files from the main folder for configuration, namely `ref.fa` (reference genome assembly), `telo.list` list of telomeres, `telo2x.fa` fasta file with each telomere motif repeated 2x, `hmm-motifs.tsv` tsv file with motifs for building HMM
  * Chromosomal end notation: proximal motifs are here rpresented as seprate "chromosomal ends" with 30 added to chromosome number. Letter `s` after chromosome number means left end of a chromosome (start) and letter `e` means right end
  * Normally, file `reads.fa` with all reads in fasta format needs to be placed to each folder. These are omitted here for space regions but these files cna be easily created from fastq files provided in ENA. For Illumina, the two reads in a read pair are distinguised by suffixes `_1` and `_2`. However, some downatream analyses can be executed from other provided files.
  * `reads-telo.blast` file obtained by snakemake by blasting reads to motifs
  * `reads-filtered.fa` reads containing a match of one of the motifs to be used form HMM analysis
  * Analysis results, namely `hmm-longest-summary.tsv` for RCA and `hmm-counts-extended-summary.tsv` for rna-seq.


```
# running RNA_seq analyses:
snakemake -s telo-ends.smk -p -c 1 --directory rna-seq all_motifs
snakemake -s telo-ends.smk -p -c 1 --directory rna-seq2 all_motifs
snakemake -s telo-ends.smk -p -c 1 --directory rna-seqN2 all_motifs

# running RCA analyses - getting read matches
snakemake -s telo-ends.smk -p -c 1 --directory rca1 rca_base
snakemake -s telo-ends.smk -p -c 1 --directory rca2 rca_base
snakemake -s telo-ends.smk -p -c 1 --directory rca3 rca_base
snakemake -s telo-ends.smk -p -c 1 --directory rca4 rca_base

# normalized RCA read counts and genome coverage
# requires placing reads.fa file into each folder
# created from fastq files submitted to ENA
snakemake -s telo-ends.smk -p -c 1 --directory rca1 rca
snakemake -s telo-ends.smk -p -c 1 --directory rca2 rca
snakemake -s telo-ends.smk -p -c 1 --directory rca3 rca
snakemake -s telo-ends.smk -p -c 1 --directory rca4 rca
```

Software needed:
* faSize, faSomeRecords from [Kent utilities](https://github.com/ucscGenomeBrowser/kent)
* blastn (blast 2.12.0 used)
* Perl (v5.38.2 used)
* bedtools (v2.31.1 used)
* rep_hmm.pl from this repository
* minimap2 (v. 2.24-r1122 used)
* samtools (v. 1.13 used)