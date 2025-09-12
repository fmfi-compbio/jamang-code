# jamang-code

Code used for *Jaminaea angkorensis* analysis from our article currently available on biorxiv:

Brona Brejova, Viktoria Hodorova, Hana Lichancova, Askar Gafurov, Dominik Bujna, Filip Brazdovic, Filip Cervenak, Tomas Petrik, Eva Hegedusova, Michaela Forgacova Jakubkova, Martina Nebohacova, Lubomir Tomaska, Matthias Sipiczki, Tomas Vinar, Jozef Nosek: **Noncanonical chromosomal-end-specific telomeric arrays in naturally telomerase-negative yeasts**. [bioRxiv 2025.09.07.674783](https://doi.org/10.1101/2025.09.07.674783)

Data can be visualized in our [genome browser](http://genome.compbio.fmph.uniba.sk/).

Each folder contains both data and scripts for a particular analysis as well as a README file describing the file and other requirements. All analyses were run on Ubuntu Linux.


## Underlying data

* [assemblies](./assemblies/): all four genomes sequenced and assembled in our study
* [proteomes](./proteomes/): predicted protein sequences of encoded by these genomes

Sequencing reads and annotated assemblies can be also found at the European Nucleotide Archive (ENA) under projects:
* [PRJEB35162](https://www.ebi.ac.uk/ena/browser/view/PRJEB35162) *Jaminaea angkorensis* (including BAL-31 digestion experiments)
* [PRJEB77303](https://www.ebi.ac.uk/ena/browser/view/PRJEB77303) *Jaminaea pallidilutea*
* [PRJEB77304](https://www.ebi.ac.uk/ena/browser/view/PRJEB77304) *Parajaminaea phylloscopi*
* [PRJEB77302](https://www.ebi.ac.uk/ena/browser/view/PRJEB77302) *Sympodiomycopsis kandeliae*


## Telomeric motif analyses 

* [motif-7mers](./motif-7mers): analysis of 7-mers in telomeric motifs
* [motif-hmm](./motif-hmm): profile HMMs of telomeric motifs used on long reads
* [telo-violin-plots](./telo-violin-plots): violin plots of telomere lengths

## Other analyses

* [bal31-dots](./bal31-dots): chromosome plots for BAL-31 experiments
* [chrom-paint](./chrom-paint): chromosome paint figure where chromosomes of *J. angkorensis* are painted by nucleotide alignments with the three target genomes
* [tandem-repeats](./tandem-repeats): a scatter plot of tandem repeats in five genomes
