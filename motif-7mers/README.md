This folder contains files needed for analysis of 7-mers in telomeric motifs

* `all-telo.fa` the list of all motifs, the main input for the analysis
* `clusters.py` a Python script that runs various analyses on the input motifs (containment distance, clustering at various thresholds)
* `motif-groups.py` a Python script for producing the figure
* `fw-c-7-sel-groups.tsv` output of `clusters.py`
* `groups-manual.tsv` reformatted version of `fw-c-7-sel-groups.tsv`
* `groups-manual-matrix12.csv` manually created layout of the figure


Commands used:
```bash
# producing fw-c-7-sel-groups.tsv
python clusters.py --containment all-telo.fa fw-c-7-sel --dist_list "0.1,0.7,0.8,1" --threshold 0.7

# reformatting fw-c-7-sel-groups.tsv to groups-manual.tsv
perl -lane 'next if $.==1; $n++; $n--; $F[0]=~s/...(..).-chr(\d+[LR])-(.*)$/$1$2$3/ or die $F[0]; $F[0]=~s/d$//; print join("\t", @F, int($n/8), $n%8); $n++' fw-c-7-sel-groups.tsv > groups-manual.tsv

# producing the figure
python3 ./motif-groups.py
```

Python 3.12.3 with libraries:
* numpy==1.26.4
* pandas==2.1.4+dfsg
* matplotlib==3.6.3
* seaborn==0.13.2
* SciPy==1.11.4
* scikit-learn==1.4.1.post1
