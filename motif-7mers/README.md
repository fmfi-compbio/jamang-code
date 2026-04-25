This folder contains files needed for analysis of 7-mers in telomeric motifs and blast analysis of these motifs

* `all-telo.fa` the list of all motifs, the main input for the analysis
* `clusters.py` a Python script that runs various analyses on the input motifs (containment distance, clustering at various thresholds)
* `motif-groups.py` a Python script for producing the figure
* `fw-c-7-sel-groups.tsv` output of `clusters.py`
* `groups-manual.tsv` reformatted version of `fw-c-7-sel-groups.tsv`
* `groups-manual-matrix12.csv` manually created layout of the figure


Commands used:
```bash
# producing fw-c-7-sel-groups.tsv
python clusters.py --add_identity --containment all-telo.fa fw-c-7-sel --dist_list "0,0.1,0.7,0.8,1" --threshold 0.7

# reformatting fw-c-7-sel-groups.tsv to groups-manual.tsv
perl -lane 'next if $.==1; $n++; $n--; $F[0]=~s/...(..).-chr(\d+[LR])-(.*)$/$1$2$3/ or die $F[0]; $F[0]=~s/d$//; print join("\t", @F, int($n/8), $n%8); $n++' fw-c-7-sel-groups.tsv > groups-manual.tsv

# producing the figure
python3 ./motif-groups.py


# blast analysis
# first double each telomeric sequence
perl -lne 'if(/^>/) { $_ .= "-2x"; } else { $_ .= $_; } print $_' all-telo.fa > all-telo2x.fa
# prepare blast db
makeblastdb -in all-telo2x.fa -dbtype nucl -out blast.tmp
# blast telo against itself
blastn -db blast.tmp -query all-telo2x.fa -task blastn -word_size 8 -outfmt "6 qaccver qlen qstart qend saccver slen sstart send nident pident bitscore evalue" -evalue 1e-3 -num_alignments 10000 -num_threads 4 > blast.tmp.out
# convert to BED-like tabular format
perl -lane '$F[2]--; ($s,$e)=@F[6,7]; $str=($s<=$e)?"+":"-"; if($str eq "-") { ($s,$e)=($e,$s); } $s--; print join("\t", $F[8], $str, @F[0,1,2,3,4,5], $s, $e, @F[9,10,11]);'  blast.tmp.out | sort -k3,3 -k1gr > blast.tmp.tsv
# remove self-aln
perl -lane 'print if $F[2] ne $F[6]' blast.tmp.tsv | sort -k12gr > all-telo2x.blast.tsv
# no cross-species alignments found
perl -lane 'print if substr($F[2],0,6) ne substr($F[6],0,6)' all-telo2x.blast.tsv
# no negative strand alignment found
perl -lane 'print if $F[1] eq "-"' all-telo2x.blast.tsv
# rm temporary files
rm blast.tmp.n* blast.tmp.out blast.tmp.tsv all-telo2x.fa
```

* Python 3.12.3 with libraries:
  * numpy==1.26.4
  * pandas==2.1.4+dfsg
  * matplotlib==3.6.3
  * seaborn==0.13.2
  * SciPy==1.11.4
  * scikit-learn==1.4.1.post1
* Perl (v5.38.2 used)
* blastn (blast 2.12.0 used)
