HMMBIN = workflow.basedir + "/../motif-hmm/hmm/"
MINIMAP_OPT = "-t 4"


MIN_MATCH_LENGTH = 30
MIN_MATCH_IDENTITY = 95

# read the list of motifs
text_file = open("telo.list", "r")
motifs = text_file.read().splitlines()
motif_results_counts = expand("{motif}/hmm.count", motif = motifs)
motif_results_longest = expand("{motif}/hmm.longest", motif = motifs)
motif_fa = expand("{motif}/reads.fa", motif = motifs)
#print(motif_results_counts)

rule all_motifs:
     input: "hmm-counts-extended-summary.tsv"

rule rca_base:
     input: "hmm-longest.tsv"

rule rca:
     input: "hmm-longest-summary.tsv", "reads-zerocov.txt"

# chromosome sizes
rule sizes:
    input:
        "{assembly}.fa"
    output:
        "{assembly}.sizes"
    shell:
        """
        faSize -detailed {input} > {output}
        """

# chromosome sizes
rule fa_count:
    input:
        "{reads}.fa"
    output:
        "{reads}.count"
    shell:
        """
        grep -c '>' {input} > {output}
        """


rule blast_reads:
  input: reads="reads.fa", telo="telo2x.fa"
  output: "reads-telo.blast"
  shell:
    """
    makeblastdb -in {input.telo} -dbtype nucl -out {output}.tmp_db
    blastn -task blastn-short -db {output}.tmp_db  -query {input.reads} -dust no  -soft_masking false -word_size 10 -evalue 1e-5 -outfmt "6 qaccver qlen qstart qend saccver slen sstart send nident pident bitscore evalue"  > {output}.tmp.blastOrig
    perl -lane '$F[2]--; ($s,$e)=@F[6,7]; $str=($s<=$e)?"+":"-"; if($str eq "-") {{ ($s,$e)=($e,$s); }} $s--; print join("\\t", $F[8], $str, @F[0,1,2,3,4,5], $s, $e, @F[9,10,11]);'  {output}.tmp.blastOrig | sort -k3,3 -k1gr > {output}
    rm {output}.tmp_db.n* {output}.tmp.blastOrig
    """

rule filter_blast:
     input: "reads-telo.blast"
     output: "reads-telo-filtered.blast"
     shell:
        """
        # keep only matches with at least min. length and identity
        perl -lane 'next if $F[5]-$F[4]<{MIN_MATCH_LENGTH} || $F[10]<{MIN_MATCH_IDENTITY}; print' {input} > {output}.tmp
	# remove _[12] suffix, if any, print read, chr end, score
	# sort by read and score
	perl -lane '$F[2]=~s/_[12]$//; print join(" ", $F[2], $F[6], $F[11])' {output}.tmp | sort -k1,1 -k3gr > {output}.tmp2
	# get for each read only the best chr end
	perl -lane 'next if $.>1 && $F[0] eq $o; print join(" ", $F[0], $F[1]); $o=$F[0];' {output}.tmp2 > {output}.tmp3
	# from blast output only get lines that agree with read-chr end
	# the lines do not need to satisfy other criteria
	perl -lane 'BEGIN {{ open $in, "<", "{output}.tmp3" or die; while($k=<$in>) {{ chomp $k; $a{{$k}}=1; }} }} $F[2]=~s/_[12]$//; $k="$F[2] $F[6]"; next unless exists $a{{$k}}; print' < {input} > {output}
	rm {output}.tmp {output}.tmp2 {output}.tmp3
        """

rule prefilter_reads:
     input: blast="reads-telo-filtered.blast", fa="reads.fa"
     output: "reads-filtered.fa"
     shell:
        """
	perl -lane 'print $F[2]' {input.blast} | sort -u > {output}.tmp
	faSomeRecords {input.fa} {output}.tmp {output}
	rm {output}.tmp
        """


rule merge_intervals_motif:
   input: blast="reads-telo-filtered.blast", sizes="reads-filtered.sizes"
   output: "{motif}/reads-telo.bed"
   shell:
     """
     mkdir -p {wildcards.motif}
     # get only our motif, reformat to bed file
     perl -lane 'next unless $F[6] eq "{wildcards.motif}2x"; print join("\t", @F[2,4,5,6,11,1])' {input.blast} | sort -k1,1 -k2g > {output}.tmp.bed
     # merge at distance 500bp
     bedtools merge -s -d 500 -i {output}.tmp.bed -c 4,5,6 -o distinct,max,distinct > {output}.tmp2.bed
     # add 500bp to both ends
     bedtools slop -b 500 -i {output}.tmp2.bed -g {input.sizes} > {output}
     rm {output}.tmp.bed {output}.tmp2.bed 
     """

rule fai:
   input: fa="{name}.fa"
   output: fai="{name}.fa.fai"
   shell:
     """
     samtools faidx {input}
     """

rule get_reads_motif:
   input: fa="reads-filtered.fa", fai="reads-filtered.fa.fai", bed="{motif}/reads-telo.bed"
   output: "{motif}/reads.fa"
   shell:
     """
     bedtools getfasta -fi {input.fa} -bed {input.bed} -s | perl -lne 's/\\(\\+\\)$/:FW/; s/\\(\\-\\)$/:RC/; print'> {output}
     """

rule run_hmm:
    input: fa="{motif}/reads.fa", tsv="hmm-motifs.tsv"
    output: "{motif}/hmm.stats", "{motif}/hmm.states", "{motif}/hmm.bed", txt="{motif}/hmm-motifs.txt"
    shell:
      """
      grep {wildcards.motif} {input.tsv} | perl -lane 'print $F[0], " ", $F[2]' > {output.txt}
      {HMMBIN}/rep_hmm.pl {wildcards.motif}/hmm-motifs.txt {HMMBIN}/rep_hmm.prob {input.fa}  {wildcards.motif}/hmm
      """

rule hmm_counts:
    input: "{motif}/hmm.stats"
    output: "{motif}/hmm.count"
    shell:
      """
      perl -lane 'if($F[1] eq "bg") {{ $r=$F[2]-$F[3]-$F[4]; $r=0 if $r<0; print join("\\t", $F[0], "bg-middle", $r); print join("\\t", $F[0], "bg-start", $F[3]); print join("\t", $F[0], "bg-end", $F[4]); }} else {{ my $s=0; foreach $x (@F[2..$#F]) {{ $s+=$x; }} print join("\\t", @F[0,1], $s) }}' {input} > {output}
      """

rule hmm_longest:
     input: bed="{motif}/hmm.bed", motif="{motif}/hmm-motifs.txt"
     output: "{motif}/hmm.longest"
     shell:
       """
       # compute motif length
       perl -lane 'print length($F[1]); END {{ die $. unless $.==1; }}' {input.motif} > {output}.tmp2
       # add record length and motif length, separate read name and other info,
       # write only telo, sort by read and length
       perl -lane '$L='`cat {output}.tmp2`'; $F[0]=~s/:(\\d+-\\d+:(FW|RC))$/\\t$1/ or die; print join("\\t", @F, $F[2]-$F[1], $L) if $F[3] eq "telo"' {input.bed} | sort -k1,1 -k6gr > {output}.tmp       
       # get the longest 
       perl -lane 'next if $.>1 && $F[0] eq $o; $o=$F[0]; print' {output}.tmp > {output}
       rm {output}.tmp {output}.tmp2
       """

rule summartize_counts:
     input: motif_results_counts
     output: "hmm-counts.tsv"
     shell:
       """
       # get lines with state stats from all chr hmm results
       grep --with-filename . {input} | perl -lne 's/\\/hmm.count:/\\t/; print' > {output}.tmp
       # put the stats on a single line with read id
       # 0:chr_end 1: read_id 2:line_id 3:telo_matches 4:bg-start 5:bg-middle 6:bg-end
       perl -lane 'die unless @F==4; $k="$F[0]\\t$F[1]"; $a{{$k}}{{$F[2]}}=$F[3]; END {{ $n=0; foreach $k (sort keys %a) {{ $n++; printf "%s\\ts%d", $k, $n; foreach $k2 (qw/telo bg-start bg-middle bg-end/) {{ $v=$a{{$k}}{{$k2}}+0; printf "\\t%d", $v; }} print ""; }} }} ' {output}.tmp > {output}.tmp2
       # remove reads with 0 matches
       perl -lane 'print if $F[3]>0' {output}.tmp2 > {output}
       rm {output}.tmp {output}.tmp2
       """

rule summarize_longest:
     input: motif_results_longest
     output: "hmm-longest.tsv"
     shell:
       """
       # get lines with longest from all chr hmm results
       grep --with-filename . {input} | perl -lne 's/\\/hmm.longest:/\\t/ or die; print' > {output}
       """

rule summarize_longest_summary:
     input: tsv="hmm-longest.tsv", r_count="reads.count", t_list="telo.list"
     output: "hmm-longest-summary.tsv"
     shell:
       """
       # get only lines with telo length >= 1.5 motif length
       # count for each chrom
       perl -lane 'print $F[0] if $F[6]>=1.5*$F[7]' {input.tsv} | sort | uniq -c > {output}.tmp
       perl -lane 'BEGIN {{ $N='`cat {input.r_count}`'; }} print join("\\t", $F[1], $F[0], $F[0]*1e6/$N)' {output}.tmp > {output}.tmp2
       # sort telo list
       sort {input.t_list} > {output}.tmp3
       # join with telo list
       join -a 1 {output}.tmp3 {output}.tmp2 | perl -lane 'if(@F==1) {{ push @F, 0, 0; }} print join("\\t", @F)' > {output}
       rm {output}.tmp {output}.tmp2
       """



rule summarize_fasta:
     input: motif_fa
     output: "motifs.fa"
     shell:
       """
       grep . {input} | perl -lne 's/\\/reads.fa:/\\t/; print' > {output}.tmp
       perl -lane 'if($F[1]=~/^>/) {{ $F[1]=">" . $F[0] . ":" . substr($F[1],1); }} print $F[1];' {output}.tmp > {output}
       rm {output}.tmp
       """

rule get_before_fasta:
     input: tsv="hmm-counts.tsv", fa="motifs.fa", fai="motifs.fa.fai"
     output: "motifs_before.fa"
     shell:
       """
       perl -lane '$s="$F[0]:$F[1]"; next if $F[4]+5<30; print join("\\t", $s, 0, $F[4]+5,$F[2],0,"-")' {input.tsv} > {output}.tmp.bed
       bedtools getfasta -nameOnly -s -fi {input.fa} -bed {output}.tmp.bed | perl -lne 's/\\(-\\)$//; print' > {output}
       rm {output}.tmp.bed
       """

rule get_after_fasta:
     input: tsv="hmm-counts.tsv", fa="motifs.fa", fai="motifs.fa.fai"
     output: "motifs_after.fa"
     shell:
       """
       perl -lane '$s="$F[0]:$F[1]"; $l=$F[3]+$F[4]+$F[5]+$F[6]; next if $F[6]+5<30; print join("\\t", $s, $l-$F[6]-5, $l,$F[2],0,"+")' {input.tsv} > {output}.tmp.bed
       bedtools getfasta -nameOnly -s -fi {input.fa} -bed {output}.tmp.bed | perl -lne 's/\\(\\+\\)$//; print' > {output}
       rm {output}.tmp.bed
       """



rule aln_genome:
     input: db="ref.fa", reads="motifs_{name}.fa"
     output: "motifs_{name}-MAP-ref.paf"
     shell:
        """
        minimap2 -t 4 -c -x map-ont -k 10 -w 5 -f 0 --secondary=no {input.db} {input.reads} > {output}
        """


rule paf_view:
    input:
        "{aln}.paf"
    output:
        "{aln}.paf.view"
    shell:
        """
        perl -lane '$id=sprintf("%.1f", $F[9]*100/$F[10]); print join("\\t", @F[9,4,0..3,5..8],$id)' {input} > {output}
        """

# select the first match in each sequence listed as the first
# (second arg. of minimap, here probably reads)
rule aln_view_first:
    input:
        "{aln}.view"
    output:
        "{aln}.first.view"
    shell:
        """
        sort -k3,3 -k5g {input} | perl -lane 'print if $.==1 || $F[2] ne $o; $o=$F[2]' > {output}
        """




rule filter_ref:
     input: "motifs_{name}-MAP-ref.paf.first.view"
     output: "motifs_{name}-MAP-ref.paf.sel.tsv",
     shell:
       """
       # all genome matches, with start in the genome
       perl -lane 'if($F[1] eq "+") {{ $p=$F[8]; }} else {{ $p=$F[9]; }} print join("\t", $F[2], "genome", $F[6] . ":" . $p . $F[1])' {input} | sort -k1,1 > {output}
       """


# {side} is before or after
rule combine_after:
     input:
       ref="motifs_{side}-MAP-ref.paf.sel.tsv",
     output:
       "motifs_{side}.tsv"
     wildcard_constraints:
        side=r"[a-z]+"
     shell:
       """
       cp {input.ref} {output}
       """

rule combine_all:
     input:
        orig="hmm-counts.tsv",
        before="motifs_before.tsv",
        after="motifs_after.tsv"
     output:
        "hmm-counts-extended.tsv"
     shell:
       """
       sort -k3,3 {input.orig} > {output}.tmp1
       sort -k1,1 {input.before} | perl -lane 'print join("\\t", @F[0..2])' > {output}.tmpB
       sort -k1,1 {input.after} | perl -lane 'print join("\\t", @F[0..2])' > {output}.tmpA
       join -a 1 -1 3 {output}.tmp1 {output}.tmpB | perl -lane 'if(scalar(@F)!=9) {{ die unless scalar(@F)==7; push @F, ("unknown", "-"); }} print join("\t", @F)' > {output}.tmp2
       join -a 1 {output}.tmp2 {output}.tmpA | perl -lane 'if(scalar(@F)!=11) {{ die unless scalar(@F)==9; push @F, ("unknown", "-"); }} print join("\t", @F)' > {output}
       rm {output}.tmp1 {output}.tmp2 {output}.tmpA {output}.tmpB       
       """

rule extended_summary:
     input: tsv="hmm-counts-extended.tsv", t_list="telo.list"
     output: "hmm-counts-extended-summary.tsv"
     shell:
       """
       # count read pairs per motif
       perl -lane '$F[2]=~s/_[12]:\\d+-\\d+:(FW|RC)$//; print "$F[1] $F[2]"' {input.tsv} | sort -u | perl -lane 'print $F[0]' | sort | uniq -c > {output}.tmp
       # join with motif table
       sort {input.t_list} > {output}.tmp2
       join -a 1 -2 2 {output}.tmp2 {output}.tmp | perl -lane 'if(@F==1) {{ push @F, 0; }} print join("\\t", @F)' > {output}
       rm {output}.tmp {output}.tmp2
       """


# nanopore aligned by minimap 
rule Nanopore_minimap:
    input:
         ref="ref.fa", reads="reads.fa"
    output:
         bam="reads.bam", bai="reads.bam.bai"
    shell:
        """
        minimap2 -a -x map-ont --secondary=no {MINIMAP_OPT} {input.ref} {input.reads} | samtools view -S -b - | samtools sort - -o {output.bam}
        samtools index {output.bam}
        """

# create coverage from bam file (for illumina / nanopore / rnaseq coverage)
rule bam_to_bedgraph:
    input:
         bam="{aln}.bam"
    output:
         cov="{aln}.bedgraph"
    shell:
        """
      	bedtools genomecov -ibam {input.bam} -bga -split > {output.cov}
	"""

rule bedgraph_to_lowcov:
    input:
         bedgraph="{name}.bedgraph"
    output:
         r"{name}-low{mincov,[0-9]+}.bed"
    shell:
        """
        # get regions with coverage lower than {wildcards.mincov}
        perl -lane 'print if $F[3]<{wildcards.mincov}' {input} > {output}
        """

rule total_zero:
    input:
         bedgraph="{name}.bedgraph"
    output:
         "{name}-zerocov.txt"
    shell:
        """
        # get sum of length of regions with 0 cov
        perl -lane '$s+=$F[2]-$F[1] if $F[3]==0; END {{ print $s; }}' {input} > {output}
        """
