HMMBIN := ./hmm

RUN := $(wildcard $(DIR)/*/reads.fa)
HMMCOUNT := $(patsubst %reads.fa,%hmm.count,$(RUN))

.SECONDARY:
# prevents deletion of intermediate targets in chained rules

.DELETE_ON_ERROR:
# delete targets if a rule fails

SHELL = /bin/bash
.SHELLFLAGS = -o pipefail -c
# pipeline fail if any fails

$(DIR)/%/hmm.stats : $(DIR)/%/reads.fa $(DIR)/hmm-motifs.tsv
	grep $* $(DIR)/hmm-motifs.tsv | perl -lane 'print $$F[0], " ", $$F[2]' > $(DIR)/$*/hmm-motifs.txt
	$(HMMBIN)/rep_hmm.pl $(DIR)/$*/hmm-motifs.txt $(HMMBIN)/rep_hmm.prob $(DIR)/$*/reads.fa $(DIR)/$*/hmm

# sum of counts of match states for each motif
$(DIR)/%/hmm.count : $(DIR)/%/hmm.stats
	perl -lane 'if($$F[1] eq "bg") { $$r=$$F[2]-$$F[3]-$$F[4]; $$r=0 if $$r<0; print join("\t", $$F[0], "bg-middle", $$r); print join("\t", $$F[0], "bg-start", $$F[3]); print join("\t", $$F[0], "bg-end", $$F[4]); } else { my $$s=0; foreach $$x (@F[2..$$#F]) { $$s+=$$x; } print join("\t", @F[0,1], $$s) }' $< > $@

$(DIR)/hmm.count: $(HMMCOUNT)
	grep --with-filename . $^ | perl -lne 's/^[^\/]+\///; s/\/hmm.count:/\t/; print' | sort > $@

$(DIR)/hmm.count.summary: $(DIR)/hmm.count
	$(HMMBIN)/hmm_summary.py < $< > $@

