default: all

include datasets.mk

BIN_PATH=~/Projects/coati/build/release/src/coati-sample

all:

.PHONY: all default

results/phase_type_data.csv: 
	( cd data/filtered_cds && \
	Rscript --vanilla ../../scripts/create_phase_type_data.R ) > $@

results/params/%.final.json: data/filtered_cds/%/.script_done
	Rscript --vanilla phase_im.R data/filtered_cds/$* results/params $(BIN_PATH)
	cp $$(ls -At results/params/$*.*.final.json | head -n 1) $@

results/estimated_params.json: $(addprefix results/params/,$(addsuffix .final.json,$(DATASETS)))
	jq -s '[.[][-1]]' $^ > $@

results/estimated_params.csv: results/estimated_params.json
	Rscript --vanilla scripts/make_params_file.R < $< > $@

results/params/%.samples.rds: results/params/%.final.json
	cp $$(ls -At results/params/$*.*.samples.rds | head -n 1) $@

final_rds: $(addprefix results/params/,$(addsuffix .samples.rds,$(DATASETS)))

.PHONY: final_rds

results/best_aln/%/.script_done: results/params/%.final.json results/params/%.samples.rds
	Rscript --vanilla phase_best_aln.R results/params/$*.samples.rds results/params/$*.final.json
	touch $@

results/best_aln_sum/%.csv: results/best_aln/%/.script_done
	(cd results/best_aln && Rscript --vanilla ../../scripts/summarize_aln.R $*) > $@

best_aln_sum: $(addprefix results/best_aln_sum/,$(addsuffix .csv,$(DATASETS)))

.PHONY: best_aln_sum