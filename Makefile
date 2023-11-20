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

results/sampled_aln/%/.script_done: results/params/%.final.json results/params/%.samples.rds
	Rscript --vanilla phase_sample_aln.R results/params/$*.samples.rds results/params/$*.final.json
	touch $@

results/best_aln_sum/%.csv: results/best_aln/%/.script_done
	(cd results/best_aln && Rscript --vanilla ../../scripts/summarize_aln.R $*) > $@

results/sampled_aln_sum/%.csv: results/sampled_aln/%/.script_done
	(cd results/sampled_aln && Rscript --vanilla ../../scripts/summarize_aln.R $*) > $@

results/sampled_aln_dnds/%.csv: results/sampled_aln/%/.script_done
	(cd results/sampled_aln && Rscript --vanilla ../../scripts/summarize_dnds.R $*) > $@

best_aln_sum: $(addprefix results/best_aln_sum/,$(addsuffix .csv,$(DATASETS)))

sampled_aln_sum: $(addprefix results/sampled_aln_sum/,$(addsuffix .csv,$(DATASETS)))

sampled_aln_dnds: $(addprefix results/sampled_aln_dnds/,$(addsuffix .csv,$(DATASETS)))

.PHONY: best_aln_sum sampled_aln_sum sampled_aln_dnds

results/observed_indels.csv.gz: $(shell find results/sampled_aln_sum -type f -name '*.csv')
	Rscript --vanilla scripts/make_obs_indels.R | gzip > $@

results/observed_dnds.csv.gz: $(shell find results/sampled_aln_dnds -type f -name '*.csv')
	Rscript --vanilla scripts/make_obs_dnds.R | gzip > $@

results/obs_indels_summary.csv: results/observed_indels.csv.gz
	Rscript --vanilla scripts/make_obs_indels_summary.R $< > $@
