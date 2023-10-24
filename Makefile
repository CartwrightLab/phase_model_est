include datasets.mk

BIN_PATH=~/Projects/coati/build/release/src/coati-sample

default: all

all:

results/params/%.final.json: data/filtered_cds/%/.script_done
	Rscript --vanilla phase_im.R data/filtered_cds/$* results/params $(BIN_PATH)
	cp $$(ls -At results/params/$*.*.final.json | head -n 1) $@


results/estimated_params.json: $(addprefix results/params/,$(addsuffix .final.json,$(DATASETS)))
	jq -s '[.[][-1]]' $^ > $@

results/estimated_params.csv: results/estimated_params.json
	Rscript --vanilla scripts/make_params_file.R < $< > $@

results/phase_type_data.csv: 
	( cd data/filtered_cds && \
	Rscript --vanilla ../../scripts/create_phase_type_data.R ) > $@
