BIN_PATH=~/Projects/coati/build/release/src/coati-sample

default: all

all:

results/params/%.final.json: data/filtered_cds/%/.script_done
	Rscript --vanilla phase_im.R data/filtered_cds/$* results/params $(BIN_PATH)
	cp $$(ls -At results/params/$*.*.final.json | head -n 1) $@
