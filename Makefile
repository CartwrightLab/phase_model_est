BIN_PATH=~/Projects/coati/build/release/src/coati-sample

default: all

all:

.PHONY: results/params/02_MusRat

results/params/02_MusRat: results/params/02_MusRat.final.json


results/params/02_MusRat.final.json:
	Rscript --vanilla phase_im.R data/filtered_cds/02_MusRat/ results/params $(BIN_PATH)
	cp $$(ls -At results/params/02_MusRat.*.final.json | head -n 1) $@
