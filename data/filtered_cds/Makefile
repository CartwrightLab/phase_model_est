default: all

all:

data/filtered_cds/02_MusRat/.script_done:
	Rscript --vanilla ../../scripts/filter_cds.R 02_MusRat/ ../../raw_data/ZouZhang/02_MusRat_aligned_cds/ '*.fasta'
	touch $@

