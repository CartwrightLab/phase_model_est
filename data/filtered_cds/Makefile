# '*.fasta'
datasets1=01_FcaCaf 02_MusRat 04_Gallus 08_Arab 09_SarMon 12_fugu
# '*_NTaln.fasta'
datasets2=06_Nematode 10_ants 11_malaria 13_Phytophthora 14_Fusarium 15_Oryza \
16_Solanum
# '*.fa'
datasets3=05_Droso
# '*-1_highqual.codon.mfa'
datasets4=07_yeast
# '*.fasta'
datasets5=ATGC001 ATGC044 ATGC073 ATGC112 ATGC143 ATGC171 ATGC192 \
ATGC242 ATGC279 ATGC352 ATGC389 ATGC013 ATGC045 ATGC075 ATGC117 ATGC145 \
ATGC173 ATGC201 ATGC243 ATGC301 ATGC361 ATGC393 ATGC015 ATGC050 ATGC078 \
ATGC120 ATGC146 ATGC178 ATGC202 ATGC244 ATGC309 ATGC369 ATGC420 ATGC016 \
ATGC056 ATGC095 ATGC127 ATGC152 ATGC181 ATGC223 ATGC246 ATGC330 ATGC371 \
ATGC428 ATGC021 ATGC058 ATGC099 ATGC134 ATGC155 ATGC182 ATGC232 ATGC258 \
ATGC332 ATGC372 ATGC431 ATGC022 ATGC069 ATGC110 ATGC136 ATGC159 ATGC186 \
ATGC235 ATGC260 ATGC339 ATGC374 ATGC028 ATGC071 ATGC111 ATGC138 ATGC169 \
ATGC188 ATGC239 ATGC277 ATGC341 ATGC376

datasets=$(datasets1) $(datasets2) $(datasets3) $(datasets4) $(datasets5)

default: all

all: $(datasets)

.PHONY: $(datasets)

$(datasets) : % : %/.script_done

$(addsuffix /.script_done,$(datasets1)) : %/.script_done : ../../raw_data/ZouZhang/%_aligned_cds
	mkdir -p $*
	Rscript --vanilla ../../scripts/filter_cds.R $*/ ../../raw_data/ZouZhang/$*_aligned_cds/ '*.fasta'
	touch $@

$(addsuffix /.script_done,$(datasets2)) : %/.script_done : ../../raw_data/ZouZhang/%_aligned_cds
	mkdir -p $*
	Rscript --vanilla ../../scripts/filter_cds.R $*/ ../../raw_data/ZouZhang/$*_aligned_cds/ '*_NTaln.fasta'
	touch $@

$(addsuffix /.script_done,$(datasets3)) : %/.script_done : ../../raw_data/ZouZhang/%_aligned_cds
	mkdir -p $*
	Rscript --vanilla ../../scripts/filter_cds.R $*/ ../../raw_data/ZouZhang/$*_aligned_cds/ '*.fa'
	touch $@

$(addsuffix /.script_done,$(datasets4)) : %/.script_done : ../../raw_data/ZouZhang/%_aligned_cds
	mkdir -p $*
	Rscript --vanilla ../../scripts/filter_cds.R $*/ ../../raw_data/ZouZhang/$*_aligned_cds/ '*-1_highqual.codon.mfa'
	touch $@

$(addsuffix /.script_done,$(datasets5)) : %/.script_done : ../../raw_data/ZouZhang/%
	mkdir -p $*
	Rscript --vanilla ../../scripts/filter_cds.R $*/ ../../raw_data/ZouZhang/$*/ '*.fasta'
	touch $@
