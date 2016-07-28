# JULY 2016
# chrisw
# preprocess data files for use with Paradigm

CNV_FILE = cnv.tab
EXP_FILE = mrna.tab
PATHWAY_FILE = pathway.tab

SUFFIX_1 = _commonCols
SUFFIX_2 = _commonRows

LIB_DIR = lib
DATA_DIR = data

TARGETS = cnv_normalizedSampleNames.tab mrna_normalizedSampleNames.tab sampleNames.tab mrna_meanshifted.tab

all: normalizeNames sampleNames.tab 

# mean-shift normalization for expression data
mrna_meanshifted.tab: $(EXP_FILE)
	Rscript $(LIB_DIR)/normalize.r $< mean_shift row \
	> 1.tmp ;
	\
	mv 1.tmp $@ ;
	\
	rm -f 1.tmp ;
	\

# list of sample names for each input file
sampleNames.tab: cnv_normalizedSampleNames.tab$(SUFFIX_2)$(SUFFIX_1) mrna_normalizedSampleNames.tab$(SUFFIX_2)$(SUFFIX_1)
	rm -f 1.tmp ;
	touch 1.tmp ;
	\
	for file in $? ; do \
	  echo $${file} \
		> 2.tmp ; \
		\
		head -n 1 $${file} \
		| cut -f 2- \
		| transpose.pl \
		>> 2.tmp ; \
		\
		paste 1.tmp 2.tmp \
		> 3.tmp ; \
		\
		mv 3.tmp 1.tmp ; \
		\
		rm -f 2.tmp 3.tmp ; \
		\
	done ;
	\
	cut -f 2- 1.tmp \
	> $@ ;
	\
	rm -f 1.tmp ;
	\

# synchronize column and row names for each input file
normalizeNames: cnv_normalizedSampleNames.tab mrna_normalizedSampleNames.tab
	python $(LIB_DIR)/commonRowOrder.py \
		-v \
		--output_suffix $(SUFFIX_2) \
		cnv_normalizedSampleNames.tab mrna_normalizedSampleNames.tab ;
	\
	python $(LIB_DIR)/commonColumnOrder.py \
		-v \
		-s $(SUFFIX_1) \
		cnv_normalizedSampleNames.tab$(SUFFIX_2) mrna_normalizedSampleNames.tab$(SUFFIX_2) ;
	\

# normalize sample names
mrna_normalizedSampleNames.tab: mrna_meanshifted.tab
	head -n 1 $< \
	| transpose.pl \
	| sed -e 's,^prad_wcdt/,,' \
		-e 's/-BL$$//' \
		-e 's/-Pro/Pro/' \
		-e 's/-BL/-/' \
	| transpose.pl \
	> 1.tmp ;
	\
	tail -n +2 $< \
	| cat.pl 1.tmp - \
	> 2.tmp ;
	\
	mv 2.tmp $@ ;
	\
	rm -f 1.tmp 2.tmp ;
	\

# normalize sample names
cnv_normalizedSampleNames.tab:
	head -n 1 $(CNV_FILE)\
	| transpose.pl \
	| sed -e 's/-BL$$//' \
		-e 's/-Pro/Pro/' \
		-e 's/-BL/-/' \
		-e 's/^DTB-095$$/DTB-095-1/' \
	| transpose.pl \
	> 1.tmp ;
	\
	tail -n +2 $(CNV_FILE) \
	| cat.pl 1.tmp - \
	> 2.tmp ;
	\
	mv 2.tmp $@ ;
	\
	rm -f 1.tmp 2.tmp ;
	\

clean_all: clean_tmp clean_targets
	rm -f $(wildcard *$(SUFFIX_1)) ;
	\
	rm -f $(wildcard *$(SUFFIX_2)) ;
	\

clean_tmp:
	rm -f $(wildcard *.tmp)

clean_targets:
	rm -f $(TARGETS)

