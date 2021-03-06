# build_shazam.mak
# JULY 2016
# chrisw
# build shazam plots for a pair of contrasting sample groups
# This makefile uses shazam scripts available from:
# 1) https://github.com/UCSC-MedBook/external-tools/tree/master/shazam
# 2) https://github.com/ucscCancer/pathmark-scripts

# need to access non-shazam-specific scripts
# filterRowNames.py
LIB_DIR=../lib

# directory of pid_xxx_pathway.tab files
PATHWAYS_DIR=pathways

# tsv describing pathways in PATHWAYS_DIR with col_1:PID , col_2:description , col_3:source
PATHWAY_LIST_FILE=pathway_pids.tab

# combined paradigm output of real samples and null samples
COMBINED_IPL_DIR=..
COMBINED_IPL_FILE=combined_paradigm_matrices.tab

# suffix part of pathway file names in PATHWAY_DIR
PATHWAY_FILE_SUFFIX=_pathway.tab

# list of pathway ids, usually in the form of PID_xxxx
PATHWAY_IDS=$(shell ls -1 $(PATHWAYS_DIR) | grep "$(PATHWAY_FILE_SUFFIX)" | sed -e 's/$(PATHWAY_FILE_SUFFIX)$$//')

# list of HUGO symbols used to retrieve paradigm features from IPL files
HUGO_SYMBOLS=../data/hugo_symbols.tab

# directory to store pathway-specific IPL files
PATHWAY_IPLS_DIR=pathway_ipls

# contrast file with col_1 sample IDs and col_2 with sample group name
CONTRAST_FILE=all_samples.tab

# sample group names from CONTRAST_FILE to use for contrast
SAMPLE_GROUP_1="real"
# to compare to nulls, use SAMPLE_GROUP_2="Null"
#SAMPLE_GROUP_2="Small Cell"
SAMPLE_GROUP_2="null"

HTML_DIR=$(SAMPLE_GROUP_1)_$(SAMPLE_GROUP_2)_HTML

#: MAKE EVERYTHING !
all: all_samples.tab build_shazams.done

#: build the pathway shazam plots and html as well as the summary html page
build_shazams.done: filter_pathway_ipls.done
	mkdir $(HTML_DIR) ;
	\
	bash build_shazam_htmls.sh \
		$(PATHWAY_IPLS_DIR) \
		$(PATHWAYS_DIR) \
		$(PATHWAY_LIST_FILE) \
		$(CONTRAST_FILE) \
		$(SAMPLE_GROUP_1) \
		$(SAMPLE_GROUP_2) \
		index.html ;
	\
	mv index.html html/. ;
	\
	mv html/* $(HTML_DIR)/. ;
	\
	rm -rf html ;
	\
	date >> $@ ;
	\

#: filter down the HUGO-only combined IPL file down to pathway-specific IPL files
filter_pathway_ipls.done: hugo_only_combined_ipls.done
	\
	mkdir -p $(PATHWAY_IPLS_DIR) ;
	\
	for pathway_id in $(PATHWAY_IDS) ; do \
		echo $${pathway_id} ; \
		\
		cat $(PATHWAYS_DIR)/$${pathway_id}$(PATHWAY_FILE_SUFFIX) \
		| awk 'BEGIN{FS="\t";OFS=FS}{if (NF < 3) {print;};}END{};' \
		| cut -f 2 \
		> 1.tmp ; \
		\
		python $(LIB_DIR)/filterRowNames.py \
			-f 1.tmp \
			-s "_transpose_$${pathway_id}.out" \
			$(COMBINED_IPL_DIR)/$(COMBINED_IPL_FILE)_HUGO_ONLY ; \
		\
		mv $(COMBINED_IPL_DIR)/$(COMBINED_IPL_FILE)_HUGO_ONLY_transpose_$${pathway_id}.out $(PATHWAY_IPLS_DIR)/. ; \
		\
		rm -f 1.tmp ; \
		\
	done ;
	\
	rm -f 1.tmp ;
	\
	date >> $@ ;
	\

#: get all sample names, including nulls
all_samples.tab: hugo_only_combined_ipls.done
	head -n 1 $(COMBINED_IPL_DIR)/$(COMBINED_IPL_FILE)_HUGO_ONLY \
	| fmt -1 \
	| awk 'BEGIN{FS="\t"; OFS=FS;}{input=$$0; if ($$1~/^na_iter_/){category="null";} else{category="real";} print input FS category;}END{}' \
	> 1.tmp ;
	\
	tail -n +2 1.tmp \
	> 2.tmp ;
	\
	mv 2.tmp $@ ;
	\
	rm -f 1.tmp 2.tmp ;
	\

#: filter the merged IPL file down to just HUGO symbol features
hugo_only_combined_ipls.done:
	python $(LIB_DIR)/filterRowNames.py \
		-f $(HUGO_SYMBOLS) \
		-s "_HUGO_ONLY" \
		$(COMBINED_IPL_DIR)/$(COMBINED_IPL_FILE) ; \
	\
	date >> $@ ;
	\

clean:
	rm -f all_samples.tab index.html hugo_only_combined_ipls.done filter_pathway_ipls.done build_shazams.done copy_index_html_to_dir.done ;
	\
	rm -rf $(HTML_DIR) html $(PATHWAY_IPLS_DIR) ;
	\
	rm -f $(wildcard *.tmp) ;
	\
	rm -f $(wildcard $(COMBINED_IPL_DIR)/$(COMBINED_IPL_FILE)_HUGO_ONLY*) ;
	\
