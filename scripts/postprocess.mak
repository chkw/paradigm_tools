# JULY 2016
# chrisw, ynewton
# basic postprocessing of Paradigm output

PARADIGM_REALS_FILE = paradigm.reals.tab
PARADIGM_NULLS_FILE = paradigm.nulls.tab
PATHWAY_FILE = pathway.tab
LIB_DIR = lib
TEMP_DIR = temp

FILTERED_RESULTS_FILE=filtered_results.tab
CONSTITUENT_PATHWAY_GENE_SETS_FILE=contituent_pathway_gene_sets.tab

TARGETS= combined_paradigm_matrices.tab compute_zscores modulation_plot.pdf paradigm.reals.modulation_filtered_3.tab

test:

all: $(TARGETS) cleanup

# combine paradigm nulls and reals matrices into one:
# python2.7 combine_matrices.py --in_dir temp/ --out_file temp_out.tab --intersect TRUE
combined_paradigm_matrices.tab: $(PARADIGM_REALS_FILE) $(PARADIGM_NULLS_FILE)
	rm -rf $(TEMP_DIR) ;
	\
	mkdir -p $(TEMP_DIR) ;
	\
	cp $(PARADIGM_REALS_FILE) $(TEMP_DIR) ;
	\
	cp $(PARADIGM_NULLS_FILE) $(TEMP_DIR) ;
	\
	python2.7 $(LIB_DIR)/combine_matrices.py \
		--in_dir $(TEMP_DIR) \
		--out_file $@ \
		--intersect TRUE ;
	\

# convert paradigm outputs to z-score matrices:
# python2.7 compute_z_matrices.py --in_paradigm paradigm_output_superpathway.tab --sample_specific_null_pull false --out_tcga yn_tcga.tab --out_null yn_null.tab
compute_zscores: combined_paradigm_matrices.tab
	python $(LIB_DIR)/compute_z_matrices.py \
		--in_paradigm $< \
		--sample_specific_null_pull true \
		--out_tcga reals.z \
		--out_null nulls.z ;
	\

# produce modulation plot:
# Rscript produce_modulation_plot.r input_tcga_zmatrix.tab input_null_zmatrix.tab superpathway.tab 0.03 0 10 0.5 output.pdf
modulation_plot.pdf: reals.z nulls.z
	Rscript $(LIB_DIR)/produce_modulation_plot.r $? $(PATHWAY_FILE) 0.03 0 10 0.5 1.tmp ;
	\
	mv 1.tmp $@ ;
	\
	rm -f 1.tmp ;
	\

# filter paradigm based on the modulation plot results (filter sd may need to change for another run)
# filter sd = 3 for paradigm run 20160721
paradigm.reals.modulation_filtered_3.tab: $(PARADIGM_REALS_FILE) reals.z
	Rscript $(LIB_DIR)/filter_paradigm_based_on_sd.r $? 0.03 3 $@;
	\

clean_all: cleanup
	rm -f reals.z nulls.z ;
	\
	rm -f $(TARGETS) ;
	\

cleanup:
	rm -rf $(TEMP_DIR) ;
	\
	rm -f $(wildcard *.tmp) ;
	\

