#!/usr/bin/Rscript

#YN 20130123

#usage, options and doc goes here
argspec <- c("filter_paradigm_based_on_sd.r - Josh's way of modulating paradigm output.
Usage:
    filter_paradigm_based_on_sd.r input_paradigm_data.tab input_tcga_zmatrix.tab 0.03 2.5 filtered_paradigm.tab
Example:
	Rscript filter_paradigm_based_on_sd.r input_paradigm_data.tab input_tcga_zmatrix.tab 0.03 2.5 filtered_paradigm.tab
Options:
	input matrix (annotated by row and column names)
	output file is specified through redirect character >")
	
read_matrix <- function(in_file){
	header <- strsplit(readLines(con=in_file, n=1), "\t")[[1]]
	cl.cols<- 1:length(header) > 1
	data_matrix.df <- read.delim(in_file, header=TRUE, row.names=NULL, stringsAsFactors=FALSE, na.strings="NA", check.names=FALSE)
	data_matrix <- as.matrix(data_matrix.df[,cl.cols])
	rownames(data_matrix) <- data_matrix.df[,1]
	return(data_matrix)
}

write_matrix <- function(data_matrix, file_name){
	header <- append(c("gene"), colnames(data_matrix))
	write.table(t(header), file_name, quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	write.table(data_matrix, file_name, quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE, append=TRUE)
}

significant_in_enough_samples <- function(zrow, std, sample_percent_cutoff){
	num_outside_std <- length(zrow[abs(zrow) > std])
	proportion_outside_std <- num_outside_std / length(zrow)
	if(proportion_outside_std >= sample_percent_cutoff)
	{
		return(1)
	}
	else
	{
		return(0)
	}
}

main <- function(argv) {  
	if(length(argv) == 1){
		if(argv==c('--help')){ 
			write(argspec, stderr());
			q();
		}
	}
		
	if(!(length(argv) == 5)){
		write("ERROR: invalid number of arguments is specified", stderr());
		q();
	}
	
	#store command line arguments in variables:
	input_paradigm_file <- argv[1]
	input_tcga_zmatrix_file <- argv[2]
	sample_percent_cutoff <- as.numeric(argv[3])
	std <- as.numeric(argv[4])
	output_file <- argv[5]
	
	#input_paradigm_file <- "/Users/ynewton/school/ucsc/projects/stuart_lab/prostate_cancer/paradigm_output_superpathway.tab"
	#input_tcga_zmatrix_file <- "/Users/ynewton/school/ucsc/projects/stuart_lab/prostate_cancer/python_tcga_nonsample_specific.tab"	
	#output_file <- "/Users/ynewton/school/ucsc/projects/stuart_lab/prostate_cancer/reduced_paradigm_output_3.5sd.tab"	
	#sample_percent_cutoff <- 0.03
	#std <- 3.5
	
	#read the input file(s):
	data_matrix <- read_matrix(input_paradigm_file)
	tcga_z_matrix <- read_matrix(input_tcga_zmatrix_file)
	tcga <- data_matrix[,substr(colnames(data_matrix), 1, 8) != "na_iter_"]	
	genes <- rownames(tcga)

	gene_positions <- apply(tcga_z_matrix, 1, significant_in_enough_samples, std=std, sample_percent_cutoff=sample_percent_cutoff)
	selected_genes <- genes[which(gene_positions == 1)]
	
	result_matrix <- tcga[rownames(tcga) %in% selected_genes,]
	write_matrix(result_matrix, output_file)
}

main(commandArgs(TRUE))
