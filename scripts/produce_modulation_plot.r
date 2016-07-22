#!/usr/bin/Rscript

#YN 20130123

#usage, options and doc goes here
argspec <- c("produce_modulation_plot.r - Josh's way of modulating paradigm output.
Usage:
    produce_modulation_plot.r input_tcga_zmatrix.tab input_null_zmatrix.tab superpathway.tab 0.03 0 10 0.5 output.pdf
Example:
	Rscript produce_modulation_plot.r input_tcga_zmatrix.tab input_null_zmatrix.tab superpathway.tab 0.03 0 10 0.5 output.pdf
Options:
	input tcga z-score matrix (annotated by row and column names)
	input null z-score matrix (annotated by row and column names)
	pathway file
	percent of samples in which sd has to be greater than the specified value
	lower bound for sd on the plot
	upper bound for sd on the plot
	sd step
	output file")
	
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

compute_binomial_p <- function(index, na_list, real_list, total_genes){
	nas <- na_list[index]
	reals <- real_list[index]
	p <- nas/total_genes
	prob <- dbinom(reals, total_genes, p)
	#prob <- choose(total_genes, reals) * (p^reals) * ((1-p)^(total_genes-reals))
	
	result <- -1*log10(prob)
	if(is.infinite(result)){
		result <- 0
	}
	#print(paste("nas = ", nas, sep=""))
	#print(paste("total_genes = ", total_genes, sep=""))
	#print(paste("p = ", p, sep=""))
	#print(paste("reals = ", reals, sep=""))
	#print(paste("prob = ", prob, sep=""))
	#print(paste("result = ", result, sep=""))
	#print("")	
	return(result)
}

compute_binomial_p_lgamma <- function(index, na_list, real_list, total_genes){
	nas <- na_list[index]
	reals <- real_list[index]
	p <- nas/total_genes
	
	m <- reals
	n <- total_genes-reals
	q <- 1 - p
	
	tmp1 <- lgamma(m + n + 1.0)
	tmp1 <- tmp1 - lgamma(n + 1.0) + lgamma(m + 1.0)
	tmp1 <- tmp1 + m*log(p) + n*log(q)
	result <- exp(tmp1)
	
	#result <- -1*log10(prob)
	#if(is.infinite(result)){
	#	result <- 0
	#}
	return(result)
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
	sink('/dev/null') 
	
	if(length(argv) == 1){
		if(argv==c('--help')){ 
			write(argspec, stderr());
			q();
		}
	}
		
	if(!(length(argv) == 8)){
		write("ERROR: invalid number of arguments is specified", stderr());
		q();
	}
	
	#store command line arguments in variables:
	input_tcga_zmatrix_file <- argv[1]
	input_null_zmatrix_file <- argv[2]
	pathway_file <- argv[3]
	sample_percent_cutoff <- as.numeric(argv[4])
	std_lower_bound <- as.numeric(argv[5])
	std_upper_bound <- as.numeric(argv[6])
	std_step <- as.numeric(argv[7])
	out_pdf_file <- argv[8]
	
	#out_pdf_file <- "/Users/ynewton/school/ucsc/projects/stuart_lab/prostate_cancer/modulation_sample_specific.pdf"
	#protein_list_file <- "/Users/ynewton/school/ucsc/projects/stuart_lab/prostate_cancer/protein_coding_genes.tab"
	#sample_percent_cutoff <- 0.03
	#std_step <- 0.5
	#std_lower_bound <- 0
	#std_upper_bound <- 10
	#input_tcga_zmatrix_file <- "/Users/ynewton/school/ucsc/projects/stuart_lab/prostate_cancer/python_tcga_sample_specific.tab"
	#input_null_zmatrix_file <- "/Users/ynewton/school/ucsc/projects/stuart_lab/prostate_cancer/python_null_sample_specific.tab"
	#pathway_file <- "/Users/ynewton/school/ucsc/projects/stuart_lab/prostate_cancer/superpathway.tab"
	
	#read the input file(s):
	tcga_z_matrix_all <- read_matrix(input_tcga_zmatrix_file)
	null_z_matrix_all <- read_matrix(input_null_zmatrix_file)
	
	pathway <- as.matrix(read.table(file=pathway_file, header = FALSE, sep = "\t", quote = "", fill=TRUE))
	pathway_proteins <- pathway[pathway[,1] == "protein",]
	genes <- pathway_proteins[,2]

	total_genes <- length(genes)
	tcga_z_matrix <- tcga_z_matrix_all[rownames(tcga_z_matrix_all) %in% genes,]
	null_z_matrix <- null_z_matrix_all[rownames(null_z_matrix_all) %in% genes,]
	
	sds <- seq(std_lower_bound,std_upper_bound,std_step)
	num_genes_real <- c()
	for(std in sds)
	{
		#print(std)
		significant_genes <- apply(tcga_z_matrix, 1, significant_in_enough_samples, std=std, sample_percent_cutoff=sample_percent_cutoff)
		num_genes_real <- c(num_genes_real, sum(unlist(significant_genes)))
	}
	num_genes_na <- c()
	for(std in sds)
	{
		#print(std)
		significant_genes <- apply(null_z_matrix, 1, significant_in_enough_samples, std=std, sample_percent_cutoff=sample_percent_cutoff)
		num_genes_na <- c(num_genes_na, sum(unlist(significant_genes)))
	}	
		
	#figure out probability of drawing that many reals
	#log_probs <- unlist(lapply(c(1:length(sds)), compute_binomial_p, na_list=num_genes_na, real_list=num_genes_real, total_genes=length(genes)))
	
	ml = log((num_genes_real+1)/(num_genes_na+1))
	z = ((num_genes_real+1)-(num_genes_na+1)) / sqrt((num_genes_na+1)*(1-(num_genes_na+1)/length(genes)))
	
	suppressPackageStartupMessages(library(calibrate))
	pdf(out_pdf_file, bg="white")
	par(mfrow=c(1,1))
	max_plot_value <- max(max(num_genes_real), max(num_genes_na))
	scaling_factor1 = max_plot_value/max(z)
	scaling_factor2 = max_plot_value/max(ml)
	plot(x=sds, y=num_genes_real, pch=20, cex=0.7, type="o", col="blue", xlab="Standard deviations", ylab=paste("Num genes with at least ", sample_percent_cutoff, " samples with significant score", sep=""), main="Modulation curves") #, ylim=c(0,max_plot_value+max_plot_value/5))
	axis(1, at=min(sds):max(sds))	
	lines(x=sds, y=num_genes_na, pch=24, cex=0.7, type="o", col="green")
	lines(x=sds, y=scaling_factor1*z, pch=22, cex=0.7, type="o", col="red")
	#lines(x=sds, y=scaling_factor2*ml, pch=18, cex=0.7, type="o", col="purple")
	#legend("topright", c("Reals", "Nulls", "Significance z-score", "ML"), cex=0.5, col=c("blue", "green", "red", "purple"), lty=1)
	legend("topright", c("Reals", "Nulls", "Significance z-score"), cex=0.5, col=c("blue", "green", "red"), lty=1)	
	try(textxy(sds, scaling_factor1*z, labs=sds));
	#try(textxy(sds, scaling_factor2*ml, labs=sds));	
	dev.off()
}

main(commandArgs(TRUE))
