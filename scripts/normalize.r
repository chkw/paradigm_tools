#!/usr/bin/Rscript

#Yulia Newton, last updated 20121217 v.3

#usage, options and doc goes here
argspec <- c("normalize.r - takes any flat file and normalizes the rows or the columns using various normalizations (median_shift, mean_shift, t_statistic (z-score), exp_fit, normal_fit, weibull_0.5_fit, weibull_1_fit, weibull_1.5_fit, weibull_5_fit). Requires a single header line and a single cloumn of annotation.
Usage:
    normalize.r input.tab norm_type norm_by > output.tab
Example:
	Rscript normalize.r test_matrix.tab median_shift column > output.tab
	Rscript normalize.r test_matrix.tab mean_shift row normals.tab > output.tab
Options:
	input matrix (annotated by row and column names)
	normalization type; available options:
		median_shift - shifts all values by the median or the row/column if no normals are specified, otherwise shifts by the median of normals
		mean_shift - shifts all values by the mean or the row/column if no normals are specified, otherwise shifts by the mean of normals
		z_score - converts all values to z-scores; if normals are specified then converts to z-scores within normal and non-normal classes separately 
		exponential_fit - (only by column) ranks data and quantile normalizes to exponential CDF
		normal_fit - (only by column) ranks data and transforms normal CDF
		weibull_0.5_fit - (only by column) ranks data and transforms Weibull CDF with scale parameter = 1 and shape parameter = 0.5
		weibull_1_fit - (only by column) ranks data and transforms Weibull CDF with scale parameter = 1 and shape parameter = 1
		weibull_1.5_fit - (only by column) ranks data and transforms Weibull CDF with scale parameter = 1 and shape parameter = 1.5
		weibull_5_fit - (only by column) ranks data and transforms Weibull CDF with scale parameter = 1 and shape parameter = 5
		log_2 - log base 2 transform
		log_10 - log base 10 transform
		exp - exponentiate the data values
		beta_0.5_0.5_fit - (only by column) ranks data and quantile normalizes to Beta(alpha=0.5, beta=0.5) 
		beta_0.05_0.05_fit - (only by column) ranks data and quantile normalizes to Beta(alpha=0.05, beta=0.05) 
		beta_2_2_fit - (only by column) ranks data and quantile normalizes to Beta(alpha=2, beta=2)
		right_shift - shifts each gene so the minimum value is 0 instead of a negative (shifts to the right)
		uniform_0_1_fit - (only by column) ranks data and transforms Uniform(0,1) CDF
		dif_by_normals - similar to t-statistic but performed for each tumor sample (done for each gene on a single tumor measurement against distribution of normal samples); (x-mean_normal)/sqrt(var_all_tumor/num_tumor_samples + var_all_normals/num_normal_samples)
		center - take the number line of values and make the center of it to be new 0 (e.g. 0 though 12 values will become -6 through 6)
	normalization by:
		row
		column
	normals_file is an optional parameter which contains either a list of column headers from the input matrix, which should be considered as normals, or a matrix of normal samples
	output file is specified through redirect character >")

read_matrix <- function(in_file){
	header <- strsplit(readLines(con=in_file, n=1), "\t")[[1]]
	cl.cols<- 1:length(header) > 1
	data_matrix.df <- read.delim(in_file, header=TRUE, row.names=NULL, stringsAsFactors=FALSE, na.strings="NA", check.names=FALSE)
	data_matrix <- as.matrix(data_matrix.df[,cl.cols])
	rownames(data_matrix) <- data_matrix.df[,1]
	return(data_matrix)
}

write_matrix <- function(data_matrix){
	header <- append(c("genes"), colnames(data_matrix))
	if(is.null(colnames(data_matrix)) || length(colnames(data_matrix))==1)
	{
		header <- c("genes", "norm_score")
	}	
	write.table(t(header), stdout(), quote=FALSE, sep="\t", row.names=FALSE, col.names=FALSE)
	write.table(data_matrix, stdout(), quote=FALSE, sep="\t", row.names=TRUE, col.names=FALSE)	
}

read_normals <- function(in_file){
	#return(as.matrix(read.table(in_file, header=FALSE, sep="", as.is = TRUE))[, 1])
	return(as.matrix(read.delim(in_file, header=FALSE, sep="", as.is = TRUE, quote="")))
}

normalize <- function(data_matrix, norm_type, normals_list, tumors_list){
	if(norm_type == 'MEDIAN_SHIFT'){
		return(shift(data_matrix, 'MEDIAN', normals_list, tumors_list))
	}
	else if(norm_type == 'MEAN_SHIFT'){
		return(shift(data_matrix, 'MEAN', normals_list, tumors_list))
	}
	else if(norm_type == 'Z_SCORE'){
		return(compute_z_score(data_matrix, normals_list, tumors_list))
	}
	else if(norm_type == 'EXPONENTIAL_FIT'){
		return(fit_distribution(data_matrix, 'EXPONENTIAL'))
	}
	else if(norm_type == 'BETA_0.5_0.5_FIT'){
		return(fit_distribution(data_matrix, 'BETA_0.5_0.5'))
	}	
	else if(norm_type == 'BETA_0.05_0.05_FIT'){
		return(fit_distribution(data_matrix, 'BETA_0.05_0.05'))
	}		
	else if(norm_type == 'BETA_2_2_FIT'){
		return(fit_distribution(data_matrix, 'BETA_2_2'))
	}		
	else if(norm_type == 'NORMAL_FIT'){
		return(fit_distribution(data_matrix, 'NORMAL'))
	}
	else if(norm_type == 'WEIBULL_0.5_FIT'){
		return(fit_distribution(data_matrix, 'WEIBULL_0.5'))
	}
	else if(norm_type == 'WEIBULL_1_FIT'){
		return(fit_distribution(data_matrix, 'WEIBULL_1'))
	}
	else if(norm_type == 'WEIBULL_1.5_FIT'){
		return(fit_distribution(data_matrix, 'WEIBULL_1.5'))
	}
	else if(norm_type == 'WEIBULL_5_FIT'){
		return(fit_distribution(data_matrix, 'WEIBULL_5'))
	}else if(norm_type == 'UNIFORM_0_1_FIT'){
		return(fit_distribution(data_matrix, 'UNIFORM_0_1'))
	}else if(norm_type == 'LOG_2'){
		return(transform(data_matrix, 'LOG_2'))
	}else if(norm_type == 'LOG_10'){
		return(transform(data_matrix, 'LOG_10'))
	}else if(norm_type == 'EXP'){
		return(transform(data_matrix, 'EXP'))
	}else if(norm_type == 'RIGHT_SHIFT'){
		return(right_shift_matrix(data_matrix))
	}else if(norm_type == 'DIF_BY_NORMALS'){
		return(dif_by_normals(data_matrix, normals_list, tumors_list))
	}else if(norm_type == 'CENTER'){
		return(center(data_matrix))
	}else{
		write("ERROR: unknown normalization type", stderr());
		q();		
	}
}

center <- function(mtrx){
	return(t(apply(mtrx,1,function(x){shift <- (max(x) - min(x))/2; middle <- max(x) - shift; return(x-middle);})))
}

right_shift_matrix <- function(data_matrix){
	return(t(apply(data_matrix, 1, right_sift_vector)))
}

right_sift_vector <- function(vect){
	min_vect <- min(vect)
	#write(min_vect, stderr());
	if(is.na(min_vect)){
		return(vect)
	}
	
	if(min_vect < 0){
		return(vect+(-1)*min_vect)
	}else{
		return(vect)
	}
}

transform <- function(data_matrix, transformation_type){
	if(transformation_type=="LOG_2" || transformation_type=="LOG_10")
	{
		tranform_type = "LOG"
		if(transformation_type=="LOG_2"){
			transform_base = 2;
		}else if(transformation_type=="LOG_10"){
			transform_base = 10;
		}else{
			write("ERROR: unknown log transform type/base detected", stderr());
			q();		
		}
		return(t(apply(data_matrix+1, 1, log, base=transform_base)))
	}else if(transformation_type=="EXP"){
		return(t(apply(data_matrix+1, 1, exp)))
	}else{
		write("ERROR: unknown transform type detected", stderr());
		q();	
	}
}

shift <- function(data_matrix, shift_type, normals_list, tumors_list){
	return(t(apply(data_matrix, 1, shift_normalize_row, norm_type=shift_type, normals_list=normals_list, tumors_list=tumors_list)))
}

shift_normalize_row <- function(data_row, norm_type, normals_list, tumors_list){
	if(length(normals_list) == 0){	#no normals are specified
		if(norm_type == 'MEDIAN'){
			row_stat <- median(data_row)
		}
		else if(norm_type == 'MEAN'){
			row_stat <- mean(data_row)
		}
		return(unlist(lapply(data_row, function(x){return(x - row_stat);})))
	}
	else{	#normals are specified
		normal_values <- data_row[normals_list]
		tumor_columns <- data_row[tumors_list]
		
		if(norm_type == 'MEDIAN'){
			row_stat <- median(normal_values)
		}
		else if(norm_type == 'MEAN'){
			row_stat <- mean(normal_values)
		}
		return(unlist(lapply(tumor_columns, function(x){return(x - row_stat);})))
	}
}

dif_by_normals <- function(data_matrix, normals_list, tumors_list){
	return(t(apply(data_matrix, 1, t_stat_normalize_row, normals_list=normals_list, tumors_list=tumors_list)))
}

t_stat_normalize_row <- function(data_row, normals_list, tumors_list){
	if(length(normals_list) == 0){	#no normals are specified
		return(data_row)
	}
	else{	#normals are specified
		#(x-mean_normal)/sqrt(var_all_tumor/num_tumor_samples + var_all_normals/num_normal_samples)
		normal_values <- data_row[normals_list]
		normal_mean <- mean(normal_values)
		normal_sd <- sd(normal_values)
		normal_var <- normal_sd^2
		
		tumor_values <- data_row[tumors_list]
		tumor_sd <- sd(tumor_values)
		tumor_var <- tumor_sd^2
				
		normalized_values <- unlist(lapply(tumor_values, function(x){return((x - normal_mean)/sqrt(normal_var/length(normals_list) + tumor_var/length(tumors_list)));}))
		return(normalized_values)
	}
}

compute_z_score <- function(data_matrix, normals_list, tumors_list){
	return(t(apply(data_matrix, 1, z_score_normalize_row, normals_list=normals_list, tumors_list=tumors_list)))
}

z_score_normalize_row <- function(data_row, normals_list, tumors_list){
	if(length(normals_list) == 0){	#no normals are specified
		row_mean <- mean(data_row)
		row_sd <- sd(data_row)
		return(unlist(lapply(data_row, function(x){return((x - row_mean)/row_sd);})))
	}
	else{	#normals are specified
		normal_values <- data_row[normals_list]
		normal_mean <- mean(normal_values)
		normal_sd <- sd(normal_values)
		normalized_normals <- unlist(lapply(normal_values, function(x){return((x - normal_mean)/normal_sd);}))
		
		tumor_values <- data_row[tumors_list]
		normalized_tumors <- unlist(lapply(tumor_values, function(x){return((x - normal_mean)/normal_sd);}))
		
		return(append(normalized_normals, normalized_tumors))
	}
}

rankNA <- function(col){		#originally written by Dan Carlin
	col[!is.na(col)]<-(rank(col[!is.na(col)])/sum(!is.na(col)))-(1/sum(!is.na(col)))
	return(col)
}

fit_distribution <- function(data_matrix, dist){
	if(dist == 'EXPONENTIAL'){
		ranked_data_matrix <- as.matrix(apply(data_matrix,1,rankNA))	#idea by Dan Carlin
		return(as.matrix(apply(ranked_data_matrix, 1, qexp)))
	}
	else if(dist == 'NORMAL'){
		ranked_data_matrix <- as.matrix(apply(data_matrix,1,rankNA))
		return(as.matrix(apply(ranked_data_matrix, 1, qnorm, mean=0, sd=1)))
	}
	else if(dist == 'WEIBULL_0.5'){
		ranked_data_matrix <- as.matrix(apply(data_matrix,1,rankNA))
		returnas.matrix((apply(ranked_data_matrix, 1, qweibull, scale=1, shape=0.5)))
	}
	else if(dist == 'WEIBULL_1'){
		ranked_data_matrix <- as.matrix(apply(data_matrix,1,rankNA))
		return(as.matrix(apply(ranked_data_matrix, 1, qweibull, scale=1, shape=1)))
	}
	else if(dist == 'WEIBULL_1.5'){
		ranked_data_matrix <- as.matrix(apply(data_matrix,1,rankNA))
		return(as.matrix(apply(ranked_data_matrix, 1, qweibull, scale=1, shape=1.5)))
	}
	else if(dist == 'WEIBULL_5'){
		ranked_data_matrix <- as.matrix(apply(data_matrix,1,rankNA))
		return(as.matrix(apply(ranked_data_matrix, 1, qweibull, scale=1, shape=5)))
	}
	else if(dist == 'BETA_0.5_0.5'){
		ranked_data_matrix <- as.matrix(apply(data_matrix,1,rankNA))
		return(as.matrix(apply(ranked_data_matrix, 1, qbeta, shape1=0.5, shape2=0.5)))
	}	
	else if(dist == 'BETA_0.05_0.05'){
		ranked_data_matrix <- as.matrix(apply(data_matrix,1,rankNA))
		return(as.matrix(apply(ranked_data_matrix, 1, qbeta, shape1=0.05, shape2=0.05)))
	}		
	else if(dist == 'BETA_2_2'){
		ranked_data_matrix <- as.matrix(apply(data_matrix,1,rankNA))
		return(as.matrix(apply(ranked_data_matrix, 1, qbeta, shape1=2, shape2=2)))
	}
	else if(dist == 'UNIFORM_0_1'){
		ranked_data_matrix <- as.matrix(apply(data_matrix,1,rankNA))
		return(as.matrix(apply(ranked_data_matrix, 1, qunif, min=0, max=1)))
	}	
}

main <- function(argv) {  
	#determine if correct number of arguments are specified and if normals are specified
	with_normals = FALSE
	
	if(length(argv) == 1){
		if(argv==c('--help')){ 
			write(argspec, stderr());
			q();
		}
	}
		
	if(!(length(argv) == 3 || length(argv) == 4)){
		write("ERROR: invalid number of arguments is specified", stderr());
		q();
	}
	
	if(length(argv) == 4){
		with_normals = TRUE
		normals_file <- argv[4]
	}
	
	#store command line arguments in variables:
	input_file <- argv[1]
	norm_type <- toupper(argv[2])
	norm_by <- toupper(argv[3])
	
	#input_file <- "/Users/ynewton/school/ucsc/projects/stuart_lab/data_normalization/test_matrix.tab"
	#norm_type <- "MEAN_SHIFT"
	#norm_by <- "ROW"
	#normals_file <- "/Users/ynewton/school/ucsc/projects/stuart_lab/data_normalization/test_matrix2.tab"
	#normals_file2 <- "/Users/ynewton/school/ucsc/projects/stuart_lab/data_normalization/normals.tab"
	
	#read the input file(s):
	data_matrix <- read_matrix(input_file)
	
	if(ncol(data_matrix) == 1){		#if the matrix is a single signature or sample (need for later to make sure we output as a column vector)
		single_sig <- TRUE
	}else{
		single_sig <- FALSE
	}
	
	if(with_normals){
		normals <- read_normals(normals_file)
		if(length(colnames(normals)) == 1){
			normals_indices <- which(colnames(data_matrix) %in% normals)
			tumor_indices <- which(!(colnames(data_matrix) %in% normals))	
		}else{
			normals_numeric <- as.matrix(normals[2:length(normals[,1]),2:length(normals[1,])])
			normals_numeric <- apply(normals_numeric, 2, as.numeric)
			rownames(normals_numeric) <- normals[,1][2:length(normals[,1])]
			colnames(normals_numeric) <- normals[1,][2:length(normals[1,])]
			
			#select only the intersection of the rows between the two matrices:
			if(ncol(normals_numeric) > 1){
				normals_numeric <- as.matrix(normals_numeric[rownames(normals_numeric) %in% rownames(data_matrix),])
				data_matrix <- data_matrix[rownames(data_matrix) %in% rownames(normals_numeric),]
				normals_numeric <- normals_numeric[order(rownames(normals_numeric)),]
				data_matrix <- data_matrix[order(rownames(data_matrix)),]
			}else{
				normals_numeric <- normals_numeric[,1]
				normals_numeric <- normals_numeric[names(normals_numeric) %in% rownames(data_matrix)]
				data_matrix <- data_matrix[rownames(data_matrix) %in% names(normals_numeric),]
				normals_numeric <- as.matrix(normals_numeric[order(names(normals_numeric))])
				data_matrix <- data_matrix[order(rownames(data_matrix)),]
			}
			
			combined_matrix <- cbind(data_matrix, normals_numeric)
			tumor_indices <- c(1:length(data_matrix[1,]))
			normals_indices <- c(length(tumor_indices)+1:length(normals_numeric[1,]))
			data_matrix <- combined_matrix	
		}
	}else{
		normals_indices <- c()
		tumor_indices <- c()
	}
	
	#if normalize by columns then transpose the matrix:
	if(norm_by == 'COLUMN'){
		data_matrix <- as.matrix(t(data_matrix))
	}
	
	#normalize:
	data_matrix <- normalize(data_matrix, norm_type, normals_indices, tumor_indices)
	
	#if normalize by columns then transpose the matrix again since we normalized the transposed matrix by row:
	if(norm_by == 'COLUMN' && !single_sig){
		data_matrix <- as.matrix(t(data_matrix))
	}	
	
	write_matrix(data_matrix)
}

main(commandArgs(TRUE))
