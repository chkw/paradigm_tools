#Yulia Newton
#python2.7 compute_z_matrices.py --in_paradigm paradigm_output_superpathway.tab --sample_specific_null_pull false --out_tcga yn_tcga.tab --out_null yn_null.tab

import optparse, sys, numpy, math, time
from numpy import array
parser = optparse.OptionParser()
parser.add_option("--in_paradigm", dest="in_paradigm", action="store", default="", help="")
parser.add_option("--sample_specific_null_pull", dest="sample_specific_null_pull", action="store", default="", help="")
parser.add_option("--out_tcga", dest="out_tcga", action="store", default="", help="")
parser.add_option("--out_null", dest="out_null", action="store", default="", help="")
opts, args = parser.parse_args()

start_time = time.time()
	
#process input arguments:
in_paradigm_file = opts.in_paradigm
sample_specific_null_pull = opts.sample_specific_null_pull.upper()
out_tcga_file = opts.out_tcga
out_null_file = opts.out_null
	
all_null_mean = 0
all_null_sd = 0
in_paradigm = open(in_paradigm_file, 'r')
out_tcga = open(out_tcga_file, 'w')
out_null = open(out_null_file, 'w')
na_cols = []
tcga_cols = []
tcga_sample_to_na_mapping = {}
tcga_col_to_sample_mapping = {}
null_col_to_sample_mapping = {}
line_num = 0
for line in in_paradigm:
	line_elems = line.strip().split("\t")
	data_elems = line_elems[1:]
	if line_num == 0:	#header line
		na_headers = []
		tcga_headers = []
		for i in range(len(data_elems)):
			elem = data_elems[i]
			if elem.startswith("na_iter_"):		#null column
				na_headers.append(elem)
				
				na_cols.append(i)
				tcga_col_to_sample_mapping[i] = elem
				
				tcga_sample_id = elem[10:]
				if not(tcga_sample_id in tcga_sample_to_na_mapping):
					tcga_sample_to_na_mapping[tcga_sample_id] = []
				tcga_sample_to_na_mapping[tcga_sample_id].append(i)
				
				null_col_to_sample_mapping [i] = elem
			else:	#real column
				tcga_headers.append(elem)
				
				tcga_cols.append(i)
				tcga_col_to_sample_mapping[i] = elem
		
		#output headers to files:
		print >> out_tcga, "gene\t"+"\t".join(tcga_headers)
		print >> out_null, "gene\t"+"\t".join(na_headers)
	else:	#data
		print line_elems[0]
		
		#compute real zscore matrix:
		tcga_z_row = []
		sample_index = 0
		for tcga_sample_col_index in tcga_cols:
			#print "tcga_sample_col_index = "+str(tcga_sample_col_index)
			tcga_val = float(data_elems[tcga_sample_col_index])
			#print "tcga_val = "+str(tcga_val)
			na_vals = []
			if sample_specific_null_pull == "FALSE" or sample_specific_null_pull == "NO":					
				if sample_index == 0:
					for col_index in na_cols:
						na_vals.append(float(data_elems[col_index]))
				
					na_vals_array = array(na_vals)
					na_mean = na_vals_array.mean()
					#na_mean = sum(na_vals) / len(na_vals)
					#print "na_mean = "+str(na_mean)
					#na_sd = math.sqrt(sum([math.pow(x-na_mean, 2) for x in na_vals]) / (len(na_vals)-1))
					na_sd = na_vals_array.std(ddof=1)
					#print "na_sd = "+str(na_sd)
					
					all_null_mean = na_mean
					#print "all_null_mean = "+str(all_null_mean)		
					all_null_sd = na_sd
					#print "all_null_sd = "+str(all_null_sd)		
				else:
					na_mean = all_null_mean
					na_sd = all_null_sd
										
			else:
				tcga_sample_id = tcga_col_to_sample_mapping[tcga_sample_col_index]
				#print "tcga_sample_id = "+str(tcga_sample_id)
				sample_specific_na_cols = tcga_sample_to_na_mapping[tcga_sample_id]
				#print "sample_specific_na_cols = "+str(sample_specific_na_cols)
				for col_index in sample_specific_na_cols:
					na_vals.append(float(data_elems[col_index]))
			
				#print "na_vals = "+str(na_vals)
			
				na_vals_array = array(na_vals)
				na_mean = na_vals_array.mean()
				#na_mean = sum(na_vals) / len(na_vals)
				#print "na_mean = "+str(na_mean)
				#na_sd2 = math.sqrt(sum([math.pow(x-na_mean, 2) for x in na_vals]) / (len(na_vals)-1))	#unbiased estimator
				na_sd = na_vals_array.std(ddof=1)
				#print "na_sd = "+str(na_sd)
			
			if not(na_sd == 0):
				tcga_z_score = (tcga_val - na_mean) / na_sd
			else:
				tcga_z_score = 0
			#print "tcga_z_score = "+str(tcga_z_score)
			tcga_z_row.append(tcga_z_score)
			
			sample_index += 1	
		
		#print the tcga row:
		print >> out_tcga, line_elems[0]+"\t"+"\t".join([str(x) for x in tcga_z_row])
		
		#compute null zscore matrix:
		null_z_row = []
		sample_index = 0
		for null_sample_col_index in na_cols:
			#print "null_sample_col_index = "+str(null_sample_col_index)
			null_val = float(data_elems[null_sample_col_index])
			#print "null_val = "+str(null_val)
			na_vals = []
			if sample_specific_null_pull == "FALSE" or sample_specific_null_pull == "NO":
				#for col_index in na_cols:
				#	#if not(col_index == null_sample_col_index):	#don't take self?
				#	na_vals.append(float(data_elems[col_index]))

				na_mean = all_null_mean
				na_sd = all_null_sd

			else:
				#print "null_sample_col_index = "+str(null_sample_col_index)
				tcga_sample_id = null_col_to_sample_mapping[null_sample_col_index][10:]
				#print "null_sample_id = "+str(null_col_to_sample_mapping[null_sample_col_index])
				#print "tcga_sample_id = "+str(tcga_sample_id)
				sample_specific_na_cols = tcga_sample_to_na_mapping[tcga_sample_id]
				#print "sample_specific_na_cols = "+str(sample_specific_na_cols)
				for col_index in sample_specific_na_cols:
					na_vals.append(float(data_elems[col_index]))
			
				na_vals_array = array(na_vals)
				na_mean = na_vals_array.mean()
				#na_mean = sum(na_vals) / len(na_vals)
				#print "na_mean = "+str(na_mean)
				#na_sd = math.sqrt(sum([math.pow(x-na_mean, 2) for x in na_vals]) / (len(na_vals)-1))
				na_sd = na_vals_array.std(ddof=1)
				#print "na_sd = "+str(na_sd)

			if not(na_sd == 0):
				tcga_z_score = (tcga_val - na_mean) / na_sd
			else:
				tcga_z_score = 0
			#print "tcga_z_score = "+str(tcga_z_score)
			null_z_row.append(tcga_z_score)

			sample_index += 1
			
		#print the null row:
		print >> out_null, line_elems[0]+"\t"+"\t".join([str(x) for x in null_z_row])
	
	line_num += 1

in_paradigm.close()
out_tcga.close()
out_null.close()

print >> sys.stdout, time.time() - start_time, "seconds"
