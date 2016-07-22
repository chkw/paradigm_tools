#Yulia Newton
#python2.7 combine_matrices.py --in_dir temp/ --out_file temp_out.tab --intersect TRUE

import optparse, sys, os, glob
from numpy import *
from scipy.special import erfc

parser = optparse.OptionParser()
parser.add_option("--in_dir", dest="in_dir", action="store", default="", help="")
parser.add_option("--out_file", dest="out_file", action="store", default="", help="")
parser.add_option("--intersect", dest="intersect", action="store", default="TRUE", help="")
opts, args = parser.parse_args()

#process input arguments:
in_dir = opts.in_dir
if not(in_dir[len(in_dir)-1] == "/"):
	in_dir = in_dir + "/"
out_file = opts.out_file
intersect_flag = opts.intersect.upper()
if intersect_flag == "TRUE" or intersect_flag == "YES":
	intersect = True
else:
	intersect = False

#iterate over all files:
all_sigs = {}
all_headers = {}
#for file_name in glob.glob(os.path.join(in_dir, "*.tab")):
for file_name in glob.glob(os.path.join(in_dir, "*")):
	input_file_name_parts = file_name.split("/")
	input_file_name = input_file_name_parts[len(input_file_name_parts)-1]
	all_sigs[input_file_name] = {}
	
	input_file = open(file_name, 'r')
	line_count = 0
	line_size = 0
	for line in input_file:
		line_elems = line.replace("\n","").split("\t")
		if line_count == 0:
			all_headers[input_file_name] = "\t".join(line_elems[1:])
			line_size = len(line_elems[1:])
		else:
			file_line = "\t".join(line_elems[1:])
			file_line = file_line.replace("\t\t","\tNA\t")
			file_line_elems = line_elems[1:]
			if len(file_line_elems) < line_size:
				len_diff = line_size - len(file_line_elems)
				for ld in range(len_diff):
					file_line = file_line + "\tNA"
			all_sigs[input_file_name][line_elems[0]] = file_line
			
		line_count += 1
	input_file.close()

'''print "all_headers:"
print all_headers
print ""'''

all_sig_names = all_sigs.keys()

'''print "all_sig_names:"
print all_sig_names
print ""

print "all_sigs[all_sig_names[0]].keys():"
print all_sigs[all_sig_names[0]].keys()'''

genes = set()
for sig in all_sig_names:
	#print sig+":"
	#print set(all_sigs[sig].keys())
	genes.update(set(all_sigs[sig].keys()))

#genes = all_sigs[all_sig_names[0]].keys()

output_file = open(out_file, "w")
new_header = ""
sig_num = 1
for sig in all_sig_names:
	if sig_num > 1:
		new_header = new_header + "\t" + all_headers[sig]
	else:
		new_header = new_header + all_headers[sig]
	sig_num += 1
print >> output_file, "id\t"+new_header
for gene in genes:
	output_line = gene

	if intersect:
		present_in_all_sigs = True
		for sig in all_sig_names:
			if gene in all_sigs[sig]:
				output_line = output_line + "\t" + all_sigs[sig][gene]
			else:
				h = all_headers[sig].split("\t")
				for i in range(len(h)):			
					output_line = output_line + "\t"
				present_in_all_sigs = False

		if present_in_all_sigs:
			print >> output_file, output_line

	else:	#union
		for sig in all_sig_names:
			if gene in all_sigs[sig]:
				output_line = output_line + "\t" + all_sigs[sig][gene]
			else:
				h = all_headers[sig].split("\t")
				for i in range(len(h)):
					output_line = output_line + "\tNA"

		print >> output_file, output_line	

output_file.close()
