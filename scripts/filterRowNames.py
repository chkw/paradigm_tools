# filterRowNames.py
# Restrict matrix data to include only specified values in column 1.
# APR2016	chrisw

# imports
from optparse import OptionParser
import sys
import os
from sets import Set

# global vars


# methods and functions

def getOptions():
	"parse options"
	parser = OptionParser(usage="%prog [options] [data matrix file names]", description="Restrict matrix data to include only specified values in column 1 of a row names file.")
	parser.add_option("-v", action="store_true", default=False, dest="verbose", help="Switch for verbose mode.")
	parser.add_option("-s", action="store", default="_filtered_rows", type="string", dest="fileNameSuffix", help="file name suffix for output")
	parser.add_option("-f", action="store", default=None, type="string", dest="rowNamesFile", help="file to get row names, one per line")

 	(options, args) = parser.parse_args()

	return (options, args)

def log(msg, die=False):
	if (verbose | die):
		sys.stderr.write(msg)
	if die:
		sys.exit(1)

#:####################################

def main():
	global verbose
	(options, args) = getOptions()
	verbose = options.verbose
	log('options:\t%s\n' % (str(options)))
	log('args:\t%s\n' % (str(args)))

	fileNameSuffix = options.fileNameSuffix
	rowNamesFile = options.rowNamesFile

	# get row names
	rowNames = Set()
	rowFileObj = open(rowNamesFile, "r")
	for line in rowFileObj:
		line = line.rstrip()
		fields = line.split("\t", 1)
		rowNames.add(fields[0])
	rowFileObj.close()
	log("rowNames:%s\n" % (str(len(rowNames))))

	# process files
	for filePath in args:
		(dirName, dataFileName) = os.path.split(filePath)
		log("processing %s\n" % (filePath))

		outFilePath = filePath + fileNameSuffix
		outFileObj = open(outFilePath, "w")

		inFileObj = open(filePath, "r")

		count = 0
		headerWritten = False
		for line in inFileObj:
			if not headerWritten:
				outFileObj.write(line)
				headerWritten = True
				continue
			fields = line.split("\t", 1)
			if fields[0] in rowNames:
				count += 1
				outFileObj.write(line)

		log("%s: %s rows\n" % (filePath, str(count)))
		outFileObj.close()
		inFileObj.close()

# main program section
if __name__ == "__main__":
	main()
