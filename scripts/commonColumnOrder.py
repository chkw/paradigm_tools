# commonColumnOrder.py
# Restrict matrix data to include only common columns, outputed in the same order.
# APR2016	chrisw

# imports
from optparse import OptionParser
import sys
import csv
import os
from sets import Set
from collections import deque

# global vars


# methods and functions

def getOptions():
	"parse options"
	parser = OptionParser(usage="%prog [options] [data matrix file names]", description="Restrict matrix data to include only common columns, outputed in the same order.")
	parser.add_option("-v", action="store_true", default=False, dest="verbose", help="Switch for verbose mode.")
	parser.add_option("-s", action="store", default="_common_cols_ordered", type="string", dest="fileNameSuffix", help="file name suffix for output")

 	(options, args) = parser.parse_args()

	return (options, args)

def log(msg, die=False):
	if (verbose | die):
		sys.stderr.write(msg)
	if die:
		sys.exit(1)

def readFileLines(filename, strip=True):
	fileLines = []
	file = open(filename, 'r')
	for line in file.readlines():
		if strip:
			line = line.rstrip("\r\n")
		fileLines.append(line)
	file.close()
	return fileLines

def readTsv(fileLines, d="\t"):
	reader = csv.DictReader(fileLines, delimiter=d)
	return reader

#:####################################

def main():
	global verbose
	(options, args) = getOptions()
	verbose = options.verbose
	log('options:\t%s\n' % (str(options)))
	log('args:\t%s\n' % (str(args)))

	fileNameSuffix = options.fileNameSuffix

	# get column names
	colNameSets = []
	firstColDeque = deque()
	for filePath in args:
		(dirName, dataFileName) = os.path.split(filePath)
		sys.stderr.write("getting columns from %s\n" % (filePath))

		fileObj = open(filePath, "r")
		line = fileObj.readline()
		fileObj.close()

		line = line.rstrip()
		fields = line.split("\t")
		fieldsDeque = deque(fields)
		firstColDeque.append(fieldsDeque.popleft())
		fieldsSet = Set(fieldsDeque)

		colNameSets.append(fieldsSet)

		log("%s has %s columns.\n" % (dataFileName, str(len(fieldsSet))))

	# get common columns
	commonColSet = Set(colNameSets[0])
	for colSet in colNameSets:
		commonColSet = commonColSet.intersection(colSet)
	commonCols = list(commonColSet)

	log("%s columns in common\n" % (str(len(commonCols))))
	log("first column names: %s\n" % (str(firstColDeque)))

 	myDialect = csv.register_dialect("myDialect", csv.excel_tab, lineterminator="\n")

	# process files
	for filePath in args:
		(dirName, dataFileName) = os.path.split(filePath)
		sys.stderr.write("processing %s\n" % (filePath))

		fileObj = open(filePath, "r")
		lines = readFileLines(filePath)
		fileObj.close()

		reader = readTsv(lines)

		fileFieldNames = list(commonCols)
		fileFieldNames.append(firstColDeque.popleft())
		fileFieldNames.reverse()

		with open(filePath + fileNameSuffix, 'w') as csvfile:
			writer = csv.DictWriter(csvfile, delimiter="\t", lineterminator="\n", fieldnames=fileFieldNames, extrasaction="ignore")
			writer.writeheader()

			processedLines = 1

			for row in reader:
				writer.writerow(row)

				processedLines += 1
				if (processedLines % 1000) == 0:
					sys.stderr.write("%s lines\n" % (str(processedLines)))

# main program section
if __name__ == "__main__":
	main()
