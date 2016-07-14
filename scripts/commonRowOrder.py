# commonRowOrder.py
# Restrict matrix data to include only common rows, outputed in the same order.
# JULY 2016	chrisw

# imports
from optparse import OptionParser
import sys

# global vars


# methods and functions

def getOptions():
	"parse options"
	parser = OptionParser(usage="%prog [options] [data matrix file names]", description="Restrict matrix data to include only common rows, outputed in the same order.")
	parser.add_option("-v", action="store_true", default=False, dest="verbose", help="Switch for verbose mode.")
	parser.add_option("--output_suffix", action="store", default="_common_rows_ordered", type="string", dest="fileNameSuffix", help="file name suffix for output")

 	(options, args) = parser.parse_args()

	return (options, args)

def log(msg, die=False):
	if (verbose | die):
		sys.stderr.write(msg)
	if die:
		sys.exit(1)

def getRowName(line, columnIndex=0, delimiter="\t"):
	rowName = None
	fields = line.split(delimiter)
	rowName = fields[columnIndex]
	return rowName

def outputOrderedRows(filePath, commonRowNames, outputFilePath, columnIndex=0, delimiter="\t"):
	'''
	assumes first line is a header lines
	'''
	headerRow = None

	# get lines keyed by row name
	fileObj = open(filePath, "r")
	linesDict = {}
	try:
		for line in fileObj.readlines():
			if headerRow == None:
				headerRow = line
			else:
				rowName = getRowName(line, columnIndex=columnIndex, delimiter=delimiter)
				linesDict[rowName] = line
	finally:
		fileObj.close()

	# output lines ordered by row names
	linesWritten = 0
	outputF = open(outputFilePath, "w")
	try:
		# write header
		outputF.write(headerRow)
		linesWritten += 1

		# write lines
		for rowName in commonRowNames:
			outputF.write(linesDict[rowName])
			linesWritten += 1
	finally:
		outputF.close()

	log("%s lines written to %s\n" % (str(linesWritten), outputFilePath))
	return None

#:####################################

def main():
	global verbose
	(options, args) = getOptions()
	verbose = options.verbose
	log('options:\t%s\n' % (str(options)))
	log('args:\t%s\n' % (str(args)))

	fileNameSuffix = options.fileNameSuffix

	# get row names
	rowNameSets = []
	for filePath in args:
		sys.stderr.write("getting row names from %s\n" % (filePath))

		rowNames = set()
		lineCount = 0
		fileObj = open(filePath, "r")
		try:
			for line in fileObj.readlines():
				if lineCount > 0:
					rowName = getRowName(line)
					rowNames.add(rowName)
				lineCount += 1
		finally:
			fileObj.close()

		rowNameSets.append(rowNames)

		log("%s has %s row names.\n" % ((filePath), str(len(rowNames))))

	# commmon row names
	commonRowNames = rowNameSets.pop()
	for rowNameSet in rowNameSets:
		commonRowNames = commonRowNames.intersection(rowNameSet)

	log("common row names: %s\n" % (str(len(commonRowNames))))

	commonRowNames = list(commonRowNames)

	for filePath in args:
		sys.stderr.write("ordering rows for %s\n" % (filePath))
		outputFilePath = filePath + options.fileNameSuffix
		outputOrderedRows(filePath, commonRowNames, outputFilePath)

# main program section
if __name__ == "__main__":
	main()
