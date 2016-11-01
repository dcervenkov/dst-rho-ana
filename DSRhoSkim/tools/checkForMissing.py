#!/usr/bin/env python

import sys
import os
import fnmatch
import re

def stringSplitByNumbers(x):
	r = re.compile('(\d+)')
	l = r.split(x)
	return [int(y) if y.isdigit() else y for y in l]

if len(sys.argv) < 2:
	print "ERROR: Missing argument(s)"
	print "USAGE: " + sys.argv[0] + " DIR(S)..."
	sys.exit(1)

for dir in sys.argv[1:]:
	print dir + ":"
	errors = 0
	indexFiles = fnmatch.filter(os.listdir(dir), "*.index")
	indexFiles = sorted(indexFiles, key = stringSplitByNumbers)
	for stream in range(20):
		indexFilesSingleStream = fnmatch.filter(indexFiles, "*s" + str(stream) + ".*")
		for exp in range(70):
			indexFilesSingleStreamSingleExp = fnmatch.filter(indexFilesSingleStream, "*exp" + str(exp) + "_*")
			prevRunStart = 0
			prevRunEnd = 0
			for indexFile in indexFilesSingleStreamSingleExp:
				prefix, exp, runs, stream = indexFile.split("_")
				runStart, runEnd = runs.split("-")
				runStart = int(runStart.lstrip('r'))
				runEnd = int(runEnd.lstrip('r'))
				if runStart != (prevRunEnd + 1) and prevRunEnd != 0:
					if runStart < (prevRunEnd + 1):
						print "ERROR  ",
					else:
						print "MISSING",
					print ": runStart = {:d} \tprevRunEnd = {:d} \tfile: {:s}".format(runStart, prevRunEnd, indexFile)
					errors += 1
				prevRunStart = runStart
				prevRunEnd = runEnd
	if errors:
		print str(errors) + " problems"
	else:
		print "All OK"
	
