#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""

import gzip, argparse
#usage: processcovfiles.py -i filename.coverage
#notmfefiles = glob.glob('*.notMFEwindows.*.coverage.gz')
#mfefiles = glob.glob('*.MFEwindows.*.coverage.gz')

parser = argparse.ArgumentParser(description='Convert coverage files to count files.')
parser.add_argument('-i', '--input', help='Coverage file name', type=str, required=True)
parser.add_argument('-o', '--outdir', help='Output directory name', type=str, required=False)
args = parser.parse_args()
filename = args.input
outdir = args.outdir
if outdir == None:
	outdir = '.'
elif outdir[-1] == '/':
	outdir = outdir[0:-1]

def importfile(filename):
	mydict = dict()
	for line in gzip.open(filename):
		line = line.strip().split()
		if len(line) == 8:
			xsome, startpos, stoppos, featureid, qv, strand, mfeindex, sirnacounts = line
			sirnacounts = int(sirnacounts)
		elif len(line) == 9:
			xsome, startpos, stoppos, featureid, qv, strand, cigar, mfeindex, sirnacounts = line
			sirnacounts = int(sirnacounts)
		elif len(line) == 6:
			xsome, startpos, stoppos, featureid, mfeindex, sirnacounts = line
			sirnacounts = int(sirnacounts)
		else:
			print('Line does not contain 8 or 9 fields. Please correct file: ' + filename)
			print(line)
			quit()
		if xsome in mydict.keys():
			if featureid in mydict[xsome].keys():
				if (startpos, stoppos) in mydict[xsome][featureid].keys():
					mydict[xsome][featureid][(startpos, stoppos)][0] += 1
					mydict[xsome][featureid][(startpos, stoppos)][1] = mydict[xsome][featureid][(startpos, stoppos)][1] + sirnacounts
				else:
					mydict[xsome][featureid][(startpos, stoppos)] = [1, sirnacounts]
			else:
				mydict[xsome][featureid] = dict()
				mydict[xsome][featureid][(startpos, stoppos)] = [1, sirnacounts]
		else:
			mydict[xsome] = dict()
			mydict[xsome][featureid] = dict()
			mydict[xsome][featureid][(startpos, stoppos)] = [1, sirnacounts]
	return mydict

def calcsirnacountspernt(mydict):
	for xsome in mydict.keys():
		for featureid in mydict[xsome].keys():
			length = 0
			sirnacounts = 0
			for mytuple in mydict[xsome][featureid].keys():
				tlength, tsirnacounts = mydict[xsome][featureid][mytuple]
				length = length + tlength
				sirnacounts = sirnacounts + tsirnacounts
			mydict[xsome][featureid] = [length, sirnacounts, float(sirnacounts)/float(length)]
	return mydict

def writecsv(mydict, filename):
	fout = open(outdir + '/' + filename.split('/')[-1].replace('coverage.gz', 'counts'), 'w')
	fout.write('Chromosome' + '\t' + 'FeatureID' + '\t' + 'Length' + '\t' + 'siRNA_Counts' + '\t' + 'CountsPerNT' + '\n')
	for xsome in sorted(mydict.keys()):
		for featureid in sorted(mydict[xsome].keys()):
			length, sirnacounts, countspernt = mydict[xsome][featureid]
			fout.write(xsome + '\t' + featureid + '\t' + str(length) + '\t' + str(sirnacounts) + '\t' + str(countspernt) + '\n')
	fout.close()

mydict = importfile(filename)
mydict = calcsirnacountspernt(mydict)
writecsv(mydict, filename)
