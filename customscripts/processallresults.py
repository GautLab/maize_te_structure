#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""

import gzip, argparse
#usage: processallresults.py -i filename.txt
#notmfefiles = glob.glob('*.notMFEwindows.*.coverage.gz')
#mfefiles = glob.glob('*.MFEwindows.*.coverage.gz')

parser = argparse.ArgumentParser(description='Convert binomial test results to csv.')
parser.add_argument('-i', '--input', help='Binomial test output filename', type=str, required=True)
args = parser.parse_args()
filename = args.input


with open(filename, 'rb') as fin:
    data = fin.read().split('\n')

with open(filename.replace('.txt','.csv'),'w') as fout:
    fout.write('Feature Type\tsiRNA Size\tBinomial Test(p-value)\tLibrary\tDescription\n')
    for i in range(0,len(data)):
        if len(data[i]) != 0:
            if i%3 == 0:
                line = data[i]
                feature, null, null, siRNAlib, null, null, null, null, null = line.split('.')
                null, null, lib, null, siRNAsize, null, null, null = siRNAlib.split('_')
            elif i%3 == 2:
                pvalue = float(data[i])
                fout.write(feature + '\t' + siRNAsize + '\t' + str(pvalue) + '\t' + lib + '\n')
