#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""
import sys, os

try:
    sys.argv[1]
except:
    print('Invalid input file or input file missing')
    quit()
else:
    if sys.argv[1] == '':
        filename = sys.argv[1]
    else:
        filename = sys.argv[1]

##pass args filename
windowSize = 110

def createWindows(contig, windowSize):
    seqSize = len(contig[1])
    lastWinIdx = seqSize - windowSize + 1
    fout = open('windows.fa', 'w')
    for i in range(0, lastWinIdx):
        appendString = '_' + str(i+1)
        fout.write(contig[0] + appendString + '\n')
        fout.write(contig[1][i:i+windowSize] + '\n')
        #print i
    fout.close()
    return 0

def contigParser(filename):
    fin = open(filename, 'r')
    header = fin.readline().replace('\n', '')
    faSeq = ''
    line = fin.readline().replace('\n', '')
    while line != '':
        faSeq = faSeq + line
        line = fin.readline().replace('\n', '')
        if len(line) != 0 and line[0] == '>':
            print('Multiple headers in this file. Please split prior to running this script')
            quit()
    return [header, faSeq]

if __name__ == '__main__':
    #filename = '10_100000809-100001283.fa'
    contig = contigParser(filename)

    outdirExt = filename.split('.')
    outdir = filename.replace('.' + outdirExt[len(outdirExt)-1], '')
    try:
        os.stat(outdir)
    except:
        os.mkdir(outdir)

    os.chdir(outdir)
    createWindows(contig, windowSize)
