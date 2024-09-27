#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""

import os, glob

myfastalist = glob.glob('G*.txt')
mysequencelist = list()

for myfastafile in myfastalist:
    mysequenceset = set()
    for line in open(myfastafile):
        mychar = line[0]
        myntset = {'A','C','T','G','N'}
        if mychar in myntset:
            line = line.strip()
            mylist = line.split()
            mysequenceset.add(mylist[0])
    mysequencelist.append(mysequenceset)

myfiledict = dict()
for line in open('joblist.txt'):
    filename = line.strip()
    myfiledict[filename] = dict()
    for length in open(filename):
        mylist = length.strip().split()
        if mylist[1] in myfiledict[filename].keys():
            myfiledict[filename][mylist[1]].append(mylist[0])
        else:
            myfiledict[filename][mylist[1]] = list()
            myfiledict[filename][mylist[1]].append(mylist[0])

mylengthslist = ['21', '22', '24']
for filename in myfiledict.keys():
    for sirnalen in mylengthslist:
        fout = open(filename.replace('sorted.lengths','siRNA.' + sirnalen + '.csv'),'w')
        fout.write('ReadID\n')
        for myread in myfiledict[filename][sirnalen]:
            fout.write(myread + '\n')
        fout.close()

def createFasta(myset,filename):
    counter = 1
    filename = filename.split('.')
    myfilename = ''
    for i in range(0,len(filename)-1):
        myfilename = myfilename + filename[i]
    fout = open(myfilename + '.fasta', 'w')
    for i in myset:
        myheader = '>' + i + '_' + str(counter)
        mySequence = i
        fout.write(myheader + '\n')
        fout.write(mySequence + '\n')
    fout.close()

def parsejoblistlenghts(filename):
    return 0
