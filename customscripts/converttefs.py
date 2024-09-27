#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""
import glob, os

'''
directory structure
./*.bed
./contigs/1_283345471-283399472
./contigs/tes/1_283345471-283399472.fasta
ignore when no _ split('.')[2] if find('_') != -1, then parse
from bed B73v4.TE.filtered_helitron.bed
mkdir helitron
1       27211   30897   helitron::1:27210-30897 42      +       3687M
folder = split()[0]_[1]-[2]
mv folder helitron/[3].replace(':','_')
'''

###required for running following subroutines
mybedfiles = glob.glob('*.bed')
myfeaturetypes = list()
mycleanbedfiles = list()
for i in mybedfiles:
    mylongstr = i.split('.')
    mystr = mylongstr[2]
    if mystr.find('_') != -1:
        myfeaturetypes.append(mystr.replace('filtered_',''))
        mycleanbedfiles.append(mylongstr[0] + '.' + mylongstr[1] + '.' + mylongstr[2] + '.bed')

###rename directories to standard naming convention
counter = 0
for bedfile in mycleanbedfiles:
    myfeaturetypedir = './contigs/' + myfeaturetypes[counter]
    if not os.path.exists(myfeaturetypedir):
        os.mkdir(myfeaturetypedir)
    for myline in open(bedfile):
        myline = myline.strip().split()
        myfeaturedir = myline[0] + '_' + myline[1] + '-' + myline[2]
        mynewdir = myline[3].replace(':','_')
        if os.path.exists('./contigs/' + myfeaturedir):
            os.rename('./contigs/' + myfeaturedir, './contigs/' + myfeaturetypes[counter] + '/' + mynewdir)
    counter += 1

###rename original fa files to proper naming convention
counter = 0
for bedfile in mycleanbedfiles:
    myfeaturetypedir = './contigs/tes/' + myfeaturetypes[counter]
    if not os.path.exists(myfeaturetypedir):
        os.mkdir(myfeaturetypedir)
    for myline in open(bedfile):
        myline = myline.strip().split()
        myfeaturefa = myline[0] + '_' + myline[1] + '-' + myline[2] + '.fa'
        mynewfa = myline[3].replace(':','_') + '.fa'
        if os.path.exists('./contigs/tes/' + myfeaturefa):
            os.rename('./contigs/tes/' + myfeaturefa, './contigs/tes/' + myfeaturetypes[counter] + '/' + mynewfa)
    counter += 1

#rename headers in all window output files
counter = 0
for bedfile in mycleanbedfiles:
    myfeaturetypedir = './contigs/' + myfeaturetypes[counter]
    #mywindowdirs = glob.glob(myfeaturetypedir + '/*/')
    for myline in open(bedfile):
        myline = myline.strip().split()
        myfeaturedir = myline[3].replace(':', '_')
        myworkingdir = './contigs/' + myfeaturetypes[0] + '/' + myfeaturedir
        if os.path.exists(myworkingdir):
            ###parse windows.fa
            fin = open(myworkingdir + '/windows.fa')
            myfile = fin.readlines()
            fin.close()
            listlen = len(myfile)
            for listindex in range(0,listlen):
                #replace >4:13818177-13818315_24 with >miRNA::4:13818177-13818315_24
                if listindex%2 == 0:
                    myfile[listindex] = '>' +  myfeaturetypes[counter] + '::' + myfile[listindex][1:-1] + '\n'
            fout = open(myworkingdir + '/windows.fa', 'w')
            fout.writelines(myfile)
            fout.close()
            ###parse windows.out
            fin = open(myworkingdir + '/windows.out')
            myfile = fin.readlines()
            fin.close()
            listlen = len(myfile)
            for listindex in range(0, listlen):
                # replace >4:13818177-13818315_24 with >miRNA::4:13818177-13818315_24
                if listindex % 3 == 0:
                    myfile[listindex] = '>' + myfeaturetypes[counter] + '::' + myfile[listindex][1:-1] + '\n'
            fout = open(myworkingdir + '/windows.out', 'w')
            fout.writelines(myfile)
            fout.close()
            ###parse windows.out.tsv
            fin = open(myworkingdir + '/windows.out.tsv')
            myfile = fin.readlines()
            fin.close()
            listlen = len(myfile)
            for listindex in range(0, listlen):
                # replace >4:13818177-13818315_24 with >miRNA::4:13818177-13818315_24
                myfile[listindex] = myfeaturetypes[counter] + '::' + myfile[listindex][0:-1] + '\n'
            fout = open(myworkingdir + '/windows.out.tsv', 'w')
            fout.writelines(myfile)
            fout.close()
    counter += 1
