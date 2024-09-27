#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 00:06:48 2014

@author: Edwin
SEE LINE 49 for specific fasta filenames to be written
"""
import sys, os
#os.chdir("/home/edwin/Downloads/new pacb asm/")
#inputfile = "C:\\Users\\Edwin\\Dropbox\\research\\ranz\\New folder\\vectorbase\\vectorbaseex.fa"
#inputfile = "/home/edwin/Downloads/new pacb asm/all_tigs.fa"

try:
    sys.argv[1]
except:
    print('Invalid input file or input file missing')
    quit()
else:
    if sys.argv[1] == '':
        inputfile = sys.argv[1]
    else:
        inputfile = sys.argv[1]

try:
    sys.argv[2]
except:
    print('Invalid delimeter or delimeter missing')
    quit()
else:
    if sys.argv[2] == '':
        delimeter = sys.argv[2]
    else:
	delimeter = sys.argv[2]
        
outdir = "contigs"

try:
    os.stat(outdir)
except:
    os.mkdir(outdir)
        
fin = open(inputfile)
foutjob = open('joblist.txt','w')
#delimeter = '|'
totallines = sum(1 for line in fin)
fin = open(inputfile)
print 'opened input file'
seq_read = ''
seq_name = ''
last_pos = 0
count = 0
#print('begin for loop')
while count < totallines:
    #print('for loop executed')
    last_pos = fin.tell()
    line = fin.readline()
    count = count + 1
    if line[0:1] == '>':
        #print 'found seq_name ' + line
        seq_name = line.replace('\n','')
        ###change replace statements here depending on fasta file
        outputfile = seq_name.split(delimeter)[0].replace('>','').replace(' ','_').replace(':','_').replace('/','_') + '.fasta'
        last_pos = fin.tell()
        line = fin.readline().replace('\n','')
        count = count + 1
        #print 'begin while loop on seq ' + line
        while line[0:1] != '>' and line[0:1] != '':
            seq_read = seq_read + line
            last_pos = fin.tell()
            line = fin.readline().replace('\n','')
            count = count + 1
            #print 'while loop next at last pos ' + str(last_pos)
        #print 'while loop ended'
        if line[0:1] == '>' or line[0:1] == '':
            print 'begin writing file ' + outputfile
            fout = open(outdir + '/' + outputfile, 'w')
            #fout = open('contigs//' + outputfile, 'w')
            fout.write(str(seq_name) + '\n')
            fout.write(str(seq_read) + '\n')
            fout.close
            #print outputfile
            foutjob.write('./contigs/' + outputfile + '\n')
            #print len(seq_read.replace("\n",""))
            seq_name = ''
            seq_read = ''
            count = count - 1
            fin.seek(last_pos)
            #print 'file written'
fin.close()
fout.close()
foutjob.close()
