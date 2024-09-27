#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""

####find nonmfes = mfes -
import sys,os,glob
import multiprocessing as mp

try:
    sys.argv[1]
except:
    print('Invalid input file or input file missing')
    quit()
else:
    if sys.argv[1] == '':
        mfecountfile = sys.argv[1]
    else:
        mfecountfile = sys.argv[1]

#countfiles = glob.glob('*.counts')
#nonmfecountfiles = glob.glob('*.not*.counts')
#mfecountfiles = set(countfiles) - set(nonmfecountfiles)
mybeddict = dict()

def writemfecountstocsv(mymfedict, filename):
    fout = open(filename.replace('.counts','.hits'),'w')
    print('Writing to: ' + filename.replace('.counts','.hits'))
    fout.write('Chromosome' + '\t' + 'SequenceID' + '\t' + 'siRNA_hits' + '\t' + 'MFE Length' + '\t' + 'Sequence_Length' + '\n')
    for xsome in mymfedict:
        for seqid in mymfedict[xsome]:
            mysum = 0
            for i in mymfedict[xsome][seqid]:
                mysum = mysum + i
            pos = seqid.split(':')[3].split('-')
            length = int(pos[1])-int(pos[0])
            mfelen = 0
            for mfelength in mybeddict[xsome][seqid]:
                mfelen = mfelen + mfelength
            fout.write(xsome + '\t' + seqid + '\t' + str(mysum) + '\t' + str(mfelen) + '\t' + str(length) + '\n')
    fout.close()

#single threaded
#bedfiles = glob.glob('*.bed')
#notbedfiles = glob.glob('*.not.bed')
#sortedbedfiles = glob.glob('*.sorted.bed')

def returnBedDict(bedfiles):
    xsomeset = set()
    mybeddict = dict()
    for bedfile in bedfiles:
        if os.stat(bedfile).st_size > 0:
            for line in open(bedfile):
                line = line.strip().split()
                # print '1st loop: ' + bedfile
                xsomeset.add(line[0])
    for xsome in xsomeset:
        mybeddict[xsome] = dict()
    for bedfile in bedfiles:
        if os.stat(bedfile).st_size > 0:
            for line in open(bedfile):
                line = line.strip().split()
                # print '2nd loop: ' + bedfile
                xsome = line[0]
                seqid = line[3]
                startpos = int(line[1])
                stoppos = int(line[2])
                if seqid in mybeddict[xsome].keys():
                    mybeddict[xsome][seqid].append(stoppos - startpos)
                else:
                    mybeddict[xsome][seqid] = [stoppos - startpos]
    return mybeddict

    # mybeddict = returnBedDict(sortedbedfiles)

#single thread
def parseCountFiles(mfecountfiles):
    for mfefile in mfecountfiles:
        mymfedict = dict()
        mymfexsomeset = set()
        # mynonmfeline[2] + '::' + mynonmfeline[0] + ':' + mynonmfeline[3] + ':' + mynonmfeline[4]
        for mfeline in open(mfefile):
            mymfexsomeset.add(mfeline.strip().split()[0])
        for xsome in mymfexsomeset:
            mymfedict[xsome] = dict()
        for mfeline in open(mfefile):
            mymfelist = mfeline.strip().split()
            xsome = mymfelist[0]
            seqid = mymfelist[3]
            count = int(mymfelist[-1])
            if seqid in mymfedict[xsome].keys():
                mymfedict[xsome][seqid].append(count)
            else:
                mymfedict[xsome][seqid] = [count]
        writemfecountstocsv(mymfedict,mfefile)

#parseCountFiles(mfecountfiles)
#multithreaded
def mpParseCountFiles(mfefile):
    mymfedict = dict()
    mymfexsomeset = set()
    # mynonmfeline[2] + '::' + mynonmfeline[0] + ':' + mynonmfeline[3] + ':' + mynonmfeline[4]
    if os.stat(mfefile).st_size > 0:
        for mfeline in open(mfefile):
            #print mfeline
            mymfexsomeset.add(mfeline.strip().split()[0])
        for xsome in mymfexsomeset:
            mymfedict[xsome] = dict()
        for mfeline in open(mfefile):
            mymfelist = mfeline.strip().split()
            #print mymfelist
            xsome = mymfelist[0]
            seqid = mymfelist[3]
            count = int(mymfelist[5])
            if seqid in mymfedict[xsome].keys():
                mymfedict[xsome][seqid].append(count)
            else:
                mymfedict[xsome][seqid] = [count]
        #print str(len(mymfedict))
        writemfecountstocsv(mymfedict,mfefile)

if __name__ == '__main__':
    mfecountfiles = glob.glob('*.counts')
    sortedbedfiles = glob.glob('*.bed')
    mybeddict = returnBedDict(sortedbedfiles)
    print len(mfecountfiles)
    print len(mybeddict)
    print mfecountfile
    mpParseCountFiles(mfecountfile)
    #threads = 128
    #pool = mp.Pool(threads)
    #pool.map(mpParseCountFiles, mfecountfiles)
