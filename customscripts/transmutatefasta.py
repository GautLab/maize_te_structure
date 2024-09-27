#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""
import glob, os, random, argparse
import multiprocessing as mp

parser = argparse.ArgumentParser(description='Create csv file of statistics of MFE calculations for windows.')
parser.add_argument('-p', '--permutations', help='Total number of random permuations to create. default value is 5', type=int, required=False)
parser.add_argument('-r', '--randomizationtype', help='Type of randomization to run; 0 sequence shuffle, 1 preserve gc content, 2 preserve only length. default value is 0', type=int, required=False)
parser.add_argument('-n', '--numberofrandom', help='Total n number random samples, for subsampling. default value is 0', type=int, required=False)
parser.add_argument('-t', '--threads', help='The number of threads to run. Default = 1', type=int, required=False)
parser.add_argument('-f', '--fileextension', help='The file extension of the fasta files. Default = fa.', type=str, required=False)
parser.add_argument('-o', '--outdir', help='The prefix of the directory to output randomized fasta files to. Default = randomfas', type=str, required=False)
args = parser.parse_args()

permutations = args.permutations
if permutations == None:
    permutations = 5
randtype = args.randomizationtype
if randtype == None:
    randtype = 0
numberofrandom = args.numberofrandom
if numberofrandom == None:
    numberofrandom = 0
threads = args.threads
if threads == None:
    threads = 1
fileextension = args.fileextension
if fileextension == None:
    fileextension = 'fa'
myoutdir = args.outdir
if myoutdir == None:
    myoutdir = 'randomfas'

def returnRandomSet(filelist, number):
    windowfiles = list()
    randomset = set()
    for windowfile in filelist:
        windowfile = windowfile.strip()
        windowfiles.append(windowfile)
    if number == 0:
        return set(windowfiles)
    else:
        while len(randomset) != number:
            randomset.add(random.choice(windowfiles))
    return randomset

def randomshufflestring(mystr):
    mynewstr = ''
    mynewstrlist = random.sample(mystr, len(mystr))
    for i in mynewstrlist:
        mynewstr = mynewstr + i
    return mynewstr

def randomshufflestringbygc(mystr):
    gccontent = calculategccontent(mystr)
    [numas, numcs, numts, numgs] = countbasecompo(mystr)
    strlen = len(mystr)
    numofgcs = numgs + numcs
    numofats = strlen - numofgcs
    newgs = random.randint(1, numofgcs)
    newcs = numofgcs - newgs
    newas = random.randint(1, numofats)
    newts = numofats - newas
    mynewstr = ''
    mynewstrlist = random.sample(mystr, strlen)
    for i in mynewstrlist:
        mynewstr = mynewstr + i
    return mynewstr

###to test
def countbasecompo(mystr):
    numas = 0
    numcs = 0
    numts = 0
    numgs = 0
    for i in mystr:
        if i == 'A':
            numas += 1
        elif i == 'C':
            numcs += 1
        elif i == 'T':
            numts += 1
        elif i == 'G':
            numgs += 1
    return [numas, numcs, numts, numgs]

def calculategccontent(mystr):
    [numas, numcs, numts, numgs] = countbasecompo(mystr)
    gccontent = float(numcs + numgs)/float(numas + numts + numgs + numcs)
    return gccontent

def testtwostrings(oristr, newstr):
    return countbasecompo(oristr) == countbasecompo(newstr)

def processrandomfastas(myrandomsubset, finaloutdir):
    for myfastafile in myrandomsubset:  # this can be mp
        fin = open(myfastafile)
        header = fin.readline().strip()
        sequence = fin.readline().strip()
        randsequence = randomshufflestring(sequence)
        if not testtwostrings(sequence, randsequence):
            print('Error: Sequences not equal. Most likely due to illegal characters in fasta file: ' + myfastafile)
            quit()
        else:
            fout = open(finaloutdir + '/' + myfastafile.replace('.' + fileextension, '_rand.' + fileextension), 'w')
            fout.write(header + '\n')
            fout.write(randsequence + '\n')
            fout.close()

def mpprocessrandomfastas(job_args):
    myfastafile = job_args[0]
    finaloutdir = job_args[1]
    fin = open(myfastafile)
    header = fin.readline().strip()
    sequence = fin.readline().strip()
    randsequence = randomshufflestring(sequence)
    if not testtwostrings(sequence, randsequence):
        print('Error: Sequences not equal. Most likely due to illegal characters in fasta file: ' + myfastafile)
        quit()
    else:
        fout = open(finaloutdir + '/' + myfastafile.replace('.' + fileextension, '_rand.' + fileextension), 'w')
        fout.write(header + '\n')
        fout.write(randsequence + '\n')
        fout.close()

if __name__ == '__main__':
    myfastafiles = glob.glob('*.' + fileextension)
    myrandomsubset = returnRandomSet(myfastafiles, numberofrandom)
    if len(myrandomsubset) > 0:
        for permutation in range(1, permutations + 1):
            finaloutdir = myoutdir + str(permutation)
            if not os.path.exists('./' + finaloutdir):
                os.mkdir(finaloutdir)
            if threads == 1:
                processrandomfastas(myrandomsubset, finaloutdir)
            elif threads > 1:
                pool = mp.Pool(threads)
                finaloutdirs = [finaloutdir] * len(myrandomsubset)
                job_args = [(item_myrandomsubset, finaloutdirs[j]) for j, item_myrandomsubset in enumerate(myrandomsubset)]
                mylist = pool.map(mpprocessrandomfastas, job_args)
            else:
                print('Invalid # of threads set. Please use a positive integer for threads')
                quit(1)
    else:
        print('No fasta files found!')
        quit(1)
