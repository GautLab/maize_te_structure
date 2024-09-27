#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""
import argparse, glob
import multiprocessing as mp

parser = argparse.ArgumentParser(description='Create"s a bed file from MFE regions found in a list of windows.out.csv')
parser.add_argument('-i', '--input', help='A file containing a list of windows.out.tsv files. Default = all directories', type=str, required=False)
parser.add_argument('-m', '--minMFE', help='The cutoff for deciding minMFE. Default = -40.0', type=float, required=False)
parser.add_argument('-w', '--window', help='The size of each window. Default = 110', type=int, required=False)
parser.add_argument('-s', '--stepsize', help='Size of steps between windows. Default = 1', type=int, required=False)
parser.add_argument('-o', '--output', help='Name of bed file to output to. Default: Prefix of feature type', type=str, required=False)
parser.add_argument('-t', '--threads', help='The number of threads to run. Default = 1', type=int, required=False)
#parser.add_argument('--multiple', help='Use this option if multiple positions exist in headers', type=bool, required=False, action='store_true')
args = parser.parse_args()

windowoutfile = args.input
if windowoutfile == None:
    windowsoutfilelist = glob.glob('*/windows.out.tsv')
else:
    fin = open(windowoutfile)
    windowsoutfilelist = fin.readlines()
    fin.close()
    for i in range(0,len(windowsoutfilelist)):
        windowsoutfilelist[i] = windowsoutfilelist[i].strip()
minMFE = args.minMFE
if minMFE == None:
    minMFE = -40.0
windowsize = args.window
if windowsize == None:
    windowsize = 110
stepsize = args.stepsize
if stepsize == None:
    stepsize = 1
outputfile = args.output
if outputfile == None:
    if windowoutfile == None:
        outputfile = windowsoutfilelist[0].split('__')[0] + '_MFE.bed'
threads = args.threads
if threads == None:
    threads = 1

def returnWindowMFEList(file):
    mymfelist = list()
    fin = open(file,'r')
    windowlines = fin.readlines()
    fin.close()
    if len(windowlines) > 0:
        for i in range(0,len(windowlines)):
            windowlines[i] = windowlines[i].strip().split('\t')
        for i in range(0, len(windowlines)):
            if len(windowlines[i]) < 2:
                print('Error in input window, line # ' + str(i) + '. Window: ' + windowlines[i][0])
                quit(1)
            else:
                header = windowlines[i][0]
                mywindow = int(header.split('_')[-1])
                mymfe = float(windowlines[i][1])
                mymfelist.append([mywindow, header, mymfe])
    return sorted(mymfelist,key=lambda l:l[0])

#no mfe
#./lnc_RNA__7_182207619-182208001/windows.out
#has mfe
#./lnc_RNA__10_142389865-142390471/windows.out

#fix this
#this is for single core processing
def processgenomefeatures(mywindowslist):
    filename = ''
    myposmasterlist = list()
    for line in mywindowslist:
        #print line
        filename = line.strip()
        #print filename
        mytemplist = returnmfepos(returnWindowMFEList(filename), minMFE)
        if mytemplist != 1:
            myposmasterlist.append(mytemplist)
    prefix = filename[2:len(filename)].split('__')[0].split('/')[-1]
    writebedfile(myposmasterlist, outputfile)

#multicore processing
def mpprocessgenomefeatures(filename):
    filename = filename.strip()
    mytemplist = list()
    try:
        mytemplist = returnmfepos(returnWindowMFEList(filename), minMFE)
        return mytemplist
    except:
        print(filename)

def returnmfepos(mfelist, minmfe):
    mfepos = list()
    try:
        header = mfelist[0][1][0:len(mfelist[0][1])-2]
    except:
        print(mfelist[0][0])
    try:
        initialpos = int(header.split('::')[1].split(':')[1].split('-')[0])
    except:
        print(header)
    minmfeexists = False
    mfeposlist = list()
    for i in range(0,len(mfelist)):
        if mfelist[i][2] <= minmfe:
            mfepos.append(i+1+initialpos)
            minmfeexists = True
    if minmfeexists:
        #print 'minMFE found!'
        mfeposlist = returnmfewindows(mfepos)
    else:
        #print "No minMFE found for minMFE=" + str(minmfe)
        return ''
    #return [header, mfeposlist, mfepos]
    return [header, mfeposlist]

def checkcontinuity(mfeposlist):
    mycontmfeposlist = list()
    mfestartpos = mfeposlist[0]
    mydtype = type(mfestartpos)
    #print(mydtype)
    if mydtype == int:
        for i in range(1, len(mfeposlist)):
            mydiff = (mfeposlist[i-1] - mfeposlist[i])
            if mydiff != -1 or (i == len(mfeposlist)-1 and mydiff == -1):
                if i < len(mfeposlist)-1:
                    mfeendpos = (mfeposlist[i - 1]) #here we return the last index of continuity
                else:
                    mfeendpos = (mfeposlist[i])
                mycontmfeposlist.append([mfestartpos, mfeendpos + windowsize-1])
                mfestartpos = mfeposlist[i]
                #print(i)
    elif mydtype == list:
        return mfeposlist
    else:
        print('Illegal datatype in position lists. Type: ' + str(mydtype))
        print(str(mfeposlist))
        quit(1)
    return mycontmfeposlist

def returnmfewindows(mfeposlist):
    mycontmfeposlist = checkcontinuity(mfeposlist)
    overlappedmfewindows = list()
    myovlflagexists = False
    for i in range(1, len(mycontmfeposlist)):
        pos = mycontmfeposlist[i]
        prevpos = mycontmfeposlist[i-1]
        maxstartpos = max(prevpos[0], pos[0])
        minendpos = min(prevpos[1], pos[1])
        if(maxstartpos - minendpos) <= 0:
            overlappedmfewindows.append([prevpos[0], pos[1]])
            myovlflagexists = True
    if myovlflagexists:
        #print(myovlflagexists)
        #print(mycontmfeposlist)
        #print(overlappedmfewindows)
        return returnmfewindows(overlappedmfewindows)
    else:
        #print(myovlflagexists)
        #print(mycontmfeposlist)
        #print(overlappedmfewindows)
        return mycontmfeposlist

def writebedfile(mymfewindows,outputfile):
    fout = open(outputfile,'w')
    for i in mymfewindows:
        if len(i) != 1:
            header = i[0]
            myposlist = i[1]
            myxsome = header.split('::')[1].split(':')[0]
            for j in myposlist:
                startpos = j[0]
                stoppos = j[1]
                fout.write(myxsome + '\t' + str(startpos) + '\t' + str(stoppos) + '\t' + header + '\n')
    fout.close()

if __name__ == '__main__':
    #threads = 8
    #print(threads)
    #print(len(windowsoutfilelist)) #prints total # of files
    #print(windowsoutfilelist[0]) # prints first file
    #print(windowsoutfilelist[-1]) # prints last file
    #fout = open(outputfile,'w')
    #fout.write('')
    #fout.close()
    if threads > 1:
        pool = mp.Pool(threads)
        myposmasterlist = pool.map(mpprocessgenomefeatures, windowsoutfilelist)
        myposmasterlist = filter(None, myposmasterlist)
        writebedfile(myposmasterlist, outputfile)
    elif threads == 1:
        processgenomefeatures(windowsoutfilelist)
    else:
        print('Illegal number of threads!. Please choose a positive integer value of 1 or greater.')
        quit(1)
