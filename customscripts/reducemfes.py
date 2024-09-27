#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""
import os,random, glob, argparse
import multiprocessing as mp

parser = argparse.ArgumentParser(description='Create csv file of statistics of MFE calculations for windows.')
parser.add_argument('-r', '--random', help='Total number of windows.out files to extract statistics from. default value is 0', type=int, required=False)
parser.add_argument('-m', '--minMFE', help='The cutoff for deciding minMFE. Default = -40.0', type=float, required=False)
parser.add_argument('-w', '--window', help='Length of window size. default value is 110', type=int, required=False)
parser.add_argument('-f', '--feature', help='Feature type to be used as a prefix. Default takes the root of the directory. i.e. rRNA__test_0-100, rRNA is the feature type.', type=str, required=False)
parser.add_argument('-t', '--threads', help='The number of threads to run. Default = 1', type=int, required=False)
args = parser.parse_args()

myrandsubset = args.random
if myrandsubset == None:
    myrandsubset = 0
minMFE = args.minMFE
if minMFE == None:
    minMFE = -40.0
windowsize = args.window
if windowsize == None:
    windowsize = 110
featuretype = args.feature
threads = args.threads
if threads == None:
    threads = 1

def writeTSVfile(mylist, filename):
    outfile = filename + '.tsv'
    if os.path.exists(outfile):
        os.remove(outfile) #remove old files
    if len(mylist) > 0:
        fout = open(outfile, 'w')
        for i in range(0,len(mylist)):
            mystring = ''
            for j in range(0,len(mylist[i])):
                mystring = mystring + str(mylist[i][j]) + '\t'
            #print(mystring[0:-1] + '\n')
            fout.write(mystring[0:-1] + '\n')
        fout.close()

def returnRandomSet(filelist,number):
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

def calculateMinMFE(mfelist):
    return(min(mfelist, key=lambda item: item[1]))

def calculateVariance(mfelist,mean):
    mysum = 0.0
    variance = 0.0
    for mfe in mfelist:
        mysum = (mfe[1] - mean)*(mfe[1] - mean) + mysum
    variance = mysum/float(len(mfelist))
    return variance

def calculateMean(mfelist):
    mysum = 0.0
    mymean = 0.0
    for mfe in mfelist:
        mysum = mysum + mfe[1]
    mymean = mysum/float(len(mfelist))
    return mymean

def calculatePCTThreshMFE(mfelist,thresholdval):
    mfecount = 0
    mfepctcov = 0.0
    for mfe in mfelist:
        if mfe[1] <= thresholdval:
            mfecount += 1
    mfepctcov = mfecount/float(len(mfelist))
    return [mfepctcov,mfecount]

def returnWindowMFEList(file):
    windowLine = " "
    mymfelist = list()
    fin = open(file,'r')
    while windowLine != "":
        windowLine = fin.readline().strip().replace('>', '')
        header = windowLine
        if header != '':
            fin.readline()
            windowLine = fin.readline().strip()
            myfloatlist = windowLine.split(' ')
            myfloatlist = myfloatlist[len(myfloatlist) - 1].replace('(', '').replace(')', '')
            myfloat = float(myfloatlist)
            mymfelist.append([header, myfloat])
    return mymfelist

def returnMin(mfelist):
    tempList = min(mfelist, key=lambda item: item[1])  # returns a list [header, mfe]
    return tempList

def returnWindowMFEStats(myset):
    mymasterlist = list()
    for file in myset:
        windowline = " "
        mymfelist = returnWindowMFEList(file)
        writeTSVfile(mymfelist,file)
        #call plotter here
        #plotMFEs(mylist)
        if len(mymfelist) == 0:
            print('Warning! Error in file: ' + file + ' with length: ' + str(len(mymfelist)))
        else:
            tempList = min(mymfelist, key=lambda item: item[1]) #returns a list [header, mfe]
            mymean = calculateMean(mymfelist)
            tempList.append(mymean)
            tempList.append(calculateVariance(mymfelist,mymean))
            mfepctcov = calculatePCTThreshMFE(mymfelist,minMFE)
            #mymfepos = returmmfepos(mymfelist,minMFE)
            tempList.append(mfepctcov[1])
            tempList.append(mfepctcov[0])
            #templist = [headerid, minMFE, meanMFE, MFEvariance, MFEcountbelowthresh, %MFEbelowthresh]
            mymasterlist.append(tempList)
    return mymasterlist

def mpreturnWindowMFEStats(file):
    mymfelist = returnWindowMFEList(file)
    writeTSVfile(mymfelist, file)
    #call plotter here
    #plotMFEs(mylist)
    if len(mymfelist) == 0:
        print('Warning! Error in file: ' + file + 'with length: ' + str(len(mymfelist)))
    else:
        tempList = min(mymfelist, key=lambda item: item[1]) #returns a list [header, mfe]
        mymean = calculateMean(mymfelist)
        tempList.append(mymean)
        tempList.append(calculateVariance(mymfelist,mymean))
        mfepctcov = calculatePCTThreshMFE(mymfelist,minMFE)
        #mymfepos = returmmfepos(mymfelist,minMFE)
        tempList.append(mfepctcov[1])
        tempList.append(mfepctcov[0])
        #templist = [headerid, minMFE, meanMFE, MFEvariance, MFEcountbelowthresh, %MFEbelowthresh]
        return tempList

def writeCSVListH(mylist, filename, headers):
    fout = open(filename, 'w')
    fout.write(str(headers[0]))
    for j in range(1,len(headers)):
        fout.write(',' + str(headers[j]))
    fout.write('\n')
    for list in mylist:
        fout.write(str(list[0]))
        for i in range(1,len(list)):
            fout.write(',' + str(list[i]))
        fout.write('\n')
    fout.close()
    return 0

if __name__ == '__main__':
    mywindowsdirs = glob.glob('*/windows.out')
    myheaders = ['header', 'minMFE', 'meanMFE', 'MFEvariance', 'MFEcountltthresh', 'PctMFEltthresh']
    if (len(mywindowsdirs) > 0):
        if featuretype == None:
            featuretype = mywindowsdirs[0].split('/')[0].split('__')[0]
    else:
        print('No windows.out files exist in subdirectories. Please rerun RNAfold')
        quit(1)
    mycsvfile = 'my_' + featuretype + '_mfestats.csv'
    mymfefileset = returnRandomSet(mywindowsdirs, myrandsubset)
    if threads > 1:
        pool = mp.Pool(threads)
        mymfelistmp = pool.map(mpreturnWindowMFEStats, mymfefileset)
        mymfelist = filter(None,mymfelistmp)
    elif threads == 1:
        mymfelist = returnWindowMFEStats(mymfefileset)
    else:
        print('Invalid # of threads set. Please use a positive integer for threads')
        quit(1)
    writeCSVListH(mymfelist, mycsvfile, myheaders)
