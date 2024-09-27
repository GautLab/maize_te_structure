#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""

import gzip, argparse
#usage: processmfebeds.py -i annotation.bed -m mfe.bed

parser = argparse.ArgumentParser(description='Creates nonMFE bed files from MFE and annotation bed files.')
parser.add_argument('-i', '--input', help='Annotation bed file name', type=str, required=True)
parser.add_argument('-m', '--mfe', help='MFE region bed file name', type=str, required=True)
parser.add_argument('-o', '--outdir', help='Output directory name', type=str, required=False)
args = parser.parse_args()
filename = args.input
mfefilename = args.mfe
outdir = args.outdir
if outdir == None:
    outdir = '.'
elif outdir[-1] == '/':
    outdir = outdir[0:-1]

def readannotation(filename):
    mydict = dict()
    myxsomeset = set()
    fin = open(filename)
    mylines = fin.readlines()
    for i in range(0, len(mylines)):
        mylines[i] = mylines[i].strip().split()
    mylines = sorted(mylines, key = lambda x: (x[0], x[3], x[1], x[2]))
    for i in range(0, len(mylines)):
        myxsomeset.add(mylines[i][0])
    for mykey in myxsomeset:
        mydict[mykey] = dict()
    for i in range(0, len(mylines)):
        mykey = mylines[i][0]
        myfeatureid = mylines[i][3]
        if myfeatureid in mydict[mykey].keys():
            mydict[mykey][myfeatureid].append([int(mylines[i][1]), int(mylines[i][2])])
        else:
            mydict[mykey][myfeatureid] = [[int(mylines[i][1]), int(mylines[i][2])]]
        #mydict[mylines[i][0]].append([int(mylines[i][1]), int(mylines[i][2]), mylines[i][3]])
    for mykey in mydict.keys():
        for featureid in mydict[mykey].keys():
            mydict[mykey][featureid] = sorted(mydict[mykey][featureid], key = lambda x: (x[0], x[1]))
    return mydict

def mergeoverlaps(featureregions):
    newregions = list()
    if len(featureregions) > 1:
        firststartpos = featureregions[0][0]
        firststoppos = featureregions[0][1]
        startpos = featureregions[0][0]
        stoppos = featureregions[0][1]
        laststartpos = featureregions[-1][0]
        laststoppos = featureregions[-1][1]
        for i in range(0, len(featureregions)-1):
            nextstartpos = featureregions[i+1][0]
            nextstoppos = featureregions[i+1][1]
            mydiff = nextstartpos - stoppos
            if mydiff < 0:
                #print('execute special condition mydiff < 0 @ ' + str(i))
                #print(str(i) + ':' + str(startpos) + ':' + str(stoppos))
                stoppos = nextstoppos
                #print(str(i) + ':' + str(startpos) + ':' + str(stoppos))
            else:
                #print('execute special condition mydiff >= 0 @ ' + str(i))
                #print(str(i) + ':' + str(startpos) + ':' + str(stoppos))
                newregions.append([startpos, stoppos])
                startpos = nextstartpos
                stoppos = nextstoppos
                #print(str(i) + ':' + str(startpos) + ':' + str(stoppos))
        if startpos - firststartpos == 0:
            newregions.append([startpos, stoppos])
        if len(newregions) == 0:
            # check if all regions overlap and append them here
            #print('execute no length in newregions')
            #print(str(i) + ':' + str(startpos) + ':' + str(stoppos))
            newregions.append([startpos, stoppos])
        if laststartpos - newregions[-1][1] > 0:
            if startpos == laststartpos:
                #print('execute laststartpos - lastnewregion > 0 and startpos == laststartpos')
                #print(str(i) + ':' + str(startpos) + ':' + str(newregions[-1][1]))
                #print(str(i) + ':' + str(laststartpos) + ':' + str(newregions[-1][1]))
                # add final missing region
                newregions.append([laststartpos, laststoppos])
            else:
                #print('execute laststartpos - lastnewregion > 0 and startpos != laststartpos')
                #print(str(i) + ':' + str(startpos) + ':' + str(newregions[-1][1]))
                #print(str(i) + ':' + str(laststartpos) + ':' + str(newregions[-1][1]))
                # add final missing region
                newregions.append([startpos, laststoppos])
    else:
        newregions = featureregions[:]
    return newregions

def comparefeaturemferegion(mfefeature, annotationfeature):
    # merge mfe regions if they overlap
    mfefeature = mergeoverlaps(mfefeature)
    nonmfefeature = list()
    annotationstart = annotationfeature[0][0]
    annotationstop = annotationfeature[0][1]
    mfewindows = len(mfefeature)
    if mfewindows == 0:
        nonmfefeature.append([annotationstart, annotationstop])
        return nonmfefeature
    elif mfewindows > 0:
        mfestart = mfefeature[0][0]
        mfestop = mfefeature[0][1]
        if mfewindows == 1:
            if mfestart == annotationstart and mfestop == annotationstop:
                return nonmfefeature
            elif mfestart == annotationstart:
                # mfe starts at beginning of annotated feature
                nonmfestart = mfestop
                nonmfestop = annotationstop
                nonmfefeature.append([nonmfestart, nonmfestop])
            elif mfestop == annotationstop:
                # mfe ends at end of annotated feature
                nonmfestart = annotationstart
                nonmfestop = mfestart
                nonmfefeature.append([nonmfestart, nonmfestop])
            else:
                # mfe is inside of annotated feature. append twice
                nonmfestart = annotationstart
                nonmfestop = mfestart
                nonmfefeature.append([nonmfestart, nonmfestop])
                nonmfestart = mfestop
                nonmfestop = annotationstop
                nonmfefeature.append([nonmfestart, nonmfestop])
        else:
            for i in range(0, mfewindows-1):
                mfestart = mfefeature[i][0]
                mfestop = mfefeature[i][1]
                nextmfestart = mfefeature[i+1][0]
                nextmfestop = mfefeature[i+1][1]
                if mfestart == annotationstart:
                    # mfe starts at beginning of annotated feature
                    #print('mfestart = annotationstart @ i = ' + str(i))
                    nonmfestart = mfestop
                    nonmfestop = nextmfestart
                    #nonmfefeature.append([nonmfestart, nonmfestop, annotationfeature[0][2], annotationfeature[0][3]])
                    #print('start: ' + str(nonmfestart) + ' stop: ' + str(nonmfestop))
                else:
                    # mfe is in the middle of the annotation
                    #print('mfestart != annotationstart @ i = ' + str(i))
                    if i == 0:
                        # base case where mfe is not at beginning
                        # takes care of first mfe region
                        nonmfestart = annotationstart
                        nonmfestop = mfestart
                        nonmfefeature.append([nonmfestart, nonmfestop])
                        if mfewindows == 2:
                            # case where i = 0 is the second to last
                            #print('mfewindows = 2 and mfe regions do not start at beginning of annotation')
                            nonmfestart = mfestop
                            nonmfestop = nextmfestart
                            nonmfefeature.append([nonmfestart, nonmfestop])
                            #print('start: ' + str(nonmfestart) + ' stop: ' + str(nonmfestop))
                    else:
                        # case where i > 0
                        #print('mfestart != annotationstart  and i != 0 @ i = ' + str(i))
                        if i < mfewindows-1:
                            # case where i isn't second to last
                            #print('execute case where i < ' + str(len(mfefeature)-1) + ' and i is ' + str(i))
                            nonmfestart = mfefeature[i-1][1]
                            nonmfestop = mfestart
                            nonmfefeature.append([nonmfestart, nonmfestop])
                            #print('start: ' + str(nonmfestart) + ' stop: ' + str(nonmfestop))
                            #nonmfestart = mfestop
                            #nonmfestop = nextmfestart
                            #nonmfefeature.append([nonmfestart, nonmfestop, annotationfeature[0][2], annotationfeature[0][3]])
                            #print('start: ' + str(nonmfestart) + ' stop: ' + str(nonmfestop))
                            if i+1 == mfewindows-1:
                                # case where i is second to last? we need the last case
                                #print('execute case where i == ' + str(len(mfefeature)-1) + ' and i is ' + str(i))
                                nonmfestart = mfestop
                                nonmfestop = nextmfestart
                                nonmfefeature.append([nonmfestart, nonmfestop])
                                #print('start: ' + str(nonmfestart) + ' stop: ' + str(nonmfestop))
            if mfefeature[-1][1] != annotationstop:
                # take care of the last mfe region where it doesn't end at the annotation end
                nonmfestart = nextmfestop
                nonmfestop = annotationstop
                nonmfefeature.append([nonmfestart, nonmfestop])
                #print('start: ' + str(nonmfestart) + ' stop: ' + str(nonmfestop))
    return nonmfefeature

def generatenonmfes(mymfes, myannotation):
    mynonmfes = dict()
    for xsome in myannotation.keys():
        mynonmfes[xsome] = dict()
        for featureid in myannotation[xsome].keys():
            if xsome in mymfes.keys():
                if featureid in mymfes[xsome].keys():
                    #do something
                    annotationfeature = myannotation[xsome][featureid]
                    mfefeature = mymfes[xsome][featureid]
                    mynonmfes[xsome][featureid] = comparefeaturemferegion(mfefeature, annotationfeature)
                else:
                    #featureid does not exist in mymfes[xsome]
                    mynonmfes[xsome][featureid] = myannotation[xsome][featureid][:]
            #mynonmfes[myxsome][featureid] = sorted(mynonmfes[myxsome][featureid], key = lambda x: (x[0], x[1]))
            else:
                mynonmfes[xsome] = myannotation[xsome].copy()
    return mynonmfes

def writecsv(mydict, filename):
    filename = filename.split('/')[-1]
    fout = open(outdir + '/' + filename, 'w')
    for key in sorted(mydict.keys()):
        for featureid in sorted(mydict[key].keys()):
            if len(mydict[key][featureid]) > 0:
                for i in range(0, len(mydict[key][featureid])):
                    try:
                        mystring = key + '\t' + str(mydict[key][featureid][i][0]) + '\t' + str(mydict[key][featureid][i][1]) + '\t' + featureid + '\n'
                    except:
                        print(key + ': ' + str(i))
                        print(mydict[key][featureid][i])
                    else:
                        fout.write(mystring)

#mfefilename = 'MFEbeds/lnc_RNA.MFEwindows.sorted.bed'
#filename = 'sortedbeds/lnc_RNA.sorted.bed'
mymfes = readannotation(mfefilename)
myannotation = readannotation(filename)
mynonmfes = generatenonmfes(mymfes, myannotation)
writecsv(mynonmfes, mfefilename.replace('MFE.sorted', 'nonMFE'))
#mfefeature = mymfes['1']['lnc_RNA::1:1030350-1031090']
#annotationfeature = myannotation['1']['lnc_RNA::1:1030350-1031090']
#mergeoverlaps(annotationfeature)
#mergeoverlaps(mfefeature)
#comparefeaturemferegion(mfefeature, annotationfeature)
