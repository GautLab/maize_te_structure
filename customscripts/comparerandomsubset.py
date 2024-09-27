#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""
import pandas as pd

def returncleanlist(fn):
    fin = open(fn)
    randomlist = fin.readlines()
    fin.close()
    for i in range(0,len(randomlist)):
        randomlist[i] = randomlist[i].strip().split(',')
    for i in range(0,len(randomlist)):
        mystr = ''
        for j in randomlist[i][0].split('_')[0:-1]:
            mystr = mystr + j + '_'
        mystr = mystr[0:-1]
        randomlist[i][0] = mystr
    return randomlist

def writecsvfile(myposmasterlist,fn):
    fout = open(fn,'w')
    for myposlist in myposmasterlist:
        mystr = ''
        for j in myposlist:
            mystr = mystr + j + ','
    	fout.write(mystr[0:-1] + '\n')
    fout.close()

allsololtrs = pd.read_csv(listofsololtrfiles[0], header=None)
randsltr1 = pd.read_csv(listofsololtrfiles[1], header=None)
randsltr2 = pd.read_csv(listofsololtrfiles[2], header=None)
randsltr3 = pd.read_csv(listofsololtrfiles[3], header=None)
randsltr4 = pd.read_csv(listofsololtrfiles[4], header=None)
randsltr5 = pd.read_csv(listofsololtrfiles[5], header=None)
randsltrsubset = pd.read_csv(listofsololtrfiles[6], header=None)

allsololtrs.columns = ['Header', 'minMFE', 'meanMFE', 'MFEvariance', 'MFEcountltthresh', 'PctMFEltthresh']
randsltr1.columns = ['Header', 'minMFE', 'meanMFE', 'MFEvariance', 'MFEcountltthresh', 'PctMFEltthresh']
randsltr2.columns = ['Header', 'minMFE', 'meanMFE', 'MFEvariance', 'MFEcountltthresh', 'PctMFEltthresh']
randsltr3.columns = ['Header', 'minMFE', 'meanMFE', 'MFEvariance', 'MFEcountltthresh', 'PctMFEltthresh']
randsltr4.columns = ['Header', 'minMFE', 'meanMFE', 'MFEvariance', 'MFEcountltthresh', 'PctMFEltthresh']
randsltr5.columns = ['Header', 'minMFE', 'meanMFE', 'MFEvariance', 'MFEcountltthresh', 'PctMFEltthresh']
randsltrsubset.columns = ['Header']
realsltr1 = pd.merge(allsololtrs, randsltrsubset, how='inner', on='Header')

#pd.merge(randsltr1[['Header','minMFE']],randsltr2[['Header','minMFE']])
wholepd = pd.merge(pd.merge(pd.merge(pd.merge(pd.merge(realsltr1[['Header','minMFE']],randsltr1[['Header','minMFE']], on='Header', how='outer'),randsltr2[['Header','minMFE']], on='Header', how='outer'),randsltr3[['Header','minMFE']], on='Header', how='outer'),randsltr4[['Header','minMFE']], on='Header', how='outer'),randsltr5[['Header','minMFE']], on='Header', how='outer')
