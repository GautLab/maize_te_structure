#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""
import os,random, glob, argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description='Graphs a single or list of csv files containing MFE calculated statistics.')
parser.add_argument('-i', '--input', help='A file containing a list of csv files', type=str, required=True)
parser.add_argument('-g', '--graph', help='The type of graph to plot. Chosen from headers of csv file. Default = minMFE', type=str, required=False)
parser.add_argument('-o', '--output', help='Name of png file to output to. Default: Prefix of input filename.png', type=str, required=False)
args = parser.parse_args()

csvlistfilename = args.input
mygraphtype = args.graph
if mygraphtype == None:
    mygraphtype = 'minMFE'
outputfile = args.output
if outputfile == None:
    outputfile = csvlistfilename.replace('.txt','.png')

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

def graphColumn(mycsvhfile,columnname):
    mypddf = pd.read_csv(mycsvhfile)
    sns.set(font_scale=2)
    sns.set_style('ticks')
    sns.set_palette('colorblind')
    sns.distplot(mypddf[columnname])
    sns.despine()
    plt.show()
    return 0

def saveGraphColumn(mycsvhfile,columnname,filename):
    plt.clf()
    mypddf = pd.read_csv(mycsvhfile)
    sns.set(font_scale=2)
    sns.set_style('ticks')
    sns.set_palette('colorblind')
    x = mypddf[columnname]
    snsplt = sns.distplot(x)
    plt.axvline(x.mean(), color='k', linestyle='dashed', linewidth=1)
    sns.despine()
    myfig = snsplt.get_figure()
    myfig.savefig(filename)
    plt.clf()
    return 0

def saveGraphsColumn(mycsvhfilelist,columnname,filename,mypalette):
    plt.clf()
    sns.set(font_scale=2)
    sns.set_style('ticks')
    sns.set_palette(mypalette)
    mypddflist = list()
    for file in mycsvhfilelist:
        mypd = pd.read_csv(file)
        mypddflist.append(mypd)
    for i in range(0,len(mypddflist)):
        snsplt = sns.distplot(mypddflist[i][columnname], label=mycsvhfilelist[i].replace('my','').replace('mfes.csv',''))
        #x = mypddflist[i][columnname]
        #plt.axvline(x.mean(), color='k', linewidth=1)
    sns.despine()
    plt.legend(loc=2, fontsize='xx-small')
    myfig = snsplt.get_figure()
    myfig.savefig(filename)
    plt.clf()
    return 0

if __name__ == '__main__':
    mycsvlist = list()
    for csvfile in open(csvlistfilename):
        mycsvlist.append(csvfile.strip())
    saveGraphsColumn(mycsvlist, mygraphtype, outputfile, 'colorblind')
    '''
    mycsvfiles = glob.glob('*.csv')
    brandonslist = ['mygenemfes.csv', 'mymiRNAmfes.csv', 'mygff3filesmfes.csv', 'mytemfes.csv']
    saveGraphsColumn(brandonslist,'minMFE','BrandonsGraph.png','colorblind')
    mytelist = ['mysolo_LTRmfes.csv', 'myLTR_retrotransposonmfes.csv', 'myLINE_elementmfes.csv', 'myterminal_inverted_repeat_elementmfes.csv', 'myhelitronmfes.csv', 'mySINE_elementmfes.csv']
    '''
