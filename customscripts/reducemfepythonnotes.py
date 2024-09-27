#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""
import os,random, glob, argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#filename = 'windows.out.tsv'
#mydict = dict()
#for line in open(filename,'r'):
#    line = line.replace('\n', '').split()
#    if line[0] in mydict.keys():
#        mydict[line[0]].append(line[1])
#    else:
#        mydict[line[0]] = [line[1]]

parser = argparse.ArgumentParser(description='Create csv file of statistics of MFE calculations for windows.')
parser.add_argument('-r', '--random', help='Total number of windows.out files to extract statistics from. default value is 0', type=int, required=False)
args = parser.parse_args()

myrandsubset = args.random
if myrandsubset == None:
    myrandsubset = 0

windowsize = 110

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

'''
###FINISH ME!!!
def plotMFEsList(myMFElist):
    #create pd dataframe
    mymfepd = pd.read_csv()
    return 0

###FINISH ME!!!
def plotMFEsFile(MFEfilename,header):
    #create pd dataframe
    mymfedf = pd.read_csv(MFEfilename)
    #begin sns
    sns.set(font_scale=2)
    sns.set_style('ticks')
    sns.set_palette('colorblind')
    sns.despine()
    return 0
'''

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

#no mfe
#./lnc_RNA__7_182207619-182208001/windows.out
#has mfe
#./lnc_RNA__10_142389865-142390471/windows.out

# execution order:
# myposmasterlist = list()
# myposmasterlist.append(returnmfepos(returnWindowMFEList(filename),-40.0))

def processgenomefeatures(windowslist):
    filename = ''
    myposmasterlist = list()
    for line in windowslist:
        filename = line.strip()
        mytemplist = returnmfepos(returnWindowMFEList(filename), -40.0)
        if mytemplist != 1:
            myposmasterlist.append(mytemplist)
    prefix = filename[2:len(filename)].split('__')[0]
    writebedfile(myposmasterlist,prefix)
    return [myposmasterlist,prefix]

def returnmfepos(mfelist, minmfe):
    mfepos = list()
    header = mfelist[0][0][0:len(mfelist[0][0])-2]
    initialpos = int(header.split('::')[1].split(':')[1].split('-')[0])
    minmfeexists = False
    mfeposlist = list()
    for i in range(0,len(mfelist)):
        if mfelist[i][1] <= minmfe:
            mfepos.append(i+1+initialpos)
            minmfeexists = True
    if minmfeexists:
        #print 'minMFE found!'
        mfeposlist = checkcontinuity(mfepos)
    else:
        #print "No minMFE found for minMFE=" + str(minmfe)
        return 1
    #return [header, mfeposlist, mfepos]
    return [header, mfeposlist]

def checkcontinuity(mfelist):
    pos = startpos = mfelist[0]
    myposlist = list()
    for i in mfelist:
        if i > pos:
            stoppos = pos-1
            myposlist.append([startpos,stoppos+windowsize-1])
            startpos = i
            #print str(i) +  '>' + str(pos)
            pos = i+1
        else:
            pos += 1
        if(i == mfelist[len(mfelist)-1]):
            stoppos = mfelist[len(mfelist)-1]+windowsize-1
            myposlist.append([startpos, stoppos])
    return myposlist

def writebedfile(myposmasterlist,prefix):
    fout = open(prefix + '.bed','w')
    for i in myposmasterlist:
        if len(i) != 1:
            header = i[0]
            myposlist = i[1]
            myxsome = header.split('::')[1].split(':')[0]
            for j in myposlist:
                startpos = j[0]
                stoppos = j[1]
                fout.write(myxsome + '\t' + str(startpos) + '\t' + str(stoppos) + '\t' + header + '\t42\t+\t' + str(stoppos-startpos) + 'M\n')
    fout.close()

def returnWindowMFEStats(myset):
    mymasterlist = list()
    counter = 0
    minmfe = -40
    for file in myset:
        windowline = " "
        mymfelist = list()
        tempList = list()
        fin = open(file, 'r')
        counter += 1
        #print(counter)
        while windowline != "":
            windowline = fin.readline().strip().replace('>', '')
            header = windowline
            if header != '':
                myseq = fin.readline()
                windowline = fin.readline().strip()
                floatline = windowline.split(' ')
                floatline = floatline[len(floatline) - 1].replace('(', '').replace(')', '')
                myfloat = float(floatline)
                mymfelist.append([header, myfloat])
        #call plotter here
        #plotMFEs(mylist)
        tempList = min(mymfelist, key=lambda item: item[1]) #returns a list [header, mfe]
        mymean = calculateMean(mymfelist)
        tempList.append(mymean)
        tempList.append(calculateVariance(mymfelist,mymean))
        mfepctcov = calculatePCTThreshMFE(mymfelist,minmfe)
        #mymfepos = returmmfepos(mymfelist,minmfe)
        tempList.append(mfepctcov[1])
        tempList.append(mfepctcov[0])
        #templist = [headerid, minMFE, meanMFE, MFEvariance, MFEcountbelowthresh, %MFEbelowthresh]
        mymasterlist.append(tempList)
    return mymasterlist

'''
def writeWindowSet(myset, filename):
    fout = open(filename, 'w')
    for file in myset:
        fout.write(file + '\n')
    fout.close()

def writeCSVList(mylist, filename):
    fout = open(filename, 'w')
    for list in mylist:
        fout.write(str(list[0]))
        for i in range(1,len(list)):
            fout.write(',' + str(list[i]))
        fout.write('\n')
    fout.close()
    return 0
'''

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

'''

mycsvfiles = glob.glob('*.csv')
brandonslist = ['mygenemfes.csv', 'mymiRNAmfes.csv', 'mygff3filesmfes.csv', 'mytemfes.csv']
saveGraphsColumn(brandonslist,'minMFE','BrandonsGraph.png','colorblind')
mytelist = ['mysolo_LTRmfes.csv', 'myLTR_retrotransposonmfes.csv', 'myLINE_elementmfes.csv', 'myterminal_inverted_repeat_elementmfes.csv', 'myhelitronmfes.csv', 'mySINE_elementmfes.csv']

'''

if __name__ == '__main__':
    #mywindowsdirs = ['//networkshare//maize//rnafold//categories//gene//contigs','//networkshare//maize//rnafold//categories//pre_miRNA//contigs','//networkshare//maize//rnafold//categories//gff3files//contigs','//networkshare//maize//rnafold//categories//rRNA//contigs','//networkshare//maize//rnafold//categories//lnc_RNA//contigs','//networkshare//maize//rnafold//categories//snoRNA//contigs','//networkshare//maize//rnafold//categories//miRNA//contigs','//networkshare//maize//rnafold//categories//snRNA//contigs','//networkshare//maize//rnafold//categories//mRNA//contigs','//networkshare//maize//rnafold//categories//SRP_RNA//contigs','//networkshare//maize//rnafold//categories//ncRNA_gene//contigs','//networkshare//maize//rnafold//categories//tRNA//contigs','//networkshare//maize//rnafold//alexdata//sirevirus//contigs','//networkshare//maize//rnafold//alexdata//LTR//contigs']
    #mywindowsdirs = list()
    mywindowsdirs = glob.glob('*/windows.out')
    myheaders = ['header', 'minMFE', 'meanMFE', 'MFEvariance', 'MFEcountltthresh', 'PctMFEltthresh']
    if (len(mywindowsdirs) > 0):
        mytype = mywindowsdirs[0].split('/')[0].split('__')[0]
    else:
        print('No windows.out files exist in subdirectories. Please rerun RNAfold')
        quit(1)
    mycsvfile = 'my_' + mytype + '_mfes.csv'
    mymfefileset = returnRandomSet(mywindowsdirs, myrandsubset)
    mymfelist = returnWindowMFEStats(mymfefileset)
    writeCSVListH(mymfelist, mycsvfile, myheaders)
    for i in range(1,len(myheaders)):
        mypngfile = 'my' + mytype + 'calc' + myheaders[i] + 's.png'
        saveGraphColumn(mycsvfile,myheaders[i],mypngfile)

    myfilelist = list()
    for line in open('windowoutfiles.txt'):
        myfilelist.append(line.strip())
    mytemplist = processgenomefeatures(myfilelist)

    ###end for loop

mycsvdir = '//networkshare//maize//rnafold//categories//myoutputs'
os.chdir(mycsvdir)
mycsvlist = glob.glob('*.csv')
myheaders = ['minMFE', 'meanMFE', 'MFEvariance', 'MFEcountltthresh', 'PctMFEltthresh']
for csv in mycsvlist:
    for header in myheaders:
        mypngfilename = csv.replace('.csv','') + '_' + header + '.png'
        saveGraphColumn(csv,header,mypngfilename)

header = ['header', 'minMFE', 'meanMFE', 'MFEvariance', 'MFEcountltthresh', 'PctMFEltthresh']
mymirnacsvhfile = 'mymiRNAcalcmfes.csv'
mymirnapngfile = 'mytecalcminmfes.png'
os.chdir('//networkshare//maize//rnafold//categories//genes//contigs//miRNA')
mymirnafileset = returnRandomSet('windowoutfiles.txt',0)
#writeWindowSet(mymirnaset,'allmymiRNA.txt')
#myFile = mymirnaset.pop()
mymirnaset = returnWindowMFEStats(mymirnafileset)
writeCSVListH(mymirnaset, mymirnacsvhfile,header)
#graphColumn(mymirnacsvhfile,'minMFE')
saveGraphColumn(mymirnacsvhfile,'minMFE',mymirnapngfile)



header = ['header', 'minMFE', 'meanMFE', 'MFEvariance', 'MFEcountltthresh', 'PctMFEltthresh']
mygenecsvhfile = 'mygenecalcmfes.csv'
mygenepngfile = 'mygenecalcminmfes.png'
os.chdir('//networkshare//maize//rnafold//categories//genes//contigs')
mygenefileset = returnRandomSet('windowoutfiles.txt',0)
#writeWindowSet(mygeneset,'myrandom250genes.txt')
#myFile = mygeneset.pop()
mygenelist = returnWindowMFEStats(mygenefileset)
writeCSVList(mygenelist, mygenecsvhfile, header)
#graphColumn(mygenecsvhfile,'minMFE')
saveGraphColumn(mygenecsvhfile,'minMFE',mygenepngfile)

header = ['header', 'minMFE', 'meanMFE', 'MFEvariance', 'MFEcountltthresh', 'PctMFEltthresh']
mytecsvhfile = 'mytecalcmfes.csv'
mytepngfile = 'mytecalcminmfes.png'
os.chdir('//networkshare//maize//rnafold//TEs//contigs')
mytefileset = returnRandomSet('windowoutfiles.txt',0)
#writeWindowSet(mytefileset,'myrandom250tes.txt')
mytelist = returnWindowMFEStats(mytefileset)
writeCSVList(mytelist, mytecsvhfile, header)
#graphColumn(mytecsvhfile,'MFEcountltthresh')
saveGraphColumn(mytecsvhfile,'minMFE',mytepngfile)


##########
import os,random
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
os.chdir('//networkshare//maize//rnafold')
mymirnacsvhfile = 'mymiRNAcalcmfes.csv'
mygenecsvhfile = 'mygenecalcmfes.csv'
mytecsvhfile = 'mytecalcmfes.csv'

mygenedf = pd.read_csv(mygenecsvhfile)
mytedf = pd.read_csv(mytecsvhfile)
mymirnadf = pd.read_csv(mymirnacsvhfile)

mycol = 'PctMFEltthresh'
sns.set(font_scale=2)
sns.set_style('ticks')
#sns.set_palette('colorblind')
sns.distplot(mygenedf[mycol], label='Gene')
sns.distplot(mytedf[mycol], label='TE')
sns.distplot(mymirnadf[mycol], label='miRNA')
sns.despine()
plt.legend()
plt.show()
#mean of all counts per type
#binomial test for the porportion
###################################test code for 250 in windows
##tes
myteset = set()
os.chdir('TEs')
for line in open('myrandom250tes.txt'):
    myteset.add(line.strip())

mymastertelist = list()
counter = 0
for file in myteset:
    mymastertelist.append(returnWindowMFEList(file))
    counter +=1
    print(counter)

##genes
mygeneset = set()
#os.chdir('genes')
for line in open('myrandom250genes.txt'):
    mygeneset.add(line.strip())

counter = 0
mymastergenelist = list()
for file in mygeneset:
    mymastergenelist.append(returnWindowMFEList(file))
    counter += 1
    print(counter)

###################################test code for windows.out
windowline = " "
mylist = list()
filename = ''
fin = open(filename,'r')
while windowline != "":
    header = fin.readline().strip().replace('>', '')
    if header != '':
        myseq = fin.readline()
        myline = fin.readline().strip()
        myline = myline.split(' ')
        myline = myline[len(myline) - 1].replace('(', '').replace(')', '')
        myfloat = float(myline)
        mylist.append([header, myfloat])

#SRP_RNA::1:229766138-229766442
###################################individual code
mymasterlist = list()
counter = 0
for line in open('windowoutfiles.txt','r'):
    line = line.strip()
    fin = open(line,'r')
    wline = " "
    mylist = list()
    #headermin = ''
    #mymin = 10000
    while wline != "":
        header = fin.readline().strip().replace('>','')
        if header != '':
            myseq = fin.readline()
            myline = fin.readline().strip()
            myline = myline.split(' ')
            myline = myline[len(myline)-1].replace('(','').replace(')','')
            myfloat = float(myline)
            mylist.append([header,myfloat])
    counter += 1
    #mymasterlist.append(mylist)
    mymasterlist.append(min(mylist, key=lambda item: item[1]))


for key in mydict:
    mylist.append([int(key.split('_')[3]), key, min(mydict[key])])

fout = open('windows.out.red.tsv','w')
for val in mylist:
    fout.write(str(val[0]) + '\t' + val[1] + '\t' + val[2] + '\n')

fout.close()

'''
sns.distplot(mygenesmfedf['minMFE'])
Out[8]: <matplotlib.axes._subplots.AxesSubplot at 0x7fd12ecb6590>
sns.set(font_scale=2)
sns.set_style('ticks')
sns.set_palette('colorblind')
sns.distplot(mygenesmfedf['minMFE'])
Out[10]: <matplotlib.axes._subplots.AxesSubplot at 0x7fd12ed11d50>
sns.despine()
sns.distplot(mytesmfedf['minMFE'])
Out[12]: <matplotlib.axes._subplots.AxesSubplot at 0x7fd12ed11d50>
import scipy.stats as stats
stats.shapiro(mytesmfedf['minMFE'])
Out[14]: (0.9755634665489197, 0.00026618794072419405)
stats.shapiro(mygenesmfedf['minMFE'])
Out[15]: (0.9856363534927368, 0.012953350320458412)
stats.f_oneway(mygenesmfedf['minMFE'],mytesmfedf['minMFE'])
Out[16]: F_onewayResult(statistic=19.088020256506073, pvalue=1.5197866219841495e-05)
stats.f_oneway(mygenesmfedf['meanMFE'],mytesmfedf['meanMFE'])
Out[17]: F_onewayResult(statistic=7.9871007785979877, pvalue=0.004900527102706647)
sns.set(font_scale=2)
sns.set_style('ticks')
sns.set_palette('colorblind')
sns.distplot(mygenesmfedf['meanMFE'])
Out[19]: <matplotlib.axes._subplots.AxesSubplot at 0x7fd12e4e1e50>
sns.distplot(mytesmfedf['meanMFE'])
Out[20]: <matplotlib.axes._subplots.AxesSubplot at 0x7fd12e4e1e50>
stats.shapiro(mygenesmfedf['meanMFE'])
Out[21]: (0.8915553092956543, 2.090320971870341e-12)
stats.shapiro(mytesmfedf['meanMFE'])
Out[22]: (0.8611882328987122, 3.063789320986374e-14)
stats.kstest(mygenesmfedf['meanMFE'],'norm')
Out[23]: KstestResult(statistic=1.0, pvalue=0.0)
stats.kstest(mytesmfedf['meanMFE'],'norm')
Out[24]: KstestResult(statistic=1.0, pvalue=0.0)
stats.kstest(mytesmfedf['minMFE'],'norm')
Out[25]: KstestResult(statistic=1.0, pvalue=0.0)
stats.kstest(mygenesmfedf['minMFE'],'norm')
Out[26]: KstestResult(statistic=1.0, pvalue=0.0)
stats.kstest(mytesmfedf['meanMFE'],mygenesmfedf['meanMFE'])
'''



#how to process transcripts
line = ''
mylist = line.split(':')
header =  mylist[0][0:-2].replace('>','')
positions = list()
for i in mylist[1:]:
    mypos = i[0:-2]
    mypos = mypos.split('-')
    startpos = int(mypos[0])
    stoppos = int(mypos[1])
    positions.append([startpos,stoppos])

# for processing MASiVE gff
# header = [0:-2]
# grep -a header Zmay_v4_chr_ALL.fasta.clean.MASiVE.gff | awk '{ print $1"\t"$4"\t"$5"\t"$7 }'
# cut -f2 -d"/" windowoutfiles.txt
mybedset = set()

len(mybedset)

mygffset-mybedset
mybedset-mygffset

mygffdict = dict()
mybeddict = dict()
for key in mygffset:
    mygffdict[key] = []

for key in mybedset:
    mybeddict[key] = []

for line in open(mygfffile):
    mylist = line.strip().split()
    xsome = mylist[0]
    startpos = int(mylist[2])
    stoppos = int(mylist[3])

mygffdict
{'10': [], 'Pt': [], '1': [], 'Mt': [], '3': [], '2': [], '5': [], '4': [], '7': [], '6': [], '9': [], '8': []}
mygfffile = ''
mybedfile = ''

for line in open(mygfffile):
    mylist = line.strip().split()
    xsome = mylist[0]
    startpos = int(mylist[3])
    stoppos = int(mylist[4])

for line in open(mygfffile):
    mylist = line.strip().split()
    startpos = int(mylist[3])
    stoppos = int(mylist[4])
    xsome = mylist[0]
    mygffdict[xsome].append([startpos,stoppos])

for line in open(mybedfile):
    mylist = line.strip().split()
    xsome = mylist[0]
    startpos = int(mylist[3])
    startpos = int(mylist[2])
    stoppos = int(mylist[3])
    mybeddict[xsome].append([startpos,stoppos])


####find nonmfes = mfes -
import os,glob
import multiprocessing as mp

os.chdir('/networkshare/maize/rnafold/sRNA_reads')

countfiles = glob.glob('*.counts')
nonmfecountfiles = glob.glob('*.not*.counts')
mfecountfiles = set(countfiles) - set(nonmfecountfiles)
#nonmfecountfiles = set(nonmfecountfiles)

#mfefile = 'snoRNA.sorted_maizeAGPv4_SRR032091_cleana_trimm_lt28_Zea_mays.AGPv4.dna.toplevel_gtq10_mapped.21.sorted.col4.coverage.counts'
#nonmfefile = 'snoRNA.not_maizeAGPv4_SRR032091_cleana_trimm_lt28_Zea_mays.AGPv4.dna.toplevel_gtq10_mapped.21.sorted.col4.coverage.counts'

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
bedfiles = glob.glob('*.bed')
notbedfiles = glob.glob('*.not.bed')
sortedbedfiles = glob.glob('*.sorted.bed')

def returnBedDict(bedfiles):
    xsomeset = set()
    mybeddict = dict()
    for bedfile in bedfiles:
        if os.stat(bedfile).st_size > 0:
            for line in open(bedfile):
                line = line.strip().split()
                #print '1st loop: ' + bedfile
                xsomeset.add(line[0])
    for xsome in xsomeset:
        mybeddict[xsome] = dict()
    for bedfile in bedfiles:
        if os.stat(bedfile).st_size > 0:
            for line in open(bedfile):
                line = line.strip().split()
                #print '2nd loop: ' + bedfile
                xsome = line[0]
                seqid = line[3]
                startpos = int(line[1])
                stoppos = int(line[2])
                if seqid in mybeddict[xsome].keys():
                    mybeddict[xsome][seqid].append(stoppos-startpos)
                else:
                    mybeddict[xsome][seqid] = [stoppos-startpos]
    return mybeddict

mybeddict = returnBedDict(sortedbedfiles)

#single thread
def parseCountFiles(mfecountfiles):
    for mfefile in mfecountfiles:
        mymfedict = dict()2[0]
            seqid = mymfelist[3]
            count = int(mymfelist[8])
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
    for mfeline in open(mfefile):
        mymfexsomeset.add(mfeline.strip().split()[0])
    for xsome in mymfexsomeset:
        mymfedict[xsome] = dict()
    for mfeline in open(mfefile):
        mymfelist = mfeline.strip().split()
        xsome = mymfelist[0]
        seqid = mymfelist[3]
        count = int(mymfelist[5])
        if seqid in mymfedict[xsome].keys():
            mymfedict[xsome][seqid].append(count)
        else:
            mymfedict[xsome][seqid] = [count]
    print str(len(mymfedict))
    writemfecountstocsv(mymfedict,mfefile)

pool = mp.Pool(11)
pool.map(mpParseCountFiles, mfecountfiles)

###
#here we will convert gff3 files to bed files with proper naming conventions used by bedtools
#awk '{ print $1"\t"$4"\t"$5"\t"$3"::"$1":"$4-1"-"$5  }' ${i} > $(basename ${i} .gff3).annotation.bed

#use below as a template for extracting the id and then create a 2nd column of the bedtools new id
#head -n 10 Zea_mays.AGPv4.39.chr_miRNA.gff3 | awk '{ print $9 }' | cut -f1 -d";" | cut -f2 -d":"




writemfecountstocsv(mynonmfedict,nonmfefile)

glob.glob("*.hits")
'snoRNA.sorted_maizeAGPv4_SRR032091_cleana_trimm_lt28_Zea_mays.AGPv4.dna.toplevel_gtq10_mapped.21.sorted.col4.coverage.hits'
'snoRNA.not_maizeAGPv4_SRR032091_cleana_trimm_lt28_Zea_mays.AGPv4.dna.toplevel_gtq10_mapped.21.sorted.col4.coverage.hits'
####

sorted(mybeddict[xsome]) # sort by startpos
sorted(mygffdict[xsome]) # sort by startpos




###build algo for detecting overlap
def detectoverlap(alpha, beta, delta, gamma):
    mymin = min(alpha, delta)
    mymax = max(beta, gamma)
    firstlen = beta - alpha
    secondlen = gamma - delta
    if ((mymax-mymin) <= (firstlen + secondlen)):
        return 1
    else:
        return 0
