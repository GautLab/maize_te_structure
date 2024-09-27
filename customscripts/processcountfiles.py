#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""

import gzip, argparse, glob
import pandas as pd
from os import path
#usage: processcovfiles.py -i filename.coverage

parser = argparse.ArgumentParser(description='Process count files and output large csv.')
parser.add_argument('-i', '--input', help='Annotation file name', type=str, required=True)
parser.add_argument('-c', '--countfiles', help='Directory containing count files', type=str, required=False)
parser.add_argument('-o', '--outdir', help='Output directory name', type=str, required=False)
args = parser.parse_args()
annotationfile = args.input
countfiles = args.countfiles
if countfiles == None:
	countfiles = 'countfiles'
outdir = args.outdir
if outdir == None:
	outdir = '.'
elif outdir[-1] == '/':
	outdir = outdir[0:-1]

#annotationfiles = glob.glob('MFEwindowregions/sortedbeds/*.bed')
#annotationfile = annotationfiles[7]
mypd = pd.read_csv(annotationfile, sep='\t', header=None)
mypd = mypd[[0,1,2,3,5]]
mypd.columns = ['Chromosome', 'StartPos', 'StopPos', 'FeatureID', 'Strand']

myfeatureids = set(mypd['FeatureID'])
featureid = annotationfile.split('/')[-1].split('.')[0]
mfe='MFE'
nonmfe='non' + mfe

libraries = list()
for i in open('libraries.txt'):
	libraries.append(i.strip())

sirnalens = list()
for i in open('sirnalengths.txt'):
	sirnalens.append(i.strip())

#mfefilenames = glob.glob(countfiles + '/' + featureid + '_MFE_*.counts')
mfecounter = 0
nonmfecounter = 0
for library in libraries:
	for sirnalen in sirnalens:
		#library = libraries[15]
		#sirnalen = sirnalens[2]
		mfefile = countfiles + '/' + featureid + '_' + mfe + '_maizeAGPv4.' + library + '.ca.trimm_' + sirnalen + '.ucounts.mapped.counts'
		#mfefile = countfiles + '/' + 'lnc_RNA_nonMFE_maizeAGPv4.SRX483603.ca.trimm_24.ucounts.mapped.sorted.bed.counts'
		nonmfefile = countfiles + '/' + featureid + '_' + nonmfe + '_maizeAGPv4.' + library + '.ca.trimm_' + sirnalen + '.ucounts.mapped.counts'
		myprefix = library + '_' + sirnalen
		# Do first library
		if path.exists(mfefile):
			mfecounter += 1
			mymfelibpd = pd.read_csv(mfefile, sep='\t')
			if mfecounter == 1:
				mymfelibpd = mymfelibpd[['FeatureID', 'Length', 'siRNA_Counts', 'CountsPerNT']]
				mymfelibpd.rename(columns = {'Length' : 'MFE_Region_Length'}, inplace = True)
			else:
				mymfelibpd = mymfelibpd[['FeatureID', 'siRNA_Counts', 'CountsPerNT']]
			if len(set(mymfelibpd['FeatureID']) - myfeatureids) > 0:
				# do something here if featureids don't match with annotation
				print('FeatureIDs do not match in ' + mfefile)
		else:
			# make an empty pandas df
			mymfelibpd = pd.DataFrame(columns = ['FeatureID', 'siRNA_Counts', 'CountsPerNT'])
		mymfelibpd.rename(columns = {'siRNA_Counts' : 'MFE_' + myprefix + '_siRNA_Species_Counts'}, inplace = True)
		mymfelibpd.rename(columns = {'CountsPerNT' : 'MFE_' + myprefix + '_siRNA_Species_CountsPerNT'}, inplace = True)
		mypd = pd.merge(mypd, mymfelibpd, how='left', on=['FeatureID','FeatureID'])
		if path.exists(nonmfefile):
			nonmfecounter += 1
			mynonmfelibpd = pd.read_csv(nonmfefile, sep='\t')
			if nonmfecounter == 1:
				mynonmfelibpd = mynonmfelibpd[['FeatureID', 'Length', 'siRNA_Counts', 'CountsPerNT']]
				mynonmfelibpd.rename(columns = {'Length' : 'nonMFE_Region_Length'}, inplace = True)
			else:
				mynonmfelibpd = mynonmfelibpd[['FeatureID', 'siRNA_Counts', 'CountsPerNT']]
			if len(set(mynonmfelibpd['FeatureID']) - myfeatureids) > 0:
				# do something here if featureids don't match with annotation
				print('FeatureIDs do not match in ' + nonmfefile)
		else:
			# make an empty pandas df
			mynonmfelibpd = pd.DataFrame(columns = ['FeatureID', 'siRNA_Counts', 'CountsPerNT'])
		mynonmfelibpd.rename(columns = {'siRNA_Counts' : 'nonMFE_' + myprefix + '_siRNA_Species_Counts'}, inplace = True)
		mynonmfelibpd.rename(columns = {'CountsPerNT' : 'nonMFE_' + myprefix + '_siRNA_Species_CountsPerNT'}, inplace = True)
		mypd = pd.merge(mypd, mynonmfelibpd, how='left', on=['FeatureID','FeatureID'])

mypd = mypd.fillna(0)

mypd.to_csv(annotationfile.split('/')[-1].replace('.sorted.bed','.siRNA_count_stats.csv'), sep='\t', header=True)
