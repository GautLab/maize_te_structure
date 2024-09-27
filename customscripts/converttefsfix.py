#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""

import glob, os, sys
import multiprocessing as mp


#fix error in headers where dupes exist: i.e. LINE_element::LINE_element::10:119228570-119229550_1
feature = 'solo_LTR'
startdir = './contigs/'

def fixheaders(feature, startdir):
  mydirs = glob.glob(startdir + feature + '/*')
  #fasta files
  for mydir in mydirs:
      fixfile(mydir, '.fa', 2)
  for mydir in mydirs:
      fixfile(mydir, '.out', 3)
  for mydir in mydirs:
      fixfile(mydir, '.out.tsv', 1)

def fixfile(mydir, extension, mylineindex):
      if os.path.exists(mydir + '/windows.' + extension):
          fin = open(mydir + '/windows.' + extension)
          myfile = fin.readlines()
          fin.close()
          listlen = len(myfile)
          for listindex in range(0,listlen):
              if listindex%mylineindex == 0:
                  myfile[listindex] = myfile[listindex].replace(feature + '::' + feature, feature)
          fout = open(mydir + '/windows.' + extension,'w')
          fout.writelines(myfile)
          fout.close()
