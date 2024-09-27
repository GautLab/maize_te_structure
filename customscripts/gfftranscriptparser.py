#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
@author: Edwin
"""
import sys, os, math
import subprocess as sp
import multiprocessing as mp

try:
    sys.argv[1]
except:
    print('Invalid input gff3 file or gff3 file missing')
    quit()
else:
    if sys.argv[1] == '':
        gff3fname = sys.argv[1]
    else:
        gff3fname = sys.argv[1]

##pass args gff3 filename

try:
    sys.argv[2]
except:
    print('Invalid input file or input file missing')
    quit()
else:
    if sys.argv[2] == '':
        genome = sys.argv[2]
    else:
        genome = sys.argv[2]

try:
    sys.argv[3]
except:
    print('Invalid number of threads')
    quit()
else:
    if sys.argv[3] == '':
        threads = sys.argv[3]
    else:
        threads = sys.argv[3]

try:
    sys.argv[4]
except:
    print('Invalid or missing output directory')
    quit()
else:
    if sys.argv[3] == '':
        outdir = sys.argv[4]
    else:
        outdir = sys.argv[4]

mytranscriptsdict = dict()
mytrlendict = dict()
#gff3fname = 'Zea_mays.AGPv4.39.chr_exon.gff3'
#genome = 'Zea_mays.AGPv4.dna.toplevel.fa'
#threads = 11
#outdir = 'contigs/'
#created by using bioawk -cfastx '{print $name, length($seq)}' transcripts.fasta
sortfile = 'contigs/transcripts.lengths'
outfile = 'contigs/transcripts.sorted.lengths'

def parsegff3(gff3fname):
    mygff3dict = dict()
    counter = 0
    for line in open(gff3fname):
        counter += 1
        if line[0] != '#':
            myline = line.strip().split()
            xsome = myline[0]
            db = myline[1]
            type = myline[2]
            pos1 = int(myline[3])
            pos2 = int(myline[4])
            direction = myline[6]
            transcriptID = myline[8].split(';')[0].split(':')[1]
            #if direction == '+':
            startpos = min(pos1, pos2)
            stoppos = max(pos1, pos2)
            #else:
                #startpos = max(pos1, pos2)
                #stoppos = min(pos1, pos2)
            if transcriptID in mygff3dict.keys():
                if direction == mygff3dict[transcriptID][3]:
                    if xsome == mygff3dict[transcriptID][0]:
                        mygff3dict[transcriptID][4].append([startpos, stoppos])
                    else:
                        print('Error chromosome of exons for transcriptID:' + transcriptID + ' do not match')
                        break
                else:
                    print('Error in directions of exons for transcriptID:' + transcriptID + ' do not match')
                    break
            else:
                mygff3dict[transcriptID] = [xsome, db, type, direction, [[startpos, stoppos]]]
        if counter % 10000 == 0:
            print 'Exons processed: ' + str("{:,}".format(counter))
    print str('{:,}'.format(counter)) + ' Total exon"s processed'
    print str('{:,}'.format(len(mygff3dict))) + ' Total transcript"s processed'
    return mygff3dict

def aggregatetranscript(gff3fname, mygff3dict):
    counter = 0
    for line in open(gff3fname):
        counter += 1
        if line[0] != '#':
            myline = line.strip().split()
            xsome = myline[0]
            db = myline[1]
            type = myline[2]
            pos1 = int(myline[3])
            pos2 = int(myline[4])
            direction = myline[6]
            transcriptID = myline[8].split(';')[0].split(':')[1]
            #if direction == '+':
            startpos = min(pos1, pos2)
            stoppos = max(pos1, pos2)
            #else:
                #startpos = max(pos1, pos2)
                #stoppos = min(pos1, pos2)
            if transcriptID in mygff3dict.keys():
                olddirection = mygff3dict[transcriptID][3]
                if olddirection == direction:
                    if xsome == mygff3dict[transcriptID][0]:
                        mygff3dict[transcriptID][4].append([startpos, stoppos])
                    else:
                        print('Error chromosome of UTR for transcriptID:' + transcriptID + ' do not match')
                        break
                else:
                    print('Error in directions of UTR for transcriptID:' + transcriptID + ' do not match')
                    break
            else:
                print('Error in transcriptID:' + transcriptID + '. Does not exist in reference dictionary')
                break
        if counter % 10000 == 0:
            print 'UTR"s processed: ' + str('{:,}'.format(counter))
    print str('{:,}'.format(counter)) + ' Total UTR"s processed'
    print str('{:,}'.format(len(mygff3dict))) + ' Total transcript"s processed'
    return mygff3dict

def returnfasta(key,mygfflist,genome):
    myheader = key
    mysequence = ''
    for i in mygfflist[4]:
        myfalist = sp.check_output(['faidx', genome, mygfflist[0] + ':' + str(i[0]) + '-' + str(i[1])]).split()
        myheader = myheader + '_' + mygfflist[0] + ':' + str(i[0]) + '-' + str(i[1])
        for j in range(1,len(myfalist)):
            mysequence = mysequence + myfalist[j]
    if mygfflist[3] == '-':
        mysequence = reversecomplement(mysequence)
    return [myheader,mysequence]

def reversecomplement(seq):
    seq.upper()
    rev_seq = ''
    for i in seq:
        if i == 'A':
            rev_seq = rev_seq + 'T'
        elif i == 'T':
            rev_seq = rev_seq + 'A'
        elif i == 'C':
            rev_seq = rev_seq + 'G'
        elif i == 'G':
            rev_seq = rev_seq + 'C'
        else:
            rev_seq = rev_seq + i
    rev_seq = rev_seq[::-1]
    return rev_seq

def checkoverlap(mylist):
    for i in range(1,len(mylist)):
        myvar1 = mylist[i-1][0]
        myvar2 = mylist[i-1][1]
        myvar3 = mylist[i][0]
        myvar4 = mylist[i][1]
        if((max(myvar2,myvar4) - min(myvar1,myvar3)) <= (myvar2 - myvar1 + myvar4 - myvar3)):
            return 1
            print str(myvar1) + ':' + str(myvar2) + ',' + str(myvar3) + ':' + str(myvar4)
    return 0

def sorttranscripts(mytranscriptdict):
    counter = 1
    for key in mytranscriptdict:
        #direction = mytranscriptdict[key][3]
        #if direction == '+':
            #sorted(mytranscriptdict[key][4], key=lambda l:l[0]) #sort list in ascending based on index 0
        #else:
            #sorted(mytranscriptdict[key][4], key=lambda l:l[1], reverse=True) #sort list in descending based on index 1
        mylist = mytranscriptsdict[key][4]
        sorted(mylist, key=lambda l: l[0])
        if checkoverlap(mylist) != 1:
            mytranscriptsdict[key][4] = mylist
            counter += 1
        else:
            print 'Error! overlap detected in key:' + key
            quit(1)
        if counter % 10000 == 0:
            print 'Transcript"s processed: ' + str('{:,}'.format(counter))
        print str('{:,}'.format(counter)) + ' Total transcripts"s processed'
    return mytranscriptdict

def writefasta(key, genome, outdir):
    mygfflist = mytranscriptsdict[key]
    myfasta = returnfasta(key, mygfflist, genome)
    fout = open(outdir + key + '.fasta', 'w')
    fout.write('>' + myfasta[0] + '\n')
    fout.write(myfasta[1] + '\n')
    fout.close()
    return

###testing
gff3fname = 'Zea_mays.AGPv4.39.chr_exon.gff3'
mytranscriptsdict = parsegff3(gff3fname)

#only run if above are CDS, else utrs are not needed for exons since they include UTRs
gff35prutrfname = 'Zea_mays.AGPv4.39.chr_five_prime_UTR.gff3'
gff33prutrfname = 'Zea_mays.AGPv4.39.chr_three_prime_UTR.gff3'
mytranscriptsdict = aggregatetranscript(gff35prutrfname,mytranscriptsdict)
mytranscriptsdict = aggregatetranscript(gff33prutrfname,mytranscriptsdict)

counter = 0
mykeys = list()
for i in mytranscriptsdict:
    mykeys.append(i)
    counter += 1
    if counter > 6:
            break

mylist = mytranscriptsdict[mykeys[2]][4]
mylist = mytranscriptsdict[mykeys[3]][4]
mylistr = mytranscriptsdict[mykeys[0]][4]
mylistr = mytranscriptsdict[mykeys[1]][4]
mylistr = mytranscriptsdict[mykeys[4]][4]
mylistr = mytranscriptsdict[mykeys[5]][4]
mylistr = mytranscriptsdict[mykeys[6]][4]

j = 0
for i in range(0,len(mylistr)):
    if mylistr[i][0] < j:
            print 'error ' + str(i) + ' not > ' + str(j)
    j = i


#sorting
#sorted(mylistr, key=lambda l:l[0]) #sort 2d list ascending based on key | [[KEYi],[INTi]]
#sorted(mylistr, key=lambda l:l[0], reverse=True) #sort 2d list descending based on key | [[KEYi],[INTi]]
#sorted(mylistr, key=lambda l:l[1]) #sort 2d list ascending based on key | [[INTi],[KEYi]]
#sorted(mylistr, key=lambda l:l[1], reverse=True) #sort 2d list descending based on key | [[INTi],[KEYi]]

####My Main!!!
if __name__ == '__main__':
    #gff3fname = 'Zea_mays.AGPv4.39.chr_exon.gff3'
    #genome = 'Zea_mays.AGPv4.dna.toplevel.fa'
    #outdir = 'contigs/'
    #threads = 11
    ###fix because it fails when dir exists
    if(os.system('mkdir -p' + outdir)):
        print 'Cannot create ' + outdir + ' directory.'
        quit()
    else:
        print 'Directory' + outdir + ' created.'
    mytranscriptsdict = parsegff3(gff3fname)
    mytranscriptsdict = sorttranscripts(mytranscriptsdict)
#key = mykeys[0]
    pool = mp.Pool(processes=threads)
    for key in mytranscriptsdict:
        pool.apply_async(writefasta, args=(key, genome, outdir))
    pool.close()
    pool.join()

#when extracting sequences, sort and extract right strand
def intsort(mylist):
    mymaxbase10power = int(math.ceil(math.log10(max(max(mylist)))))
    mybase10intcounts = list()
    for base in range(0,mymaxbase10power):
        #mybase10intcounts.append([0] * 10) #create 2d list/array
        templist = list()
        for index in range(0,10):
            templist.append([])
        mybase10intcounts.append(templist)
    for i in mylist:
        #myintcounts[i] += 1
        #myint%(10**8)/(10**7) #where 8 = mybase10 and 7 = mybase10-1
        myint = 0
        #for base in range(mymaxbase10power*-1,1):
        for base in range(0,mymaxbase10power):
            if base == 0:
                mybase10intcounts[base][myint%(10 ** base+1)]
            else:
                mybase10intcounts[base][myint%(10 ** base+1) / (10 ** base+1 )]  # where 8 = mybase10 and 7 = mybase10-1

### parsing through lengths and sorting
def writetrlengths(outfile):
    fout = open(outfile,'w')
    for tr in mytrlengthslist:
        twochars = tr[0][0:2]
        if twochars == 'Zm' or twochars == 'GR':
            fout.write(tr[0] + '_' + tr[1] + '.fasta\t' + str(tr[2]) + '\n')
        elif twochars == 'EN':
            fout.write(tr[0] + '-' + tr[1] + '.fasta\t' + str(tr[2]) + '\n')
        elif twochars == 'zm':
            fout.write(tr[0] + tr[1] + '.fasta\t' + str(tr[2]) + '\n')
    fout.close()
#testid = zma-mir396g
outfile = 'contigs/transcripts.sorted.lengths'
sortfile = 'contigs/transcripts.lengths'
mytrlendict = dict()

def createsortedtranscripts(outfile,sortfile):
    for line in open(sortfile):
        line = line.strip()
        twochars = line[0:2]
        if twochars == 'Zm' or twochars == 'GR':
            line = line.split('_')
            key = line[0]
            transcriptid = line[1]
            trlength = int(line[len(line) - 1].split()[1])
        elif twochars == 'EN':
            line = line.split('-')
            key = line[0]
            transcriptid = line[1].split('_')[0]
            trlength = int(line[len(line) - 1].split()[1])
        elif twochars == 'zm':
            line = line.split('_')
            key = line[0][0:len(line[0]) - 1]
            transcriptid = line[0][len(line[0]) - 1]
            trlength = int(line[len(line) - 1].split()[1])
        if key in mytrlendict.keys():
                mytrlendict[key].append([transcriptid,trlength])
        else:
                mytrlendict[key] = [[transcriptid,trlength]]

len(mytrlendict)
for key in mytrlendict:
    mytrlendict[key] = sorted(mytrlendict[key], key=lambda l:l[1], reverse=True)

counter = 0
mykeys = list()
for key in mytrlendict:
    mykeys.append(key)
    counter += 1
    if counter > 6:
            break

for key in mykeys:
    print mytrlendict[key]

mytrlengthslist = list()
for key in mytrlendict:
    mytrlengthslist.append([key,mytrlendict[key][0][0],mytrlendict[key][0][1]])

#single threaded
writetrlengths(outfile)
#multithreaded
pool = mp.Pool(processes=threads)
for key in mytranscriptsdict:
    pool.apply_async(writetrlengths, args=(outfile))
pool.close()
pool.join()


#fout = open(outfile,'w')
#for tr in mytrlengthslist:
#    fout.write(tr[0] + '_' + tr[1] + '.fasta\t' + str(tr[2]) + '\n')

#fout.close()

