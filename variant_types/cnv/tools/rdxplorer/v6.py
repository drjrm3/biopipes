#!/usr/bin/env python

import os
import os.path
import globals as glob
import file_utils as fu
import sys
import shutil
import AnalyzeSequence as ans


def execute(com, debug=False):
    if debug==True:
        print (com)
    os.system(com)

""" Extract One Chromosome """
def pileupChromsExtractOne(pileup, outdir, chromOfInterest, debug=False):
    #chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20","21","22", "X", "Y"]
    com = 'awk \'{ if($1=="' + chromOfInterest + '" || $1=="chr'+chromOfInterest+'" )  printf "%s\\t%s\\n", $2,$4}\'' + ' ' + pileup +' > ' +  outdir + '/' + 'chr' + str(chromOfInterest) + '.txt'
    execute(com, debug)

""" Extract All Chromosomes in the BAM file """
def pileupChromsExtractMany(pileup, outdir, debug=False):

    com = 'awk \'{ if($1=="1" || $1=="chr1" )  printf "%s\\t%s\\n", $2,$4  > "' +  outdir + '/' + 'chr1.txt"  \n'\
    ' else if($1=="2" || $1=="chr2" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr2.txt"  \n'\
    ' else if($1=="3" || $1=="chr3" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr3.txt"  \n'\
    ' else if($1=="4" || $1=="chr4" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr4.txt"  \n'\
    ' else if($1=="5" || $1=="chr5" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr5.txt"  \n'\
    ' else if($1=="6" || $1=="chr6" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr6.txt"  \n'\
    ' else if($1=="7" || $1=="chr7" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr7.txt"  \n'\
    ' else if($1=="8" || $1=="chr8" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr8.txt"  \n'\
    ' else if($1=="9" || $1=="chr9" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr9.txt"  \n'\
    ' else if($1=="10" || $1=="chr10" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr10.txt"  \n'\
    ' else if($1=="11" || $1=="chr11" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr11.txt"  \n'\
    ' else if($1=="12" || $1=="chr12" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr12.txt"  \n'\
    ' else if($1=="13" || $1=="chr13" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr13.txt"  \n'\
    ' else if($1=="14" || $1=="chr14" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr14.txt"  \n'\
    ' else if($1=="15" || $1=="chr15" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr15.txt"  \n'\
    ' else if($1=="16" || $1=="chr16" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr16.txt"  \n'\
    ' else if($1=="17" || $1=="chr17" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr17.txt"  \n'\
    ' else if($1=="18" || $1=="chr18" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr18.txt"  \n'\
    ' else if($1=="19" || $1=="chr19" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr19.txt"  \n'\
    ' else if($1=="20" || $1=="chr20" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr20.txt"  \n'\
    ' else if($1=="21" || $1=="chr21" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr21.txt"  \n'\
    ' else if($1=="22" || $1=="chr22" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chr22.txt"  \n'\
    ' else if($1=="X" || $1=="chrX" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chrX.txt"  \n'\
    ' else if($1=="Y" || $1=="chrY" )  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chrY.txt"  \n'\
    ' else  printf "%s\\t%s\\n", $2,$4 > "' +  outdir + '/' + 'chrO.txt"  }\' ' + pileup


    execute(com, debug)
    fu.delete(filename=outdir + '/' + 'chrO.txt')







""" Generate list of chromosomes """
def generateList(gender='M'):

    lastChrom = 23
    if str(gender) == 'M':
        lastChrom = 24

    chromosomes = []
    for i in range(1, (lastChrom+1)):
        chr = "chr"
        if i == 23:
            chr = chr+'X'
        elif i == 24:
            chr = chr+'Y'
        else:
            chr = chr+str(i)
        chromosomes.append(chr)

    return chromosomes

#########################################################################################################
#########################################################################################################

###     Call methods below from outside of this script

#########################################################################################################
#########################################################################################################


""" Extract One or All Chromosomes, as specified """
""" This method is used outside of the script """
def pileupChromsExtract(pileup, outdir, debug=False, chromOfInterest="All"):
    if chromOfInterest.lower()=='All'.lower():
        pileupChromsExtractMany(pileup=pileup, outdir=outdir, debug=debug)
    else:
        pileupChromsExtractOne(pileup=pileup, outdir=outdir, chromOfInterest=chromOfInterest, debug=debug)


""" Generate list of chromosomes as specified """
""" If one chromosome is specified, returns list with one element """
""" This method is used outside of the script """
def generateListOfChrom(gender='M', chromOfInterest='All'):
    chromosomes=[]
    if chromOfInterest.lower()=='All'.lower():
        chromosomes = generateList(gender=gender)
    else:
        chromosomes.append('chr'+str(chromOfInterest))

    return chromosomes
