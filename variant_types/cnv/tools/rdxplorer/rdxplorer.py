#!/usr/bin/env python

import pileup as plp
import v6
import globals as glob
import file_utils as fu
import AnalyzeSequence as ans
import os
import sys
import rpy2.robjects as ro
import rpy2.robjects.numpy2ri
import rdxplorer_api as rdxp

if int(len(sys.argv)) < 12:
    rdxp.usage()

else:
    path2bam=sys.argv[1]
    reference=sys.argv[2]
    wrkgdir=sys.argv[3]
    chromOfInterest=sys.argv[4]

    gender=sys.argv[5]
    hg=sys.argv[6]
    winSize=sys.argv[7]
    baseCopy=sys.argv[8]
    filter=sys.argv[9]
    sumWithZero=sys.argv[10]
    debug=sys.argv[11]
    delete=sys.argv[12]

    debug=fu.str2bool(debug)
    delete=fu.str2bool(delete)
    sumWithZero=fu.str2bool(sumWithZero)
    baseCopy=int(baseCopy)
    winSize=int(winSize)
    filter=int(filter)

    if rdxp.complainAndBail() == True:
        if debug==True:
            print("The following arguments have been accepted:")
            a=0
            for arg in sys.argv:
                if a==0:
                    print("Program: " + arg)
                elif a==1:
                    print("Bam file name: " + arg)
                else:
                    print ('\t' + arg)
                a=a+1

        accepted_chromosomes = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20","21","22", "X", "Y"]
        if fu.find_first_index(accepted_chromosomes, chromOfInterest) < 0:
            chromOfInterest = 'All'


        rdxp.init(path2bam, reference, wrkgdir, chromOfInterest, gender, hg, winSize, baseCopy, filter, sumWithZero, debug, delete)
        
