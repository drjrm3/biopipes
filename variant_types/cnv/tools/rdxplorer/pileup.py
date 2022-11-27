#!/usr/bin/env python

import v6
import globals as glob
import file_utils as fu
import AnalyzeSequence as ans
import os
import numpy as np




#############################################################################################


def qpileup(path2bam, reference, chromOfInterest='All', qual=20, wrkgdir='.', gender='M', debug=False, delete=True):

    print ("Quality " + str(qual) + " pileup is being generated")

    if fu.isExist(path2bam) == False:
        print ('File does not exists: ' + path2bam + ' ' )
    elif fu.isExist(reference) == False:
        print ('File does not exists: ' + reference + ' ' )
    else:

        q_mapped_sam = wrkgdir+'/'+'q' + str(qual)+ '.mapped.sam'
        q_pileup = wrkgdir+'/'+ 'q' + str(qual)+ '.pileup'
        q_d10_pileup  = wrkgdir+'/'+ 'q' + str(qual)+ '.d10.pileup'


        com = ''


        if qual==0:
            if chromOfInterest.lower()=='All'.lower():
                com = glob.SAMTOOLS_PATH + '/samtools view -F 0x400 ' + path2bam + ' | awk \'$3!="*" && $5==' + str(qual) + '\'  > ' + q_mapped_sam
            else:
                com = glob.SAMTOOLS_PATH + '/samtools view -F 0x400 ' + path2bam + ' | awk \'$3!="*" && $5==' + str(qual) + '\' | awk \'$3=="' + chromOfInterest + '" || $3=="chr' + chromOfInterest + '"\'  > ' + q_mapped_sam
        else:
            if chromOfInterest.lower()=='All'.lower():
                com = glob.SAMTOOLS_PATH + '/samtools view -F 0x400 ' + path2bam + ' | awk \'$3!="*" && $5>=' + str(qual) + '\'  > ' + q_mapped_sam
            else:
                com = glob.SAMTOOLS_PATH + '/samtools view -F 0x400 ' + path2bam + ' | awk \'$3!="*" && $5>=' + str(qual)  + '\' | awk \'$3=="' + chromOfInterest + '" || $3=="chr' + chromOfInterest + '"\'  > ' + q_mapped_sam



        v6.execute(com, debug)

        if fu.isExist(q_mapped_sam):
            com = 'wc -l ' + q_mapped_sam + ' > ' + q_mapped_sam + '.log'
            v6.execute(com, debug)
            com = glob.SAMTOOLS_PATH + '/samtools pileup -Sf ' + reference + '  ' + q_mapped_sam + ' > ' + q_pileup
            v6.execute(com, debug)
        else:
            print ('File does not exists: ' + q_mapped_sam + ' ' )


        if fu.isExist(q_pileup):
            com = 'wc -l ' + q_pileup + ' > ' + q_pileup + '.log'
            v6.execute(com, debug)
        else:
            print ('File does not exists: ' + q_pileup + ' ' )

        if delete==True:
            fu.delete(q_mapped_sam)
            fu.delete(q_d10_pileup)



#############################################################################################


def depth(wrkgdir='.', chromOfInterest='All', hg='hg19', winSize=100, gender='M', sumWithZero=True, debug=False, delete=True):

    chromosomes = v6.generateListOfChrom(gender=gender, chromOfInterest=chromOfInterest)
    print('\t'.join(chromosomes))

    pileups = ["q0", "q20"]

    if sumWithZero==False:
        pileups = ["q20"]

    for p in pileups:

        v6.pileupChromsExtract(wrkgdir +'/'+ str(p)+'.pileup', wrkgdir, debug, chromOfInterest)
        print("# Created " + str(p)+'.pileup')

        for chrom in chromosomes:
            filename = wrkgdir+'/'+ chrom +'.txt'
            if (fu.isExist(filename) and fu.fileSize(filename) > 0):
                ans.base4qual(chrom=chrom.replace("chr", ""), pos_depth=filename, base_hdffile=wrkgdir+'/'+ chrom+str(p)+'base.txt', win_hdffile=wrkgdir+'/'+ chrom+str(p)+'win' + str(winSize) + '.txt', hg=hg, winSize=winSize, compress=False, sep="\t", debug=debug)
                print(chrom +str(p)+'base.txt  generated' )
                if delete==True:
                    fu.delete(filename)
                    print("Deleting " + wrkgdir+'/'+ chrom +'.txt' )
            else:
                print ("File " + filename + ' does not exist or zero size')

        if delete==True:
            fu.delete(wrkgdir +'/'+ str(p)+'.pileup')




def mergeWin(wrkgdir='.', chromOfInterest='All', winSize=100, gender='M', debug=False, delete=True):
    chromosomes=v6.generateListOfChrom(gender=gender, chromOfInterest=chromOfInterest)
    #print(chromosomes)
    pileups = ["q0", "q20"]
    for chrom in chromosomes:
        print(str(chrom))
        win_hdffile1=wrkgdir+'/'+ chrom+str(pileups[0])+'win' + str(winSize) + '.txt'
        win_hdffile2=wrkgdir+'/'+ chrom+str(pileups[1])+'win' + str(winSize) + '.txt'
        win_hdffile=wrkgdir+'/'+ chrom+'Sum'+str(pileups[0])+str(pileups[1])+'win'+str(winSize)+'.txt'


        if (fu.isExist(win_hdffile1) and fu.fileSize(win_hdffile1) > 0) and (fu.isExist(win_hdffile2) and fu.fileSize(win_hdffile2) > 0):
            w0=np.array(fu.read_one_float_col(filename=win_hdffile1))
            w20=np.array(fu.read_one_float_col(filename=win_hdffile2))
            w=w0+w20
            fu.save2txt(read_data=w, txtfile=win_hdffile, compress=False, debug=debug)

        elif (fu.isExist(win_hdffile2) and fu.fileSize(win_hdffile2) > 0):
            w20=np.array(fu.read_one_float_col(filename=win_hdffile2))
            w=w20
            fu.save2txt(read_data=w, txtfile=win_hdffile, compress=False, debug=debug)

        if delete==True:
            print('File' + str(win_hdffile1) + ' was kept' )
            print('File' + str(win_hdffile2) + ' was kept' )
            #fu.delete(win_hdffile1)
            #fu.delete(win_hdffile2)
