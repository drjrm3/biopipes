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




""" This mode is Not Yet Implemented """
def usage_cli():
    print('\n')
    print('...................RDXplorer arguments ...................')
    print('\n')
    print('\t"-b",  "--bam" - input BAM file.  Default: None')
    print('\t"-r", "--reference - reference fasta file.  Default: None')
    print('\t"-w", "--wrkgdir"  - Output directory.  Default: directory where the bam file is')
    print('\t"-g", "--gender"  - Gender Default: M')
    print('\t"--hg"  - Human Genome build.  Default: HG19')
    print('\t"--winSize" - window size.  Default: 100')
    print('\t"--baseCopy" - Base copy number.  Default: 2')
    print('\t"-f", "--filter" - Filter Summary Table.  Default: 10 win sizes, 1000 bp for 100bp Win')
    print('\t"-z" - Sum q0 and q20 depth.  Default: True')
    print('\t"-s" - Silent.  Default: True( Print output) ')
    print('\t"-k" - Keep.  Default: True( Delete tmp files) ')
    print('\t"-h", " - Help. Print this message ')
    print('\n')
    print('\tFirst two arguments are required, third is highly recommended' )
    print('\n')
    print('.......................................................................')

""" This mode is Not Yet Implemented """
def usage():
    print('\n')
    print('...................RDXplorer arguments ...................')
    print('\n')
    print('\t - input BAM file. Type: String. Default: None')
    print('\t - reference fasta file. Type: String. Default: None')
    print('\t - Output directory. Type: String. Default: directory where the bam file is')
    print('\t - Gender. Type: String.  Default: M')
    print('\t - Human Genome build. Type: String.  Default: HG19')
    print('\t - Window size. Type: Integer.  Default: 100')
    print('\t - Base copy number. Type: Integer.  Default: 2')
    print('\t - Filter Summary Table.  Default: 10 win sizes, 1000 bp for 100bp Win')
    print('\t - Sum q0 and q20 depth while calculating windows average. Type: Boolead. Default: True. If "False", only q20 is used')
    print('\t - Debug mode. Type: Boolean.  Default: True( Print output) ')
    print('\t - Delete temporary files. Type: Boolean. Default: True( Delete tmp files) ')
    print('\n')
    print('.......................................................................')



def complainAndBail():
    if glob.APP_HOME == '' or glob.APP_HOME is None:
        print ("APP_HOME is not set. Please configure the globals.py")
        return False
    elif glob.SAMTOOLS_PATH == '' or glob.SAMTOOLS_PATH is None:
        print ("SAMTOOLS_PATH is not set. Please configure the globals.py")
        return False
    elif fu.isExist(glob.SAMTOOLS_PATH+'/'+'samtools') == False:
        print ("Samtools is essential. Please install it and/or configure globals.py")
        return False
    else:
        return True

def init(path2bam, reference, wrkgdir, chromOfInterest='All', gender='M', hg='hg19', winSize=100, baseCopy=2, filter=10, sumWithZero=True, debug=True, delete=True):

    if path2bam is None:
        usage()
        sys.exit(2)
    elif reference is None:
        usage()
        sys.exit(2)
    elif wrkgdir is None:
        wrkgdir = os.path.dirname(path2bam)
    else:
        if wrkgdir != '.':
            fu.mkdirp(wrkgdir)

        #generate q0 pileup if sumWithZero is True

        if sumWithZero==True:
            plp.qpileup(path2bam, reference=reference, chromOfInterest=chromOfInterest, qual=0, wrkgdir=wrkgdir,  gender=gender, debug=debug, delete=delete)
        plp.qpileup(path2bam, reference=reference, chromOfInterest=chromOfInterest, qual=20, wrkgdir=wrkgdir,  gender=gender, debug=debug, delete=delete)
        plp.depth(wrkgdir=wrkgdir, chromOfInterest=chromOfInterest, hg=hg, winSize=winSize, gender=gender, sumWithZero=sumWithZero, debug=debug, delete=delete)

        #merge, if sumWithZero is True
        if sumWithZero==True:
            plp.mergeWin(wrkgdir=wrkgdir, chromOfInterest=chromOfInterest, winSize=winSize, gender=gender, debug=debug, delete=delete)


        ro.r('source("' + glob.APP_HOME + '/rlib/ewtMain.R")')
        runEwtMain = ro.r(['runEwtMain'])
        chromosomes = v6.generateListOfChrom(gender=gender, chromOfInterest=chromOfInterest)
        pileups = ["q0", "q20"]


        for chrom in chromosomes:
            countwinfile=chrom+'Sum'+str(pileups[0])+str(pileups[1])+'win'+str(winSize)+'.txt'
            if sumWithZero == False:
                countwinfile=chrom+str(pileups[1])+'win' + str(winSize) + '.txt'
            if gender=='M' and (chrom=='Y' or chrom=='X'):
                baseCopy = 1
            if fu.isExist(wrkgdir + '/' + countwinfile) and fu.fileSize(wrkgdir + '/' + countwinfile) > 0:
                print(str(chrom) + " " + countwinfile)
                #runEwtMain(apphome=glob.APP_HOME, wrkgdir=wrkgdir, countwinfile=countwinfile, chr=chrom, hg=hg, gender=gender, filter=filter, baseCopy=baseCopy, dirsep = "/")
                runEwtMain(apphome=glob.APP_HOME, wrkgdir=wrkgdir, countwinfile=countwinfile, chr=chrom, hg=hg, win=winSize, gender=gender, filter=filter, baseCopy=baseCopy, dirsep = "/")
                com='gzip -f ' + wrkgdir + '/' + chrom + '.nzd'
                v6.execute(com)
                if delete==True:
                    #fu.delete(wrkgdir + '/' + countwinfile)
                    print('File' + str(countwinfile) + ' was kept' )
                    fu.delete(wrkgdir + '/' + chrom + '.gcc')
                else:
                    com='gzip -f ' + wrkgdir + '/' + chrom + '.gcc'
                    v6.execute(com)
