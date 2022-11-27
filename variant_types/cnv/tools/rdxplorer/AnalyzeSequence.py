#!/usr/bin/env python
import numpy
import os
import globals as glob
import file_utils as fu

chromSize18 = [

    #HG18
    247249719, #chr1
    242951149, #chr2
    199501827, #chr3
    191273063, #chr4
    180857866, #chr5
    170899992, #chr6
    158821424, #chr7
    146274826, #chr8
    140273252, #chr9
    135374737, #chr10
    134452384, #chr11
    132349534, #chr12
    114142980, #chr13
    106368585, #chr14
    100338915, #chr15
    88827254, #chr16
    78774742, #chr17
    76117153, #chr18
    63811651, #chr19
    62435964, #chr20
    46944323, #chr21
    49691432, #chr22
    154913754, #chrX
    57772954   #chrY

    ]


chromSize19 = [

    #HG19
    249250621,#1
    243199373,#2
    198022430,#3
    191154276,#4
    180915260,#5
    171115067,#6
    159138663,#7
    146364022,#8
    141213431,#9
    135534747,#10
    135006516,#11
    133851895,#12
    115169878,#13
    107349540,#14
    102531392,#15
    90354753,#16
    81195210,#17
    78077248,#18
    59128983,#19
    63025520,#20
    48129895,#21
    51304566,#22
    155270560,#23
    59373566#24

    ]

def initChrom(chrom, hg='hg19'):
    """ chrom must be in form 1 - 22, X or Y"""
    c=0
    if chrom == "X":
        c = 23
    elif chrom == "Y":
        c = 24
    else:
        c = int(chrom)

    base = 0
    if hg=='hg19':
        base = chromSize19[c-1]
    else:
        base = chromSize18[c-1]

    print (' Base for chrom ' + str (chrom) + ' for ' +hg)

    return numpy.zeros((base,), dtype=numpy.int)


def base4qual(chrom, pos_depth, base_hdffile, win_hdffile, hg='hg19', winSize=100, compress=True, sep="\t", debug=False):

    """ pos_depth  is a file with 2 columns - genome position and depth of coverage"""

    b = initChrom(chrom, hg)
    if debug == True:
      print("b " + str(len(b)) )
    fh = open(pos_depth, "r")
    for line in fh:
        line = line.strip('\r\n')
        fields = line.strip().split(sep)
        position = int(fields[0])
        depth = int(fields[1])
        if position <=len(b):
            b[position-1] = depth

    print (' Base calculated')

    ####
    win=[]
    numOfWin=len(b)/winSize
    for i in range(0,numOfWin):
        sum=0
        for j in range(0, winSize):
            sum=sum+b[i*winSize+j]
        win.append(sum)

    print (' Window calculated')
    ####
    #fu.save2txt(read_data=b, txtfile=base_hdffile, compress=False, debug=debug)
    fu.save2txt(read_data=win, txtfile=win_hdffile, compress=False, debug=debug)

