#!/bin/bash

##### these parameters need to be changed


path2bam='/arch//tmp/Example.bam'
reference='/work/bio/references/b37/hg19.fa'
wrkgdir='/arch/tmp/Example'
chromOfInterest='17'



#### these parameters may be kept default. See readme file for explanation

gender='F'
hg='hg19'
winSize=100
baseCopy=2
filter=10
sumWithZero=True
debug=True
delete=True

BEFORE=`date '+%s'`

python rdxplorer.py ${path2bam} ${reference} ${wrkgdir} ${chromOfInterest} ${gender} ${hg} ${winSize} ${baseCopy} ${filter} ${sumWithZero} ${debug} ${delete}

AFTER=`date '+%s'`
TIME=$(($AFTER - $BEFORE))
echo "Doing this took $TIME seconds "
echo "DONE."
