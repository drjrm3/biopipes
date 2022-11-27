#!/bin/bash

##### these parameters need to be changed

path2bam='/share/ngs_data/mssm/WG'
reference='/work/bio/resources/human_g1k_v37.fasta'
wrkgdir='/share/ngs_data/mssm/WG/rdx'
chromOfInterest='All'



#### these parameters may be kept default. See readme file for explanation

gender='M'
hg='hg19'
winSize=100
baseCopy=2
filter=10
sumWithZero=True
debug=True
delete=True

for index in AC724_102_WG_child  AC724_101_WG_mother  AC724_201_WG_father
do
 BEFORE=`date '+%s'`
 bam=${path2bam}/${index}/${index}_rmdup.bam
 wrkg=${wrkgdir}/${index}
 echo "Processing $bam"
 python rdxplorer.py ${bam} ${reference} ${wrkg} ${chromOfInterest} ${gender} ${hg} ${winSize} ${baseCopy} ${filter} ${sumWithZero} ${debug} ${delete}

 AFTER=`date '+%s'`
 TIME=$(($AFTER - $BEFORE))
 echo "Doing this took $TIME seconds "
 echo "DONE $index"
done
