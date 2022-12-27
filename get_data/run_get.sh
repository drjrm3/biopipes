#!/usr/bin/env bash

BIOPROJECT=PRJEB35986
REGEX=ERR

esearch -db sra -query $BIOPROJECT | efetch -format runinfo > $BIOPROJECT.txt

cut -d ',' -f 1 $BIOPROJECT.txt | grep $REGEX | while read ACCESSION; do
    #echo "Copying accession $ACCESSION."
    echo "fastq-dump --split-3 $ACCESSION -O $BIOPROJECT" #> $ACCESSION.fastq
done | parallel -j 4
