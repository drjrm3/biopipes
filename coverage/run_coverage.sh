#!/usr/bin/env bash

time ./coverage.sh /nvme/SRS290927.bam create_beds/hg38_chrlens.txt 100 SRS290927.coverage

awk '{print $1, $2, $3, $4, $4/($3-$2)}' SRS290927.bedcov | tr ' ' '\t' > SRS290927.bedcov2

samtools coverage  -r chr6 /nvme/SRS290927.bam

# meandepth
grep chr6 SRS290927.bedcov2 | awk '{sum += $NF}END{print sum/NR}'

# coverage 
grep chr6 SRS290927.bedcov2 | awk '{if($4 != 0)sum += 1}END{print sum/NR}'
