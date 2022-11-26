#!/usr/bin/env bash

time ./coverage.sh -i /scratch/genomics/coverage/SRS290927.bam -o /scratch/genomics/coverage/SRS290927.cov #-s 4 -p 4

#samtools coverage  -r chr6 /nvme/SRS290927.bam

# meandepth
#grep chr6 SRS290927.bedcov2 | awk '{sum += $NF}END{print sum/NR}'

# coverage 
#grep chr6 SRS290927.bedcov2 | awk '{if($4 != 0)sum += 1}END{print sum/NR}'
