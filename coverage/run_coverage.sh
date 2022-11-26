#!/usr/bin/env bash

#python3 coverage.py -b /nvme/SRS290927.bam create_beds/hg38_chrlens.txt 100 SRS290927.coverage
python3 coverage.py \
    --input-bamfile /nvme/SRS290927.bam \
    --output-coverage-file SRS290927.coverage \
    --proces 2
    #--chrom chr1
    #--chr-lengths create_beds/hg38_chrlens.txt

#samtools coverage  -r chr6 /nvme/SRS290927.bam

# meandepth
#grep chr6 SRS290927.bedcov2 | awk '{sum += $NF}END{print sum/NR}'

# coverage 
#grep chr6 SRS290927.bedcov2 | awk '{if($4 != 0)sum += 1}END{print sum/NR}'
