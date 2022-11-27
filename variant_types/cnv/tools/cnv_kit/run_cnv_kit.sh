#!/usr/bin/env bash

#pip3 install cnvkit

#cnvkit.py --version

#cnvkit.py batch -h

cnvkit.py batch \
    -m wgs \
    -n \
    --fasta /scratch/ref/reffa/hg38/hg38.fa \
    /nvme/SRS290927.bam \

