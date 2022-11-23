#!/usr/bin/env bash
# TODO: Move this to python?

TMPDIR=$(mktemp -d)
echo $TMPDIR
set -e

NPROCS=16
NJOBS=128

BAM=$1
CHRLENS=$2
BIN=$3
OUT=$4

echo -n "[*] Generating bed file ... "
./create_beds/create_bed.py -i $CHRLENS -o out_${BIN}.bed -b $BIN
echo "Done"

NLINES_TOTAL=$(wc -l out_${BIN}.bed | awk '{print $1}')
NLINES_PER_BED=$(( $NLINES_TOTAL / $NJOBS ))

split -l $NLINES_PER_BED -d -a6 out_${BIN}.bed $TMPDIR/bed_split.

for BED in $TMPDIR/bed_split.*; do
    IDX=$(echo $BED | awk -F'.' '{print $(NF)}')
    echo "samtools bedcov $BED $BAM > $TMPDIR/tmp.$IDX.bedcov"
done > cmds.txt
parallel -j $NPROCS < cmds.txt
cat $TMPDIR/tmp.*.bedcov > $(basename $BAM .bam).bedcov

rm -rf $TMPDIR

