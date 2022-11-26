#!/usr/bin/env bash

# Default values.
BAMFILE=""
COVFILE=""
STHREADS=2
PROCS=$(cat /proc/cpuinfo | grep proc | wc -l | awk '{print $1/2}')
CHROM_FILTERS="chrUn _alt _random"

# Usage.
usage() {
cat << EOF
Usage: coverage.sh -i | --input-bam  BAMFILE
                   -o | --output-cov COVFILE
                  [-s | --samtools-threads STHREADS (default: $STHREADS)] 
                  [-p | --processes PROCS (default: $PROCS)]"
                  [-f | --filters CHROM_FILTERS (default: "$CHROM_FILTERS")]
                  [-h | --help]
EOF
}

# Get parse arguments and validate them.
PARSED_ARGUMENTS=$(getopt -a \
    -n coverage.sh \
    -o i:o:s:p:h \
    --long input-bam:,output-coverage:,samtools-threads:,processes:,help \
    -- "$@")
VALID_ARGUMENTS=$?
if [ "$VALID_ARGUMENTS" != "0" ]; then
    usage && exit 1
fi


# Read arguments.
eval set -- "$PARSED_ARGUMENTS"
while [[ ! -z $1 ]]; do
    case "$1" in
        -i | --input-bam)        BAMFILE="$2"        ; shift 2 ;;
        -o | --output-cov)       COVFILE="$2"        ; shift 2 ;;
        -s | --samtools-threads) STHREADS="$2"       ; shift 2 ;; 
        -f | --filters)          CHROM_FILTERS="$2" ; shift 2 ;; 
        -p | --processes)        PROCS="$2"          ; shift   ;; 
        -h | --help) usage && exit 1 ;;
        --) shift; break ;;
        *) echo "Unexpected option: $1 - this should not happen."
           usage ;;
    esac
done

# Check inputs.
[ -f $BAMFILE ] || { >&2 echo "ERROR: Input BAMFILE not given." && usage && exit 1; }
[ $COVFILE ]    || { >&2 echo "ERROR: COVFILE not given."       && usage && exit 1; }

# Process chroms.
while read CHROM; do
    FILTER=""
    for CHROM_FILTER in $CHROM_FILTERS; do
        if [[ "$CHROM" == *"$CHROM_FILTER"* ]]; then FILTER="TRUE"; fi
    done
    if [[ $FILTER ]]; then continue; fi

    echo "samtools depth -@ $STHREADS -r $CHROM $BAMFILE > $COVFILE.$CHROM.tmp"
done < <(samtools view -H $BAMFILE | grep '@SQ' | awk '{print $2}' | sed 's/SN://g') > cmds.txt


