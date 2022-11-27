#!/usr/bin/env bash

# Default values.
BAMFILE=""
BEDFILE=""
COVFILE=""
PROCS=$(cat /proc/cpuinfo | grep proc | wc -l | awk '{print $1/2}')
CHROM_FILTERS="chrUn _alt _random"

# Usage.
usage() {
cat << EOF
Usage: coverage.sh -i | --input-bam  BAMFILE
                   -b | --bed-file   BEDFILE
                   -o | --output-cov COVFILE
                  [-p | --processes PROCS (default: $PROCS)]"
                  [-f | --filters CHROM_FILTERS (default: "$CHROM_FILTERS")]
                  [-h | --help]
EOF
}

# Get parse arguments and validate them.
PARSED_ARGUMENTS=$(getopt -a \
    -n coverage.sh \
    -o i:b:o:p:h \
    --long input-bam:,bed-file:,output-coverage:,processes:,help \
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
        -b | --bed-file)         BEDFILE="$2"        ; shift 2 ;;
        -o | --output-cov)       COVFILE="$2"        ; shift 2 ;;
        -f | --filters)          CHROM_FILTERS="$2"  ; shift 2 ;; 
        -p | --processes)        PROCS="$2"          ; shift 2 ;; 
        -h | --help) usage && exit 1 ;;
        --) shift; break ;;
        *) echo "Unexpected option: $1 - this should not happen."
           usage ;;
    esac
done

# Check inputs.
[ -f $BAMFILE ] || { >&2 echo "ERROR: Input BAMFILE not given." && usage && exit 1; }
[ -f $BEDFILE ] || { >&2 echo "ERROR: Input BEDFILE not given." && usage && exit 1; }
[ $COVFILE ]    || { >&2 echo "ERROR: COVFILE not given."       && usage && exit 1; }

# Process chroms.
samtools bedcov $BEDFILE $BAMFILE > $COVFILE

