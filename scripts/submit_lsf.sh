#!/usr/bin/env bash

MEMORY=32000
LOCAL_CORES=20

PROFILE="lsf"
LOG_DIR=logs/
JOB_NAME="snakemake_master_process."$(date --iso-8601='minutes')

echo "Have you remembered doing conda activate basic_bioinfo and pipenv shell?"

mkdir -p $LOG_DIR

bsub -R "select[mem>$MEMORY] rusage[mem=$MEMORY] span[hosts=1]" \
    -n "$LOCAL_CORES" \
    -M "$MEMORY" \
    -o "$LOG_DIR"/"$JOB_NAME".o \
    -e "$LOG_DIR"/"$JOB_NAME".e \
    -J "$JOB_NAME" \
      snakemake --local-cores "$LOCAL_CORES" --profile "$PROFILE" --stats "$LOG_DIR"/snakemake_stats "$@"

exit 0
