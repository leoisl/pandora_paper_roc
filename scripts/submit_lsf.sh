#!/usr/bin/env bash

MEMORY=1000000
LOCAL_CORES=110
cluster="bigmem"

PROFILE="lsf"
LOG_DIR=logs/
JOB_NAME="snakemake_master_process."$(date --iso-8601='minutes')

mkdir -p $LOG_DIR

bsub -R "select[mem>$MEMORY] rusage[mem=$MEMORY] span[hosts=1]" \
    -n "$LOCAL_CORES" \
    -M "$MEMORY" \
    -o "$LOG_DIR"/"$JOB_NAME".o \
    -e "$LOG_DIR"/"$JOB_NAME".e \
    -J "$JOB_NAME" \
    -P "$cluster" \
      snakemake --local-cores "$LOCAL_CORES" --profile "$PROFILE" --stats "$LOG_DIR"/snakemake_stats "$@"

exit 0
