#!/bin/bash

# --- Configuration ---
DUMP_DIR="/path/to/dump/dir"    # directory containing dump files
DUMP_GLOB="dump.*"              # adjust to match your dump file naming
OUTPUT_DIR="$DUMP_DIR/voxel_output"
FINAL_H5="$DUMP_DIR/trajectory.h5"
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

N=$(ls "$DUMP_DIR"/$DUMP_GLOB | wc -l)
ARRAY_MAX=$((N - 1))

echo "Found $N dump files in $DUMP_DIR, submitting array 0-${ARRAY_MAX}"

ARRAY_JOB=$(sbatch --parsable \
    --array=0-${ARRAY_MAX}%20 \
    --export=DUMP_DIR=$DUMP_DIR,DUMP_GLOB=$DUMP_GLOB,OUTPUT_DIR=$OUTPUT_DIR,SCRIPT_DIR=$SCRIPT_DIR \
    "$SCRIPT_DIR/voxel_analysis.slurm")

echo "Array job ID: $ARRAY_JOB"

MERGE_JOB=$(sbatch --parsable \
    --dependency=afterok:$ARRAY_JOB \
    --export=OUTPUT_DIR=$OUTPUT_DIR,FINAL_H5=$FINAL_H5,SCRIPT_DIR=$SCRIPT_DIR \
    "$SCRIPT_DIR/merge_h5.slurm")

echo "Merge job ID: $MERGE_JOB (runs after array completes)"
