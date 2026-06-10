#!/bin/bash

# =============================================================================
# Dielectric constant pipeline — step 1 submission
#
# Submits a SLURM job array (one task per dump file) then chains the
# aggregate step to run automatically when all array tasks finish.
#
# Usage:
#   bash submit_pipeline.sh
#
# Edit the Configuration section below before running.
# =============================================================================

# --- Configuration ---
DUMP_DIR="/path/to/dump/dir"      # directory containing LAMMPS dump files
DUMP_GLOB="dump.*.lammpstrj"      # glob matching your dump file names
OUTPUT_DIR="$DUMP_DIR/dipole_output"

CUTOFF=1.2     # O-H bond cutoff in Angstroms
TYPE_O=1       # LAMMPS atom type for oxygen
TYPE_H=2       # LAMMPS atom type for hydrogen

VENV_PATH="/home1/lkyamamo/venv/struc_analysis"

# Maximum simultaneous array tasks (throttle to avoid overwhelming the scheduler)
MAX_SIMULTANEOUS=50

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# ---------------------

N=$(ls "$DUMP_DIR"/$DUMP_GLOB 2>/dev/null | wc -l)

if [ "$N" -eq 0 ]; then
    echo "ERROR: No files matching '$DUMP_GLOB' found in '$DUMP_DIR'" >&2
    exit 1
fi

ARRAY_MAX=$((N - 1))
echo "Found $N dump files in $DUMP_DIR"
echo "Submitting array 0-${ARRAY_MAX} (max ${MAX_SIMULTANEOUS} simultaneous)"

mkdir -p "$SCRIPT_DIR/logs"
mkdir -p "$OUTPUT_DIR"

ARRAY_JOB=$(sbatch --parsable \
    --array=0-${ARRAY_MAX}%${MAX_SIMULTANEOUS} \
    --export=DUMP_DIR="$DUMP_DIR",DUMP_GLOB="$DUMP_GLOB",OUTPUT_DIR="$OUTPUT_DIR",SCRIPT_DIR="$SCRIPT_DIR",CUTOFF="$CUTOFF",TYPE_O="$TYPE_O",TYPE_H="$TYPE_H",VENV_PATH="$VENV_PATH" \
    "$SCRIPT_DIR/1.array.slurm")

if [ -z "$ARRAY_JOB" ]; then
    echo "ERROR: Array job submission failed" >&2
    exit 1
fi

echo "Array job ID: $ARRAY_JOB"

MERGE_JOB=$(sbatch --parsable \
    --dependency=afterok:$ARRAY_JOB \
    --export=OUTPUT_DIR="$OUTPUT_DIR",SCRIPT_DIR="$SCRIPT_DIR",VENV_PATH="$VENV_PATH" \
    "$SCRIPT_DIR/1.aggregate.slurm")

echo "Aggregate job ID: $MERGE_JOB (runs after array completes)"
echo ""
echo "Monitor with: squeue -j $ARRAY_JOB,$MERGE_JOB"
echo "Final output: $OUTPUT_DIR/dipole_output.txt"
echo "Warnings log: $OUTPUT_DIR/dipole_warnings.log"
