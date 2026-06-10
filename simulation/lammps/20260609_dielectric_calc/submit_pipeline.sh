#!/bin/bash

# =============================================================================
# Dielectric constant pipeline — full pipeline submission
#
# Step 1 : SLURM job array  — one task per dump file, computes per-frame dipole
# Step 1b: Aggregate job    — combines per-frame outputs into dipole_output.txt
# Step 2 : dipoleStd job    — computes dielectric constant from dipole_output.txt
#
# All three steps are chained automatically via SLURM dependencies.
#
# Usage:
#   bash submit_pipeline.sh
#
# Edit the Configuration section below before running.
# =============================================================================

# --- Configuration ---
DUMP_DIR="../dumps"              # directory containing LAMMPS dump files
DUMP_GLOB="dielectric.*.custom" # glob matching your dump file names
OUTPUT_DIR="$DUMP_DIR/dipole_output"

CUTOFF=1.2     # O-H bond cutoff in Angstroms
TYPE_O=1       # LAMMPS atom type for oxygen
TYPE_H=2       # LAMMPS atom type for hydrogen

# Step 2 (2.dipoleStd.py) parameters
TEMPERATURE=298.0   # simulation temperature in Kelvin
LA=37.2514          # box dimension a in Angstroms
LB=37.2514          # box dimension b in Angstroms
LC=37.2514          # box dimension c in Angstroms

VENV_PATH="/home1/lkyamamo/venv/struc_analysis"

# Maximum simultaneous array tasks (throttle to avoid overwhelming the scheduler)
MAX_SIMULTANEOUS=64

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
# ---------------------

N=$(ls "$DUMP_DIR"/$DUMP_GLOB 2>/dev/null | wc -l)

if [ "$N" -eq 0 ]; then
    echo "ERROR: No files matching '$DUMP_GLOB' found in '$DUMP_DIR'" >&2
    exit 1
fi

mkdir -p "$SCRIPT_DIR/logs"
mkdir -p "$OUTPUT_DIR"

# --- Build list of unprocessed dump files ---
DUMP_PREFIX="${DUMP_GLOB%%\**}"
DUMP_SUFFIX="${DUMP_GLOB##*\*}"
PENDING_FILE="$SCRIPT_DIR/pending_dumps.txt"
> "$PENDING_FILE"

skipped=0
for DUMP in $(ls "$DUMP_DIR"/$DUMP_GLOB | sort); do
    BASENAME=$(basename "$DUMP")
    TIMESTEP="${BASENAME#$DUMP_PREFIX}"
    TIMESTEP="${TIMESTEP%$DUMP_SUFFIX}"
    if [ -f "$OUTPUT_DIR/dipole_${TIMESTEP}.txt" ]; then
        skipped=$((skipped + 1))
    else
        echo "$DUMP" >> "$PENDING_FILE"
    fi
done

N_PENDING=$(wc -l < "$PENDING_FILE")
echo "Found $N dump files total — $skipped already processed, $N_PENDING pending"

if [ "$N_PENDING" -eq 0 ]; then
    echo "All frames already processed. Nothing to submit."
    echo "Run 1.aggregate.py directly if you need to regenerate dipole_output.txt."
    exit 0
fi

ARRAY_MAX=$((N_PENDING - 1))
echo "Submitting array 0-${ARRAY_MAX} (max ${MAX_SIMULTANEOUS} simultaneous)"

ARRAY_JOB=$(sbatch --parsable \
    --array=0-${ARRAY_MAX}%${MAX_SIMULTANEOUS} \
    --export=PENDING_FILE="$PENDING_FILE",DUMP_GLOB="$DUMP_GLOB",OUTPUT_DIR="$OUTPUT_DIR",SCRIPT_DIR="$SCRIPT_DIR",CUTOFF="$CUTOFF",TYPE_O="$TYPE_O",TYPE_H="$TYPE_H",VENV_PATH="$VENV_PATH" \
    "$SCRIPT_DIR/1.array.slurm")

if [ -z "$ARRAY_JOB" ]; then
    echo "ERROR: Array job submission failed" >&2
    exit 1
fi

echo "Array job ID: $ARRAY_JOB"

MERGE_JOB=$(sbatch --parsable \
    --dependency=afterany:$ARRAY_JOB \
    --export=OUTPUT_DIR="$OUTPUT_DIR",SCRIPT_DIR="$SCRIPT_DIR",VENV_PATH="$VENV_PATH" \
    "$SCRIPT_DIR/1.aggregate.slurm")

echo "Aggregate job ID: $MERGE_JOB (runs after array completes)"

DIPOLE_STD_JOB=$(sbatch --parsable \
    --dependency=afterok:$MERGE_JOB \
    --export=OUTPUT_DIR="$OUTPUT_DIR",SCRIPT_DIR="$SCRIPT_DIR",VENV_PATH="$VENV_PATH",TEMPERATURE="$TEMPERATURE",LA="$LA",LB="$LB",LC="$LC" \
    "$SCRIPT_DIR/2.dipoleStd.slurm")

echo "dipoleStd job ID: $DIPOLE_STD_JOB (runs after aggregate completes)"
echo ""
echo "Monitor with: squeue -j $ARRAY_JOB,$MERGE_JOB,$DIPOLE_STD_JOB"
echo "Dipole output : $OUTPUT_DIR/dipole_output.txt"
echo "Warnings log  : $OUTPUT_DIR/dipole_warnings.log"
echo "Dielectric CSV: $OUTPUT_DIR/dipole_output_timestep_data_binned.csv"
