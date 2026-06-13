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
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DUMP_DIR="../dumps"              # directory containing LAMMPS dump files
DUMP_GLOB="dielectric.*.custom" # glob matching your dump file names
OUTPUT_DIR="$SCRIPT_DIR/dipole_output"
START_TIMESTEP=0
END_TIMESTEP=60000000
DUMP_EVERY=10



CUTOFF=1.2     # O-H bond cutoff in Angstroms
TYPE_O=1       # LAMMPS atom type for oxygen
TYPE_H=2       # LAMMPS atom type for hydrogen

# Step 2 (2.dipoleStd.py) parameters
TEMPERATURE=298.0     # simulation temperature in Kelvin
LA=37.2514            # box dimension a in Angstroms
LB=37.2514            # box dimension b in Angstroms
LC=37.2514            # box dimension c in Angstroms
AVERAGING_METHOD=binned  # windowed | cumulative | hybrid | binned

VENV_PATH="/home1/lkyamamo/venv/struc_analysis"

# Maximum simultaneous array tasks (throttle to avoid overwhelming the scheduler)
MAX_SIMULTANEOUS=64


# ---------------------

# check to see if the dump files exist
if [ ! -f "$DUMP_DIR/${DUMP_GLOB/\*/$START_TIMESTEP}" ]; then
    echo "ERROR: No files matching '${DUMP_GLOB/\*/$START_TIMESTEP}' found in '$DUMP_DIR'" >&2
    exit 1
fi

mkdir -p "$SCRIPT_DIR/logs"

PENDING_FILE="$SCRIPT_DIR/pending_dumps.txt"
OUT_GLOB="dipole_*.txt"
> "$PENDING_FILE"

# see if the output directory exists. if exists, assume processed files exist
if [ -d "$OUTPUT_DIR" ]; then
    echo "Existing processed frames"
    skipped=0
    toprocess=0
    missing=0
    for ((i=START_TIMESTEP; i<=END_TIMESTEP; i+=DUMP_EVERY)); do
        if [ ! -f "$OUTPUT_DIR/${OUT_GLOB/\*/$i}" ]; then
            TEMP="${DUMP_GLOB/\*/$i}"
            if [ -f "$DUMP_DIR/$TEMP" ]; then
                echo "$DUMP_DIR/$TEMP" >> "$PENDING_FILE"
                toprocess=$((toprocess + 1))
            else
                echo "Missing timestep: $i"
                missing=$((missing + 1))
            fi
        else
            skipped=$((skipped + 1))
        fi
    done

    if [ $skipped == $((($END_TIMESTEP - $START_TIMESTEP)/$DUMP_EVERY)+1)]; then
        echo "All frames already processed"
        echo "Run 1.aggregate.py directly if you need to regenerate dipole_output.txt."
        exit 0
    fi

    echo "Skipped: $skipped"
    echo "Pending: $toprocess"
    echo "Missing: $missing"
# no ouptut dir therefore cannot have any processed frames
else
    mkdir -p "$OUTPUT_DIR"
    echo "No previous frames processed"
    toprocess=0
    missing=0
    for ((i=START_TIMESTEP; i<=END_TIMESTEP; i+=DUMP_EVERY)); do
        TEMP="${DUMP_GLOB/\*/$i}"
        if [ -f "$DUMP_DIR/$TEMP" ]; then
            echo "$DUMP_DIR/$TEMP" >> "$PENDING_FILE"
            toprocess=$((toprocess + 1))
        else                
            echo "Missing timestep: $i"
            missing=$((missing + 1))
        fi
    done
    echo "Skipped: 0"
    echo "Pending: $toprocess"
    echo "Missing: $missing"
fi




ARRAY_MAX=$((toprocess - 1))
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
    --export=OUTPUT_DIR="$OUTPUT_DIR",SCRIPT_DIR="$SCRIPT_DIR",VENV_PATH="$VENV_PATH",TEMPERATURE="$TEMPERATURE",LA="$LA",LB="$LB",LC="$LC",AVERAGING_METHOD="$AVERAGING_METHOD" \
    "$SCRIPT_DIR/2.dipoleStd.slurm")

echo "dipoleStd job ID: $DIPOLE_STD_JOB (runs after aggregate completes)"
echo ""
echo "Monitor with: squeue -j $ARRAY_JOB,$MERGE_JOB,$DIPOLE_STD_JOB"
echo "Dipole output : $OUTPUT_DIR/dipole_output.txt"
echo "Warnings log  : $OUTPUT_DIR/dipole_warnings.log"
echo "Dielectric CSV: $OUTPUT_DIR/dipole_output_timestep_data_${AVERAGING_METHOD}.csv"
