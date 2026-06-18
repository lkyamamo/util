#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3800M
#SBATCH --time=08:00:00
#SBATCH --output=logs/calc_%j.out

# Dielectric constant pipeline — multi-trajectory MPI.
#
# One MPI rank per trajectory file. Each rank reads its assigned file
# from line 0 to EOF with no seeking. Rank count is set automatically
# from the number of files found in DUMP_DIR.
#
# Usage (from this directory on the HPC login node):
#   sbatch 0.submit.slurm
#
# Edit the Configuration section below before running.

set -euo pipefail

SCRIPT_DIR="$SLURM_SUBMIT_DIR"
cd "$SCRIPT_DIR"

# --- Configuration ---
DUMP_DIR="../dumps"
DUMP_GLOB="dielectric.*.custom"
OUTPUT_DIR="$SCRIPT_DIR/dipole_output"

DUMP_EVERY=10
CUTOFF=1.2
TYPE_O=1
TYPE_H=2

TEMPERATURE=303.0
AVERAGING_METHOD=cumulative

VENV_PATH="/home1/lkyamamo/venv/struc_analysis"

# --- Setup ---
NUM_FILES=$(ls "$DUMP_DIR"/$DUMP_GLOB 2>/dev/null | wc -l)
if [ "$NUM_FILES" -eq 0 ]; then
    echo "ERROR: No files matching $DUMP_GLOB found in $DUMP_DIR"
    exit 1
fi

mkdir -p "$SCRIPT_DIR/logs" "$OUTPUT_DIR/ranks"

module purge
module load usc

if [ -n "$VENV_PATH" ]; then
    source "$VENV_PATH/bin/activate"
fi

echo "Job ${SLURM_JOB_ID} starting at $(date)"
echo "Trajectory files: ${NUM_FILES}  Ranks: ${NUM_FILES}"

srun --ntasks="$NUM_FILES" python "$SCRIPT_DIR/1.calc_mpi.py" \
    --dump-dir        "$DUMP_DIR" \
    --dump-glob       "$DUMP_GLOB" \
    --dump-every      "$DUMP_EVERY" \
    --cutoff          "$CUTOFF" \
    --type-o          "$TYPE_O" \
    --type-h          "$TYPE_H" \
    --output-dir      "$OUTPUT_DIR"

echo "Concatenating ${NUM_FILES} rank files..."
DIPOLE_FILE="$OUTPUT_DIR/dipole_output.txt"
> "$DIPOLE_FILE"
for ((r=0; r<NUM_FILES; r++)); do
    cat "$OUTPUT_DIR/ranks/dipole_rank_${r}.txt" >> "$DIPOLE_FILE"
done
echo "  Wrote: $DIPOLE_FILE ($(wc -l < "$DIPOLE_FILE") lines)"

echo "Running dipole_std..."
python "$SCRIPT_DIR/2.dipole_std.py" \
    "$DIPOLE_FILE" \
    "$TEMPERATURE" \
    "$AVERAGING_METHOD"

echo "Job ${SLURM_JOB_ID} finished at $(date)"
