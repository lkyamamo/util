#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --ntasks=64
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3800M
#SBATCH --time=08:00:00
#SBATCH --output=logs/calc_%j.out

# Dielectric constant pipeline (MPI).
#
# Each rank processes an equal consecutive share of frames, appending
# results directly to ranks/dipole_rank_N.txt. Rank 0 concatenates all
# rank files and runs dipole_std after a barrier.
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
START_TIMESTEP=0
END_TIMESTEP=60000000
DUMP_EVERY=10

CUTOFF=1.2
TYPE_O=1
TYPE_H=2

TEMPERATURE=303.0
LA=37.2514
LB=37.2514
LC=37.2514
AVERAGING_METHOD=cumulative

VENV_PATH="/home1/lkyamamo/venv/struc_analysis"

# --- Setup ---
NUM_FRAMES=$(( (END_TIMESTEP - START_TIMESTEP) / DUMP_EVERY + 1 ))

mkdir -p "$SCRIPT_DIR/logs" "$OUTPUT_DIR/ranks"

module purge
module load usc

if [ -n "$VENV_PATH" ]; then
    source "$VENV_PATH/bin/activate"
fi

echo "Job ${SLURM_JOB_ID} starting at $(date)"
echo "Frames: ${NUM_FRAMES}  Ranks: ${SLURM_NTASKS}"

srun python "$SCRIPT_DIR/1.calc_mpi.py" \
    --num-frames       "$NUM_FRAMES" \
    --start-timestep   "$START_TIMESTEP" \
    --dump-every       "$DUMP_EVERY" \
    --dump-dir         "$DUMP_DIR" \
    --dump-glob        "$DUMP_GLOB" \
    --cutoff           "$CUTOFF" \
    --type-o           "$TYPE_O" \
    --type-h           "$TYPE_H" \
    --output-dir       "$OUTPUT_DIR"

echo "Concatenating ${SLURM_NTASKS} rank files..."
DIPOLE_FILE="$OUTPUT_DIR/dipole_output.txt"
> "$DIPOLE_FILE"
for ((r=0; r<SLURM_NTASKS; r++)); do
    cat "$OUTPUT_DIR/ranks/dipole_rank_${r}.txt" >> "$DIPOLE_FILE"
done
echo "  Wrote: $DIPOLE_FILE ($(wc -l < "$DIPOLE_FILE") lines)"

echo "Running dipole_std..."
python "$SCRIPT_DIR/2.dipole_std.py" \
    "$DIPOLE_FILE" \
    "$TEMPERATURE" \
    "$LA" "$LB" "$LC" \
    "$AVERAGING_METHOD"

echo "Job ${SLURM_JOB_ID} finished at $(date)"
