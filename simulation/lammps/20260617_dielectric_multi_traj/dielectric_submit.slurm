#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3800M
#SBATCH --time=08:00:00
#SBATCH --output=logs/calc_%j.out

# Dielectric constant pipeline — multi-trajectory MPI, pipeline-driven variant.
#
# Same one-rank-per-trajectory-file design as 0.submit.slurm in this
# directory, but configuration comes from environment variables (set via
# sbatch --export by submit_dielectric_pipeline.sh) instead of a hardcoded
# block at the top — mirrors analysis/distributions/20260608_GrNrBaSqw/distribution_submit.slurm.
#
# 0.submit.slurm itself is left untouched for manual/standalone runs.

set -euo pipefail

SCRIPT_DIR="$SLURM_SUBMIT_DIR"
cd "$SCRIPT_DIR"

# --- Configuration — set these via sbatch --export before submitting ---
DUMP_DIR="${DUMP_DIR:-$SCRIPT_DIR/dumps}"
DUMP_GLOB="${DUMP_GLOB:-dielectric.*.custom}"
OUTPUT_DIR="$SCRIPT_DIR/dipole_output"

DUMP_EVERY="${DUMP_EVERY:-10}"
CUTOFF="${CUTOFF:-1.2}"
TYPE_O="${TYPE_O:-1}"
TYPE_H="${TYPE_H:-2}"

TEMPERATURE="${TEMPERATURE:?TEMPERATURE must be set}"
LA="${LA:?LA must be set}"
LB="${LB:?LB must be set}"
LC="${LC:?LC must be set}"
AVERAGING_METHOD="${AVERAGING_METHOD:-cumulative}"

VENV_PATH="${VENV_PATH:-/home1/lkyamamo/venv/struc_analysis}"

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

DIPOLE_FILE="$OUTPUT_DIR/dipole_output.txt"

if [ -f "$DIPOLE_FILE" ]; then
    echo "dipole_output.txt already exists — skipping concatenation and validation."
else
    echo "Verifying all ranks completed..."
    missing=0
    for ((r=0; r<NUM_FILES; r++)); do
        if [ ! -f "$OUTPUT_DIR/ranks/dipole_rank_${r}.done" ]; then
            echo "  ERROR: rank ${r} did not complete (missing .done flag)" >&2
            missing=$((missing + 1))
        fi
    done
    if [ "$missing" -gt 0 ]; then
        echo "ERROR: ${missing} rank(s) failed — dipole_output.txt will not be written." >&2
        exit 1
    fi

    echo "Concatenating ${NUM_FILES} rank files..."
    > "$DIPOLE_FILE"
    for ((r=0; r<NUM_FILES; r++)); do
        cat "$OUTPUT_DIR/ranks/dipole_rank_${r}.txt" >> "$DIPOLE_FILE"
    done
    echo "  Wrote: $DIPOLE_FILE ($(wc -l < "$DIPOLE_FILE") lines)"
fi

echo "Running dipole_std..."
# Teed to a fixed filename (rather than relying on job stdout capture) so
# aggregate_dielectric_vs_temperature.py has one reliable path to parse
# regardless of sbatch vs --interactive execution.
python "$SCRIPT_DIR/2.dipole_std.py" \
    "$DIPOLE_FILE" \
    "$TEMPERATURE" \
    "$LA" "$LB" "$LC" \
    "$AVERAGING_METHOD" \
    | tee "$OUTPUT_DIR/summary.txt"

echo "Job ${SLURM_JOB_ID} finished at $(date)"
