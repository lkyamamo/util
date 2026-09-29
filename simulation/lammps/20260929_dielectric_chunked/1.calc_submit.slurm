#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
# One rank per task, and the rank count IS the task count (see NRANKS below).
# Unlike the sweep pipeline — whose driver passes --ntasks on the sbatch command
# line — nothing here overrides the header, so this is the only place the rank
# count is set. It must not exceed the number of dumps/dielectric.*.custom files,
# and it must not change between a run and its resume: 2.calc_mpi.py assigns
# files to ranks by rank count, records the count it started with, and refuses to
# resume under a different one.
#SBATCH --ntasks=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=08:00:00
#SBATCH --job-name=dielectric-calc
# logs/ has to exist before the job starts — SLURM opens this file itself, well
# before the mkdir below runs. "mkdir -p logs" once, in this directory.
#SBATCH --output=logs/calc_%j.out

# Stage 2 of 2 — dipole series and dielectric constant from the chunked dumps.
#
# The manual equivalent of the sweep pipeline's dielectric calc stage
# (dielectric_submit.slurm in 20260617_dielectric_multi_traj, driven by
# submit_temperature_sweep.sh): every knob that pipeline exports is a hardcoded
# value in the Configuration block below instead.
#
#   2.calc_mpi.py    MPI over the chunk files. Each rank takes a contiguous
#                    block of dumps/dielectric.N.custom, reads each one
#                    sequentially, and appends one dipole line per frame to
#                    dipole_output/ranks/dipole_rank_<r>.txt. Checkpointed by
#                    line count, so a killed job resumes where it stopped.
#   3.dipole_std.py  Concatenated dipole series -> <M>, <M^2>, the variance and
#                    the dielectric constant, plus a per-timestep CSV.
#
# Usage (from this directory on the HPC, after 0.production_submit.slurm has
# finished — or at least after some chunks have):
#   sbatch 1.calc_submit.slurm
#
# Resuming: submit it again with the SAME --ntasks. To start over instead, delete
# dipole_output/ entirely (not just dipole_output.txt) — and note that with
# DELETE_PROCESSED_DUMPS=1 the dumps needed to redo it are gone.
#
# Edit the Configuration section below before running.

set -euo pipefail

# $SLURM_SUBMIT_DIR under sbatch; the script's own directory when run by hand
# with bash inside an salloc.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

# --- Configuration ---
# The chunk files 0.production_submit.slurm wrote, in this directory. Absolute,
# so the manifest 2.calc_mpi.py pins the file list in holds absolute paths.
DUMP_DIR="$SCRIPT_DIR/dumps"
DUMP_GLOB="dielectric.*.custom"
OUTPUT_DIR="$SCRIPT_DIR/dipole_output"

# Must match 0.production_submit.slurm's DUMP_EVERY and CHARGE_H. There is no
# driver here to keep the two stages consistent (see README.md).
DUMP_EVERY=10
CHARGE_H=0.406988

# O-H bond detection: cutoff in Angstrom, and the LAMMPS atom types.
CUTOFF=1.2
TYPE_O=1
TYPE_H=2

# KELVIN — must match 0.production_submit.slurm's TARGET_TEMP. Together with the
# box edges it sets 3.dipole_std.py's prefactor, so a wrong value here silently
# scales the dielectric constant.
TEMPERATURE=303.15

# Box edges in Angstrom. 37.2514 is the cube ICE_CUBIC replicated 6x6x6 deforms
# to at 1.0 g/cc; read them off the production run's box bounds if unsure.
LA=37.2514
LB=37.2514
LC=37.2514

# windowed | cumulative | hybrid | binned — how 3.dipole_std.py averages along
# the series for its per-timestep CSV. The final printed constant is over all
# frames either way.
AVERAGING_METHOD=cumulative

# Delete each chunk file once its frames have been reduced to dipole lines and
# those lines are fsynced. IRREVERSIBLE — a deleted dump means re-running
# 0.production_submit.slurm. 0 here (the pipeline defaults it to 1 because a
# sweep's production runs write far more trajectory than disk can hold); set 1
# when this run's dumps are too big to keep.
DELETE_PROCESSED_DUMPS=0

# Needs numpy, scipy and mpi4py. Empty uses whatever python is on PATH.
VENV_PATH="/home1/lkyamamo/venv/struc_analysis"

# --- Setup ---
# The manifest, once 2.calc_mpi.py has written one, is the list this run was
# launched against. Counting the directory instead would shrink as
# DELETE_PROCESSED_DUMPS removes finished files, so a resume could see fewer
# files than ranks and fail a check it created itself.
MANIFEST="$OUTPUT_DIR/ranks/dumps_manifest.txt"
if [ -s "$MANIFEST" ]; then
    NUM_FILES=$(grep -c . "$MANIFEST")
    echo "Trajectory count from $MANIFEST (the list this run started against)."
else
    NUM_FILES=$(ls "$DUMP_DIR"/$DUMP_GLOB 2>/dev/null | wc -l)
    if [ "$NUM_FILES" -eq 0 ]; then
        echo "ERROR: No files matching $DUMP_GLOB found in $DUMP_DIR" >&2
        echo "       Has 0.production_submit.slurm run?" >&2
        exit 1
    fi
fi

# One rank per task allocated. Fewer ranks than files is fine — each rank then
# works through its block of files in sequence — but more is not.
NUM_RANKS="${SLURM_NTASKS:-1}"
if [ "$NUM_RANKS" -gt "$NUM_FILES" ]; then
    echo "ERROR: $NUM_RANKS ranks requested, but there are only $NUM_FILES trajectory files." >&2
    echo "       Lower #SBATCH --ntasks in this script to at most $NUM_FILES." >&2
    exit 1
fi

mkdir -p logs "$OUTPUT_DIR/ranks"

module purge
module load usc

if [ -n "$VENV_PATH" ]; then
    source "$VENV_PATH/bin/activate"
fi

echo "Job ${SLURM_JOB_ID:-<interactive>} starting at $(date)"
echo "Trajectory files: ${NUM_FILES}  Ranks: ${NUM_RANKS}"

calc_args=()
if [ "$DELETE_PROCESSED_DUMPS" = "1" ]; then
    calc_args+=(--delete-processed-dumps)
    echo "Processed trajectory files WILL BE DELETED as each one is reduced."
fi

srun --ntasks="$NUM_RANKS" python "$SCRIPT_DIR/2.calc_mpi.py" \
    "${calc_args[@]+"${calc_args[@]}"}" \
    --dump-dir        "$DUMP_DIR" \
    --dump-glob       "$DUMP_GLOB" \
    --dump-every      "$DUMP_EVERY" \
    --cutoff          "$CUTOFF" \
    --type-o          "$TYPE_O" \
    --type-h          "$TYPE_H" \
    --charge-h        "$CHARGE_H" \
    --output-dir      "$OUTPUT_DIR"

DIPOLE_FILE="$OUTPUT_DIR/dipole_output.txt"

if [ -f "$DIPOLE_FILE" ]; then
    echo "dipole_output.txt already exists — skipping concatenation and validation."
else
    echo "Verifying all ranks completed..."
    missing=0
    for ((r=0; r<NUM_RANKS; r++)); do
        if [ ! -f "$OUTPUT_DIR/ranks/dipole_rank_${r}.done" ]; then
            echo "  ERROR: rank ${r} did not complete (missing .done flag)" >&2
            missing=$((missing + 1))
        fi
    done
    if [ "$missing" -gt 0 ]; then
        echo "ERROR: ${missing} rank(s) failed — dipole_output.txt will not be written." >&2
        exit 1
    fi

    # Built under a temporary name and moved into place only once every rank file
    # has been appended. Writing straight to $DIPOLE_FILE meant a job killed
    # mid-concatenation (walltime, node failure) left a PARTIAL dipole_output.txt
    # behind — and the existence check above would then skip both the validation
    # and the rebuild on the next run, handing 3.dipole_std.py a truncated series
    # and a silently wrong dielectric constant. The rename is atomic within a
    # filesystem, so the file either does not exist or is complete.
    echo "Concatenating ${NUM_RANKS} rank files..."
    tmp_dipole="$DIPOLE_FILE.partial.$$"
    trap 'rm -f "$tmp_dipole"' EXIT
    > "$tmp_dipole"
    for ((r=0; r<NUM_RANKS; r++)); do
        cat "$OUTPUT_DIR/ranks/dipole_rank_${r}.txt" >> "$tmp_dipole"
    done
    mv "$tmp_dipole" "$DIPOLE_FILE"
    trap - EXIT
    echo "  Wrote: $DIPOLE_FILE ($(wc -l < "$DIPOLE_FILE") lines)"
fi

echo "Running dipole_std..."
# Teed to a fixed filename as well as the job log, so the numbers have one
# reliable path whether this ran under sbatch or by hand.
python "$SCRIPT_DIR/3.dipole_std.py" \
    "$DIPOLE_FILE" \
    "$TEMPERATURE" \
    "$LA" "$LB" "$LC" \
    "$AVERAGING_METHOD" \
    | tee "$OUTPUT_DIR/summary.txt"

echo "Job ${SLURM_JOB_ID:-<interactive>} finished at $(date)"
