#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --constraint=epyc-7513
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --time=08:00:00
#SBATCH --job-name=dielectric-production
# logs/ has to exist before the job starts — SLURM opens this file itself, well
# before the mkdir below runs. "mkdir -p logs" once, in this directory.
#SBATCH --output=logs/production_%j.out

# Stage 1 of 2 — chunked NVT dielectric production for a SINGLE temperature.
#
# The manual equivalent of the sweep pipeline's dielectric production stage
# (jobs/slurm/lammps_dielectric_submit.slurm + that pipeline's
# dielectric-production.input): same restartable chunk loop, same by-hand dipole
# cross-check, but one hardcoded temperature and configuration block instead of
# a driver exporting environment variables, and no cascade in front of it — the
# thermalized structure is assumed to exist already (THERMALIZED_STRUCTURE).
#
# Writes N_TIMES trajectory files, dumps/dielectric.{1..N_TIMES}.custom, which
# stage 2 (1.calc_submit.slurm) reduces to a dipole series and a dielectric
# constant.
#
# Usage (from this directory on the HPC):
#   sbatch 0.production_submit.slurm
#
# Resuming: just submit it again. chunks_done.txt records which chunks finished,
# and the run picks up at the first one that did not (see "Resume point" below).
# To start over instead, delete chunks_done.txt, restart.*, dumps/ and
# dipole_by_hand.txt from this directory first.
#
# Edit the Configuration section below before running.

set -euo pipefail
ulimit -s unlimited

# $SLURM_SUBMIT_DIR under sbatch; the script's own directory when run by hand
# with bash inside an salloc.
SCRIPT_DIR="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
cd "$SCRIPT_DIR"

# --- Configuration ---
# The already-thermalized, already-at-density structure this production run
# starts from (symlinked in as start.data). A cascade stop
# (cascade/thermalized_T<C>.data) or any earlier thermalization's final.data.
THERMALIZED_STRUCTURE="../thermalized.data"

# Symlinked in as OH.usc, which dielectric-production.input's pair_coeff names.
POTENTIAL_FILE="/scratch1/lkyamamo/finalized-modifications/potentials/20260729_OH_B_theta_mod_v3.usc"

# Must be built with pair_style usc.
LMP_BIN="/home1/lkyamamo/executables/lammps/lmp_mpi_shock_2019"

INPUT_FILE="dielectric-production.input"

# KELVIN, unlike the sweep pipeline's Celsius lists — nothing converts for you
# here. Keep it equal to TEMPERATURE in 1.calc_submit.slurm, which feeds
# 3.dipole_std.py's prefactor.
TARGET_TEMP=303.15

# The chunking. N_TIMES chunk files of NVT_LENGTH steps each; DUMP_EVERY is the
# frame cadence within a chunk. 128 x 468750 x 0.25 fs = 15 ns total.
N_TIMES=128
NVT_LENGTH=468750
DUMP_EVERY=10

# Hydrogen partial charge (e) for the by-hand LAMMPS dipole. Oxygen's is
# derived from charge neutrality (-2x this). Keep it equal to CHARGE_H in
# 1.calc_submit.slurm, or the two dipole series stop being comparable.
CHARGE_H=0.406988

# --- Setup ---
if [ ! -f "$THERMALIZED_STRUCTURE" ]; then
    echo "ERROR: THERMALIZED_STRUCTURE does not exist: $THERMALIZED_STRUCTURE" >&2
    echo "       This stage does not thermalize anything — point it at a finished" >&2
    echo "       structure equilibrated at ${TARGET_TEMP} K and the target density." >&2
    exit 1
fi
if [ ! -f "$POTENTIAL_FILE" ]; then
    echo "ERROR: POTENTIAL_FILE does not exist: $POTENTIAL_FILE" >&2
    exit 1
fi

mkdir -p logs dumps

# Symlinks rather than copies, under the names the input script uses, so the run
# directory records which structure and potential it used without duplicating
# either. -sfn so a resume refreshes them instead of failing on the existing link.
ln -sfn "$(realpath "$THERMALIZED_STRUCTURE")" start.data
ln -sfn "$(realpath "$POTENTIAL_FILE")" OH.usc

module purge
module load usc
module load fftw

NTASKS="${SLURM_NTASKS:-64}"

echo "Job ${SLURM_JOB_ID:-<interactive>} starting at $(date)"
echo "Run directory: $SCRIPT_DIR"
echo "start.data -> $(realpath start.data)"
echo "OH.usc     -> $(realpath OH.usc)"
echo "TARGET_TEMP=$TARGET_TEMP K  N_TIMES=$N_TIMES  NVT_LENGTH=$NVT_LENGTH  DUMP_EVERY=$DUMP_EVERY  CHARGE_H=$CHARGE_H"

# --- Resume point ---
#
# chunks_done.txt is appended by dielectric-production.input after each chunk's
# restart file is safely written, so its last line is the newest chunk whose
# restart.<n> is known-complete. Resuming from there redoes at most the chunk
# that was in flight when the job died, whose dump is rewritten from scratch.
#
# Worked out here, from this directory, rather than passed in: re-submitting
# this script is the whole resume procedure.
START_CHUNK=1
RESTART_FILE=none

if [ -s chunks_done.txt ]; then
    # Last NON-BLANK line, not simply the last line: a stray trailing newline
    # would otherwise read as "nothing completed" and throw away every finished
    # chunk by starting over.
    last_done="$(grep -v '^[[:space:]]*$' chunks_done.txt | tail -n 1 | tr -d '[:space:]')"
    if [ -n "$last_done" ] && [ -f "restart.${last_done}" ]; then
        if [ "$last_done" -ge "$N_TIMES" ] 2>/dev/null; then
            echo "All $N_TIMES chunks already complete — nothing to run."
            echo "Next: sbatch 1.calc_submit.slurm"
            date
            exit 0
        fi
        START_CHUNK=$((last_done + 1))
        RESTART_FILE="restart.${last_done}"
        echo "RESUMING: $last_done of $N_TIMES chunks done; restarting at chunk $START_CHUNK from $RESTART_FILE"
    else
        echo "WARNING: chunks_done.txt names chunk '$last_done' but restart.${last_done} is missing." >&2
        echo "         Starting from chunk 1." >&2
    fi
else
    echo "No chunks_done.txt — starting from chunk 1."
fi

# A resumed run gets its own log rather than overwriting the one that recorded
# how the previous attempt died.
LOG_NAME="log.lammps"
if [ "$START_CHUNK" -ne 1 ]; then
    LOG_NAME="log.lammps.resume${START_CHUNK}"
fi
echo "LAMMPS log: $LOG_NAME"

srun --mpi=pmix_v5 -n "$NTASKS" "$LMP_BIN" \
  -log "$LOG_NAME" \
  -var TARGET_TEMP  "$TARGET_TEMP" \
  -var N_TIMES      "$N_TIMES" \
  -var NVT_LENGTH   "$NVT_LENGTH" \
  -var DUMP_EVERY   "$DUMP_EVERY" \
  -var CHARGE_H     "$CHARGE_H" \
  -var START_CHUNK  "$START_CHUNK" \
  -var RESTART_FILE "$RESTART_FILE" \
  -in "$INPUT_FILE"

echo "Job ${SLURM_JOB_ID:-<interactive>} finished at $(date)"
echo "Next: sbatch 1.calc_submit.slurm"
