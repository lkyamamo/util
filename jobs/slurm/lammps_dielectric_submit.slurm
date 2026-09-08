#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --partition=priya
#SBATCH --time=01:00:00
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --output=STREAM_OUTPUT
#SBATCH --job-name=L_dielectric
#SBATCH --constraint=epyc-7513
#SBATCH --mail-type=all
#SBATCH --mail-user=lkyamamo@usc.edu
set -euo pipefail
ulimit -s unlimited
echo “starting simulation **************************************”
date

############################
# Paths
#
# Canonical inputs live in ../input_files/ relative to submit (never copied).
# Symlinks with the same basenames are created into the run (submit) directory so LAMMPS
# can resolve the same filenames as before while cwd is the run directory — dumps, restarts,
# log, etc., go alongside the job submission directory.
############################

RUN_DIR=$SLURM_SUBMIT_DIR
INPUT_DIR="$(cd "$RUN_DIR/../input_files" && pwd)"

# Overridable via LMP_BIN (sbatch --export), matching lammps_submit.slurm and
# lammps_cascade_submit.slurm. Unset leaves the previous hardcoded default in
# effect, so existing callers are unaffected.
export lmp="${LMP_BIN:-/home1/lkyamamo/executables/lammps/lmp_mpi_shock_2019}"

echo "Run directory (cwd for LAMMPS): $RUN_DIR"
echo "Input files symlinked from: $INPUT_DIR"

cd "$RUN_DIR"
shopt -s nullglob
for entry in "$INPUT_DIR"/*; do
  ln -sfn "$entry" "$(basename "$entry")"
done
shopt -u nullglob

mkdir -p dumps

############################
# Load modules
############################

module purge
module load usc
module load fftw

############################
# Run LAMMPS
#
# TARGET_TEMP/N_TIMES/NVT_LENGTH/DUMP_EVERY/CHARGE_H are set via sbatch
# --export by submit_dielectric_pipeline.sh and consumed by the input script
# as LAMMPS command-line variables (${TARGET_TEMP} etc.). CHARGE_H is the
# same value passed to 1.calc_mpi.py's --charge-h in stage 2 — input scripts
# that do a by-hand dipole calc (e.g. OH-dielectric-therm.input) read it from
# here rather than hardcoding it, so the two independently-computed dipoles
# stay comparable. Plain input scripts that don't reference ${CHARGE_H}
# simply ignore it.
############################

echo "Running in $(pwd)"
echo "TARGET_TEMP=$TARGET_TEMP  N_TIMES=$N_TIMES  NVT_LENGTH=$NVT_LENGTH  DUMP_EVERY=$DUMP_EVERY  CHARGE_H=$CHARGE_H"

############################
# Resume point
#
# chunks_done.txt is appended by in.input after each chunk's restart file is
# safely written, so its last line is the newest chunk whose restart.<n> is
# known-complete. Resuming from there redoes at most the chunk that was in
# flight when the job died, whose dump is rewritten from scratch.
#
# Deliberately worked out here rather than passed in by
# submit_temperature_sweep.sh: this directory is the only thing that knows how
# far the run got, so re-running "sbatch lammps_submit.slurm" by hand resumes
# just as correctly as going through the pipeline.
#
# To force a clean start, delete chunks_done.txt (and ideally the restart.*
# files) from this directory.
############################

START_CHUNK=1
RESTART_FILE=none

if [ -s chunks_done.txt ]; then
    # Last NON-BLANK line, not simply the last line: a stray trailing newline
    # would otherwise read as "nothing completed" and silently throw away every
    # finished chunk by starting over.
    last_done="$(grep -v '^[[:space:]]*$' chunks_done.txt | tail -n 1 | tr -d '[:space:]')"
    if [ -n "$last_done" ] && [ -f "restart.${last_done}" ]; then
        if [ "$last_done" -ge "$N_TIMES" ] 2>/dev/null; then
            echo "All $N_TIMES chunks already complete — nothing to run."
            date
            echo "simulation finished **************************************"
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

srun --mpi=pmix_v5 -n $SLURM_NTASKS $lmp \
  -log "$LOG_NAME" \
  -var TARGET_TEMP "$TARGET_TEMP" \
  -var N_TIMES "$N_TIMES" \
  -var NVT_LENGTH "$NVT_LENGTH" \
  -var DUMP_EVERY "$DUMP_EVERY" \
  -var CHARGE_H "$CHARGE_H" \
  -var START_CHUNK "$START_CHUNK" \
  -var RESTART_FILE "$RESTART_FILE" \
  -in in.input
date
echo “simulation finished **************************************”
echo

exit 0
