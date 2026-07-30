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

export lmp="/home1/lkyamamo/executables/lammps/lmp_mpi_shock_2019"

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
# TARGET_TEMP/N_TIMES/NVT_LENGTH/DUMP_EVERY are set via sbatch --export by
# submit_dielectric_pipeline.sh and consumed by the input script as LAMMPS
# command-line variables (${TARGET_TEMP} etc.).
############################

echo "Running in $(pwd)"
echo "TARGET_TEMP=$TARGET_TEMP  N_TIMES=$N_TIMES  NVT_LENGTH=$NVT_LENGTH  DUMP_EVERY=$DUMP_EVERY"

srun --mpi=pmix_v5 -n $SLURM_NTASKS $lmp \
  -log log.lammps \
  -var TARGET_TEMP "$TARGET_TEMP" \
  -var N_TIMES "$N_TIMES" \
  -var NVT_LENGTH "$NVT_LENGTH" \
  -var DUMP_EVERY "$DUMP_EVERY" \
  -in in.input
date
echo “simulation finished **************************************”
echo

exit 0
