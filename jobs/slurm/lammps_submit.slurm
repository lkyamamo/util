#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --partition=priya
#SBATCH --time=01:00:00
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --output=STREAM_OUTPUT
#SBATCH --job-name=L_cov
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

# The phonon build, which carries the PHONON package that OH-therm.input's
# dynamical_matrix block needs (RUN_DYNMAT=1). Override with LMP_BIN
# (submit_pipeline.sh's --lmp-bin) to go back to lmp_mpi_shock_2019 or any other
# build in that directory.
export lmp="${LMP_BIN:-/home1/lkyamamo/executables/lammps/lmp_mpi_phonon_2019}"

echo "Run directory (cwd for LAMMPS): $RUN_DIR"
echo "Input files symlinked from: $INPUT_DIR"

cd "$RUN_DIR"
shopt -s nullglob
for entry in "$INPUT_DIR"/*; do
  ln -sfn "$entry" "$(basename "$entry")"
done
shopt -u nullglob

############################
# Load modules
############################

module purge
module load usc
module load fftw

############################
# Run LAMMPS
############################

echo "Running in $(pwd)"

# Forward selected environment variables to LAMMPS as -var, but only the ones
# actually set — an unset name leaves the input script's own
# "variable NAME index <default>" in effect. That is what keeps this template
# backward compatible: with none of these exported the array is empty and the
# srun line below is byte-identical to what it always was.
lmp_vars=()
for _v in RUN_DYNMAT DYNMAT_MIN_STYLE DYNMAT_MIN_ETOL DYNMAT_MIN_FTOL \
          DYNMAT_MIN_MAXITER DYNMAT_MIN_MAXEVAL DYNMAT_DISPLACEMENT \
          DYNMAT_FILE DYNMAT_REF_FILE DYNMAT_BINARY; do
  if [ -n "${!_v:-}" ]; then
    lmp_vars+=(-var "$_v" "${!_v}")
    echo "  -var $_v ${!_v}"
  fi
done
unset _v

srun --mpi=pmix_v5 -n $SLURM_NTASKS $lmp \
  -log log.lammps \
  "${lmp_vars[@]+"${lmp_vars[@]}"}" \
  -in in.input
date
echo “simulation finished **************************************”
echo

exit 0
