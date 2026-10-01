#!/bin/bash
#BATCH -J lammps_test
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=0
#SBATCH --partition=main
#SBATCH --time=00:10:00

module purge
module load usc
module load fftw
module load openmpi

RUN_DIR=$SLURM_SUBMIT_DIR

srun --mpi=pmix_v2 -n 64 /home1/lkyamamo/executables/lammps/lmp_mpi_2019 -in ./test.lammps
LAMMPS_STATUS=$?

############################
# Ledger: record job id, exit code and script for `simledger scan`.
# Never fails the job.
############################

python3 /home1/lkyamamo/util/jobs/ledger/bin/simledger hook run --dir "$RUN_DIR" --exit-code "$LAMMPS_STATUS" \
  || echo "WARNING: ledger hook failed" >&2

