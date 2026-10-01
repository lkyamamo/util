#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --nodes=8
#SBATCH --ntasks=512
#SBATCH --partition=priya
#SBATCH --time=240:59:00
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --output=STREAM_OUTPUT
#SBATCH --job-name=gga_elastic_222
#SBATCH --constraint=epyc-7513
#SBATCH --mail-type=all
#SBATCH --mail-user=lkyamamo
ulimit -s unlimited
echo “starting simulation **************************************”
date

############################
# Paths
############################

RUN_DIR=$SLURM_SUBMIT_DIR
cp $SLURM_SUBMIT_DIR/../input_files/* $RUN_DIR

#export vsp="/home1/lkyamamo/executables/bubble_collapse_executables/vasp_gam"
export vsp="/home1/lkyamamo/executables/bubble_collapse_executables/vasp_std"

echo "Run directory: $RUN_DIR"

############################
# Load modules
############################

module purge
module load legacy/CentOS7
module load gcc/11.3.0
module load openmpi/4.1.4
module load pmix/3.2.3
module load hdf5/1.12.2
module load intel-oneapi-mkl/2022.1.0

############################
# Run VASP
############################

echo "Running in $(pwd)"

srun --mpi=pmix_v3 -n $SLURM_NTASKS $vsp
VASP_STATUS=$?
date
echo “simulation finished **************************************”
echo

############################
# Ledger: record job id, exit code and script for `simledger scan`.
# Never fails the job.
############################

/apps/spack/2406/apps/linux-rocky8-x86_64_v3/gcc-13.3.0/python-3.11.9-x74mtjf/bin/python3 /home1/lkyamamo/util/jobs/ledger/bin/simledger hook run --dir "$RUN_DIR" --exit-code "$VASP_STATUS" \
  || echo "WARNING: ledger hook failed" >&2

exit 0
