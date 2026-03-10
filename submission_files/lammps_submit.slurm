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
############################

RUN_DIR=$SLURM_SUBMIT_DIR
cp $SLURM_SUBMIT_DIR/../input_files/* $RUN_DIR

export lmp="/home1/lkyamamo/executables/lammps/lmp_mpi_2019"

echo "Run directory: $RUN_DIR"

############################
# Load modules
############################

module purge
module load usc
module load fftw

############################
# Run LAMMPS
############################

cd $RUN_DIR

echo "Running in $(pwd)"

srun --mpi=pmix_v5 -n $SLURM_NTASKS $lmp -in in.input
date
echo “simulation finished **************************************”
echo

exit 0
