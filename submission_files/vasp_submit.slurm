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
module load gcc/13.3.0 
module load openmpi/5.0.5 
module load hdf5/1.14.3 
module load intel-oneapi-mkl/2024.0.0

############################
# Run VASP
############################

echo "Running in $(pwd)"

srun --mpi=pmix_v5 -n $SLURM_NTASKS $vsp
date
echo “simulation finished **************************************”
echo
exit 0
