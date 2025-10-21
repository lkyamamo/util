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

srun --mpi=pmix_v2 -n 64 /home1/lkyamamo/lammps/7Aug2019/lmp_mpi -in ./test.lammps

