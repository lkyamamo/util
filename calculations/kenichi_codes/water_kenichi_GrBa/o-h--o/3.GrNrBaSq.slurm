#!/bin/bash
#SBATCH -J water_grba_o-h--o
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=62G
#SBATCH --partition=priya
#SBATCH --account=priyav_216
#SBATCH --time=03:00:00

module purge
module load legacy/CentOS7
module load intel-oneapi/2021.3


srun ../stat/a_o-h--o.out ../pos_vel
mv ba.dat ba_o-h--o.dat

