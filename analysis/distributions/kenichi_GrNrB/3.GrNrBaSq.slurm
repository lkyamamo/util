#!/bin/bash
#SBATCH -J Si02_glass_GrNrBaSq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=0
#SBATCH --partition=priya
#SBATCH --account=priyav_216
#SBATCH --time=04:00:00
#SBATCH --nodelist=b20-[17-18],b21-[13-14,17-18] 

cd stat
sh excutable.sh

echo "created executable"

module purge
module load legacy/CentOS7
module load intel-oneapi/2021.3

cd ../
srun stat/a.out pos_vel
cd ../..

