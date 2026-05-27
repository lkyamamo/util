#!/bin/bash

###  generate formatted file  ###
module purge
module load python
python 1.robust_lammps_to_input_xyz.py
bash 2.xyz_split.sh


##################################

module purge
module load legacy/CentOS7
module load usc gcc/9.2.0 boost
make corr 
make msd

###  correlation calulation ###

# ./corr.exe arg1 arg2 arg3 arg4

# arg1. data directory that contains separte XYZ files
# arg2. single correlation length [fs] 
# arg3. time interval between each correlations [fs]
# arg4. time step between each MD frame [fs]

# Example)
# compute correlation function for 1000 (fs), 
# start new correlation every 100 (fs) (multiple initial steps), 
# XYZ files are 1 (fs) apart
./qmd_corr.exe pos_vel 5000 10 2.5


###  MSD calulation ###

# ./msd.exe arg1 arg2 arg3 arg4 arg5

# arg1. concatenated XYZ file
# arg2. single correlation length [fs] 
# arg3. time interval between each correlations [fs]
# arg4. time step between each MD frame [fs]
# arg5. max frequency after Fourier transform [eV]

# Example)
# compute MSD etc for 1000 (fs), start new correlation every 10 (fs) (multiple initial steps), 
# MD frames are 1 (fs) apart, FT for DOS calculation upto 0.5 (ev)
./qmd_msd.exe calculation.xyz 5000 10 2.5 0.5
