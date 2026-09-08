#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --partition=priya
#SBATCH --time=24:00:00
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --output=STREAM_OUTPUT
#SBATCH --job-name=L_cascade
#SBATCH --constraint=epyc-7513
#SBATCH --mail-type=all
#SBATCH --mail-user=lkyamamo@usc.edu
set -euo pipefail
ulimit -s unlimited
echo "starting cascade **************************************"
date

############################
# Paths
#
# Same convention as lammps_submit.slurm / lammps_dielectric_submit.slurm:
# canonical inputs live in ../input_files/ (never copied), and symlinks with
# the same basenames are made into the run (submit) directory so LAMMPS
# resolves flat filenames while writing its outputs beside the submit dir.
############################

RUN_DIR=$SLURM_SUBMIT_DIR
INPUT_DIR="$(cd "$RUN_DIR/../input_files" && pwd)"

# Overridable via LMP_BIN (sbatch --export), the same convention
# lammps_submit.slurm uses and submit_temperature_sweep.sh exposes as
# --lmp-bin. The cascade is plain MD — deform, nvt, nve, write_data — so it
# needs no optional package, but the binary still has to be one built with
# pair_style usc.
export lmp="${LMP_BIN:-/home1/lkyamamo/executables/lammps/lmp_mpi_shock_2019}"

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
#
# No -var flags, unlike lammps_dielectric_submit.slurm: in.input here is
# GENERATED per sweep by generate_cascade_input.py, which bakes the box
# length, every temperature, and which stops get an NVE diffusion block
# directly into the file. The generated input is the record of what ran — read
# ../input_files/in.input to see the exact schedule.
#
# This one job produces every thermalized_T<C>.data and every
# dynamics_T<C>.lammpstrj the downstream per-temperature jobs consume, so all
# of them depend on it finishing.
############################

echo "Running in $(pwd)"
srun --mpi=pmix_v5 -n $SLURM_NTASKS $lmp -log log.lammps -in in.input

echo "--- structures written ---"
ls -1 thermalized_T*.data 2>/dev/null || echo "(none — check log.lammps)"
echo "--- diffusion trajectories written ---"
ls -1 dynamics_T*.lammpstrj 2>/dev/null || echo "(none requested)"

date
echo "cascade finished **************************************"
echo

exit 0
