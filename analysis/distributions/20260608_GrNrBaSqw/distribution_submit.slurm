#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=priya
#SBATCH --time=04:00:00
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --output=STREAM_OUTPUT
#SBATCH --job-name=analysis
#SBATCH --constraint=epyc-7513
#SBATCH --mail-type=all
#SBATCH --mail-user=lkyamamo@usc.edu
ulimit -s unlimited
echo "starting analysis **************************************"
date

############################
# Trajectory — set these before submitting
############################

# TRAJ: structural trajectory (rdf_freud.py/bad_freud.py).
# DYNAMICS_TRAJ: separate, higher-frequency trajectory (dsf.py/vdos.py) —
# resolving vibrational frequencies needs much finer time sampling than
# structural analysis does. See OH-therm.input/b-SiO-therm.input's two dump
# commands (dump.lammpstrj / dynamics.lammpstrj).
TRAJ="${TRAJ:-../OH.lammpstrj}"
DYNAMICS_TRAJ="${DYNAMICS_TRAJ:-../dynamics.lammpstrj}"
export TRAJ DYNAMICS_TRAJ
echo "Trajectory: $TRAJ"
echo "Dynamics trajectory: $DYNAMICS_TRAJ"

############################
# Run flags — set 1 to run, 0 to skip
############################

RUN_DSF="${RUN_DSF:-1}"
RUN_RDF="${RUN_RDF:-1}"
RUN_BAD="${RUN_BAD:-1}"
RUN_VDOS="${RUN_VDOS:-1}"

############################
# Paths
############################

RUN_DIR=$SLURM_SUBMIT_DIR
cd "$RUN_DIR"
echo "Running in $(pwd)"

############################
# Load modules + environment
############################

module purge
module load usc
source /home1/lkyamamo/venv/struc_analysis/bin/activate        # env with dynasor, icc_rt, freud, numpy, matplotlib

############################
# Thread count
############################

# dsf.py (dynasor/numba) uses OMP_NUM_THREADS.
# rdf_freud.py and bad_freud.py (freud/TBB) use the same variable.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "Threads: $OMP_NUM_THREADS"

# vdos.py's fft_periodogram method reads VDOS_THREADS directly (scipy.fft
# workers); its default vacf_cosine_transform method instead benefits from
# OMP_NUM_THREADS above via numpy's underlying BLAS matmul.
export VDOS_THREADS=$SLURM_CPUS_PER_TASK

############################
# Run scripts
############################

if [ "$RUN_DSF" -eq 1 ]; then
    echo "--- dsf.py ---"
    python dsf.py
fi

if [ "$RUN_RDF" -eq 1 ]; then
    echo "--- rdf_freud.py ---"
    python rdf_freud.py
fi

if [ "$RUN_BAD" -eq 1 ]; then
    echo "--- bad_freud.py ---"
    python bad_freud.py
fi

if [ "$RUN_VDOS" -eq 1 ]; then
    echo "--- vdos.py ---"
    python vdos.py
fi

date
echo "analysis finished **************************************"
echo

exit 0
