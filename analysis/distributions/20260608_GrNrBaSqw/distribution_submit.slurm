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
RUN_MSD="${RUN_MSD:-1}"

############################
# Analysis parameters
############################

# Each value below is passed to the .py scripts as an environment variable.
# Edit the text after ":-" to set one. Keys marked REQUIRED have no default —
# vdos.py/msd.py abort with a list of what is missing rather than substituting
# a number, since these determine the results. The rest fall back to the
# script's own default (see the CONFIGURATION block at the top of each .py).
# The "${VAR:-...}" form means an already-exported value wins, so the pipeline
# submitters (submit_pipeline.sh / submit_pipeline_local.sh) can drive this
# script without their settings being clobbered by the edits here.
#
# vdos.py and msd.py read the same dynamics.lammpstrj but want different
# settings — VDOS needs a short CORR_LENGTH for frequency resolution, MSD a
# long one to reach the diffusive regime — so each has its own prefixed
# variables. DYNAMICS_DT is the one and only value they share: the dt of the
# trajectory itself, which describes neither analysis in particular.

DYNAMICS_DT="${DYNAMICS_DT:-}"          # REQUIRED by vdos.py and msd.py; fs between frames

# dsf.py
DT="${DT:-}"                            # fs between dumped frames
N_FRAMES="${N_FRAMES:-}"                # max frames to read; 0 = all
STRIDE="${STRIDE:-}"                    # read every Nth frame

# vdos.py
VDOS_N_FRAMES="${VDOS_N_FRAMES:-}"
VDOS_STRIDE="${VDOS_STRIDE:-}"
VDOS_CORR_LENGTH="${VDOS_CORR_LENGTH:-}"          # REQUIRED. fs; VACF max lag. Sets the
                                                  # frequency resolution: dnu ~ 1/CORR_LENGTH
VDOS_CORR_INTERVAL="${VDOS_CORR_INTERVAL:-}"      # REQUIRED. fs; spacing between VACF origins
VDOS_MAX_FREQUENCY_EV="${VDOS_MAX_FREQUENCY_EV:-}"  # REQUIRED. eV; DOS grid upper limit
VDOS_NUM_GRIDS="${VDOS_NUM_GRIDS:-}"
VDOS_METHOD="${VDOS_METHOD:-}"                    # vacf_cosine_transform | fft_periodogram
VDOS_WINDOW="${VDOS_WINDOW:-}"
VDOS_NORMALIZATION="${VDOS_NORMALIZATION:-}"      # phonon | unit_area

# msd.py
MSD_N_FRAMES="${MSD_N_FRAMES:-}"
MSD_STRIDE="${MSD_STRIDE:-}"
MSD_CORR_LENGTH="${MSD_CORR_LENGTH:-}"            # REQUIRED. fs; max time lag. Also sets the
                                                  # diffusion fit window, so D depends on it
MSD_CORR_INTERVAL="${MSD_CORR_INTERVAL:-}"        # REQUIRED. fs; spacing between reference frames
MSD_FIT_FRACTION="${MSD_FIT_FRACTION:-}"          # REQUIRED. tail fraction used for the D fit

# Export only the ones actually set, so an empty value leaves the .py default
# in effect rather than reaching Python as an empty string.
for _var in DYNAMICS_DT DT N_FRAMES STRIDE \
            VDOS_N_FRAMES VDOS_STRIDE VDOS_CORR_LENGTH VDOS_CORR_INTERVAL \
            VDOS_MAX_FREQUENCY_EV VDOS_NUM_GRIDS VDOS_METHOD VDOS_WINDOW VDOS_NORMALIZATION \
            MSD_N_FRAMES MSD_STRIDE MSD_CORR_LENGTH MSD_CORR_INTERVAL MSD_FIT_FRACTION; do
    if [ -n "${!_var}" ]; then
        export "$_var"
        echo "  $_var=${!_var}"
    fi
done
unset _var

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

if [ "$RUN_MSD" -eq 1 ]; then
    echo "--- msd.py ---"
    python msd.py
fi

date
echo "analysis finished **************************************"
echo

exit 0
