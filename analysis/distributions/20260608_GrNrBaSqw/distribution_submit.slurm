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
# vdos_dynmat.py needs dynmat.dat, which only exists if the LAMMPS run was made
# with RUN_DYNMAT=1 — so unlike the others this defaults to off.
RUN_VDOS_DYNMAT="${RUN_VDOS_DYNMAT:-0}"

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

# dsf.py — static S(q) and dynamic S(q,w). Every key is DSF_-prefixed, like
# vdos.py's VDOS_* and msd.py's MSD_*, so nothing here can be confused with
# another analysis's setting. dsf.py errors out if the old bare names (DT,
# N_FRAMES, STRIDE, Q_MAX, ...) are set, rather than ignoring them.
#
# Note DSF_DT is dsf.py's own dt and is NOT the same key as DYNAMICS_DT above:
# vdos.py and msd.py share that one because it describes the trajectory, whereas
# this is read together with DSF_STRIDE and describes only this analysis.
DSF_DT="${DSF_DT:-}"                    # fs between consecutive dumped frames
DSF_N_FRAMES="${DSF_N_FRAMES:-}"        # frame_stop INDEX into the dump, not a count;
                                        # frames used = DSF_N_FRAMES / DSF_STRIDE
DSF_STRIDE="${DSF_STRIDE:-}"            # read every Nth frame. Skipped frames are still
                                        # parsed, so widening the span is free while
                                        # computing more frames is not

# Static S(q). Q_MAX=20 A^-1 is the range needed to Fourier transform S(q) into
# G(r) without bad termination ripples (ripple period 2*pi/Q_MAX = 0.31 A) and
# matches neutron diffraction on vitreous silica. Cost is linear in N_q x frames:
# ~10.6M q-vectors on a 43 A cell at Q_MAX=20 versus 2.3M at 12.
DSF_Q_MAX="${DSF_Q_MAX:-}"              # A^-1, static
DSF_N_Q_BINS="${DSF_N_Q_BINS:-}"        # radial q-bins after spherical averaging. Keep
                                        # <= Q_MAX / (2*pi/L) or low-q bins come back empty

# Dynamic S(q,w) — only used when COMPUTE_DYNAMIC=True in dsf.py, and NOT yet
# tuned against physics; see DYNAMIC NOTES in dsf.py before trusting the output.
# Memory is the binding constraint, not time:
#   2 x (n_pairs + 1) x N_q x (DSF_WINDOW_SIZE + 1) x 8 bytes
# which is why the dynamic run has its own q settings rather than reusing the
# static ones — 19 GB at Q_MAX_DYN=4 unpruned, versus 2373 GB at Q_MAX=20.
DSF_WINDOW_SIZE="${DSF_WINDOW_SIZE:-}"  # time lags; dnu = 1/(2 x WINDOW_SIZE x DT x STRIDE).
                                        # Must cover several periods of the slowest mode
DSF_WINDOW_STEP="${DSF_WINDOW_STEP:-}"  # frames between window origins; 1 = most averaging,
                                        # and nearly free, so leave it at 1
DSF_Q_MAX_DYN="${DSF_Q_MAX_DYN:-}"      # A^-1, dynamic; separate from DSF_Q_MAX on purpose
DSF_N_Q_BINS_DYN="${DSF_N_Q_BINS_DYN:-}"        # radial q-bins, dynamic
DSF_MAX_Q_POINTS_DYN="${DSF_MAX_Q_POINTS_DYN:-}"  # prune target; 0 = no pruning. Keeps the
                                        # low-q mesh intact and thins only above a cutoff

# Neutron-weighted S(q) columns. dsf.py delegates the weighting to dynasor,
# which weights by natural abundance — its 'H' is protium (b = -3.7406 fm) while
# rdf_freud.py treats 'H' as deuterium (+6.671 fm), opposite signs. dsf.py
# therefore REFUSES to neutron-weight an H- or D-bearing system; set this to
# 'no' there to get the unweighted partials and total instead.
DSF_NEUTRON_WEIGHTING="${DSF_NEUTRON_WEIGHTING:-}"  # yes | no

# rdf_freud.py
R_MAX="${R_MAX:-}"                      # Å; max r. Must be < half the shortest box edge
RDF_BINS="${RDF_BINS:-}"                # number of r-bins
RDF_NORMALIZATION="${RDF_NORMALIZATION:-}"        # semicolon list: unity | FZ | absolute
RDF_FUNCTIONS="${RDF_FUNCTIONS:-}"                # semicolon list: g | h | D | T
RDF_RESOLUTION_SIGMA="${RDF_RESOLUTION_SIGMA:-}"  # Å; Gaussian resolution broadening, 0 disables

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
VDOS_NORMALIZATION="${VDOS_NORMALIZATION:-}"      # phonon | unit_area (sum rule)
VDOS_WEIGHTING="${VDOS_WEIGHTING:-}"              # semicolon list: unity | coherent |
                                                  # incoherent | total (species weight)

# msd.py
MSD_N_FRAMES="${MSD_N_FRAMES:-}"
MSD_STRIDE="${MSD_STRIDE:-}"
MSD_CORR_LENGTH="${MSD_CORR_LENGTH:-}"            # REQUIRED. fs; max time lag. Also sets the
                                                  # diffusion fit window, so D depends on it
MSD_CORR_INTERVAL="${MSD_CORR_INTERVAL:-}"        # REQUIRED. fs; spacing between reference frames
MSD_FIT_FRACTION="${MSD_FIT_FRACTION:-}"          # REQUIRED. tail fraction used for the D fit

# vdos_dynmat.py — the harmonic counterpart to vdos.py. It reads the dynamical
# matrix LAMMPS wrote (not a trajectory), plus the first frame of $TRAJ for the
# per-atom element labels. DYNMAT_FILE/DYNMAT_BINARY must match what the LAMMPS
# stage used, which submit_pipeline.sh guarantees by sending both stages the same
# value.
DYNMAT_FILE="${DYNMAT_FILE:-}"                              # default dynmat.dat
DYNMAT_BINARY="${DYNMAT_BINARY:-}"                          # yes | no (default no)
VDOS_DYNMAT_MAX_FREQUENCY="${VDOS_DYNMAT_MAX_FREQUENCY:-}"  # REQUIRED. DOS grid upper
                                                            # limit, in VDOS_DYNMAT_XUNIT
VDOS_DYNMAT_XUNIT="${VDOS_DYNMAT_XUNIT:-}"                  # meV | THz | cm-1 | eV
VDOS_DYNMAT_BINS="${VDOS_DYNMAT_BINS:-}"                    # frequency grid points
VDOS_DYNMAT_SMEARING="${VDOS_DYNMAT_SMEARING:-}"            # Gaussian FWHM in XUNIT; 0 = off
VDOS_DYNMAT_MATRIX_STYLE="${VDOS_DYNMAT_MATRIX_STYLE:-}"    # regular | eskm
VDOS_DYNMAT_NORMALIZATION="${VDOS_DYNMAT_NORMALIZATION:-}"  # phonon | unit_area
VDOS_DYNMAT_WEIGHTING="${VDOS_DYNMAT_WEIGHTING:-}"        # semicolon list: unity |
                                                            # coherent | incoherent | total
VDOS_DYNMAT_PARTIAL="${VDOS_DYNMAT_PARTIAL:-}"              # yes | no
VDOS_DYNMAT_ASR="${VDOS_DYNMAT_ASR:-}"                      # none | simple
VDOS_DYNMAT_OUTPUT="${VDOS_DYNMAT_OUTPUT:-}"                # output basename
# Mode-character analysis: which motion each band is (rocking / bending /
# stretching at bridging oxygens), plus participation ratio and g(nu)/nu^2.
# Needs DYNMAT_REF_TRAJ, the coordinates LAMMPS wrote right after `minimize`.
VDOS_DYNMAT_CHARACTER="${VDOS_DYNMAT_CHARACTER:-}"          # yes | no (default no)
DYNMAT_REF_TRAJ="${DYNMAT_REF_TRAJ:-}"                      # default dynmat_ref.lammpstrj
VDOS_DYNMAT_BRIDGE_ELEMENT="${VDOS_DYNMAT_BRIDGE_ELEMENT:-}"      # default O
VDOS_DYNMAT_NEIGHBOR_ELEMENT="${VDOS_DYNMAT_NEIGHBOR_ELEMENT:-}"  # default Si
VDOS_DYNMAT_BOND_CUTOFF="${VDOS_DYNMAT_BOND_CUTOFF:-}"      # Angstrom, default 2.2

# Export only the ones actually set, so an empty value leaves the .py default
# in effect rather than reaching Python as an empty string.
for _var in DYNAMICS_DT \
            DSF_DT DSF_N_FRAMES DSF_STRIDE DSF_Q_MAX DSF_N_Q_BINS \
            DSF_WINDOW_SIZE DSF_WINDOW_STEP DSF_Q_MAX_DYN DSF_N_Q_BINS_DYN \
            DSF_MAX_Q_POINTS_DYN DSF_NEUTRON_WEIGHTING \
            R_MAX RDF_BINS RDF_NORMALIZATION RDF_FUNCTIONS RDF_RESOLUTION_SIGMA \
            VDOS_N_FRAMES VDOS_STRIDE VDOS_CORR_LENGTH VDOS_CORR_INTERVAL \
            VDOS_MAX_FREQUENCY_EV VDOS_NUM_GRIDS VDOS_METHOD VDOS_WINDOW VDOS_NORMALIZATION VDOS_WEIGHTING \
            MSD_N_FRAMES MSD_STRIDE MSD_CORR_LENGTH MSD_CORR_INTERVAL MSD_FIT_FRACTION \
            DYNMAT_FILE DYNMAT_BINARY \
            VDOS_DYNMAT_MAX_FREQUENCY VDOS_DYNMAT_XUNIT VDOS_DYNMAT_BINS \
            VDOS_DYNMAT_SMEARING VDOS_DYNMAT_MATRIX_STYLE VDOS_DYNMAT_NORMALIZATION VDOS_DYNMAT_WEIGHTING \
            VDOS_DYNMAT_PARTIAL VDOS_DYNMAT_ASR VDOS_DYNMAT_OUTPUT \
            VDOS_DYNMAT_CHARACTER DYNMAT_REF_TRAJ VDOS_DYNMAT_BRIDGE_ELEMENT \
            VDOS_DYNMAT_NEIGHBOR_ELEMENT VDOS_DYNMAT_BOND_CUTOFF; do
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

# vdos_dynmat.py diagonalizes the dynamical matrix; np.linalg.eigh is a threaded
# LAPACK call, so OMP_NUM_THREADS above is what parallelizes it.
if [ "$RUN_VDOS_DYNMAT" -eq 1 ]; then
    echo "--- vdos_dynmat.py ---"
    python vdos_dynmat.py
fi

date
echo "analysis finished **************************************"
echo

exit 0
