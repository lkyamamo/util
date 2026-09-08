#!/bin/bash
# distribution_run.sh — local runner for the analysis pipeline (all cores,
# no Slurm); mirrors distribution_submit.slurm's run flags/thread wiring.
#
# Usage:
#   ./distribution_run.sh              — use all available cores
#   ./distribution_run.sh 8            — use 8 threads

############################
# Trajectory — set these before running
############################

# TRAJ: structural trajectory (rdf_freud.py/bad_freud.py).
# DYNAMICS_TRAJ: separate, higher-frequency trajectory (dsf.py/vdos.py/msd.py) —
# resolving vibrational frequencies needs much finer time sampling than
# structural analysis does. Same defaults as distribution_submit.slurm.
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
# Thread count
############################

N_THREADS=${1:-$(nproc)}
export OMP_NUM_THREADS=$N_THREADS
export NUMBA_NUM_THREADS=${OMP_NUM_THREADS}   # numba ignores OMP_NUM_THREADS
# Thread control differs per backend, and OMP_NUM_THREADS alone does NOT reach
# all of them:
#   freud/TBB (rdf_freud.py, bad_freud.py) ignores OMP_NUM_THREADS entirely and
#     would take every core; those scripts now call freud.parallel.set_num_threads()
#     with this value. Measured: 2.65x at 4 threads (66% efficiency) but only
#     3.30x at 12 (28%), so 2-4 threads is the efficient range.
#   numba (dsf.py via dynasor) reads NUMBA_NUM_THREADS, exported below.
#   OpenBLAS/MKL (vdos.py, msd.py, vdos_dynmat.py) does read OMP_NUM_THREADS.
#   scipy.fft (vdos.py fft_periodogram) reads VDOS_THREADS, set below.
# vdos.py's fft_periodogram method reads VDOS_THREADS directly (scipy.fft
# workers); its default vacf_cosine_transform method instead benefits from
# OMP_NUM_THREADS above via numpy's underlying BLAS matmul.
export VDOS_THREADS=$N_THREADS
echo "Threads: $OMP_NUM_THREADS"

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
# dsf.py takes its frame spacing from DYNAMICS_DT above, the same key vdos.py
# and msd.py read: all three analyse dynamics.lammpstrj, so its dt is one
# number. Only the dynamic S(q,w) path uses it; static S(q) has no time axis.
DSF_N_FRAMES="${DSF_N_FRAMES:-}"        # frame_stop INDEX into the dump, not a count;
                                        # frames used = DSF_N_FRAMES / DSF_STRIDE.
                                        # BLANK or 0 = the whole trajectory
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
RDF_RESOLUTION_MODE="${RDF_RESOLUTION_MODE:-}"    # gaussian (default) | lorch
RDF_LORCH_QMAX="${RDF_LORCH_QMAX:-}"              # Å⁻¹; needed by mode=lorch
RDF_LORCH_DR="${RDF_LORCH_DR:-}"                  # Å; the parameter INSIDE M(Q), not the
                                        # resolution a paper quotes. Blank = pi/Q_MAX, the
                                        # standard Lorch choice — leave it blank unless the
                                        # paper states this parameter explicitly
RDF_ATOMS_PER_FORMULA_UNIT="${RDF_ATOMS_PER_FORMULA_UNIT:-}"  # SiO2 -> 3; needed by the
                                        # 'formula' normalization and by RDF_WRIGHT
# Wright comparison output: one extra pair of files, <date>_wright.csv/.png, holding
# T(r) Lorch-broadened and per formula unit — built to overlay directly on a published
# neutron correlation function. Needs QMAX and ATOMS_PER_FORMULA_UNIT above.
RDF_WRIGHT="${RDF_WRIGHT:-}"            # yes | no (default no)
RDF_WRIGHT_QMAX="${RDF_WRIGHT_QMAX:-}"  # Å⁻¹; the paper's Fourier truncation, e.g. 45.2

# bad_freud.py — bond angle distributions. These keys are BARE rather than
# BAD_-prefixed at the script level (the submitters map --bad-elements to
# ELEMENTS, and so on); the prefix asymmetry is historical.
#
# R_CUTOFF is both the physics and the cost knob: set each pair from the FIRST
# MINIMUM of that pair's g(r), which means running rdf_freud.py first. Cost goes
# as R^6 — measured 16x going from 4 to 6 A — because coordination grows as R^3
# and triplets as its square, so a cutoff set generously "to be safe" is
# quadratically worse than it looks.
ELEMENTS="${ELEMENTS:-}"                # SEMICOLON list, e.g. "Si;O;H". Never auto-detected
R_CUTOFF="${R_CUTOFF:-}"                # "El1-El2:value;..." keys sorted alphabetically
R_MINCUT="${R_MINCUT:-}"                # same format; excludes unphysical close contacts
TRIPLET_CUTOFFS="${TRIPLET_CUTOFFS:-}"  # pipe-separated per-triplet overrides,
                                        # label:elA-elB-elC:r_max_ab:r_min_ab:r_max_cb:r_min_cb
BAD_BINS="${BAD_BINS:-}"                # bins over 0-180 deg; 180 = 1 deg, 360 = 0.5 deg

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
VDOS_PLOT_XUNIT="${VDOS_PLOT_XUNIT:-}"            # meV | THz | cm-1 | eV (plot axis only;
                                                  # the CSV always carries all four)

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
VDOS_DYNMAT_THREADS="${VDOS_DYNMAT_THREADS:-}"              # BLAS threads for the eigensolve;
                                                            # blank = OMP_NUM_THREADS above
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
            DSF_N_FRAMES DSF_STRIDE DSF_Q_MAX DSF_N_Q_BINS \
            DSF_WINDOW_SIZE DSF_WINDOW_STEP DSF_Q_MAX_DYN DSF_N_Q_BINS_DYN \
            DSF_MAX_Q_POINTS_DYN DSF_NEUTRON_WEIGHTING \
            R_MAX RDF_BINS RDF_NORMALIZATION RDF_FUNCTIONS RDF_RESOLUTION_SIGMA \
            RDF_RESOLUTION_MODE RDF_LORCH_QMAX RDF_ATOMS_PER_FORMULA_UNIT \
            RDF_WRIGHT RDF_WRIGHT_QMAX RDF_LORCH_DR \
            ELEMENTS R_CUTOFF R_MINCUT TRIPLET_CUTOFFS BAD_BINS \
            VDOS_N_FRAMES VDOS_STRIDE VDOS_CORR_LENGTH VDOS_CORR_INTERVAL \
            VDOS_MAX_FREQUENCY_EV VDOS_NUM_GRIDS VDOS_METHOD VDOS_WINDOW VDOS_NORMALIZATION VDOS_WEIGHTING VDOS_PLOT_XUNIT \
            MSD_N_FRAMES MSD_STRIDE MSD_CORR_LENGTH MSD_CORR_INTERVAL MSD_FIT_FRACTION \
            DYNMAT_FILE DYNMAT_BINARY \
            VDOS_DYNMAT_MAX_FREQUENCY VDOS_DYNMAT_XUNIT VDOS_DYNMAT_BINS \
            VDOS_DYNMAT_SMEARING VDOS_DYNMAT_MATRIX_STYLE VDOS_DYNMAT_NORMALIZATION VDOS_DYNMAT_WEIGHTING \
            VDOS_DYNMAT_PARTIAL VDOS_DYNMAT_ASR VDOS_DYNMAT_OUTPUT VDOS_DYNMAT_THREADS \
            VDOS_DYNMAT_CHARACTER DYNMAT_REF_TRAJ VDOS_DYNMAT_BRIDGE_ELEMENT \
            VDOS_DYNMAT_NEIGHBOR_ELEMENT VDOS_DYNMAT_BOND_CUTOFF; do
    if [ -n "${!_var}" ]; then
        export "$_var"
        echo "  $_var=${!_var}"
    fi
done
unset _var

############################
# Load environment
############################

source /home1/lkyamamo/venv/struc_analysis/bin/activate

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
