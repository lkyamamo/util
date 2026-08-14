#!/bin/bash
# submit_pipeline.sh — run a LAMMPS setup+trajectory job, then distribution
# analysis (rdf/bad/dsf/vdos/msd) on its dump, via sbatch --dependency.
#
# Run from inside the LAMMPS run directory; its basename is the run id.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: submit_pipeline.sh [options]

Run from inside the LAMMPS run directory, whose name is the run id
(e.g. cd .../runs/water/N5184/0149 && /path/to/submit_pipeline.sh ...).

Requires submit_pipeline.conf next to this script (same directory), which
sets defaults for every parameter below — the script refuses to run without
it. Copy submit_pipeline.conf.example to submit_pipeline.conf and edit it
for your setup; submit_pipeline.conf itself is gitignored since it holds
your own local paths. Any flag below overrides its config value for that
one invocation.

Required:
  --input-script FILE          LAMMPS .input script (setup + trajectory generation)
  --starting-structure FILE    Starting structure .data file, symlinked as start.data
  --potential-file FILE        Potential file, symlinked into input_files/
  --analysis-parent-dir DIR    Where <run_id>_distribution_analysis/ is created

Note: --potential-link-name NAME (optional but often required in practice)
sets the symlink's name in input_files/ (e.g. "OH.usc"). It must match
whatever your in.input's pair_coeff line expects, which frequently differs
from --potential-file's own basename (e.g. a real potential file named
20260723_OH.vashishta but pair_coeff expects OH.usc). Defaults to
--potential-file's basename if omitted.

Note: --input-script's in.input must write two trajectory dumps, flat (no
subdirectory): dump.lammpstrj (read by rdf_freud.py/bad_freud.py) and a
separate, higher-frequency dynamics.lammpstrj (read by dsf.py/vdos.py/msd.py —
resolving vibrational frequencies and diffusion needs much finer time
sampling than structural analysis does). See OH-therm.input/b-SiO-therm.input
in this directory for examples of both dump commands.

Optional:
  --analysis-template-dir DIR   Distribution-analysis scripts to copy into stage 2
                                 (default: <repo>/analysis/distributions/20260608_GrNrBaSqw)
  --run-dsf {0|1}                Run dsf.py (static/dynamic structure factor) (default: 1)
  --run-rdf {0|1}                Run rdf_freud.py (radial distribution function) (default: 1)
  --run-bad {0|1}                Run bad_freud.py (bond angle distribution) (default: 1)
  --run-vdos {0|1}                Run vdos.py (vibrational density of states) (default: 0)
  --run-msd {0|1}                 Run msd.py (mean square displacement / diffusion) (default: 0)
  --run-vdos-dynmat {0|1}          Run the dynamical-matrix VDOS (default: 0). This one
                                   spans BOTH stages: it sets LAMMPS's RUN_DYNMAT=1, which
                                   appends a minimization + dynamical_matrix block to the
                                   end of in.input (acting on the final MD configuration),
                                   and then runs vdos_dynmat.py on the resulting matrix.
                                   Needs a LAMMPS build with the PHONON package — see
                                   --lmp-bin.

  --lmp-bin PATH                   LAMMPS executable (default:
                                   /home1/lkyamamo/executables/lammps/lmp_mpi_phonon_2019,
                                   set in jobs/slurm/lammps_submit.slurm). The default is
                                   the phonon build because dynamical_matrix needs the
                                   PHONON package; point this at lmp_mpi_shock_2019 to go
                                   back to the previous binary.

  Physics config — each overrides that script's own hardcoded default
  (see the CONFIGURATION block at the top of each .py) when passed; omit
  to leave the script's built-in default in effect:
    --rdf-r-max FLOAT            rdf_freud.py R_MAX (Å), e.g. 20.0
    --rdf-bins INT                rdf_freud.py BINS, e.g. 2000
    --rdf-normalization STR       rdf_freud.py RDF_NORMALIZATION, SEMICOLON-separated
                                    pair-weight conventions from 'unity', 'FZ',
                                    'absolute', e.g. "FZ;absolute;unity"
    --rdf-functions STR           rdf_freud.py RDF_FUNCTIONS, SEMICOLON-separated
                                    correlation functions from 'g', 'h', 'D', 'T',
                                    e.g. "g;h;D" (h with absolute is Soper's G_n(r))
    --rdf-resolution-sigma FLOAT  rdf_freud.py RDF_RESOLUTION_SIGMA (Å), Gaussian
                                    resolution broadening; 0 disables, default 0.1
    --bad-elements STR            bad_freud.py ELEMENTS, SEMICOLON-separated, e.g. "Si;O;H"
    --bad-r-cutoff STR             bad_freud.py R_CUTOFF, semicolon-separated pair:value
                                    entries, e.g. "H-H:2.0;H-O:1.4;O-O:2.8"
    --bad-r-mincut STR             bad_freud.py R_MINCUT, same format, e.g.
                                    "H-H:0.5;H-O:0.5;O-O:0.5"
    --bad-triplet-cutoffs STR       bad_freud.py TRIPLET_CUTOFFS, pipe-separated entries of
                                    colon-separated fields
                                    label:elA-elB-elC:r_max_ab:r_min_ab:r_max_cb:r_min_cb
                                    (leave a cutoff field empty to fall back to
                                    R_CUTOFF/R_MINCUT), e.g.
                                    "O-Si-O:O-Si-O:2.2:0.5:2.2:0.5|H-O-H:H-O-H:1.4:0.5:1.4:0.5"
    --bad-bins INT                bad_freud.py BINS, e.g. 180
    --dsf-n-frames INT             dsf.py N_FRAMES, e.g. 500
    --dsf-stride INT               dsf.py STRIDE, e.g. 1
    --dsf-neutron-weighting STR    dsf.py DSF_NEUTRON_WEIGHTING: 'yes' (default) or 'no'.
                                     'no' drops the neutron-weighted S(q) columns; required
                                     for H/D-bearing systems, which are refused outright.
    --dsf-window-size INT           dsf.py WINDOW_SIZE, e.g. 500
    --dsf-q-max FLOAT              dsf.py Q_MAX (Å⁻¹), e.g. 20.0
    --dsf-n-q-bins INT              dsf.py N_Q_BINS, e.g. 200

    --dynamics-dt FLOAT             dt of dynamics.lammpstrj in fs, e.g. 2.0 — the
                                    ONLY value dsf.py, vdos.py and msd.py share, since it
                                    describes the trajectory rather than either
                                    analysis. Every other flag below is per script.

    vdos.py — every flag below is independent of the msd.py flags that follow,
    even where they name the same quantity (VDOS wants a short correlation
    length for frequency resolution, MSD a long one for the diffusive regime):
    --vdos-n-frames INT             vdos.py N_FRAMES (max frames to read; 0 = all)
    --vdos-stride INT               vdos.py STRIDE (read every Nth frame)
    --vdos-corr-length FLOAT        vdos.py CORR_LENGTH (fs; VACF max lag /
                                    Welch-segment length), e.g. 5000
    --vdos-corr-interval FLOAT      vdos.py CORR_INTERVAL (fs; spacing between
                                    VACF reference frames / segment starts), e.g. 500
    --vdos-max-frequency-ev FLOAT   vdos.py MAX_FREQUENCY_EV (eV), e.g. 0.5
    --vdos-num-grids INT             vdos.py NUM_GRIDS (frequency grid points), e.g. 5000
    --vdos-method STR                vdos.py METHOD: 'vacf_cosine_transform' (default,
                                    matches analysis/dynamics/src/msd.cpp) or 'fft_periodogram'
    --vdos-window STR                vdos.py WINDOW: 'cosine_lag'/'none' under
                                    vacf_cosine_transform, 'hann'/'none' under fft_periodogram
    --vdos-normalization STR         vdos.py VDOS_NORMALIZATION sum rule: 'phonon'
                                       (default) or 'unit_area'
    --vdos-weighting STR             vdos.py VDOS_WEIGHTING, SEMICOLON-separated species
                                       weights from 'unity', 'coherent', 'incoherent',
                                       'total', e.g. "unity;total" (default "unity").
                                       Neutron weightings refuse H/D-bearing systems.

    msd.py — same trajectory as vdos.py, but set independently: MSD wants a long
    correlation length to reach the diffusive regime, VDOS a short one for
    frequency resolution.
    --msd-n-frames INT               msd.py N_FRAMES (max frames to read; 0 = all)
    --msd-stride INT                 msd.py STRIDE (read every Nth frame)
    --msd-corr-length FLOAT          msd.py CORR_LENGTH (fs; max time lag), e.g. 5000.
                                    Defaults to 75% of the trajectory if omitted —
                                    which also sets the diffusion fit window, so D
                                    depends on it
    --msd-corr-interval FLOAT        msd.py CORR_INTERVAL (fs; spacing between
                                    reference frames), e.g. 500
    --msd-fit-fraction FLOAT         msd.py FIT_FRACTION: fraction of the tail of the
                                    correlation window used for the diffusion-coefficient
                                    linear fit, e.g. 0.5

    Dynamical-matrix VDOS (--run-vdos-dynmat 1). The --dynmat-* flags drive the
    LAMMPS stage and reach in.input as -var; the --vdos-dynmat-* flags drive
    vdos_dynmat.py in the analysis stage. Omit any to leave the default in the
    .input file / the .py's CONFIGURATION block in effect.

    Cost before you enable this: the finite-difference loop is 6N force
    evaluations and the matrix is (3N)^2. At replicate 6 6 6 (N=5184) that is
    31104 force evaluations and a 15552x15552 matrix — ~1.9 GB, ~2.9 GB as text.
    Budget stage-1 --time accordingly and consider --dynmat-binary yes.

    --dynmat-min-style STR           min_style for the pre-dynmat minimization (default cg)
    --dynmat-min-etol FLOAT          minimize etol    (default 1.0e-12)
    --dynmat-min-ftol FLOAT          minimize ftol    (default 1.0e-12) — the one that
                                     matters most; residual forces become spurious
                                     imaginary modes in the spectrum
    --dynmat-min-maxiter INT         minimize maxiter (default 100000)
    --dynmat-min-maxeval INT         minimize maxeval (default 1000000)
    --dynmat-displacement FLOAT      finite-difference displacement in Angstrom
                                     (default 0.0001). Too small is numerical noise,
                                     too large samples anharmonicity — the spectrum
                                     should be stable across 1e-5..1e-3
    --dynmat-file NAME               matrix filename written by LAMMPS and read by
                                     vdos_dynmat.py (default dynmat.dat)
    --dynmat-binary {yes|no}         write raw float64 instead of text (default no).
                                     Halves the file and removes the text parse; both
                                     stages read this same flag, so they cannot disagree

    --vdos-dynmat-max-frequency FLOAT  REQUIRED when --run-vdos-dynmat 1. Upper limit of
                                     the DOS grid, in --vdos-dynmat-xunit's unit, e.g. 500
                                     (meV). No default: too low silently truncates
    --vdos-dynmat-xunit STR          meV | THz | cm-1 | eV (default meV). The CSV always
                                     carries all four; this sets the plot axis and the
                                     unit --vdos-dynmat-max-frequency/-smearing are in
    --vdos-dynmat-bins INT           frequency grid points (default 500)
    --vdos-dynmat-smearing FLOAT     Gaussian FWHM in xunit; 0 = plain histogram
                                     (default 0)
    --vdos-dynmat-matrix-style STR   regular | eskm (default regular) — must match the
                                     dynamical_matrix style in the .input file, since it
                                     sets the eigenvalue-to-frequency conversion
    --vdos-dynmat-normalization STR  phonon (default; integral = 3 per atom) or unit_area
    --vdos-dynmat-weighting STR      vdos_dynmat.py VDOS_DYNMAT_WEIGHTING, SEMICOLON-
                                       separated species weights from 'unity', 'coherent',
                                       'incoherent', 'total' (default "unity"). Needs
                                       --vdos-dynmat-partial yes; refuses H/D systems.
    --vdos-dynmat-partial {yes|no}   per-element partial DOS (default yes). 'no' uses
                                     eigvalsh instead of eigh, halving memory and runtime
    --vdos-dynmat-asr STR            none (default) or simple — impose the acoustic sum
                                     rule, pushing the 3 acoustic modes to exactly zero
    --vdos-dynmat-output NAME        output basename (default vdos_dynmat), giving
                                     <date>_<name>.csv and <date>_<name>.png
  None of these use commas (sbatch --export is comma-delimited and silently
  truncates any value containing one), so no quoting/encoding is needed
  beyond normal shell quoting of the whole flag value.

  --nodes N               --analysis-nodes N
  --ntasks N               --analysis-ntasks N
  --time HH:MM:SS           --analysis-time HH:MM:SS
  --job-name NAME           --analysis-job-name NAME
  --constraint NAME         --analysis-constraint NAME
  --nodelist NODES          --analysis-nodelist NODES
                           --analysis-cpus-per-task N
      Overrides for the corresponding #SBATCH directives — the left column
      overrides jobs/slurm/lammps_submit.slurm (the LAMMPS run), the right
      column overrides distribution_submit.slurm (the analysis job). Omit any
      of these to leave the template's own value in effect. --nodelist pins
      the job to specific compute node(s) (SLURM's --nodelist, e.g. "b19-05"
      or "b19-[05-07]"). --analysis-cpus-per-task has no LAMMPS-side
      counterpart; it also drives OMP_NUM_THREADS for dsf.py/rdf_freud.py/
      bad_freud.py and VDOS_THREADS for vdos.py (distribution_submit.slurm
      sets both from $SLURM_CPUS_PER_TASK).

  --skip-trajectory                Skip stage 1 (LAMMPS run) entirely and assume the
                                  trajectory already exists at run/dump.lammpstrj and
                                  run/dynamics.lammpstrj under this run directory
                                  (i.e. you already ran this pipeline, or LAMMPS
                                  directly, from here). Verified before stage 2 runs;
                                  the pipeline aborts with an error if either file is
                                  missing. --input-script/--starting-structure/
                                  --potential-file and the trajectory-creation slurm
                                  overrides (--nodes/--ntasks/--time/--job-name/
                                  --constraint/--nodelist) are ignored in this mode,
                                  and stage 2 is submitted without a --dependency
                                  (nothing to wait on).

  --force REASON                  Overwrite existing input_files/run/ (stage 1), or redo
                                  a stage-2 calculation whose output is already in
                                  <run_id>_distribution_analysis/, instead of refusing to
                                  run. REASON is required (a short explanation of why
                                  you're overwriting) and is logged, with a timestamp and
                                  the exact paths involved, to both stderr and
                                  overwrite.log in the run directory (LOG_FILE in
                                  submit_pipeline.conf).

  --clean                         Delete the whole <run_id>_distribution_analysis/
                                  directory before stage 2, rather than adding to it.
                                  Off by default. The deletion is logged to
                                  overwrite.log like --force; pass --force REASON
                                  alongside it to record why.

  --interactive                   Run both stages in the foreground via plain
                                  bash instead of submitting them with sbatch.
                                  Assumes you already have an interactive
                                  allocation (e.g. via salloc) — this does not
                                  request one. --nodes/--time/--job-name/
                                  --constraint/--nodelist (and their
                                  --analysis-* counterparts) are ignored in
                                  this mode, since no new job is submitted.
                                  --ntasks and --analysis-cpus-per-task still
                                  apply (default 64 if not given) since they
                                  drive srun -n and OMP_NUM_THREADS.

  -h, --help                    Show this help

RE-RUNNING STAGE 2
An existing <run_id>_distribution_analysis/ is added to, never replaced. The
usual way to add a calculation to a finished run is:

  submit_pipeline.sh --skip-trajectory --run-rdf 0 --run-bad 0 --run-vdos 1

which leaves the earlier rdf/bad output untouched and writes only the vdos
files. No --force is needed, because nothing existing is at risk. --force is
needed only to REDO a calculation whose output is already in the directory,
and --clean only to throw the directory away and start over. Outputs are
date-stamped (YYYYMMDD_rdfs.csv), so a forced redo on a later date lands
beside the old copy instead of replacing it; a redo on the same date replaces
it.

All terminal output from this script (not the SLURM jobs themselves) is
also appended to submit_pipeline.log in the run directory (cwd) each time
it's invoked.
EOF
}

# The values that cannot live in submit_pipeline.conf, because the conf is
# sourced with them already set and writes its own paths in terms of them.
# Everything else that used to sit here is a conf key now — see the "pipeline
# behaviour" and "pipeline paths" blocks in submit_pipeline.conf.example.
REPO_ROOT="$HOME/util"
# The run directory. This script is run from inside the LAMMPS run directory and
# its basename is the run id, so both are known before anything is read; they are
# set here, ahead of the conf, so conf values can be written relative to the run
# being submitted (e.g. INPUT_SCRIPT="$RUN_DIR/OH-therm.input"). STAGE1_DIR is
# the same path under the name the two-stage code below uses throughout.
RUN_DIR="$(pwd)"
RUN_ID="$(basename "$RUN_DIR")"
STAGE1_DIR="$RUN_DIR"

log_overwrite() {
  local msg
  msg="[$(date "+%Y-%m-%dT%H:%M:%S%z")] $1"
  echo "WARNING: $msg" >&2
  echo "$msg" >> "$LOG_FILE"
}

# The output files each calculation writes, WITHOUT the YYYYMMDD_ prefix every
# script prepends (see _dated() in each .py). Used to tell "this directory
# already holds an rdf" from "this directory holds only a vdos", so adding a
# calculation to an existing analysis directory needs no --force but redoing
# one does. Keep in sync with the OUTPUT_* constants in the .py files.
outputs_for() {
  case "$1" in
    dsf)  echo "sq.csv sq.png dsf.csv dsf.png" ;;
    rdf)  echo "rdfs.csv rdfs.png nrs.csv nrs.png" ;;
    bad)  echo "bads.csv bads.png" ;;
    vdos) echo "vdos.csv vdos.png" ;;
    msd)  echo "msd.csv msd.png" ;;
    vdos_dynmat)
      local b="${VDOS_DYNMAT_OUTPUT:-vdos_dynmat}"
      echo "$b.csv $b.png ${b}_modes.csv ${b}_character.png" ;;
  esac
}

# Input links in a reused analysis directory are usually already there from the
# earlier run. Refresh a symlink (its target can legitimately change), and
# refuse to touch a real file — a trajectory someone copied in by hand is data,
# not a link this script owns.
link_input() {
  local target="$1" linkname="$2"
  if [[ -L "$linkname" ]]; then
    ln -sfn "$target" "$linkname"
  elif [[ -e "$linkname" ]]; then
    echo "Error: $linkname already exists and is not a symlink — refusing to replace it." >&2
    echo "Move it aside, or rebuild the directory from scratch with --clean." >&2
    exit 1
  else
    ln -s "$target" "$linkname"
  fi
}

# All parameter defaults (general input, trajectory-creation slurm, analysis
# slurm, distribution code) live in submit_pipeline.conf next to this script,
# not in the script itself, so personal paths never need to be edited into
# (or committed from) tracked code. See submit_pipeline.conf.example.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="$SCRIPT_DIR/submit_pipeline.conf"
if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "Error: required config file not found: $CONFIG_FILE" >&2
  echo "Copy $SCRIPT_DIR/submit_pipeline.conf.example to $CONFIG_FILE and edit it for your setup." >&2
  exit 1
fi
# Pre-declare the per-script analysis keys as empty so that a conf predating
# any of them still works under `set -u` (an unset key would otherwise abort
# with "unbound variable"). The conf overrides whichever of these it sets.
#
# The behaviour switches and pipeline paths below are the same idea, but their
# fallbacks are real values rather than empty: they used to be assigned at the
# top of this script and are conf keys now, so these lines are what a conf
# written before they moved falls back to.
INTERACTIVE="${INTERACTIVE:-0}"
SKIP_TRAJECTORY="${SKIP_TRAJECTORY:-0}"
FORCE="${FORCE:-0}"
FORCE_REASON="${FORCE_REASON:-}"
CLEAN="${CLEAN:-0}"
LAMMPS_TEMPLATE="${LAMMPS_TEMPLATE:-$REPO_ROOT/jobs/slurm/lammps_submit.slurm}"
DUMP_FILE="${DUMP_FILE:-dump.lammpstrj}"
DYNAMICS_DUMP_FILE="${DYNAMICS_DUMP_FILE:-dynamics.lammpstrj}"
LOG_FILE="${LOG_FILE:-$RUN_DIR/overwrite.log}"
DYNAMICS_DT="${DYNAMICS_DT:-}"
RDF_R_MAX="${RDF_R_MAX:-}"
RDF_BINS_VAL="${RDF_BINS_VAL:-}"
RDF_NORMALIZATION_VAL="${RDF_NORMALIZATION_VAL:-}"
RDF_FUNCTIONS_VAL="${RDF_FUNCTIONS_VAL:-}"
RDF_RESOLUTION_SIGMA_VAL="${RDF_RESOLUTION_SIGMA_VAL:-}"
DSF_NEUTRON_WEIGHTING="${DSF_NEUTRON_WEIGHTING:-}"
VDOS_N_FRAMES="${VDOS_N_FRAMES:-}"
VDOS_STRIDE="${VDOS_STRIDE:-}"
VDOS_CORR_LENGTH="${VDOS_CORR_LENGTH:-}"
VDOS_CORR_INTERVAL="${VDOS_CORR_INTERVAL:-}"
VDOS_MAX_FREQUENCY_EV="${VDOS_MAX_FREQUENCY_EV:-}"
VDOS_NUM_GRIDS="${VDOS_NUM_GRIDS:-}"
VDOS_METHOD="${VDOS_METHOD:-}"
VDOS_WINDOW="${VDOS_WINDOW:-}"
VDOS_NORMALIZATION="${VDOS_NORMALIZATION:-}"
VDOS_WEIGHTING="${VDOS_WEIGHTING:-}"
MSD_N_FRAMES="${MSD_N_FRAMES:-}"
MSD_STRIDE="${MSD_STRIDE:-}"
MSD_CORR_LENGTH="${MSD_CORR_LENGTH:-}"
MSD_CORR_INTERVAL="${MSD_CORR_INTERVAL:-}"
MSD_FIT_FRACTION="${MSD_FIT_FRACTION:-}"
RUN_VDOS_DYNMAT="${RUN_VDOS_DYNMAT:-0}"
LMP_BIN="${LMP_BIN:-}"
DYNMAT_MIN_STYLE="${DYNMAT_MIN_STYLE:-}"
DYNMAT_MIN_ETOL="${DYNMAT_MIN_ETOL:-}"
DYNMAT_MIN_FTOL="${DYNMAT_MIN_FTOL:-}"
DYNMAT_MIN_MAXITER="${DYNMAT_MIN_MAXITER:-}"
DYNMAT_MIN_MAXEVAL="${DYNMAT_MIN_MAXEVAL:-}"
DYNMAT_DISPLACEMENT="${DYNMAT_DISPLACEMENT:-}"
DYNMAT_FILE="${DYNMAT_FILE:-}"
DYNMAT_BINARY="${DYNMAT_BINARY:-}"
VDOS_DYNMAT_MAX_FREQUENCY="${VDOS_DYNMAT_MAX_FREQUENCY:-}"
VDOS_DYNMAT_XUNIT="${VDOS_DYNMAT_XUNIT:-}"
VDOS_DYNMAT_BINS="${VDOS_DYNMAT_BINS:-}"
VDOS_DYNMAT_SMEARING="${VDOS_DYNMAT_SMEARING:-}"
VDOS_DYNMAT_MATRIX_STYLE="${VDOS_DYNMAT_MATRIX_STYLE:-}"
VDOS_DYNMAT_NORMALIZATION="${VDOS_DYNMAT_NORMALIZATION:-}"
VDOS_DYNMAT_WEIGHTING="${VDOS_DYNMAT_WEIGHTING:-}"
VDOS_DYNMAT_PARTIAL="${VDOS_DYNMAT_PARTIAL:-}"
VDOS_DYNMAT_ASR="${VDOS_DYNMAT_ASR:-}"
VDOS_DYNMAT_OUTPUT="${VDOS_DYNMAT_OUTPUT:-}"
VDOS_DYNMAT_CHARACTER="${VDOS_DYNMAT_CHARACTER:-}"
DYNMAT_REF_TRAJ="${DYNMAT_REF_TRAJ:-}"
VDOS_DYNMAT_BRIDGE_ELEMENT="${VDOS_DYNMAT_BRIDGE_ELEMENT:-}"
VDOS_DYNMAT_NEIGHBOR_ELEMENT="${VDOS_DYNMAT_NEIGHBOR_ELEMENT:-}"
VDOS_DYNMAT_BOND_CUTOFF="${VDOS_DYNMAT_BOND_CUTOFF:-}"

# shellcheck source=/dev/null
source "$CONFIG_FILE"

# VDOS_DT used to carry the dt for both scripts; it is now DYNAMICS_DT. Because
# submit_pipeline.conf is gitignored it does not travel with a pull, so a conf
# still carrying VDOS_DT would otherwise be silently ignored.
if [[ -n "${VDOS_DT:-}" ]]; then
  echo "Error: VDOS_DT is no longer used — it has been renamed to DYNAMICS_DT" >&2
  echo "(flag --dynamics-dt), the one parameter vdos.py and msd.py still share." >&2
  echo "Rename it in $CONFIG_FILE; see submit_pipeline.conf.example." >&2
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-script) INPUT_SCRIPT="$2"; shift 2 ;;
    --starting-structure) STARTING_STRUCTURE="$2"; shift 2 ;;
    --potential-file) POTENTIAL_FILE="$2"; shift 2 ;;
    --potential-link-name) POTENTIAL_LINK_NAME="$2"; shift 2 ;;
    --analysis-parent-dir) ANALYSIS_PARENT_DIR="$2"; shift 2 ;;
    --analysis-template-dir) ANALYSIS_TEMPLATE_DIR="$2"; shift 2 ;;
    --run-dsf) RUN_DSF="$2"; shift 2 ;;
    --run-rdf) RUN_RDF="$2"; shift 2 ;;
    --run-bad) RUN_BAD="$2"; shift 2 ;;
    --run-vdos) RUN_VDOS="$2"; shift 2 ;;
    --run-msd) RUN_MSD="$2"; shift 2 ;;
    --nodes) NODES="$2"; shift 2 ;;
    --ntasks) NTASKS="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --job-name) JOB_NAME="$2"; shift 2 ;;
    --constraint) CONSTRAINT="$2"; shift 2 ;;
    --nodelist) NODELIST="$2"; shift 2 ;;
    --analysis-nodes) ANALYSIS_NODES="$2"; shift 2 ;;
    --analysis-ntasks) ANALYSIS_NTASKS="$2"; shift 2 ;;
    --analysis-time) ANALYSIS_TIME="$2"; shift 2 ;;
    --analysis-job-name) ANALYSIS_JOB_NAME="$2"; shift 2 ;;
    --analysis-constraint) ANALYSIS_CONSTRAINT="$2"; shift 2 ;;
    --analysis-nodelist) ANALYSIS_NODELIST="$2"; shift 2 ;;
    --analysis-cpus-per-task) ANALYSIS_CPUS_PER_TASK="$2"; shift 2 ;;
    --rdf-r-max) RDF_R_MAX="$2"; shift 2 ;;
    --rdf-bins) RDF_BINS_VAL="$2"; shift 2 ;;
    --rdf-normalization) RDF_NORMALIZATION_VAL="$2"; shift 2 ;;
    --rdf-functions) RDF_FUNCTIONS_VAL="$2"; shift 2 ;;
    --rdf-resolution-sigma) RDF_RESOLUTION_SIGMA_VAL="$2"; shift 2 ;;
    --bad-elements) BAD_ELEMENTS="$2"; shift 2 ;;
    --bad-r-cutoff) BAD_R_CUTOFF="$2"; shift 2 ;;
    --bad-r-mincut) BAD_R_MINCUT="$2"; shift 2 ;;
    --bad-triplet-cutoffs) BAD_TRIPLET_CUTOFFS="$2"; shift 2 ;;
    --bad-bins) BAD_BINS_VAL="$2"; shift 2 ;;
    --dsf-n-frames) DSF_N_FRAMES="$2"; shift 2 ;;
    --dsf-stride) DSF_STRIDE="$2"; shift 2 ;;
    --dsf-neutron-weighting) DSF_NEUTRON_WEIGHTING="$2"; shift 2 ;;
    --dsf-window-size) DSF_WINDOW_SIZE="$2"; shift 2 ;;
    --dsf-q-max) DSF_Q_MAX="$2"; shift 2 ;;
    --dsf-n-q-bins) DSF_N_Q_BINS="$2"; shift 2 ;;
    --dynamics-dt) DYNAMICS_DT="$2"; shift 2 ;;
    --vdos-n-frames) VDOS_N_FRAMES="$2"; shift 2 ;;
    --vdos-stride) VDOS_STRIDE="$2"; shift 2 ;;
    --vdos-corr-length) VDOS_CORR_LENGTH="$2"; shift 2 ;;
    --vdos-corr-interval) VDOS_CORR_INTERVAL="$2"; shift 2 ;;
    --vdos-max-frequency-ev) VDOS_MAX_FREQUENCY_EV="$2"; shift 2 ;;
    --vdos-num-grids) VDOS_NUM_GRIDS="$2"; shift 2 ;;
    --vdos-method) VDOS_METHOD="$2"; shift 2 ;;
    --vdos-window) VDOS_WINDOW="$2"; shift 2 ;;
    --vdos-normalization) VDOS_NORMALIZATION="$2"; shift 2 ;;
    --vdos-weighting) VDOS_WEIGHTING="$2"; shift 2 ;;
    --msd-n-frames) MSD_N_FRAMES="$2"; shift 2 ;;
    --msd-stride) MSD_STRIDE="$2"; shift 2 ;;
    --msd-corr-length) MSD_CORR_LENGTH="$2"; shift 2 ;;
    --msd-corr-interval) MSD_CORR_INTERVAL="$2"; shift 2 ;;
    --msd-fit-fraction) MSD_FIT_FRACTION="$2"; shift 2 ;;
    --run-vdos-dynmat) RUN_VDOS_DYNMAT="$2"; shift 2 ;;
    --lmp-bin) LMP_BIN="$2"; shift 2 ;;
    --dynmat-min-style) DYNMAT_MIN_STYLE="$2"; shift 2 ;;
    --dynmat-min-etol) DYNMAT_MIN_ETOL="$2"; shift 2 ;;
    --dynmat-min-ftol) DYNMAT_MIN_FTOL="$2"; shift 2 ;;
    --dynmat-min-maxiter) DYNMAT_MIN_MAXITER="$2"; shift 2 ;;
    --dynmat-min-maxeval) DYNMAT_MIN_MAXEVAL="$2"; shift 2 ;;
    --dynmat-displacement) DYNMAT_DISPLACEMENT="$2"; shift 2 ;;
    --dynmat-file) DYNMAT_FILE="$2"; shift 2 ;;
    --dynmat-binary) DYNMAT_BINARY="$2"; shift 2 ;;
    --vdos-dynmat-max-frequency) VDOS_DYNMAT_MAX_FREQUENCY="$2"; shift 2 ;;
    --vdos-dynmat-xunit) VDOS_DYNMAT_XUNIT="$2"; shift 2 ;;
    --vdos-dynmat-bins) VDOS_DYNMAT_BINS="$2"; shift 2 ;;
    --vdos-dynmat-smearing) VDOS_DYNMAT_SMEARING="$2"; shift 2 ;;
    --vdos-dynmat-matrix-style) VDOS_DYNMAT_MATRIX_STYLE="$2"; shift 2 ;;
    --vdos-dynmat-normalization) VDOS_DYNMAT_NORMALIZATION="$2"; shift 2 ;;
    --vdos-dynmat-weighting) VDOS_DYNMAT_WEIGHTING="$2"; shift 2 ;;
    --vdos-dynmat-partial) VDOS_DYNMAT_PARTIAL="$2"; shift 2 ;;
    --vdos-dynmat-asr) VDOS_DYNMAT_ASR="$2"; shift 2 ;;
    --vdos-dynmat-output) VDOS_DYNMAT_OUTPUT="$2"; shift 2 ;;
    --vdos-dynmat-character) VDOS_DYNMAT_CHARACTER="$2"; shift 2 ;;
    --dynmat-ref-traj) DYNMAT_REF_TRAJ="$2"; shift 2 ;;
    --vdos-dynmat-bridge-element) VDOS_DYNMAT_BRIDGE_ELEMENT="$2"; shift 2 ;;
    --vdos-dynmat-neighbor-element) VDOS_DYNMAT_NEIGHBOR_ELEMENT="$2"; shift 2 ;;
    --vdos-dynmat-bond-cutoff) VDOS_DYNMAT_BOND_CUTOFF="$2"; shift 2 ;;
    --force) FORCE="1"; FORCE_REASON="$2"; shift 2 ;;
    --clean) CLEAN="1"; shift 1 ;;
    --interactive) INTERACTIVE="1"; shift 1 ;;
    --skip-trajectory) SKIP_TRAJECTORY="1"; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

PIPELINE_LOG="$STAGE1_DIR/submit_pipeline.log"
exec > >(tee -a "$PIPELINE_LOG") 2>&1
echo "=== submit_pipeline.sh started $(date "+%Y-%m-%dT%H:%M:%S%z") — run id: $RUN_ID ==="
echo "Logging this run's output to: $PIPELINE_LOG"

missing=0
require() {
  if [[ -z "$1" ]]; then
    echo "Missing required argument: $2" >&2
    missing=1
  fi
}
if [[ "$SKIP_TRAJECTORY" != "1" ]]; then
  require "$INPUT_SCRIPT" --input-script
  require "$STARTING_STRUCTURE" --starting-structure
  require "$POTENTIAL_FILE" --potential-file
fi
require "$ANALYSIS_PARENT_DIR" --analysis-parent-dir
if [[ "$missing" -ne 0 ]]; then
  usage
  exit 1
fi

check_zero_or_one() {
  if [[ "$1" != "0" && "$1" != "1" ]]; then
    echo "Invalid value for $2: '$1' (must be 0 or 1)" >&2
    exit 1
  fi
}
check_zero_or_one "$RUN_DSF" --run-dsf
check_zero_or_one "$RUN_RDF" --run-rdf
check_zero_or_one "$RUN_BAD" --run-bad
check_zero_or_one "$RUN_VDOS" --run-vdos
check_zero_or_one "$RUN_MSD" --run-msd
check_zero_or_one "$RUN_VDOS_DYNMAT" --run-vdos-dynmat

# vdos.py and msd.py refuse to guess the parameters that determine their
# numbers. Catch a missing one here, before anything is submitted or run,
# rather than after the (long) LAMMPS stage has already completed.
check_required_analysis_params() {
  local missing=()
  if [[ "$RUN_VDOS" == "1" ]]; then
    [[ -n "$VDOS_CORR_LENGTH" ]]      || missing+=("  --vdos-corr-length / VDOS_CORR_LENGTH           fs; VACF max lag (sets frequency resolution)")
    [[ -n "$VDOS_CORR_INTERVAL" ]]    || missing+=("  --vdos-corr-interval / VDOS_CORR_INTERVAL       fs; spacing between VACF reference frames")
    [[ -n "$VDOS_MAX_FREQUENCY_EV" ]] || missing+=("  --vdos-max-frequency-ev / VDOS_MAX_FREQUENCY_EV eV; upper limit of the DOS grid")
  fi
  if [[ "$RUN_MSD" == "1" ]]; then
    [[ -n "$MSD_CORR_LENGTH" ]]   || missing+=("  --msd-corr-length / MSD_CORR_LENGTH             fs; max time lag")
    [[ -n "$MSD_CORR_INTERVAL" ]] || missing+=("  --msd-corr-interval / MSD_CORR_INTERVAL         fs; spacing between reference frames")
    [[ -n "$MSD_FIT_FRACTION" ]]  || missing+=("  --msd-fit-fraction / MSD_FIT_FRACTION           tail fraction used for the D fit")
  fi
  if [[ "$RUN_VDOS" == "1" || "$RUN_MSD" == "1" ]]; then
    [[ -n "$DYNAMICS_DT" ]] || missing+=("  --dynamics-dt / DYNAMICS_DT                     fs between dynamics.lammpstrj frames")
  fi
  # Checked here rather than in stage 2 for the same reason as the others, but it
  # matters more: the dynamical-matrix run is the expensive part of stage 1, and
  # discovering a missing grid limit afterwards would waste all of it.
  if [[ "$RUN_VDOS_DYNMAT" == "1" ]]; then
    [[ -n "$VDOS_DYNMAT_MAX_FREQUENCY" ]] || missing+=("  --vdos-dynmat-max-frequency / VDOS_DYNMAT_MAX_FREQUENCY  upper limit of the DOS grid, in VDOS_DYNMAT_XUNIT")
  fi
  if (( ${#missing[@]} )); then
    echo "Error: required analysis parameter(s) not set:" >&2
    printf '%s\n' "${missing[@]}" >&2
    echo >&2
    echo "These determine the numbers vdos.py/msd.py produce, so nothing guesses them." >&2
    echo "Set them in $CONFIG_FILE or pass the matching flag." >&2
    exit 1
  fi
}
check_required_analysis_params

if [[ "$FORCE" == "1" && -z "$FORCE_REASON" ]]; then
  echo 'Error: --force requires a non-empty reason, e.g. --force "re-running after fixing potential file"' >&2
  exit 1
fi

if [[ "$INTERACTIVE" == "1" && -z "${SLURM_JOB_ID:-}" ]]; then
  echo "Error: --interactive requires an active Slurm allocation (run 'salloc ...' first)." >&2
  echo "No \$SLURM_JOB_ID is set in this shell, so srun would fail to allocate resources —" >&2
  echo "and lammps_submit.slurm's unconditional 'exit 0' would silently mask that failure" >&2
  echo "instead of stopping the pipeline before stage 2." >&2
  exit 1
fi

if [[ "$SKIP_TRAJECTORY" == "1" ]]; then
  missing_traj=0
  required_stage1_files=("$STAGE1_DIR/run/$DUMP_FILE" "$STAGE1_DIR/run/$DYNAMICS_DUMP_FILE")
  # The dynamical matrix is a stage-1 product too, so skipping stage 1 means it
  # must already be there. Checked here so re-running only the post-processing
  # fails up front instead of leaving a dangling symlink for stage 2 to trip on.
  if [[ "$RUN_VDOS_DYNMAT" == "1" ]]; then
    required_stage1_files+=("$STAGE1_DIR/run/${DYNMAT_FILE:-dynmat.dat}")
  fi
  for f in "${required_stage1_files[@]}"; do
    if [[ ! -e "$f" ]]; then
      echo "Error: --skip-trajectory was given but $f does not exist." >&2
      missing_traj=1
    fi
  done
  if [[ "$missing_traj" -ne 0 ]]; then
    echo "Re-run without --skip-trajectory to generate the trajectory, or fix the path above." >&2
    exit 1
  fi
  echo "--skip-trajectory: found existing trajectory in $STAGE1_DIR/run — skipping stage 1."
else
  if [[ -e "$STAGE1_DIR/input_files" || -e "$STAGE1_DIR/run" ]]; then
    if [[ "$FORCE" == "1" ]]; then
      log_overwrite "--force ($FORCE_REASON): removing existing $STAGE1_DIR/input_files and/or $STAGE1_DIR/run (run id: $RUN_ID)"
      rm -rf "$STAGE1_DIR/input_files" "$STAGE1_DIR/run"
    else
      echo "Error: $STAGE1_DIR already has input_files/ or run/ — refusing to overwrite (use --force to override)." >&2
      exit 1
    fi
  fi

  ############################
  # Stage 1: LAMMPS run (setup + trajectory generation)
  ############################

  mkdir -p "$STAGE1_DIR/input_files" "$STAGE1_DIR/run"

  cp "$INPUT_SCRIPT" "$STAGE1_DIR/input_files/in.input"
  ln -s "$(realpath "$STARTING_STRUCTURE")" "$STAGE1_DIR/input_files/start.data"
  ln -s "$(realpath "$POTENTIAL_FILE")" "$STAGE1_DIR/input_files/${POTENTIAL_LINK_NAME:-$(basename "$POTENTIAL_FILE")}"

  cp "$LAMMPS_TEMPLATE" "$STAGE1_DIR/lammps_submit.slurm"
  cp "$LAMMPS_TEMPLATE" "$STAGE1_DIR/run/lammps_submit.slurm"
fi

############################
# Stage 2: distribution analysis (reads the trajectory as input; never modifies it)
############################

mkdir -p "$ANALYSIS_PARENT_DIR"
ANALYSIS_PARENT_DIR="$(cd "$ANALYSIS_PARENT_DIR" && pwd)"
STAGE2_DIR="$ANALYSIS_PARENT_DIR/${RUN_ID}_distribution_analysis"

# An existing analysis directory is ADDED TO, not replaced: the RUN_* flags pick
# which calculations run, and a run that only asks for rdf must not destroy the
# vdos output sitting next to it. Only --clean removes the directory, and only
# a calculation whose own output is already there needs --force.
if [[ -e "$STAGE2_DIR" && "$CLEAN" == "1" ]]; then
  log_overwrite "--clean${FORCE_REASON:+ ($FORCE_REASON)}: removing existing $STAGE2_DIR"
  rm -rf "$STAGE2_DIR"
fi

if [[ -e "$STAGE2_DIR" ]]; then
  echo "Adding to existing analysis directory $STAGE2_DIR (use --clean to rebuild it from scratch)."
  existing_outputs=()
  for calc in dsf rdf bad vdos msd vdos_dynmat; do
    run_var="RUN_$(echo "$calc" | tr '[:lower:]' '[:upper:]')"
    [[ "${!run_var}" == "1" ]] || continue
    for name in $(outputs_for "$calc"); do
      # The scripts date-stamp every output, so an earlier run's files are found
      # by a leading-date glob rather than by exact name.
      for f in "$STAGE2_DIR"/[0-9]*_"$name"; do
        [[ -e "$f" ]] && existing_outputs+=("$(basename "$f")")
      done
    done
  done

  if [[ ${#existing_outputs[@]} -gt 0 ]]; then
    if [[ "$FORCE" == "1" ]]; then
      log_overwrite "--force ($FORCE_REASON): re-running calculations whose output already exists in $STAGE2_DIR: ${existing_outputs[*]}"
    else
      echo "Error: $STAGE2_DIR already holds output for calculations this run would redo:" >&2
      printf '  %s\n' "${existing_outputs[@]}" >&2
      echo "Nothing else in that directory is at stake — only the calculations above are repeats." >&2
      echo "Options: --force REASON to redo them (a run on the same date overwrites the files" >&2
      echo "listed above; a run on a later date writes new YYYYMMDD_ copies alongside them);" >&2
      echo "--run-<name> 0 to drop the repeats and add only what is missing; or --clean to" >&2
      echo "delete the whole analysis directory and start over." >&2
      exit 1
    fi
  fi
fi

mkdir -p "$STAGE2_DIR"

cp "$ANALYSIS_TEMPLATE_DIR/dsf.py" \
   "$ANALYSIS_TEMPLATE_DIR/rdf_freud.py" \
   "$ANALYSIS_TEMPLATE_DIR/bad_freud.py" \
   "$ANALYSIS_TEMPLATE_DIR/vdos.py" \
   "$ANALYSIS_TEMPLATE_DIR/msd.py" \
   "$ANALYSIS_TEMPLATE_DIR/vdos_dynmat.py" \
   "$ANALYSIS_TEMPLATE_DIR/distribution_submit.slurm" \
   "$STAGE2_DIR/"
link_input "$STAGE1_DIR/run/$DUMP_FILE" "$STAGE2_DIR/$DUMP_FILE"
link_input "$STAGE1_DIR/run/$DYNAMICS_DUMP_FILE" "$STAGE2_DIR/$DYNAMICS_DUMP_FILE"
# The dynamical matrix, like the trajectories, is consumed read-only. No second
# link is needed for the atom->element mapping: vdos_dynmat.py takes that from
# the first frame of $DUMP_FILE above, which lists the same atoms in the same
# ID order the matrix rows use.
if [[ "$RUN_VDOS_DYNMAT" == "1" ]]; then
  link_input "$STAGE1_DIR/run/${DYNMAT_FILE:-dynmat.dat}" "$STAGE2_DIR/${DYNMAT_FILE:-dynmat.dat}"
  # The minimized geometry, needed only by the mode-character analysis. Linked
  # unconditionally so turning VDOS_DYNMAT_CHARACTER on later needs no re-run.
  _ref="${DYNMAT_REF_TRAJ:-dynmat_ref.lammpstrj}"
  link_input "$STAGE1_DIR/run/$_ref" "$STAGE2_DIR/$_ref"
fi

# None of these values may contain a comma — sbatch --export is
# comma-delimited and silently truncates anything after an embedded one.
export_vars="ALL,TRAJ=$DUMP_FILE,DYNAMICS_TRAJ=$DYNAMICS_DUMP_FILE,RUN_DSF=$RUN_DSF,RUN_RDF=$RUN_RDF,RUN_BAD=$RUN_BAD,RUN_VDOS=$RUN_VDOS,RUN_MSD=$RUN_MSD,RUN_VDOS_DYNMAT=$RUN_VDOS_DYNMAT"
[[ -n "$RDF_R_MAX" ]]           && export_vars+=",R_MAX=$RDF_R_MAX"
[[ -n "$RDF_BINS_VAL" ]]        && export_vars+=",RDF_BINS=$RDF_BINS_VAL"
# rdf_freud.py accepts ; as well as , for these two, since a comma here
# would be eaten by --export above.
[[ -n "$RDF_NORMALIZATION_VAL" ]] && export_vars+=",RDF_NORMALIZATION=$RDF_NORMALIZATION_VAL"
[[ -n "$RDF_FUNCTIONS_VAL" ]]   && export_vars+=",RDF_FUNCTIONS=$RDF_FUNCTIONS_VAL"
[[ -n "$RDF_RESOLUTION_SIGMA_VAL" ]] && export_vars+=",RDF_RESOLUTION_SIGMA=$RDF_RESOLUTION_SIGMA_VAL"
[[ -n "$BAD_ELEMENTS" ]]        && export_vars+=",ELEMENTS=$BAD_ELEMENTS"
[[ -n "$BAD_R_CUTOFF" ]]        && export_vars+=",R_CUTOFF=$BAD_R_CUTOFF"
[[ -n "$BAD_R_MINCUT" ]]        && export_vars+=",R_MINCUT=$BAD_R_MINCUT"
[[ -n "$BAD_TRIPLET_CUTOFFS" ]] && export_vars+=",TRIPLET_CUTOFFS=$BAD_TRIPLET_CUTOFFS"
[[ -n "$BAD_BINS_VAL" ]]        && export_vars+=",BAD_BINS=$BAD_BINS_VAL"
[[ -n "$DSF_N_FRAMES" ]]        && export_vars+=",DSF_N_FRAMES=$DSF_N_FRAMES"
[[ -n "$DSF_STRIDE" ]]          && export_vars+=",DSF_STRIDE=$DSF_STRIDE"
[[ -n "$DSF_NEUTRON_WEIGHTING" ]] && export_vars+=",DSF_NEUTRON_WEIGHTING=$DSF_NEUTRON_WEIGHTING"
[[ -n "$DSF_WINDOW_SIZE" ]]     && export_vars+=",DSF_WINDOW_SIZE=$DSF_WINDOW_SIZE"
[[ -n "$DSF_Q_MAX" ]]           && export_vars+=",DSF_Q_MAX=$DSF_Q_MAX"
[[ -n "$DSF_N_Q_BINS" ]]        && export_vars+=",DSF_N_Q_BINS=$DSF_N_Q_BINS"
# vdos.py and msd.py read the same dynamics.lammpstrj but want different
# settings, so each reads its own prefixed vars and nothing else. DYNAMICS_DT
# is the one and only value they share: the dt of the trajectory itself.
[[ -n "$DYNAMICS_DT" ]]           && export_vars+=",DYNAMICS_DT=$DYNAMICS_DT"
[[ -n "$VDOS_N_FRAMES" ]]         && export_vars+=",VDOS_N_FRAMES=$VDOS_N_FRAMES"
[[ -n "$VDOS_STRIDE" ]]           && export_vars+=",VDOS_STRIDE=$VDOS_STRIDE"
[[ -n "$VDOS_CORR_LENGTH" ]]      && export_vars+=",VDOS_CORR_LENGTH=$VDOS_CORR_LENGTH"
[[ -n "$VDOS_CORR_INTERVAL" ]]    && export_vars+=",VDOS_CORR_INTERVAL=$VDOS_CORR_INTERVAL"
[[ -n "$VDOS_MAX_FREQUENCY_EV" ]] && export_vars+=",VDOS_MAX_FREQUENCY_EV=$VDOS_MAX_FREQUENCY_EV"
[[ -n "$VDOS_NUM_GRIDS" ]]        && export_vars+=",VDOS_NUM_GRIDS=$VDOS_NUM_GRIDS"
[[ -n "$VDOS_METHOD" ]]           && export_vars+=",VDOS_METHOD=$VDOS_METHOD"
[[ -n "$VDOS_WINDOW" ]]           && export_vars+=",VDOS_WINDOW=$VDOS_WINDOW"
[[ -n "$VDOS_NORMALIZATION" ]]    && export_vars+=",VDOS_NORMALIZATION=$VDOS_NORMALIZATION"
# Semicolon-separated; a comma here would be eaten by --export above.
[[ -n "$VDOS_WEIGHTING" ]]       && export_vars+=",VDOS_WEIGHTING=$VDOS_WEIGHTING"
[[ -n "$MSD_N_FRAMES" ]]          && export_vars+=",MSD_N_FRAMES=$MSD_N_FRAMES"
[[ -n "$MSD_STRIDE" ]]            && export_vars+=",MSD_STRIDE=$MSD_STRIDE"
[[ -n "$MSD_CORR_LENGTH" ]]       && export_vars+=",MSD_CORR_LENGTH=$MSD_CORR_LENGTH"
[[ -n "$MSD_CORR_INTERVAL" ]]     && export_vars+=",MSD_CORR_INTERVAL=$MSD_CORR_INTERVAL"
[[ -n "$MSD_FIT_FRACTION" ]]      && export_vars+=",MSD_FIT_FRACTION=$MSD_FIT_FRACTION"
# vdos_dynmat.py. DYNMAT_FILE and DYNMAT_BINARY are the two names that also go to
# stage 1 (below) — one value each, so the writer and the reader cannot disagree
# about the filename or the text/binary format.
[[ -n "$DYNMAT_FILE" ]]                 && export_vars+=",DYNMAT_FILE=$DYNMAT_FILE"
[[ -n "$DYNMAT_BINARY" ]]               && export_vars+=",DYNMAT_BINARY=$DYNMAT_BINARY"
[[ -n "$VDOS_DYNMAT_MAX_FREQUENCY" ]]   && export_vars+=",VDOS_DYNMAT_MAX_FREQUENCY=$VDOS_DYNMAT_MAX_FREQUENCY"
[[ -n "$VDOS_DYNMAT_XUNIT" ]]           && export_vars+=",VDOS_DYNMAT_XUNIT=$VDOS_DYNMAT_XUNIT"
[[ -n "$VDOS_DYNMAT_BINS" ]]            && export_vars+=",VDOS_DYNMAT_BINS=$VDOS_DYNMAT_BINS"
[[ -n "$VDOS_DYNMAT_SMEARING" ]]        && export_vars+=",VDOS_DYNMAT_SMEARING=$VDOS_DYNMAT_SMEARING"
[[ -n "$VDOS_DYNMAT_MATRIX_STYLE" ]]    && export_vars+=",VDOS_DYNMAT_MATRIX_STYLE=$VDOS_DYNMAT_MATRIX_STYLE"
[[ -n "$VDOS_DYNMAT_NORMALIZATION" ]]   && export_vars+=",VDOS_DYNMAT_NORMALIZATION=$VDOS_DYNMAT_NORMALIZATION"
# Semicolon-separated; a comma here would be eaten by --export above.
[[ -n "$VDOS_DYNMAT_WEIGHTING" ]]      && export_vars+=",VDOS_DYNMAT_WEIGHTING=$VDOS_DYNMAT_WEIGHTING"
[[ -n "$VDOS_DYNMAT_PARTIAL" ]]         && export_vars+=",VDOS_DYNMAT_PARTIAL=$VDOS_DYNMAT_PARTIAL"
[[ -n "$VDOS_DYNMAT_ASR" ]]             && export_vars+=",VDOS_DYNMAT_ASR=$VDOS_DYNMAT_ASR"
[[ -n "$VDOS_DYNMAT_OUTPUT" ]]          && export_vars+=",VDOS_DYNMAT_OUTPUT=$VDOS_DYNMAT_OUTPUT"
[[ -n "$VDOS_DYNMAT_CHARACTER" ]]       && export_vars+=",VDOS_DYNMAT_CHARACTER=$VDOS_DYNMAT_CHARACTER"
[[ -n "$DYNMAT_REF_TRAJ" ]]             && export_vars+=",DYNMAT_REF_TRAJ=$DYNMAT_REF_TRAJ"
[[ -n "$VDOS_DYNMAT_BRIDGE_ELEMENT" ]]  && export_vars+=",VDOS_DYNMAT_BRIDGE_ELEMENT=$VDOS_DYNMAT_BRIDGE_ELEMENT"
[[ -n "$VDOS_DYNMAT_NEIGHBOR_ELEMENT" ]] && export_vars+=",VDOS_DYNMAT_NEIGHBOR_ELEMENT=$VDOS_DYNMAT_NEIGHBOR_ELEMENT"
[[ -n "$VDOS_DYNMAT_BOND_CUTOFF" ]]     && export_vars+=",VDOS_DYNMAT_BOND_CUTOFF=$VDOS_DYNMAT_BOND_CUTOFF"

# Stage 1's own environment. Until now stage 1 needed none — everything it used
# was baked into in.input — so this is the first --export it gets. RUN_VDOS_DYNMAT
# reaches LAMMPS under the name its input scripts use, RUN_DYNMAT; the rest pass
# through jobs/slurm/lammps_submit.slurm, which turns each set name into a -var.
stage1_export="ALL,RUN_DYNMAT=$RUN_VDOS_DYNMAT"
[[ -n "$LMP_BIN" ]]             && stage1_export+=",LMP_BIN=$LMP_BIN"
[[ -n "$DYNMAT_MIN_STYLE" ]]    && stage1_export+=",DYNMAT_MIN_STYLE=$DYNMAT_MIN_STYLE"
[[ -n "$DYNMAT_MIN_ETOL" ]]     && stage1_export+=",DYNMAT_MIN_ETOL=$DYNMAT_MIN_ETOL"
[[ -n "$DYNMAT_MIN_FTOL" ]]     && stage1_export+=",DYNMAT_MIN_FTOL=$DYNMAT_MIN_FTOL"
[[ -n "$DYNMAT_MIN_MAXITER" ]]  && stage1_export+=",DYNMAT_MIN_MAXITER=$DYNMAT_MIN_MAXITER"
[[ -n "$DYNMAT_MIN_MAXEVAL" ]]  && stage1_export+=",DYNMAT_MIN_MAXEVAL=$DYNMAT_MIN_MAXEVAL"
[[ -n "$DYNMAT_DISPLACEMENT" ]] && stage1_export+=",DYNMAT_DISPLACEMENT=$DYNMAT_DISPLACEMENT"
[[ -n "$DYNMAT_FILE" ]]         && stage1_export+=",DYNMAT_FILE=$DYNMAT_FILE"
[[ -n "$DYNMAT_BINARY" ]]       && stage1_export+=",DYNMAT_BINARY=$DYNMAT_BINARY"
# One value again: LAMMPS writes this file, vdos_dynmat.py reads it.
[[ -n "$DYNMAT_REF_TRAJ" ]]     && stage1_export+=",DYNMAT_REF_FILE=$DYNMAT_REF_TRAJ"

enabled_scripts=""
[[ "$RUN_DSF" == "1" ]] && enabled_scripts+="dsf.py "
[[ "$RUN_RDF" == "1" ]] && enabled_scripts+="rdf_freud.py "
[[ "$RUN_BAD" == "1" ]] && enabled_scripts+="bad_freud.py "
[[ "$RUN_VDOS" == "1" ]] && enabled_scripts+="vdos.py "
[[ "$RUN_MSD" == "1" ]] && enabled_scripts+="msd.py "
[[ "$RUN_VDOS_DYNMAT" == "1" ]] && enabled_scripts+="vdos_dynmat.py "
enabled_scripts="${enabled_scripts% }"

############################
# Execution: one branch runs both stages
############################

if [[ "$INTERACTIVE" == "1" ]]; then
  if [[ "$SKIP_TRAJECTORY" == "1" ]]; then
    echo "--skip-trajectory: not running LAMMPS, using existing trajectory in $STAGE1_DIR/run."
  else
    echo "Running LAMMPS interactively in $STAGE1_DIR/run ..."
    (
      export SLURM_SUBMIT_DIR="$STAGE1_DIR/run"
      export SLURM_NTASKS="${NTASKS:-64}"
      IFS=',' read -ra _kv_pairs <<< "${stage1_export#ALL,}"
      for pair in "${_kv_pairs[@]}"; do export "$pair"; done
      cd "$STAGE1_DIR/run" && bash lammps_submit.slurm
    )
    echo "LAMMPS run finished."
  fi

  echo "Running distribution analysis interactively in $STAGE2_DIR ..."
  (
    export SLURM_SUBMIT_DIR="$STAGE2_DIR"
    export SLURM_CPUS_PER_TASK="${ANALYSIS_CPUS_PER_TASK:-64}"
    IFS=',' read -ra _kv_pairs <<< "${export_vars#ALL,}"
    for pair in "${_kv_pairs[@]}"; do export "$pair"; done
    cd "$STAGE2_DIR" && bash distribution_submit.slurm
  )
  echo "Distribution analysis finished."

  cat <<SUMMARY

Pipeline completed interactively:
  LAMMPS run (setup + trajectory) : $STAGE1_DIR$([[ "$SKIP_TRAJECTORY" == "1" ]] && echo "  (skipped — existing trajectory)")
  distribution analysis           : $STAGE2_DIR
  analysis scripts enabled        : ${enabled_scripts:-none}

$(if [[ -n "$enabled_scripts" ]]; then
  echo "Reminder: any physics config (ELEMENTS, R_CUTOFF, DT, WINDOW_SIZE, CORR_LENGTH, ...) not"
  echo "passed via --rdf-*/--bad-*/--dsf-*/--vdos-*/--msd-*/--vdos-dynmat-* flags is using each script's own default —"
  echo "check $STAGE2_DIR/{${enabled_scripts// /,}} if that's not what you want."
fi)
SUMMARY
else
  JOBID1=""
  if [[ "$SKIP_TRAJECTORY" == "1" ]]; then
    echo "--skip-trajectory: not submitting a LAMMPS job, using existing trajectory in $STAGE1_DIR/run."
  else
    sbatch_args=()
    [[ -n "$NODES" ]]      && sbatch_args+=(--nodes="$NODES")
    [[ -n "$NTASKS" ]]     && sbatch_args+=(--ntasks="$NTASKS")
    [[ -n "$TIME" ]]       && sbatch_args+=(--time="$TIME")
    [[ -n "$JOB_NAME" ]]   && sbatch_args+=(--job-name="$JOB_NAME")
    [[ -n "$CONSTRAINT" ]] && sbatch_args+=(--constraint="$CONSTRAINT")
    [[ -n "$NODELIST" ]]   && sbatch_args+=(--nodelist="$NODELIST")

    echo "Submitting LAMMPS run from $STAGE1_DIR/run ..."
    JOBID1="$(cd "$STAGE1_DIR/run" && sbatch --parsable \
      "${sbatch_args[@]+"${sbatch_args[@]}"}" \
      --export="$stage1_export" \
      lammps_submit.slurm)"
    echo "  LAMMPS job id: $JOBID1"
  fi

  analysis_sbatch_args=()
  [[ -n "$ANALYSIS_NODES" ]]      && analysis_sbatch_args+=(--nodes="$ANALYSIS_NODES")
  [[ -n "$ANALYSIS_NTASKS" ]]     && analysis_sbatch_args+=(--ntasks="$ANALYSIS_NTASKS")
  [[ -n "$ANALYSIS_TIME" ]]       && analysis_sbatch_args+=(--time="$ANALYSIS_TIME")
  [[ -n "$ANALYSIS_JOB_NAME" ]]   && analysis_sbatch_args+=(--job-name="$ANALYSIS_JOB_NAME")
  [[ -n "$ANALYSIS_CONSTRAINT" ]] && analysis_sbatch_args+=(--constraint="$ANALYSIS_CONSTRAINT")
  [[ -n "$ANALYSIS_NODELIST" ]]   && analysis_sbatch_args+=(--nodelist="$ANALYSIS_NODELIST")
  [[ -n "$ANALYSIS_CPUS_PER_TASK" ]] && analysis_sbatch_args+=(--cpus-per-task="$ANALYSIS_CPUS_PER_TASK")
  # distribution_submit.slurm writes to the fixed filename STREAM_OUTPUT, which
  # sbatch truncates by default. Now that a second run can land in a directory an
  # earlier one already used, append instead — otherwise adding a calculation
  # keeps every earlier CSV but silently destroys the log explaining them.
  analysis_sbatch_args+=(--open-mode=append)

  dependency_args=()
  [[ -n "$JOBID1" ]] && dependency_args+=(--dependency=afterok:"$JOBID1")

  echo "Submitting distribution analysis from $STAGE2_DIR${JOBID1:+, dependent on job $JOBID1} ..."
  JOBID2="$(cd "$STAGE2_DIR" && sbatch --parsable \
    "${dependency_args[@]+"${dependency_args[@]}"}" \
    "${analysis_sbatch_args[@]+"${analysis_sbatch_args[@]}"}" \
    --export="$export_vars" \
    distribution_submit.slurm)"
  echo "  distribution-analysis job id: $JOBID2"

  cat <<SUMMARY

Pipeline submitted:
  LAMMPS run (setup + trajectory) : $STAGE1_DIR  ${JOBID1:+(job $JOBID1)}${JOBID1:-(skipped — existing trajectory)}
  distribution analysis           : $STAGE2_DIR  (job $JOBID2)
  analysis scripts enabled        : ${enabled_scripts:-none}

$(if [[ -n "$enabled_scripts" ]]; then
  echo "Reminder: any physics config (ELEMENTS, R_CUTOFF, DT, WINDOW_SIZE, CORR_LENGTH, ...) not"
  echo "passed via --rdf-*/--bad-*/--dsf-*/--vdos-*/--msd-*/--vdos-dynmat-* flags is using each script's own default —"
  echo "check $STAGE2_DIR/{${enabled_scripts// /,}} if that's not what you want."
fi)
SUMMARY
fi
