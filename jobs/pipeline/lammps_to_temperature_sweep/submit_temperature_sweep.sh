#!/bin/bash
# submit_temperature_sweep.sh — one thermalization cascade, then per-temperature
# diffusion and/or dielectric calculations, then aggregation of each route into
# one CSV.
#
# Stage 1 is a SINGLE LAMMPS job: deform to the target density, ramp to the
# highest requested temperature, then descend, writing a thermalized structure
# at every stop and an NVE trajectory at the stops where diffusion was asked
# for. Everything downstream fans out from it as independent jobs:
#
#   cascade ──afterok──> T<C> dielectric production ──afterok──> T<C> dipole calc ─┐
#           │                (one independent job per dielectric temperature)      ├─> aggregate
#           └──afterok──> T<C> msd ────────────────────────────────────────────────┘
#                            (one independent job per diffusion temperature)
#
# Run from inside a sweep directory; its basename is the sweep id.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: submit_temperature_sweep.sh [options]

Run from inside a sweep directory, whose name is the sweep id. Copy this
script and its conf into an empty directory and run it there:

    mkdir -p .../sweeps/water-OH-0001 && cd .../sweeps/water-OH-0001
    cp <util>/jobs/pipeline/lammps_to_temperature_sweep/submit_temperature_sweep.sh .
    cp <util>/jobs/pipeline/lammps_to_temperature_sweep/submit_temperature_sweep.conf.example \
       submit_temperature_sweep.conf
    $EDITOR submit_temperature_sweep.conf
    ./submit_temperature_sweep.sh --diffusion-temperatures "5;25;45" --dry-run

The sweep then carries the exact script and settings it was run with, next
to its own results. Running the script in place from the util checkout works
identically — nothing about the two cases differs.

Only the conf is read from this script's own directory. Everything else the
pipeline needs — generate_cascade_input.py, OH-cascade-preamble.input,
dielectric-production.input, the aggregators, and the slurm templates — is
read from the util checkout, found via REPO_ROOT:

  1. REPO_ROOT set to a literal path in the conf. The shipped conf sets
     /home1/lkyamamo/util, so a copied script finds the pipeline files with
     no further configuration. This also governs the conf's own $REPO_ROOT
     paths, so the two cannot disagree.
  2. Otherwise REPO_ROOT exported in the environment.
  3. Otherwise, if this script is still sitting in the checkout, that checkout.
  4. Otherwise $HOME/util.

     Only 1 applies with the shipped conf. Blank the conf line (or write
     REPO_ROOT="${REPO_ROOT:-}") to fall through to 2-4, which is what you
     want on a machine where the checkout is somewhere else — a laptop, or a
     worktree you are testing from.

If REPO_ROOT is wrong the script says so and exits before creating anything.

The conf sets defaults for every parameter below — the script refuses to run
without one. Any flag below overrides its config value for that one
invocation.

TEMPERATURES ARE IN CELSIUS. They are converted to Kelvin for LAMMPS and for
2.dipole_std.py's prefactor; the T<C>/ directories and the CSVs' temperature_C
column keep the Celsius value you asked for, and the CSVs also carry
temperature_K so nobody has to redo the conversion.

Required:
  --starting-structure FILE    Seed structure (e.g. ICE_CUBIC.data), symlinked
                                as start.data. Read for its atom counts by the
                                density solve; its own Masses section is NOT
                                used (see --density).
  --potential-file FILE        Potential file, symlinked into input_files/
  --analysis-parent-dir DIR    Where the per-temperature analysis dirs are made

  At least one of:
  --diffusion-temperatures LIST   SEMICOLON-separated CELSIUS list, e.g. "5;25;45"
  --dielectric-temperatures LIST  SEMICOLON-separated CELSIUS list, e.g. "25;45;65"

      Semicolons, not commas: these values reach the SLURM stages through
      sbatch --export, which is comma-delimited and truncates silently.

      Giving only one list runs only that route. The cascade still thermalizes
      and writes a structure at every temperature in the UNION of the two, so
      a diffusion-only sweep and a dielectric-only sweep produce structures the
      same way. A temperature in both lists is one cascade stop feeding both
      routes.

      Both are empty in the shipped .conf.example, so which routes run is an
      explicit per-invocation choice. If you do set one in your .conf, note
      that passing the OTHER flag does not turn it off — each flag overrides
      only its own variable, so a conf with both set would submit both routes
      no matter which single flag you pass. Pass an explicit empty string to
      suppress a route the conf sets:

          --dielectric-temperatures ""    run only the diffusion route

Cascade (stage 1) config:
  --replicate "NX NY NZ"       Supercell replication (default: "6 6 6"). This is
                                the ONLY system-size knob: the preamble
                                (OH-cascade-preamble.input) is a fixed constant —
                                always O and H, always the same masses and
                                pair_coeff — so there is no --preamble or
                                --elements flag.
  --density FLOAT              Target mass density in g/cc for the deform
                                (default: 1.0). The cube edge is solved from the
                                replicated atom count and the PREAMBLE's masses,
                                not start.data's: ICE_CUBIC.data declares H as
                                1.0 while LAMMPS integrates the preamble's
                                1.00784, which is the difference between a
                                37.2514 A and a 37.2406 A cube at 1 g/cc.
  --melt-temperature FLOAT     Optional CELSIUS melt above the highest requested
                                temperature, held and then descended from but
                                never saved (default: none — ramp straight to
                                the highest requested temperature). This is
                                OH-therm.input's 400 K step; 127 C is that same
                                temperature.
  --timestep FLOAT             LAMMPS timestep in PICOSECONDS, metal units
                                (default: 0.00025, i.e. 0.25 fs). Note the unit:
                                the MSD parameters below are in femtoseconds.
  --seed INT                   velocity create seed (default: 156467)
  --deform-length INT          Steps for the 10 K deform-to-density block
  --ramp-length INT            Steps for each temperature ramp
  --hold-length INT            Steps held at each temperature before saving it
  --nve-length INT             Steps of NVE production at each diffusion
                                temperature
  --dynamics-dump-every INT    Steps between frames in dynamics_T<C>.lammpstrj

Diffusion (MSD) config — same meaning as msd.py's own environment variables:
  --msd-corr-length FLOAT      FEMTOSECONDS; max time lag. REQUIRED when a
                                diffusion list is given.
  --msd-corr-interval FLOAT    FEMTOSECONDS; spacing between reference frames.
                                REQUIRED when a diffusion list is given.
  --msd-fit-fraction FLOAT     Tail fraction of the window used for the D fit.
                                REQUIRED when a diffusion list is given.
  --dynamics-dt FLOAT          FEMTOSECONDS between dumped frames. Normally
                                omitted: it is DERIVED as
                                --timestep x --dynamics-dump-every x 1000, which
                                is exact because this pipeline sets both of
                                those itself. Set it only to override.

  msd.py hard-fails on the three REQUIRED values rather than guessing, so this
  script checks them before submitting anything — a diffusion sweep should not
  discover a missing fit fraction after the cascade has already run.

Dielectric production (stage 2a) config — LAMMPS -var overrides:
  --dielectric-n-chunks INT       Number of chunked dielectric.N.custom files
                                  (LAMMPS ${N_TIMES}) (default: 128)
  --dielectric-chunk-length INT   NVT steps per chunk (LAMMPS ${NVT_LENGTH})
                                  (default: 468750)
  --dielectric-dump-every INT     Steps between dumped frames within a chunk
                                  (LAMMPS ${DUMP_EVERY}) (default: 10)

Dielectric calc (stage 2b) config — same meaning as 1.calc_mpi.py /
2.dipole_std.py's own flags:
  --dielectric-cutoff FLOAT        O-H bond cutoff (Angstrom) (default: 1.2)
  --dielectric-type-o INT          LAMMPS atom type for oxygen (default: 1)
  --dielectric-type-h INT          LAMMPS atom type for hydrogen (default: 2)
  --dielectric-charge-h FLOAT      Hydrogen partial charge (e) (default:
                                  0.406988). Oxygen's charge isn't a separate
                                  parameter — the formula assumes charge
                                  neutrality (O = -2x this value). Also passed
                                  into stage 2a as the LAMMPS -var CHARGE_H so
                                  an input script with a by-hand dipole
                                  cross-check uses the same value.
  --dielectric-averaging-method STR  windowed|cumulative|hybrid|binned
                                      (default: cumulative)
  --dielectric-ntasks INT          Ranks for the stage-2b MPI calc (default:
                                  empty — one rank per dielectric.N.custom
                                  file). Must not exceed the chunk-file count.
  --dielectric-mem STR              Stage 2b's sbatch --mem (default: 0)
  --la/--lb/--lc FLOAT             Box lengths (Angstrom) feeding
                                  2.dipole_std.py's prefactor. Default to the
                                  cube edge the density solve produced — the box
                                  is no longer independent of this pipeline, so
                                  these are overrides, not required inputs.

SLURM overrides (omit any to leave the template's own value in effect):
  --nodes/--ntasks/--time/--job-name/--constraint/--nodelist
      Stage 1, the cascade (jobs/slurm/lammps_cascade_submit.slurm).
  --dielectric-nodes/--dielectric-time/--dielectric-job-name/
  --dielectric-constraint/--dielectric-nodelist
      Stage 2a, dielectric production (jobs/slurm/lammps_dielectric_submit.slurm).
  --analysis-nodes/--analysis-time/--analysis-job-name/
  --analysis-constraint/--analysis-nodelist
      Stages 2b and 3, the per-temperature analysis jobs.

Other:
  --analysis-template-dir DIR      Dielectric-calc scripts to copy into stage 2b
                                    (default: <repo>/simulation/lammps/20260617_dielectric_multi_traj)
  --msd-template-dir DIR           MSD scripts to copy into stage 3
                                    (default: <repo>/analysis/distributions/20260608_GrNrBaSqw)
  --venv PATH                      Virtualenv to activate before the density
                                    solve and the aggregators, which need numpy
                                    and ase (default: from the conf). Empty
                                    uses whatever python is already on PATH.
  --lmp-bin PATH                   LAMMPS executable for both the cascade and
                                    the dielectric production stage (default:
                                    empty — each SLURM template's own default,
                                    /home1/lkyamamo/executables/lammps/lmp_mpi_shock_2019).
                                    It must be built with pair_style usc; the
                                    cascade otherwise needs no optional package.
  --potential-link-name NAME       Symlink name in input_files/ (default:
                                    OH.usc). The preamble's pair_coeff line is
                                    fixed at OH.usc, so anything else is
                                    rejected up front rather than failing inside
                                    LAMMPS after the job has queued.

  --skip-cascade                   Skip stage 1 and assume cascade/ already has
                                  the thermalized_T<C>.data and
                                  dynamics_T<C>.lammpstrj files this sweep needs.
                                  Each one is verified to exist before anything
                                  is submitted, and the downstream jobs are
                                  submitted without a --dependency (nothing to
                                  wait on).

  --dry-run                        Print every sbatch command and the dependency
                                  graph, then exit WITHOUT submitting anything
                                  and without creating any directories. Use this
                                  to check the fan-out before spending queue
                                  time.

  --stagger-seconds N              Mitigations against SLURM's "user env
                                  retrieval failed" transient launch failure
                                  (default: 5, no-ops if 0): sleeps N seconds
                                  between per-temperature submissions, and adds
                                  sbatch --begin=now+<random 0..N>seconds to
                                  every submission. The random jitter also
                                  decorrelates this invocation from OTHER
                                  concurrently-running invocations, which the
                                  sleep cannot do.

  --force REASON                  Overwrite existing cascade/, T<C>/ or analysis
                                  directories instead of refusing to run. REASON
                                  is required and is logged with a timestamp to
                                  overwrite.log next to this script.

  --interactive                   Run every stage in the foreground via plain
                                  bash instead of sbatch. Assumes you already
                                  have an allocation (salloc); this does not
                                  request one. Temperatures run sequentially.

  -h, --help                    Show this help

NOTE ON THE SINGLE CASCADE JOB: because all of the thermalization is one job,
no downstream stage can start until EVERY temperature has finished, even though
the highest-temperature structure is written early on. That is the cost of
descending through one continuous trajectory instead of annealing each
temperature independently.

All terminal output from this script (not the SLURM jobs themselves) is also
appended to submit_temperature_sweep.log in the sweep directory (cwd).
EOF
}

INTERACTIVE="0"
SKIP_CASCADE="0"
DRY_RUN="0"
STAGGER_SECONDS="5"
FORCE="0"
FORCE_REASON=""

# Where THIS copy of the script lives. The .conf is always read from here, so
# the intended workflow works with no extra flags:
#
#   cp submit_temperature_sweep.sh submit_temperature_sweep.conf <empty dir>/
#   cd <empty dir> && ./submit_temperature_sweep.sh ...
#
# Running it in place from the util checkout works the same way — in that case
# "here" just happens to be the checkout.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# The preamble's pair_coeff line names this file; --potential-link-name must
# match it or LAMMPS fails on a missing potential after the job has queued.
REQUIRED_POTENTIAL_LINK_NAME="OH.usc"

log_overwrite() {
  local msg
  msg="[$(date "+%Y-%m-%dT%H:%M:%S%z")] $1"
  echo "WARNING: $msg" >&2
  echo "$msg" >> "$LOG_FILE"
}

# Sets BEGIN_ARGS to --begin=now+<jitter>seconds, jitter uniform in
# [0, STAGGER_SECONDS]. Applied to every sbatch call, not just the
# per-temperature loop: separate concurrent invocations can't coordinate their
# submission timing with each other except through randomness.
random_begin_args() {
  BEGIN_ARGS=()
  if [[ "$STAGGER_SECONDS" -gt 0 ]]; then
    local jitter=$((RANDOM % (STAGGER_SECONDS + 1)))
    BEGIN_ARGS=(--begin="now+${jitter}seconds")
  fi
}

# Celsius -> Kelvin. bash has no float arithmetic, hence awk.
to_kelvin() { awk -v c="$1" 'BEGIN{printf "%.2f", c + 273.15}'; }

# Path-safe Celsius label. "%g" collapses 25 and 25.0 to one spelling, so the
# same temperature written two ways can't produce two directories — matching
# generate_cascade_input.py's celsius_label().
celsius_label() { awk -v c="$1" 'BEGIN{printf "%g", c}'; }

# Semicolon-separated list -> global TEMP_ARRAY of normalized Celsius labels.
parse_temp_list() {
  TEMP_ARRAY=()
  local raw="$1" flag="$2" item
  [[ -z "$raw" ]] && return 0
  local IFS=';'
  for item in $raw; do
    item="$(echo "$item" | tr -d '[:space:]')"
    [[ -z "$item" ]] && continue
    if ! [[ "$item" =~ ^-?[0-9]+\.?[0-9]*$ ]]; then
      echo "Error: $flag contains a non-numeric entry: '$item'" >&2
      if [[ "$item" == *,* ]]; then
        echo "       (commas are not separators here — use semicolons)" >&2
      fi
      exit 1
    fi
    TEMP_ARRAY+=("$(celsius_label "$item")")
  done
}

# All defaults live in submit_temperature_sweep.conf next to this script, not
# in the script itself, so personal paths never need to be committed. Every key
# is pre-declared below so an older conf still works under set -u.
CONFIG_FILE="$SCRIPT_DIR/submit_temperature_sweep.conf"
if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "Error: required config file not found: $CONFIG_FILE" >&2
  echo "This script reads its config from its own directory, so copy BOTH files together:" >&2
  echo "  cp <util>/jobs/pipeline/lammps_to_temperature_sweep/submit_temperature_sweep.sh ." >&2
  echo "  cp <util>/jobs/pipeline/lammps_to_temperature_sweep/submit_temperature_sweep.conf.example submit_temperature_sweep.conf" >&2
  echo "then edit submit_temperature_sweep.conf for this sweep." >&2
  exit 1
fi

VENV_PATH=""
STARTING_STRUCTURE=""; POTENTIAL_FILE=""; POTENTIAL_LINK_NAME=""; LMP_BIN=""
ANALYSIS_PARENT_DIR=""; ANALYSIS_TEMPLATE_DIR=""; MSD_TEMPLATE_DIR=""
DIFFUSION_TEMPERATURES=""; DIELECTRIC_TEMPERATURES=""
REPLICATE=""; DENSITY=""; MELT_TEMPERATURE=""; TIMESTEP=""; SEED=""
DEFORM_LENGTH=""; RAMP_LENGTH=""; HOLD_LENGTH=""; NVE_LENGTH=""
DYNAMICS_DUMP_EVERY=""; DYNAMICS_DT=""
MSD_CORR_LENGTH=""; MSD_CORR_INTERVAL=""; MSD_FIT_FRACTION=""
DIELECTRIC_N_CHUNKS=""; DIELECTRIC_CHUNK_LENGTH=""; DIELECTRIC_DUMP_EVERY=""
DIELECTRIC_CUTOFF=""; DIELECTRIC_TYPE_O=""; DIELECTRIC_TYPE_H=""
DIELECTRIC_CHARGE_H=""; DIELECTRIC_AVERAGING_METHOD=""; DIELECTRIC_NTASKS=""
DIELECTRIC_MEM=""; LA=""; LB=""; LC=""
NODES=""; NTASKS=""; TIME=""; JOB_NAME=""; CONSTRAINT=""; NODELIST=""
DIELECTRIC_NODES=""; DIELECTRIC_TIME=""; DIELECTRIC_JOB_NAME=""
DIELECTRIC_CONSTRAINT=""; DIELECTRIC_NODELIST=""
ANALYSIS_NODES=""; ANALYSIS_TIME=""; ANALYSIS_JOB_NAME=""
ANALYSIS_CONSTRAINT=""; ANALYSIS_NODELIST=""

############################
# Where the rest of the pipeline lives
#
# This script is meant to be COPIED, with its conf, into a sweep directory —
# so its own location says nothing about where the other pipeline files are.
# Everything except the conf is resolved from the util checkout:
#
#   1. REPO_ROOT exported in the environment, if set.
#   2. REPO_ROOT set in the conf — the reliable way to point a copied script
#      at a specific checkout.
#   3. Otherwise, if this script IS still sitting in the checkout (its
#      siblings are beside it), that checkout — so running it in place needs
#      no configuration at all.
#   4. Otherwise $HOME/util.
#
# A fallback is seeded BEFORE the conf is sourced, and the conf gets the last
# word. The shipped conf writes
#
#     REPO_ROOT="${REPO_ROOT:-}"
#
# which keeps the seeded value, so paths further down it can be written as
# "$REPO_ROOT/starting-structures/..." and land in the same checkout the
# pipeline files come from. Replacing that with a literal path overrides the
# seed, and — because it is one assignment read top to bottom — the literal is
# then what those paths use too. Either way REPO_ROOT has exactly one value,
# and the conf's paths can never disagree with the pipeline's.
############################

if [[ -z "${REPO_ROOT:-}" ]]; then
  if [[ -f "$SCRIPT_DIR/generate_cascade_input.py" ]]; then
    REPO_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
  else
    REPO_ROOT="$HOME/util"
  fi
fi

# shellcheck source=/dev/null
source "$CONFIG_FILE"

PIPELINE_DIR="$REPO_ROOT/jobs/pipeline/lammps_to_temperature_sweep"

CASCADE_TEMPLATE="$REPO_ROOT/jobs/slurm/lammps_cascade_submit.slurm"
DIELECTRIC_TEMPLATE="$REPO_ROOT/jobs/slurm/lammps_dielectric_submit.slurm"
LOG_FILE="$PIPELINE_DIR/overwrite.log"
GENERATE_SCRIPT="$PIPELINE_DIR/generate_cascade_input.py"
PREAMBLE_FILE="$PIPELINE_DIR/OH-cascade-preamble.input"
DIELECTRIC_INPUT="$PIPELINE_DIR/dielectric-production.input"
AGGREGATE_TEMPLATE="$PIPELINE_DIR/aggregate_submit.slurm"
AGGREGATE_DIELECTRIC_SCRIPT="$PIPELINE_DIR/aggregate_dielectric_vs_temperature.py"
AGGREGATE_DIFFUSION_SCRIPT="$PIPELINE_DIR/aggregate_diffusion_vs_temperature.py"

# Fail here, with the path that was wrong, rather than midway through creating
# directories — a copied script pointed at a bad REPO_ROOT should say so before
# it has done anything.
if [[ ! -f "$GENERATE_SCRIPT" ]]; then
  echo "Error: cannot find the pipeline files under REPO_ROOT=$REPO_ROOT" >&2
  echo "  expected: $GENERATE_SCRIPT" >&2
  echo "Set REPO_ROOT in $CONFIG_FILE to your util checkout." >&2
  exit 1
fi

while [[ $# -gt 0 ]]; do
  case "$1" in
    --starting-structure) STARTING_STRUCTURE="$2"; shift 2 ;;
    --potential-file) POTENTIAL_FILE="$2"; shift 2 ;;
    --potential-link-name) POTENTIAL_LINK_NAME="$2"; shift 2 ;;
    --lmp-bin) LMP_BIN="$2"; shift 2 ;;
    --venv) VENV_PATH="$2"; shift 2 ;;
    --analysis-parent-dir) ANALYSIS_PARENT_DIR="$2"; shift 2 ;;
    --analysis-template-dir) ANALYSIS_TEMPLATE_DIR="$2"; shift 2 ;;
    --msd-template-dir) MSD_TEMPLATE_DIR="$2"; shift 2 ;;
    --diffusion-temperatures) DIFFUSION_TEMPERATURES="$2"; shift 2 ;;
    --dielectric-temperatures) DIELECTRIC_TEMPERATURES="$2"; shift 2 ;;
    --replicate) REPLICATE="$2"; shift 2 ;;
    --density) DENSITY="$2"; shift 2 ;;
    --melt-temperature) MELT_TEMPERATURE="$2"; shift 2 ;;
    --timestep) TIMESTEP="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --deform-length) DEFORM_LENGTH="$2"; shift 2 ;;
    --ramp-length) RAMP_LENGTH="$2"; shift 2 ;;
    --hold-length) HOLD_LENGTH="$2"; shift 2 ;;
    --nve-length) NVE_LENGTH="$2"; shift 2 ;;
    --dynamics-dump-every) DYNAMICS_DUMP_EVERY="$2"; shift 2 ;;
    --dynamics-dt) DYNAMICS_DT="$2"; shift 2 ;;
    --msd-corr-length) MSD_CORR_LENGTH="$2"; shift 2 ;;
    --msd-corr-interval) MSD_CORR_INTERVAL="$2"; shift 2 ;;
    --msd-fit-fraction) MSD_FIT_FRACTION="$2"; shift 2 ;;
    --dielectric-n-chunks) DIELECTRIC_N_CHUNKS="$2"; shift 2 ;;
    --dielectric-chunk-length) DIELECTRIC_CHUNK_LENGTH="$2"; shift 2 ;;
    --dielectric-dump-every) DIELECTRIC_DUMP_EVERY="$2"; shift 2 ;;
    --dielectric-cutoff) DIELECTRIC_CUTOFF="$2"; shift 2 ;;
    --dielectric-type-o) DIELECTRIC_TYPE_O="$2"; shift 2 ;;
    --dielectric-type-h) DIELECTRIC_TYPE_H="$2"; shift 2 ;;
    --dielectric-charge-h) DIELECTRIC_CHARGE_H="$2"; shift 2 ;;
    --dielectric-averaging-method) DIELECTRIC_AVERAGING_METHOD="$2"; shift 2 ;;
    --dielectric-ntasks) DIELECTRIC_NTASKS="$2"; shift 2 ;;
    --dielectric-mem) DIELECTRIC_MEM="$2"; shift 2 ;;
    --la) LA="$2"; shift 2 ;;
    --lb) LB="$2"; shift 2 ;;
    --lc) LC="$2"; shift 2 ;;
    --nodes) NODES="$2"; shift 2 ;;
    --ntasks) NTASKS="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --job-name) JOB_NAME="$2"; shift 2 ;;
    --constraint) CONSTRAINT="$2"; shift 2 ;;
    --nodelist) NODELIST="$2"; shift 2 ;;
    --dielectric-nodes) DIELECTRIC_NODES="$2"; shift 2 ;;
    --dielectric-time) DIELECTRIC_TIME="$2"; shift 2 ;;
    --dielectric-job-name) DIELECTRIC_JOB_NAME="$2"; shift 2 ;;
    --dielectric-constraint) DIELECTRIC_CONSTRAINT="$2"; shift 2 ;;
    --dielectric-nodelist) DIELECTRIC_NODELIST="$2"; shift 2 ;;
    --analysis-nodes) ANALYSIS_NODES="$2"; shift 2 ;;
    --analysis-time) ANALYSIS_TIME="$2"; shift 2 ;;
    --analysis-job-name) ANALYSIS_JOB_NAME="$2"; shift 2 ;;
    --analysis-constraint) ANALYSIS_CONSTRAINT="$2"; shift 2 ;;
    --analysis-nodelist) ANALYSIS_NODELIST="$2"; shift 2 ;;
    --force) FORCE="1"; FORCE_REASON="$2"; shift 2 ;;
    --interactive) INTERACTIVE="1"; shift 1 ;;
    --skip-cascade) SKIP_CASCADE="1"; shift 1 ;;
    --dry-run) DRY_RUN="1"; shift 1 ;;
    --stagger-seconds) STAGGER_SECONDS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

# The density solve (numpy via box_size.py, ase via box_size_from_data.py) runs
# in this shell, before anything is submitted, so the venv has to be active
# here — distribution_submit.slurm and dielectric_submit.slurm activate it
# themselves inside their jobs. After the argument loop so --venv applies.
if [[ -n "$VENV_PATH" ]]; then
  # shellcheck source=/dev/null
  source "$VENV_PATH/bin/activate"
fi

SWEEP_DIR="$(pwd)"
SWEEP_ID="$(basename "$SWEEP_DIR")"
if [[ "$DRY_RUN" != "1" ]]; then
  PIPELINE_LOG="$SWEEP_DIR/submit_temperature_sweep.log"
  exec > >(tee -a "$PIPELINE_LOG") 2>&1
  echo "=== submit_temperature_sweep.sh started $(date "+%Y-%m-%dT%H:%M:%S%z") — sweep id: $SWEEP_ID ==="
  echo "Logging this run's output to: $PIPELINE_LOG"
else
  echo "=== DRY RUN — nothing will be submitted or created (sweep id: $SWEEP_ID) ==="
fi

############################
# Validation — everything that can fail should fail here, before the cascade
############################

missing=0
require() {
  if [[ -z "$1" ]]; then
    echo "Missing required argument: $2" >&2
    missing=1
  fi
}
if [[ "$SKIP_CASCADE" != "1" ]]; then
  require "$STARTING_STRUCTURE" --starting-structure
  require "$POTENTIAL_FILE" --potential-file
fi
require "$ANALYSIS_PARENT_DIR" --analysis-parent-dir
if [[ "$missing" -ne 0 ]]; then
  usage
  exit 1
fi

parse_temp_list "$DIFFUSION_TEMPERATURES" --diffusion-temperatures
DIFFUSION_TEMPS=("${TEMP_ARRAY[@]+"${TEMP_ARRAY[@]}"}")
parse_temp_list "$DIELECTRIC_TEMPERATURES" --dielectric-temperatures
DIELECTRIC_TEMPS=("${TEMP_ARRAY[@]+"${TEMP_ARRAY[@]}"}")

N_DIFFUSION="${#DIFFUSION_TEMPS[@]}"
N_DIELECTRIC="${#DIELECTRIC_TEMPS[@]}"

if [[ "$N_DIFFUSION" -eq 0 && "$N_DIELECTRIC" -eq 0 ]]; then
  echo "Error: at least one of --diffusion-temperatures / --dielectric-temperatures must be non-empty." >&2
  echo "Give one to run only that route; give both to run both." >&2
  exit 1
fi

# msd.py refuses to guess the parameters that determine its result, so a
# diffusion sweep must not get as far as the cascade without them.
if [[ "$N_DIFFUSION" -gt 0 ]]; then
  msd_missing=()
  [[ -n "$MSD_CORR_LENGTH" ]]   || msd_missing+=("  --msd-corr-length / MSD_CORR_LENGTH       fs; max time lag")
  [[ -n "$MSD_CORR_INTERVAL" ]] || msd_missing+=("  --msd-corr-interval / MSD_CORR_INTERVAL   fs; spacing between reference frames")
  [[ -n "$MSD_FIT_FRACTION" ]]  || msd_missing+=("  --msd-fit-fraction / MSD_FIT_FRACTION     tail fraction used for the D fit")
  if [[ "${#msd_missing[@]}" -gt 0 ]]; then
    echo "Error: --diffusion-temperatures was given, but these msd.py parameters are unset:" >&2
    printf '%s\n' "${msd_missing[@]}" >&2
    echo "These determine the numbers msd.py produces, so nothing guesses them." >&2
    exit 1
  fi
fi

if [[ "$FORCE" == "1" && -z "$FORCE_REASON" ]]; then
  echo 'Error: --force requires a non-empty reason, e.g. --force "re-running after fixing the potential"' >&2
  exit 1
fi

POTENTIAL_LINK_NAME="${POTENTIAL_LINK_NAME:-$REQUIRED_POTENTIAL_LINK_NAME}"
if [[ "$POTENTIAL_LINK_NAME" != "$REQUIRED_POTENTIAL_LINK_NAME" ]]; then
  echo "Error: --potential-link-name is '$POTENTIAL_LINK_NAME', but the fixed preamble" >&2
  echo "($PREAMBLE_FILE) has 'pair_coeff * * $REQUIRED_POTENTIAL_LINK_NAME O H'." >&2
  echo "LAMMPS would not find the potential. Use --potential-link-name $REQUIRED_POTENTIAL_LINK_NAME." >&2
  exit 1
fi

if [[ "$INTERACTIVE" == "1" && "$DRY_RUN" != "1" && -z "${SLURM_JOB_ID:-}" ]]; then
  echo "Error: --interactive requires an active Slurm allocation (run 'salloc ...' first)." >&2
  echo "No \$SLURM_JOB_ID is set in this shell, so srun would fail to allocate resources —" >&2
  echo "and the submit templates' unconditional 'exit 0' would mask that failure." >&2
  exit 1
fi

# No exported value may contain a comma: sbatch --export is comma-delimited and
# truncates silently, which would corrupt a path or a parameter without any error.
for _v in ANALYSIS_PARENT_DIR DIELECTRIC_AVERAGING_METHOD REPLICATE; do
  if [[ "${!_v}" == *,* ]]; then
    echo "Error: $_v contains a comma ('${!_v}'), which sbatch --export cannot carry." >&2
    exit 1
  fi
done

############################
# Derived values
############################

# DYNAMICS_DT is derived, not guessed: this pipeline sets both the timestep and
# the dump cadence itself, so the frame spacing is exact. --timestep is in
# picoseconds (LAMMPS metal units) and msd.py wants femtoseconds, hence x1000.
if [[ -z "$DYNAMICS_DT" && "$N_DIFFUSION" -gt 0 ]]; then
  DYNAMICS_DT="$(awk -v ts="$TIMESTEP" -v every="$DYNAMICS_DUMP_EVERY" \
    'BEGIN{printf "%g", ts * every * 1000}')"
  echo "Derived DYNAMICS_DT = $DYNAMICS_DT fs (--timestep $TIMESTEP ps x --dynamics-dump-every $DYNAMICS_DUMP_EVERY x 1000)"
fi

# The cube edge the cascade will deform to. Asking the generator for it (rather
# than recomputing here) keeps one implementation of the density solve.
BOX_LENGTH=""
if [[ "$N_DIELECTRIC" -gt 0 || "$SKIP_CASCADE" != "1" ]]; then
  if [[ -n "$STARTING_STRUCTURE" ]]; then
    BOX_LENGTH="$(python "$GENERATE_SCRIPT" \
      --repo-root "$REPO_ROOT" \
      --start-data "$STARTING_STRUCTURE" \
      --replicate "$REPLICATE" \
      --density "$DENSITY" \
      --out /dev/null \
      --print-box-length)"
    echo "Solved box edge: $BOX_LENGTH A ($DENSITY g/cc, replicate $REPLICATE)"
  fi
fi
# The box is no longer independent of the pipeline, so these default to the
# solved edge instead of being three required flags.
LA="${LA:-$BOX_LENGTH}"; LB="${LB:-$BOX_LENGTH}"; LC="${LC:-$BOX_LENGTH}"
if [[ "$N_DIELECTRIC" -gt 0 && ( -z "$LA" || -z "$LB" || -z "$LC" ) ]]; then
  echo "Error: the dielectric route needs box lengths, and none could be solved" >&2
  echo "(--skip-cascade without --starting-structure). Pass --la/--lb/--lc." >&2
  exit 1
fi

CASCADE_DIR="$SWEEP_DIR/cascade"
INPUT_DIR="$SWEEP_DIR/input_files"

############################
# Dry run: print the plan and stop
############################

if [[ "$DRY_RUN" == "1" ]]; then
  echo ""
  echo "Sweep directory     : $SWEEP_DIR"
  echo "Diffusion temps (C) : ${DIFFUSION_TEMPS[*]+"${DIFFUSION_TEMPS[*]}"}"
  echo "Dielectric temps (C): ${DIELECTRIC_TEMPS[*]+"${DIELECTRIC_TEMPS[*]}"}"
  echo "Box edge            : ${BOX_LENGTH:-<unsolved>} A"
  echo "DYNAMICS_DT         : ${DYNAMICS_DT:-<n/a>} fs"
  echo ""
  echo "Job graph:"
  if [[ "$SKIP_CASCADE" == "1" ]]; then
    echo "  (stage 1 skipped — reusing $CASCADE_DIR)"
    cascade_dep="none"
  else
    echo "  [J1] cascade                                  $CASCADE_DIR"
    echo "         sbatch $CASCADE_TEMPLATE"
    cascade_dep="J1"
  fi
  for T in "${DIELECTRIC_TEMPS[@]+"${DIELECTRIC_TEMPS[@]}"}"; do
    TK="$(to_kelvin "$T")"
    echo "  [P${T}] dielectric production T=${T}C (${TK} K)   after: $cascade_dep"
    echo "         start.data -> cascade/thermalized_T${T}C.data"
    echo "  [C${T}] dielectric calc       T=${T}C             after: P${T}"
    echo "         -> $ANALYSIS_PARENT_DIR/${SWEEP_ID}_T${T}C_dielectric_calc"
  done
  for T in "${DIFFUSION_TEMPS[@]+"${DIFFUSION_TEMPS[@]}"}"; do
    TK="$(to_kelvin "$T")"
    echo "  [M${T}] msd                   T=${T}C (${TK} K)   after: $cascade_dep"
    echo "         DYNAMICS_TRAJ = cascade/dynamics_T${T}C.lammpstrj"
    echo "         -> $ANALYSIS_PARENT_DIR/${SWEEP_ID}_T${T}C_msd"
  done
  [[ "$N_DIELECTRIC" -gt 0 ]] && echo "  [A1] aggregate dielectric  after: all C* -> $SWEEP_DIR/dielectric_vs_temperature.csv"
  [[ "$N_DIFFUSION"  -gt 0 ]] && echo "  [A2] aggregate diffusion   after: all M* -> $SWEEP_DIR/diffusion_vs_temperature.csv"
  echo ""
  echo "Dry run complete — nothing was submitted or created."
  exit 0
fi

############################
# Stage 1: the cascade
############################

mkdir -p "$ANALYSIS_PARENT_DIR"
ANALYSIS_PARENT_DIR="$(cd "$ANALYSIS_PARENT_DIR" && pwd)"

CASCADE_JOBID=""

if [[ "$SKIP_CASCADE" == "1" ]]; then
  # Verify every file the downstream stages will want, before submitting any of
  # them — a missing trajectory should be an error now, not a failed job later.
  cascade_missing=()
  for T in "${DIELECTRIC_TEMPS[@]+"${DIELECTRIC_TEMPS[@]}"}"; do
    [[ -f "$CASCADE_DIR/thermalized_T${T}C.data" ]] || cascade_missing+=("thermalized_T${T}C.data")
  done
  for T in "${DIFFUSION_TEMPS[@]+"${DIFFUSION_TEMPS[@]}"}"; do
    [[ -f "$CASCADE_DIR/dynamics_T${T}C.lammpstrj" ]] || cascade_missing+=("dynamics_T${T}C.lammpstrj")
  done
  if [[ "${#cascade_missing[@]}" -gt 0 ]]; then
    echo "Error: --skip-cascade was given, but $CASCADE_DIR is missing:" >&2
    printf '  %s\n' "${cascade_missing[@]}" >&2
    echo "Re-run without --skip-cascade, or fix the temperature lists." >&2
    exit 1
  fi
  echo "--skip-cascade: found every required file in $CASCADE_DIR."
else
  if [[ -e "$INPUT_DIR" || -e "$CASCADE_DIR" ]]; then
    if [[ "$FORCE" == "1" ]]; then
      log_overwrite "--force ($FORCE_REASON): removing existing $INPUT_DIR and/or $CASCADE_DIR (sweep id: $SWEEP_ID)"
      rm -rf "$INPUT_DIR" "$CASCADE_DIR"
    else
      echo "Error: $SWEEP_DIR already has input_files/ or cascade/ — refusing to overwrite (use --force to override)." >&2
      exit 1
    fi
  fi

  mkdir -p "$INPUT_DIR" "$CASCADE_DIR"

  echo ""
  echo "--- Generating the cascade input ---"
  gen_args=(
    --repo-root "$REPO_ROOT"
    --start-data "$STARTING_STRUCTURE"
    --replicate "$REPLICATE"
    --density "$DENSITY"
    --diffusion-temps-c "$DIFFUSION_TEMPERATURES"
    --dielectric-temps-c "$DIELECTRIC_TEMPERATURES"
    --timestep "$TIMESTEP"
    --seed "$SEED"
    --deform-length "$DEFORM_LENGTH"
    --ramp-length "$RAMP_LENGTH"
    --hold-length "$HOLD_LENGTH"
    --nve-length "$NVE_LENGTH"
    --dynamics-dump-every "$DYNAMICS_DUMP_EVERY"
    --out "$INPUT_DIR/in.input"
  )
  [[ -n "$MELT_TEMPERATURE" ]] && gen_args+=(--melt-temperature "$MELT_TEMPERATURE")
  python "$GENERATE_SCRIPT" "${gen_args[@]}"

  ln -s "$(realpath "$STARTING_STRUCTURE")" "$INPUT_DIR/start.data"
  ln -s "$(realpath "$POTENTIAL_FILE")" "$INPUT_DIR/$POTENTIAL_LINK_NAME"
  cp "$CASCADE_TEMPLATE" "$CASCADE_DIR/lammps_submit.slurm"

  if [[ "$INTERACTIVE" == "1" ]]; then
    echo "Running the cascade interactively in $CASCADE_DIR ..."
    (
      export SLURM_SUBMIT_DIR="$CASCADE_DIR"
      export SLURM_NTASKS="${NTASKS:-64}"
      [[ -n "$LMP_BIN" ]] && export LMP_BIN
      cd "$CASCADE_DIR" && bash lammps_submit.slurm
    )
    echo "Cascade finished."
  else
    sbatch_args=()
    [[ -n "$NODES" ]]      && sbatch_args+=(--nodes="$NODES")
    [[ -n "$NTASKS" ]]     && sbatch_args+=(--ntasks="$NTASKS")
    [[ -n "$TIME" ]]       && sbatch_args+=(--time="$TIME")
    [[ -n "$JOB_NAME" ]]   && sbatch_args+=(--job-name="$JOB_NAME")
    [[ -n "$CONSTRAINT" ]] && sbatch_args+=(--constraint="$CONSTRAINT")
    [[ -n "$NODELIST" ]]   && sbatch_args+=(--nodelist="$NODELIST")

    random_begin_args
    echo "Submitting the cascade from $CASCADE_DIR ..."
    cascade_export="ALL"
    [[ -n "$LMP_BIN" ]] && cascade_export+=",LMP_BIN=$LMP_BIN"
    CASCADE_JOBID="$(cd "$CASCADE_DIR" && sbatch --parsable \
      "${sbatch_args[@]+"${sbatch_args[@]}"}" \
      "${BEGIN_ARGS[@]+"${BEGIN_ARGS[@]}"}" \
      --export="$cascade_export" \
      lammps_submit.slurm)"
    echo "  cascade job id: $CASCADE_JOBID"
  fi
fi

cascade_dependency_args() {
  DEP_ARGS=()
  [[ -n "$CASCADE_JOBID" ]] && DEP_ARGS=(--dependency=afterok:"$CASCADE_JOBID")
}

stagger() {
  if [[ "$1" -gt 1 && "$INTERACTIVE" != "1" && "$STAGGER_SECONDS" -gt 0 ]]; then
    echo "Staggering submission by ${STAGGER_SECONDS}s ..."
    sleep "$STAGGER_SECONDS"
  fi
}

############################
# Stage 2: dielectric, one independent pair of jobs per temperature
############################

DIELECTRIC_JOBIDS=()
DIELECTRIC_ENTRIES="$SWEEP_DIR/.dielectric_entries"
: > "$DIELECTRIC_ENTRIES"

INDEX=0
for T in "${DIELECTRIC_TEMPS[@]+"${DIELECTRIC_TEMPS[@]}"}"; do
  INDEX=$((INDEX + 1))
  stagger "$INDEX"

  TK="$(to_kelvin "$T")"
  echo ""
  echo "--- Dielectric, T = ${T} C (${TK} K) ---"

  STAGE_DIR="$SWEEP_DIR/T${T}C"
  CALC_DIR="$ANALYSIS_PARENT_DIR/${SWEEP_ID}_T${T}C_dielectric_calc"
  THERMALIZED="$CASCADE_DIR/thermalized_T${T}C.data"

  for existing in "$STAGE_DIR" "$CALC_DIR"; do
    if [[ -e "$existing" ]]; then
      if [[ "$FORCE" == "1" ]]; then
        log_overwrite "--force ($FORCE_REASON): removing existing $existing (sweep id: $SWEEP_ID, T=${T}C)"
        rm -rf "$existing"
      else
        echo "Error: $existing already exists — refusing to overwrite (use --force to override)." >&2
        exit 1
      fi
    fi
  done

  mkdir -p "$STAGE_DIR/input_files" "$STAGE_DIR/run" "$CALC_DIR/logs"

  cp "$DIELECTRIC_INPUT" "$STAGE_DIR/input_files/in.input"
  # The structure this temperature's production run starts from is the one the
  # cascade wrote for it — not a hand-picked file reused across temperatures.
  ln -s "$THERMALIZED" "$STAGE_DIR/input_files/start.data"
  ln -s "$(realpath "$POTENTIAL_FILE")" "$STAGE_DIR/input_files/$POTENTIAL_LINK_NAME"
  cp "$DIELECTRIC_TEMPLATE" "$STAGE_DIR/run/lammps_submit.slurm"

  cp "$ANALYSIS_TEMPLATE_DIR/1.calc_mpi.py" \
     "$ANALYSIS_TEMPLATE_DIR/2.dipole_std.py" \
     "$ANALYSIS_TEMPLATE_DIR/dielectric_submit.slurm" \
     "$CALC_DIR/"
  ln -s "$STAGE_DIR/run/dumps" "$CALC_DIR/dumps"

  prod_export="TARGET_TEMP=$TK,"
  [[ -n "$LMP_BIN" ]] && prod_export+="LMP_BIN=$LMP_BIN,"
  prod_export+="N_TIMES=$DIELECTRIC_N_CHUNKS,NVT_LENGTH=$DIELECTRIC_CHUNK_LENGTH,DUMP_EVERY=$DIELECTRIC_DUMP_EVERY,CHARGE_H=$DIELECTRIC_CHARGE_H"
  calc_export="ALL,DUMP_DIR=$CALC_DIR/dumps,DUMP_EVERY=$DIELECTRIC_DUMP_EVERY,CUTOFF=$DIELECTRIC_CUTOFF,TYPE_O=$DIELECTRIC_TYPE_O,TYPE_H=$DIELECTRIC_TYPE_H,CHARGE_H=$DIELECTRIC_CHARGE_H,TEMPERATURE=$TK,LA=$LA,LB=$LB,LC=$LC,AVERAGING_METHOD=$DIELECTRIC_AVERAGING_METHOD,NRANKS=$DIELECTRIC_NTASKS"

  if [[ "$INTERACTIVE" == "1" ]]; then
    echo "Running dielectric production interactively in $STAGE_DIR/run ..."
    (
      export SLURM_SUBMIT_DIR="$STAGE_DIR/run"
      export SLURM_NTASKS="${NTASKS:-64}"
      export TARGET_TEMP="$TK" N_TIMES="$DIELECTRIC_N_CHUNKS" NVT_LENGTH="$DIELECTRIC_CHUNK_LENGTH" DUMP_EVERY="$DIELECTRIC_DUMP_EVERY" CHARGE_H="$DIELECTRIC_CHARGE_H"
      [[ -n "$LMP_BIN" ]] && export LMP_BIN
      cd "$STAGE_DIR/run" && bash lammps_submit.slurm
    )
    echo "Running dielectric calc interactively in $CALC_DIR ..."
    (
      export SLURM_SUBMIT_DIR="$CALC_DIR"
      IFS=',' read -ra _kv <<< "${calc_export#ALL,}"
      for pair in "${_kv[@]}"; do export "$pair"; done
      cd "$CALC_DIR" && bash dielectric_submit.slurm
    )
  else
    prod_args=()
    [[ -n "$DIELECTRIC_NODES" ]]      && prod_args+=(--nodes="$DIELECTRIC_NODES")
    [[ -n "$NTASKS" ]]                && prod_args+=(--ntasks="$NTASKS")
    [[ -n "$DIELECTRIC_TIME" ]]       && prod_args+=(--time="$DIELECTRIC_TIME")
    [[ -n "$DIELECTRIC_JOB_NAME" ]]   && prod_args+=(--job-name="${DIELECTRIC_JOB_NAME}-T${T}C")
    [[ -n "$DIELECTRIC_CONSTRAINT" ]] && prod_args+=(--constraint="$DIELECTRIC_CONSTRAINT")
    [[ -n "$DIELECTRIC_NODELIST" ]]   && prod_args+=(--nodelist="$DIELECTRIC_NODELIST")

    cascade_dependency_args
    random_begin_args
    PROD_JOBID="$(cd "$STAGE_DIR/run" && sbatch --parsable \
      "${DEP_ARGS[@]+"${DEP_ARGS[@]}"}" \
      "${prod_args[@]+"${prod_args[@]}"}" \
      "${BEGIN_ARGS[@]+"${BEGIN_ARGS[@]}"}" \
      --export="ALL,$prod_export" \
      lammps_submit.slurm)"
    echo "  production job id: $PROD_JOBID"

    calc_args=()
    [[ -n "$ANALYSIS_NODES" ]]      && calc_args+=(--nodes="$ANALYSIS_NODES")
    [[ -n "$ANALYSIS_TIME" ]]       && calc_args+=(--time="$ANALYSIS_TIME")
    [[ -n "$ANALYSIS_JOB_NAME" ]]   && calc_args+=(--job-name="${ANALYSIS_JOB_NAME}-diel-T${T}C")
    [[ -n "$ANALYSIS_CONSTRAINT" ]] && calc_args+=(--constraint="$ANALYSIS_CONSTRAINT")
    [[ -n "$ANALYSIS_NODELIST" ]]   && calc_args+=(--nodelist="$ANALYSIS_NODELIST")
    [[ -n "$DIELECTRIC_NTASKS" ]]   && calc_args+=(--ntasks="$DIELECTRIC_NTASKS")
    calc_args+=(--mem="$DIELECTRIC_MEM")

    random_begin_args
    CALC_JOBID="$(cd "$CALC_DIR" && sbatch --parsable \
      --dependency=afterok:"$PROD_JOBID" \
      "${calc_args[@]+"${calc_args[@]}"}" \
      "${BEGIN_ARGS[@]+"${BEGIN_ARGS[@]}"}" \
      --export="$calc_export" \
      dielectric_submit.slurm)"
    echo "  dielectric-calc job id: $CALC_JOBID"
    DIELECTRIC_JOBIDS+=("$CALC_JOBID")
  fi

  echo "${T}:${TK}:${CALC_DIR}/dipole_output/summary.txt" >> "$DIELECTRIC_ENTRIES"
done

############################
# Stage 3: diffusion, one independent job per temperature
############################

DIFFUSION_JOBIDS=()
DIFFUSION_ENTRIES="$SWEEP_DIR/.diffusion_entries"
: > "$DIFFUSION_ENTRIES"

INDEX=0
for T in "${DIFFUSION_TEMPS[@]+"${DIFFUSION_TEMPS[@]}"}"; do
  INDEX=$((INDEX + 1))
  stagger "$INDEX"

  TK="$(to_kelvin "$T")"
  echo ""
  echo "--- Diffusion, T = ${T} C (${TK} K) ---"

  MSD_DIR="$ANALYSIS_PARENT_DIR/${SWEEP_ID}_T${T}C_msd"
  TRAJ="$CASCADE_DIR/dynamics_T${T}C.lammpstrj"

  if [[ -e "$MSD_DIR" ]]; then
    if [[ "$FORCE" == "1" ]]; then
      log_overwrite "--force ($FORCE_REASON): removing existing $MSD_DIR (sweep id: $SWEEP_ID, T=${T}C)"
      rm -rf "$MSD_DIR"
    else
      echo "Error: $MSD_DIR already exists — refusing to overwrite (use --force to override)." >&2
      exit 1
    fi
  fi

  mkdir -p "$MSD_DIR"
  # distribution_submit.slurm is reused rather than replaced by a new runner:
  # it already has the RUN_MSD / MSD_* / DYNAMICS_TRAJ environment contract, and
  # every other RUN_* flag is switched off below.
  cp "$MSD_TEMPLATE_DIR/msd.py" \
     "$MSD_TEMPLATE_DIR/distribution_submit.slurm" \
     "$MSD_DIR/"

  msd_export="ALL,RUN_DSF=0,RUN_RDF=0,RUN_BAD=0,RUN_VDOS=0,RUN_VDOS_DYNMAT=0,RUN_MSD=1,RUN_PLOTS=0,DYNAMICS_TRAJ=$TRAJ,DYNAMICS_DT=$DYNAMICS_DT,MSD_CORR_LENGTH=$MSD_CORR_LENGTH,MSD_CORR_INTERVAL=$MSD_CORR_INTERVAL,MSD_FIT_FRACTION=$MSD_FIT_FRACTION"

  if [[ "$INTERACTIVE" == "1" ]]; then
    echo "Running msd interactively in $MSD_DIR ..."
    (
      export SLURM_SUBMIT_DIR="$MSD_DIR"
      IFS=',' read -ra _kv <<< "${msd_export#ALL,}"
      for pair in "${_kv[@]}"; do export "$pair"; done
      cd "$MSD_DIR" && bash distribution_submit.slurm
    )
  else
    msd_args=()
    [[ -n "$ANALYSIS_NODES" ]]      && msd_args+=(--nodes="$ANALYSIS_NODES")
    [[ -n "$ANALYSIS_TIME" ]]       && msd_args+=(--time="$ANALYSIS_TIME")
    [[ -n "$ANALYSIS_JOB_NAME" ]]   && msd_args+=(--job-name="${ANALYSIS_JOB_NAME}-msd-T${T}C")
    [[ -n "$ANALYSIS_CONSTRAINT" ]] && msd_args+=(--constraint="$ANALYSIS_CONSTRAINT")
    [[ -n "$ANALYSIS_NODELIST" ]]   && msd_args+=(--nodelist="$ANALYSIS_NODELIST")

    cascade_dependency_args
    random_begin_args
    MSD_JOBID="$(cd "$MSD_DIR" && sbatch --parsable \
      "${DEP_ARGS[@]+"${DEP_ARGS[@]}"}" \
      "${msd_args[@]+"${msd_args[@]}"}" \
      "${BEGIN_ARGS[@]+"${BEGIN_ARGS[@]}"}" \
      --export="$msd_export" \
      distribution_submit.slurm)"
    echo "  msd job id: $MSD_JOBID"
    DIFFUSION_JOBIDS+=("$MSD_JOBID")
  fi

  echo "${T}:${TK}:${MSD_DIR}" >> "$DIFFUSION_ENTRIES"
done

############################
# Aggregation: one job per route, fanning in
############################

DIELECTRIC_CSV="$SWEEP_DIR/dielectric_vs_temperature.csv"
DIFFUSION_CSV="$SWEEP_DIR/diffusion_vs_temperature.csv"

run_aggregate_interactive() {
  local script="$1" entries="$2" output="$3"
  local entry_args=() line
  while IFS= read -r line; do
    [[ -n "$line" ]] && entry_args+=(--entry "$line")
  done < "$entries"
  [[ "${#entry_args[@]}" -eq 0 ]] && return 0
  python "$script" "${entry_args[@]}" --output "$output"
}

submit_aggregate() {
  local script="$1" entries="$2" output="$3" name="$4"; shift 4
  local jobids=("$@")
  [[ "${#jobids[@]}" -eq 0 ]] && return 0

  local dependency="afterok"
  local jobid
  for jobid in "${jobids[@]}"; do dependency+=":$jobid"; done

  random_begin_args
  local aggregate_jobid
  aggregate_jobid="$(cd "$SWEEP_DIR" && sbatch --parsable --dependency="$dependency" \
    "${BEGIN_ARGS[@]+"${BEGIN_ARGS[@]}"}" \
    --job-name="${name}-aggregate" \
    --export="ALL,ENTRIES_FILE=$entries,OUTPUT_CSV=$output,AGGREGATE_SCRIPT=$script,VENV_PATH=$VENV_PATH" \
    "$AGGREGATE_TEMPLATE")"
  echo "  ${name} aggregation job id: $aggregate_jobid"
}

echo ""
if [[ "$INTERACTIVE" == "1" ]]; then
  echo "Aggregating ..."
  [[ "$N_DIELECTRIC" -gt 0 ]] && run_aggregate_interactive "$AGGREGATE_DIELECTRIC_SCRIPT" "$DIELECTRIC_ENTRIES" "$DIELECTRIC_CSV"
  [[ "$N_DIFFUSION"  -gt 0 ]] && run_aggregate_interactive "$AGGREGATE_DIFFUSION_SCRIPT"  "$DIFFUSION_ENTRIES"  "$DIFFUSION_CSV"

  cat <<SUMMARY

Pipeline completed interactively:
  Sweep directory      : $SWEEP_DIR
  Diffusion temps (C)  : ${DIFFUSION_TEMPS[*]+"${DIFFUSION_TEMPS[*]}"}
  Dielectric temps (C) : ${DIELECTRIC_TEMPS[*]+"${DIELECTRIC_TEMPS[*]}"}
SUMMARY
  [[ "$N_DIELECTRIC" -gt 0 ]] && echo "  Dielectric summary   : $DIELECTRIC_CSV"
  [[ "$N_DIFFUSION"  -gt 0 ]] && echo "  Diffusion summary    : $DIFFUSION_CSV"
else
  echo "Submitting aggregation jobs ..."
  submit_aggregate "$AGGREGATE_DIELECTRIC_SCRIPT" "$DIELECTRIC_ENTRIES" "$DIELECTRIC_CSV" "dielectric" \
    "${DIELECTRIC_JOBIDS[@]+"${DIELECTRIC_JOBIDS[@]}"}"
  submit_aggregate "$AGGREGATE_DIFFUSION_SCRIPT" "$DIFFUSION_ENTRIES" "$DIFFUSION_CSV" "diffusion" \
    "${DIFFUSION_JOBIDS[@]+"${DIFFUSION_JOBIDS[@]}"}"

  cat <<SUMMARY

Pipeline submitted:
  Sweep directory      : $SWEEP_DIR
  Diffusion temps (C)  : ${DIFFUSION_TEMPS[*]+"${DIFFUSION_TEMPS[*]}"}
  Dielectric temps (C) : ${DIELECTRIC_TEMPS[*]+"${DIELECTRIC_TEMPS[*]}"}
  Cascade job          : ${CASCADE_JOBID:-<skipped>}
  Dielectric-calc jobs : ${DIELECTRIC_JOBIDS[*]+"${DIELECTRIC_JOBIDS[*]}"}
  MSD jobs             : ${DIFFUSION_JOBIDS[*]+"${DIFFUSION_JOBIDS[*]}"}
SUMMARY
  [[ "$N_DIELECTRIC" -gt 0 ]] && echo "  Dielectric summary (once its aggregation finishes) : $DIELECTRIC_CSV"
  [[ "$N_DIFFUSION"  -gt 0 ]] && echo "  Diffusion summary  (once its aggregation finishes) : $DIFFUSION_CSV"
fi
