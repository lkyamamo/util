#!/bin/bash
# submit_pipeline_local.sh — run the full LAMMPS setup+trajectory job, then
# distribution analysis (rdf/bad/dsf/vdos/msd), entirely locally via mpirun/plain
# python — no Slurm, no modules, no sbatch. Self-contained: does not shell
# out to jobs/slurm/lammps_submit.slurm or distribution_submit.slurm at all;
# the commands those templates run are inlined directly below (minus the
# SBATCH/module/srun machinery, which don't apply outside a cluster).
#
# Run from inside the LAMMPS run directory; its basename is the run id.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: submit_pipeline_local.sh [options]

Run from inside the LAMMPS run directory, whose name is the run id
(e.g. cd .../runs/water/N5184/0149 && /path/to/submit_pipeline_local.sh ...).

Required:
  --input-script FILE          LAMMPS .input script (setup + trajectory generation)
  --starting-structure FILE    Starting structure .data file, symlinked as start.data
  --potential-file FILE        Potential file, symlinked into input_files/
  --analysis-parent-dir DIR    Where <run_id>_distribution_analysis/ is created
  --lmp-bin FILE               Local LAMMPS executable (no HPC default exists)
  --venv DIR                   Path to the analysis environment (dynasor/freud/numpy/
                                matplotlib installed) — a plain venv/virtualenv
                                (DIR/bin/activate is sourced) or a conda env
                                (detected via DIR/conda-meta/, activated via
                                `conda activate DIR`)

Note: --potential-link-name NAME (optional but often required in practice)
sets the symlink's name in input_files/ (e.g. "OH.usc"). It must match
whatever your in.input's pair_coeff line expects, which frequently differs
from --potential-file's own basename (e.g. a real potential file named
20260723_OH.vashishta but pair_coeff expects OH.usc). Defaults to
--potential-file's basename if omitted.

Note: --input-script's in.input must write two trajectory dumps, flat (no
subdirectory): dump.lammpstrj (read by rdf_freud.py/bad_freud.py) and a
separate, higher-frequency dynamics.lammpstrj (read by dsf.py/vdos.py/msd.py
— resolving vibrational frequencies and diffusion needs much finer time
sampling than structural analysis does). See OH-therm.input/b-SiO-therm.input
in this directory for examples of both dump commands.

Optional:
  --analysis-template-dir DIR   Distribution-analysis scripts to copy into stage 2
                                 (default: <repo>/analysis/distributions/20260608_GrNrBaSqw)
  --run-dsf {0|1}                Run dsf.py (static/dynamic structure factor) (default: 1)
  --run-rdf {0|1}                Run rdf_freud.py (radial distribution function) (default: 1)
  --run-bad {0|1}                Run bad_freud.py (bond angle distribution) (default: 1)
  --run-vdos {0|1}                Run vdos.py (vibrational density of states) (default: 1)
  --run-msd {0|1}                 Run msd.py (mean square displacement / diffusion) (default: 1)

  Physics config — same flags and formats as submit_pipeline.sh, so the same
  invocation recipe works against either script:
    --rdf-r-max FLOAT            rdf_freud.py R_MAX (Å), e.g. 20.0
    --rdf-bins INT                rdf_freud.py BINS, e.g. 2000
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
    --dsf-dt FLOAT                dsf.py DT (fs between dumped frames), e.g. 1.0
    --dsf-n-frames INT             dsf.py N_FRAMES, e.g. 500
    --dsf-stride INT               dsf.py STRIDE, e.g. 1
    --dsf-window-size INT           dsf.py WINDOW_SIZE, e.g. 500
    --dsf-q-max FLOAT              dsf.py Q_MAX (Å⁻¹), e.g. 20.0
    --dsf-n-q-bins INT              dsf.py N_Q_BINS, e.g. 200

    vdos.py reuses --dsf-n-frames/--dsf-stride above for its own N_FRAMES/
    STRIDE (same dump, same meaning — no separate flags):
    --vdos-dt FLOAT                 vdos.py TIME_UNIT (fs between dumped frames),
                                    e.g. 1.0 — falls back to --dsf-dt/DT if omitted,
                                    but set this explicitly if running vdos.py
                                    without dsf.py or with a different DT
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
    --vdos-normalization STR         vdos.py NORMALIZATION: 'phonon' (default) or 'unit_area'

    msd.py reuses --vdos-dt/--vdos-corr-length/--vdos-corr-interval above for
    its own TIME_UNIT/CORR_LENGTH/CORR_INTERVAL (same trajectory, same
    meaning, same env vars — no separate flags):
    --msd-fit-fraction FLOAT         msd.py FIT_FRACTION: fraction of the tail of the
                                    correlation window used for the diffusion-coefficient
                                    linear fit, e.g. 0.5

  --ntasks N                    MPI ranks for `mpirun -np` (default: 4)
  --analysis-cpus-per-task N    Sets OMP_NUM_THREADS for dsf.py/rdf_freud.py/
                                bad_freud.py and VDOS_THREADS for vdos.py (default: 4)

  --force REASON                  Overwrite existing input_files/run/ (stage 1) or an
                                  existing <run_id>_distribution_analysis/ (stage 2)
                                  instead of refusing to run. REASON is required (a
                                  short explanation of why you're overwriting) and is
                                  logged, with a timestamp and the exact path removed,
                                  to both stderr and
                                  jobs/pipeline/lammps_to_analysis/overwrite.log.

  -h, --help                    Show this help

All terminal output from this script is also appended to
submit_pipeline_local.log in the run directory (cwd) each time it's invoked.
EOF
}

REPO_ROOT="$HOME/util"
DUMP_FILE="dump.lammpstrj"
DYNAMICS_DUMP_FILE="dynamics.lammpstrj"
FORCE="0"
FORCE_REASON=""
LOG_FILE="$REPO_ROOT/jobs/pipeline/lammps_to_analysis/overwrite.log"

log_overwrite() {
  local msg
  msg="[$(date "+%Y-%m-%dT%H:%M:%S%z")] $1"
  echo "WARNING: $msg" >&2
  echo "$msg" >> "$LOG_FILE"
}

# Activate --venv, whether it's a plain venv/virtualenv or a conda env.
# conda-created environments always contain a conda-meta/ directory — that's
# the reliable marker (rather than assuming based on path).
#
# A plain venv's own bin/activate script is safe to just source. A conda
# env's bin/activate is NOT: 'conda activate' only works via the conda shell
# FUNCTION that 'conda init' installs into interactive shells (~/.zshrc,
# ~/.bash_profile) — it modifies the calling shell's own environment, which
# a bare subprocess can never do. This script runs as a non-interactive bash
# subshell, so it inherits $PATH (the plain `conda` executable still runs)
# but not that shell function — calling `conda activate` directly fails with
# "CondaError: Run 'conda init' before 'conda activate'". Fix: install the
# same hook conda init would, scoped to just this script, via
# `eval "$(conda shell.bash hook)"`, before calling `conda activate`.
activate_analysis_env() {
  local env_path="$1"
  if [[ -d "$env_path/conda-meta" ]]; then
    if ! command -v conda >/dev/null 2>&1; then
      echo "Error: --venv ($env_path) looks like a conda environment, but 'conda' isn't on PATH." >&2
      exit 1
    fi
    eval "$(conda shell.bash hook)"
    conda activate "$env_path"
  else
    source "$env_path/bin/activate"
  fi
}

INPUT_SCRIPT=""
STARTING_STRUCTURE=""
POTENTIAL_FILE=""
POTENTIAL_LINK_NAME=""
ANALYSIS_PARENT_DIR=""
ANALYSIS_TEMPLATE_DIR="$REPO_ROOT/analysis/distributions/20260608_GrNrBaSqw"
LMP_BIN=""
VENV=""
RUN_DSF="1"
RUN_RDF="1"
RUN_BAD="1"
RUN_VDOS="1"
RUN_MSD="1"
NTASKS=""
ANALYSIS_CPUS_PER_TASK=""

RDF_R_MAX=""
RDF_BINS_VAL=""
BAD_ELEMENTS=""
BAD_R_CUTOFF=""
BAD_R_MINCUT=""
BAD_TRIPLET_CUTOFFS=""
BAD_BINS_VAL=""
DSF_DT=""
DSF_N_FRAMES=""
DSF_STRIDE=""
DSF_WINDOW_SIZE=""
DSF_Q_MAX=""
DSF_N_Q_BINS=""
VDOS_DT=""
VDOS_CORR_LENGTH=""
VDOS_CORR_INTERVAL=""
VDOS_MAX_FREQUENCY_EV=""
VDOS_NUM_GRIDS=""
VDOS_METHOD=""
VDOS_WINDOW=""
VDOS_NORMALIZATION=""
MSD_FIT_FRACTION=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-script) INPUT_SCRIPT="$2"; shift 2 ;;
    --starting-structure) STARTING_STRUCTURE="$2"; shift 2 ;;
    --potential-file) POTENTIAL_FILE="$2"; shift 2 ;;
    --potential-link-name) POTENTIAL_LINK_NAME="$2"; shift 2 ;;
    --analysis-parent-dir) ANALYSIS_PARENT_DIR="$2"; shift 2 ;;
    --analysis-template-dir) ANALYSIS_TEMPLATE_DIR="$2"; shift 2 ;;
    --lmp-bin) LMP_BIN="$2"; shift 2 ;;
    --venv) VENV="$2"; shift 2 ;;
    --run-dsf) RUN_DSF="$2"; shift 2 ;;
    --run-rdf) RUN_RDF="$2"; shift 2 ;;
    --run-bad) RUN_BAD="$2"; shift 2 ;;
    --run-vdos) RUN_VDOS="$2"; shift 2 ;;
    --run-msd) RUN_MSD="$2"; shift 2 ;;
    --ntasks) NTASKS="$2"; shift 2 ;;
    --analysis-cpus-per-task) ANALYSIS_CPUS_PER_TASK="$2"; shift 2 ;;
    --rdf-r-max) RDF_R_MAX="$2"; shift 2 ;;
    --rdf-bins) RDF_BINS_VAL="$2"; shift 2 ;;
    --bad-elements) BAD_ELEMENTS="$2"; shift 2 ;;
    --bad-r-cutoff) BAD_R_CUTOFF="$2"; shift 2 ;;
    --bad-r-mincut) BAD_R_MINCUT="$2"; shift 2 ;;
    --bad-triplet-cutoffs) BAD_TRIPLET_CUTOFFS="$2"; shift 2 ;;
    --bad-bins) BAD_BINS_VAL="$2"; shift 2 ;;
    --dsf-dt) DSF_DT="$2"; shift 2 ;;
    --dsf-n-frames) DSF_N_FRAMES="$2"; shift 2 ;;
    --dsf-stride) DSF_STRIDE="$2"; shift 2 ;;
    --dsf-window-size) DSF_WINDOW_SIZE="$2"; shift 2 ;;
    --dsf-q-max) DSF_Q_MAX="$2"; shift 2 ;;
    --dsf-n-q-bins) DSF_N_Q_BINS="$2"; shift 2 ;;
    --vdos-dt) VDOS_DT="$2"; shift 2 ;;
    --vdos-corr-length) VDOS_CORR_LENGTH="$2"; shift 2 ;;
    --vdos-corr-interval) VDOS_CORR_INTERVAL="$2"; shift 2 ;;
    --vdos-max-frequency-ev) VDOS_MAX_FREQUENCY_EV="$2"; shift 2 ;;
    --vdos-num-grids) VDOS_NUM_GRIDS="$2"; shift 2 ;;
    --vdos-method) VDOS_METHOD="$2"; shift 2 ;;
    --vdos-window) VDOS_WINDOW="$2"; shift 2 ;;
    --vdos-normalization) VDOS_NORMALIZATION="$2"; shift 2 ;;
    --msd-fit-fraction) MSD_FIT_FRACTION="$2"; shift 2 ;;
    --force) FORCE="1"; FORCE_REASON="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

STAGE1_DIR="$(pwd)"
RUN_ID="$(basename "$STAGE1_DIR")"
PIPELINE_LOG="$STAGE1_DIR/submit_pipeline_local.log"
exec > >(tee -a "$PIPELINE_LOG") 2>&1
echo "=== submit_pipeline_local.sh started $(date "+%Y-%m-%dT%H:%M:%S%z") — run id: $RUN_ID ==="
echo "Logging this run's output to: $PIPELINE_LOG"

missing=0
require() {
  if [[ -z "$1" ]]; then
    echo "Missing required argument: $2" >&2
    missing=1
  fi
}
require "$INPUT_SCRIPT" --input-script
require "$STARTING_STRUCTURE" --starting-structure
require "$POTENTIAL_FILE" --potential-file
require "$ANALYSIS_PARENT_DIR" --analysis-parent-dir
require "$LMP_BIN" --lmp-bin
require "$VENV" --venv
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

if [[ "$FORCE" == "1" && -z "$FORCE_REASON" ]]; then
  echo 'Error: --force requires a non-empty reason, e.g. --force "re-running after fixing potential file"' >&2
  exit 1
fi

############################
# Stage 1 setup (unconditional)
############################

if [[ -e "$STAGE1_DIR/input_files" || -e "$STAGE1_DIR/run" ]]; then
  if [[ "$FORCE" == "1" ]]; then
    log_overwrite "--force ($FORCE_REASON): removing existing $STAGE1_DIR/input_files and/or $STAGE1_DIR/run (run id: $RUN_ID)"
    rm -rf "$STAGE1_DIR/input_files" "$STAGE1_DIR/run"
  else
    echo "Error: $STAGE1_DIR already has input_files/ or run/ — refusing to overwrite (use --force to override)." >&2
    exit 1
  fi
fi

mkdir -p "$STAGE1_DIR/input_files" "$STAGE1_DIR/run"
cp "$INPUT_SCRIPT" "$STAGE1_DIR/input_files/in.input"
ln -s "$(realpath "$STARTING_STRUCTURE")" "$STAGE1_DIR/input_files/start.data"
ln -s "$(realpath "$POTENTIAL_FILE")" "$STAGE1_DIR/input_files/${POTENTIAL_LINK_NAME:-$(basename "$POTENTIAL_FILE")}"

############################
# Stage 2 setup (unconditional)
############################

mkdir -p "$ANALYSIS_PARENT_DIR"
ANALYSIS_PARENT_DIR="$(cd "$ANALYSIS_PARENT_DIR" && pwd)"
STAGE2_DIR="$ANALYSIS_PARENT_DIR/${RUN_ID}_distribution_analysis"

if [[ -e "$STAGE2_DIR" ]]; then
  if [[ "$FORCE" == "1" ]]; then
    log_overwrite "--force ($FORCE_REASON): removing existing $STAGE2_DIR"
    rm -rf "$STAGE2_DIR"
  else
    echo "Error: $STAGE2_DIR already exists — refusing to overwrite (use --force to override)." >&2
    exit 1
  fi
fi

mkdir -p "$STAGE2_DIR"
cp "$ANALYSIS_TEMPLATE_DIR/dsf.py" \
   "$ANALYSIS_TEMPLATE_DIR/rdf_freud.py" \
   "$ANALYSIS_TEMPLATE_DIR/bad_freud.py" \
   "$ANALYSIS_TEMPLATE_DIR/vdos.py" \
   "$ANALYSIS_TEMPLATE_DIR/msd.py" \
   "$STAGE2_DIR/"
ln -s "$STAGE1_DIR/run/$DUMP_FILE" "$STAGE2_DIR/$DUMP_FILE"
ln -s "$STAGE1_DIR/run/$DYNAMICS_DUMP_FILE" "$STAGE2_DIR/$DYNAMICS_DUMP_FILE"

enabled_scripts=""
[[ "$RUN_DSF" == "1" ]] && enabled_scripts+="dsf.py "
[[ "$RUN_RDF" == "1" ]] && enabled_scripts+="rdf_freud.py "
[[ "$RUN_BAD" == "1" ]] && enabled_scripts+="bad_freud.py "
[[ "$RUN_VDOS" == "1" ]] && enabled_scripts+="vdos.py "
[[ "$RUN_MSD" == "1" ]] && enabled_scripts+="msd.py "
enabled_scripts="${enabled_scripts% }"

############################
# Stage 1 execution — inlined from jobs/slurm/lammps_submit.slurm
# (module purge/load and srun dropped; mpirun used instead)
############################

echo "Running LAMMPS locally in $STAGE1_DIR/run ..."
(
  cd "$STAGE1_DIR/run"
  shopt -s nullglob
  for entry in "$STAGE1_DIR/input_files"/*; do
    ln -sfn "$entry" "$(basename "$entry")"
  done
  shopt -u nullglob

  mpirun -np "${NTASKS:-4}" "$LMP_BIN" -log log.lammps -in in.input
)
echo "LAMMPS run finished."

############################
# Stage 2 execution — inlined from distribution_submit.slurm
# (module purge/load dropped; venv path is a flag, not a hardcoded HPC path)
############################

echo "Running distribution analysis locally in $STAGE2_DIR ..."
(
  cd "$STAGE2_DIR"
  activate_analysis_env "$VENV"
  export OMP_NUM_THREADS="${ANALYSIS_CPUS_PER_TASK:-4}"
  export VDOS_THREADS="${ANALYSIS_CPUS_PER_TASK:-4}"
  export TRAJ="$DUMP_FILE" DYNAMICS_TRAJ="$DYNAMICS_DUMP_FILE" RUN_DSF RUN_RDF RUN_BAD RUN_VDOS RUN_MSD
  [[ -n "$RDF_R_MAX" ]]           && export R_MAX="$RDF_R_MAX"
  [[ -n "$RDF_BINS_VAL" ]]        && export RDF_BINS="$RDF_BINS_VAL"
  [[ -n "$BAD_ELEMENTS" ]]        && export ELEMENTS="$BAD_ELEMENTS"
  [[ -n "$BAD_R_CUTOFF" ]]        && export R_CUTOFF="$BAD_R_CUTOFF"
  [[ -n "$BAD_R_MINCUT" ]]        && export R_MINCUT="$BAD_R_MINCUT"
  [[ -n "$BAD_TRIPLET_CUTOFFS" ]] && export TRIPLET_CUTOFFS="$BAD_TRIPLET_CUTOFFS"
  [[ -n "$BAD_BINS_VAL" ]]        && export BAD_BINS="$BAD_BINS_VAL"
  [[ -n "$DSF_DT" ]]              && export DT="$DSF_DT"
  [[ -n "$DSF_N_FRAMES" ]]        && export N_FRAMES="$DSF_N_FRAMES"
  [[ -n "$DSF_STRIDE" ]]          && export STRIDE="$DSF_STRIDE"
  [[ -n "$DSF_WINDOW_SIZE" ]]     && export WINDOW_SIZE="$DSF_WINDOW_SIZE"
  [[ -n "$DSF_Q_MAX" ]]           && export Q_MAX="$DSF_Q_MAX"
  [[ -n "$DSF_N_Q_BINS" ]]        && export N_Q_BINS="$DSF_N_Q_BINS"
  [[ -n "$VDOS_DT" ]]               && export TIME_UNIT="$VDOS_DT"
  [[ -n "$VDOS_CORR_LENGTH" ]]      && export CORR_LENGTH="$VDOS_CORR_LENGTH"
  [[ -n "$VDOS_CORR_INTERVAL" ]]    && export CORR_INTERVAL="$VDOS_CORR_INTERVAL"
  [[ -n "$VDOS_MAX_FREQUENCY_EV" ]] && export MAX_FREQUENCY_EV="$VDOS_MAX_FREQUENCY_EV"
  [[ -n "$VDOS_NUM_GRIDS" ]]        && export NUM_GRIDS="$VDOS_NUM_GRIDS"
  [[ -n "$VDOS_METHOD" ]]           && export VDOS_METHOD="$VDOS_METHOD"
  [[ -n "$VDOS_WINDOW" ]]           && export VDOS_WINDOW="$VDOS_WINDOW"
  [[ -n "$VDOS_NORMALIZATION" ]]    && export VDOS_NORMALIZATION="$VDOS_NORMALIZATION"
  [[ -n "$MSD_FIT_FRACTION" ]]       && export MSD_FIT_FRACTION="$MSD_FIT_FRACTION"

  if [[ "$RUN_DSF" == "1" ]]; then echo "--- dsf.py ---"; python dsf.py; fi
  if [[ "$RUN_RDF" == "1" ]]; then echo "--- rdf_freud.py ---"; python rdf_freud.py; fi
  if [[ "$RUN_BAD" == "1" ]]; then echo "--- bad_freud.py ---"; python bad_freud.py; fi
  if [[ "$RUN_VDOS" == "1" ]]; then echo "--- vdos.py ---"; python vdos.py; fi
  if [[ "$RUN_MSD" == "1" ]]; then echo "--- msd.py ---"; python msd.py; fi
)
echo "Distribution analysis finished."

cat <<SUMMARY

Pipeline completed locally:
  LAMMPS run (setup + trajectory) : $STAGE1_DIR
  distribution analysis           : $STAGE2_DIR
  analysis scripts enabled        : ${enabled_scripts:-none}

$(if [[ -n "$enabled_scripts" ]]; then
  echo "Reminder: any physics config (ELEMENTS, R_CUTOFF, DT, WINDOW_SIZE, CORR_LENGTH, ...) not"
  echo "passed via --rdf-*/--bad-*/--dsf-*/--vdos-*/--msd-* flags is using each script's own default —"
  echo "check $STAGE2_DIR/{${enabled_scripts// /,}} if that's not what you want."
fi)
SUMMARY
