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

  Physics config — each overrides that script's own hardcoded default
  (see the CONFIGURATION block at the top of each .py) when passed; omit
  to leave the script's built-in default in effect:
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

  --force REASON                  Overwrite existing input_files/run/ (stage 1) or an
                                  existing <run_id>_distribution_analysis/ (stage 2)
                                  instead of refusing to run. REASON is required (a
                                  short explanation of why you're overwriting) and is
                                  logged, with a timestamp and the exact path removed,
                                  to both stderr and
                                  jobs/pipeline/lammps_to_analysis/overwrite.log.

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

All terminal output from this script (not the SLURM jobs themselves) is
also appended to submit_pipeline.log in the run directory (cwd) each time
it's invoked.
EOF
}

INTERACTIVE="0"
SKIP_TRAJECTORY="0"
REPO_ROOT="$HOME/util"
LAMMPS_TEMPLATE="$REPO_ROOT/jobs/slurm/lammps_submit.slurm"
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
# shellcheck source=/dev/null
source "$CONFIG_FILE"

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
    --interactive) INTERACTIVE="1"; shift 1 ;;
    --skip-trajectory) SKIP_TRAJECTORY="1"; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

STAGE1_DIR="$(pwd)"
RUN_ID="$(basename "$STAGE1_DIR")"
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
  for f in "$STAGE1_DIR/run/$DUMP_FILE" "$STAGE1_DIR/run/$DYNAMICS_DUMP_FILE"; do
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
   "$ANALYSIS_TEMPLATE_DIR/distribution_submit.slurm" \
   "$STAGE2_DIR/"
ln -s "$STAGE1_DIR/run/$DUMP_FILE" "$STAGE2_DIR/$DUMP_FILE"
ln -s "$STAGE1_DIR/run/$DYNAMICS_DUMP_FILE" "$STAGE2_DIR/$DYNAMICS_DUMP_FILE"

# None of these values may contain a comma — sbatch --export is
# comma-delimited and silently truncates anything after an embedded one.
export_vars="ALL,TRAJ=$DUMP_FILE,DYNAMICS_TRAJ=$DYNAMICS_DUMP_FILE,RUN_DSF=$RUN_DSF,RUN_RDF=$RUN_RDF,RUN_BAD=$RUN_BAD,RUN_VDOS=$RUN_VDOS,RUN_MSD=$RUN_MSD"
[[ -n "$RDF_R_MAX" ]]           && export_vars+=",R_MAX=$RDF_R_MAX"
[[ -n "$RDF_BINS_VAL" ]]        && export_vars+=",RDF_BINS=$RDF_BINS_VAL"
[[ -n "$BAD_ELEMENTS" ]]        && export_vars+=",ELEMENTS=$BAD_ELEMENTS"
[[ -n "$BAD_R_CUTOFF" ]]        && export_vars+=",R_CUTOFF=$BAD_R_CUTOFF"
[[ -n "$BAD_R_MINCUT" ]]        && export_vars+=",R_MINCUT=$BAD_R_MINCUT"
[[ -n "$BAD_TRIPLET_CUTOFFS" ]] && export_vars+=",TRIPLET_CUTOFFS=$BAD_TRIPLET_CUTOFFS"
[[ -n "$BAD_BINS_VAL" ]]        && export_vars+=",BAD_BINS=$BAD_BINS_VAL"
[[ -n "$DSF_DT" ]]              && export_vars+=",DT=$DSF_DT"
[[ -n "$DSF_N_FRAMES" ]]        && export_vars+=",N_FRAMES=$DSF_N_FRAMES"
[[ -n "$DSF_STRIDE" ]]          && export_vars+=",STRIDE=$DSF_STRIDE"
[[ -n "$DSF_WINDOW_SIZE" ]]     && export_vars+=",WINDOW_SIZE=$DSF_WINDOW_SIZE"
[[ -n "$DSF_Q_MAX" ]]           && export_vars+=",Q_MAX=$DSF_Q_MAX"
[[ -n "$DSF_N_Q_BINS" ]]        && export_vars+=",N_Q_BINS=$DSF_N_Q_BINS"
[[ -n "$VDOS_DT" ]]               && export_vars+=",TIME_UNIT=$VDOS_DT"
[[ -n "$VDOS_CORR_LENGTH" ]]      && export_vars+=",CORR_LENGTH=$VDOS_CORR_LENGTH"
[[ -n "$VDOS_CORR_INTERVAL" ]]    && export_vars+=",CORR_INTERVAL=$VDOS_CORR_INTERVAL"
[[ -n "$VDOS_MAX_FREQUENCY_EV" ]] && export_vars+=",MAX_FREQUENCY_EV=$VDOS_MAX_FREQUENCY_EV"
[[ -n "$VDOS_NUM_GRIDS" ]]        && export_vars+=",NUM_GRIDS=$VDOS_NUM_GRIDS"
[[ -n "$VDOS_METHOD" ]]           && export_vars+=",VDOS_METHOD=$VDOS_METHOD"
[[ -n "$VDOS_WINDOW" ]]           && export_vars+=",VDOS_WINDOW=$VDOS_WINDOW"
[[ -n "$VDOS_NORMALIZATION" ]]    && export_vars+=",VDOS_NORMALIZATION=$VDOS_NORMALIZATION"
[[ -n "$MSD_FIT_FRACTION" ]]       && export_vars+=",MSD_FIT_FRACTION=$MSD_FIT_FRACTION"

enabled_scripts=""
[[ "$RUN_DSF" == "1" ]] && enabled_scripts+="dsf.py "
[[ "$RUN_RDF" == "1" ]] && enabled_scripts+="rdf_freud.py "
[[ "$RUN_BAD" == "1" ]] && enabled_scripts+="bad_freud.py "
[[ "$RUN_VDOS" == "1" ]] && enabled_scripts+="vdos.py "
[[ "$RUN_MSD" == "1" ]] && enabled_scripts+="msd.py "
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
  echo "passed via --rdf-*/--bad-*/--dsf-*/--vdos-*/--msd-* flags is using each script's own default —"
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
      "${sbatch_args[@]+"${sbatch_args[@]}"}" lammps_submit.slurm)"
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
  echo "passed via --rdf-*/--bad-*/--dsf-*/--vdos-*/--msd-* flags is using each script's own default —"
  echo "check $STAGE2_DIR/{${enabled_scripts// /,}} if that's not what you want."
fi)
SUMMARY
fi
