#!/bin/bash
# submit_pipeline.sh — chain a setup LAMMPS run into a trajectory-generating
# LAMMPS run plus distribution analysis (rdf/bad/dsf), via sbatch --dependency.
#
# Run from inside the stage-1 (setup) run directory; its basename is the run id.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: submit_pipeline.sh [options]

Run from inside the stage-1 (setup) run directory, whose name is the run id
(e.g. cd .../runs/water/N5184/0149 && /path/to/submit_pipeline.sh ...).

Required:
  --setup-input-script FILE        LAMMPS .input template for stage 1 (setup run)
  --setup-starting-structure FILE  Starting structure .data file, symlinked as start.data
  --setup-potential-file FILE      Potential file, symlinked by basename into both stages'
                                    input_files/ (stage 2 uses the same potential as stage 1)
  --analysis-parent-dir DIR        Where <run_id>_distribution_analysis/ is created
  --trajectory-input-script FILE   LAMMPS .input template for stage 2 (trajectory run)

Note: --setup-input-script's in.input must end with exactly one
`write_data final.data` (the standardized stage-1 output name) — this is
what gets symlinked as stage-2's start.data. See OH-therm.input/
b-SiO-therm.input in this directory for examples.

Note: --trajectory-input-script's in.input must write its trajectory dump,
flat (no subdirectory), to dump.lammpstrj (the standardized stage-2 output
name) — this is what the copied dsf.py/rdf_freud.py/bad_freud.py read.

Optional:
  --analysis-template-dir DIR   Distribution-analysis scripts to copy into stage 2
                                 (default: <repo>/analysis/distributions/20260608_GrNrBaSqw)
  --run-dsf {0|1}                Run dsf.py (static/dynamic structure factor) (default: 1)
  --run-rdf {0|1}                Run rdf_freud.py (radial distribution function) (default: 1)
  --run-bad {0|1}                Run bad_freud.py (bond angle distribution) (default: 1)

  --setup-nodes N               --trajectory-nodes N
  --setup-ntasks N               --trajectory-ntasks N
  --setup-time HH:MM:SS           --trajectory-time HH:MM:SS
  --setup-job-name NAME           --trajectory-job-name NAME
  --setup-constraint NAME         --trajectory-constraint NAME
      Per-stage overrides for the corresponding #SBATCH directives in
      jobs/slurm/lammps_submit.slurm. Omit any of these to leave the
      template's own value in effect.

  -h, --help                    Show this help
EOF
}

REPO_ROOT="$HOME/util"
LAMMPS_TEMPLATE="$REPO_ROOT/jobs/slurm/lammps_submit.slurm"
STAGE1_OUTPUT="final.data"
DUMP_FILE="dump.lammpstrj"

SETUP_INPUT_SCRIPT=""
SETUP_STARTING_STRUCTURE=""
SETUP_POTENTIAL_FILE=""
ANALYSIS_PARENT_DIR=""
TRAJECTORY_INPUT_SCRIPT=""
ANALYSIS_TEMPLATE_DIR="$REPO_ROOT/analysis/distributions/20260608_GrNrBaSqw"
RUN_DSF="1"
RUN_RDF="1"
RUN_BAD="1"
SETUP_NODES=""
SETUP_NTASKS=""
SETUP_TIME=""
SETUP_JOB_NAME=""
SETUP_CONSTRAINT=""
TRAJECTORY_NODES=""
TRAJECTORY_NTASKS=""
TRAJECTORY_TIME=""
TRAJECTORY_JOB_NAME=""
TRAJECTORY_CONSTRAINT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --setup-input-script) SETUP_INPUT_SCRIPT="$2"; shift 2 ;;
    --setup-starting-structure) SETUP_STARTING_STRUCTURE="$2"; shift 2 ;;
    --setup-potential-file) SETUP_POTENTIAL_FILE="$2"; shift 2 ;;
    --analysis-parent-dir) ANALYSIS_PARENT_DIR="$2"; shift 2 ;;
    --trajectory-input-script) TRAJECTORY_INPUT_SCRIPT="$2"; shift 2 ;;
    --analysis-template-dir) ANALYSIS_TEMPLATE_DIR="$2"; shift 2 ;;
    --run-dsf) RUN_DSF="$2"; shift 2 ;;
    --run-rdf) RUN_RDF="$2"; shift 2 ;;
    --run-bad) RUN_BAD="$2"; shift 2 ;;
    --setup-nodes) SETUP_NODES="$2"; shift 2 ;;
    --setup-ntasks) SETUP_NTASKS="$2"; shift 2 ;;
    --setup-time) SETUP_TIME="$2"; shift 2 ;;
    --setup-job-name) SETUP_JOB_NAME="$2"; shift 2 ;;
    --setup-constraint) SETUP_CONSTRAINT="$2"; shift 2 ;;
    --trajectory-nodes) TRAJECTORY_NODES="$2"; shift 2 ;;
    --trajectory-ntasks) TRAJECTORY_NTASKS="$2"; shift 2 ;;
    --trajectory-time) TRAJECTORY_TIME="$2"; shift 2 ;;
    --trajectory-job-name) TRAJECTORY_JOB_NAME="$2"; shift 2 ;;
    --trajectory-constraint) TRAJECTORY_CONSTRAINT="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

missing=0
require() {
  if [[ -z "$1" ]]; then
    echo "Missing required argument: $2" >&2
    missing=1
  fi
}
require "$SETUP_INPUT_SCRIPT" --setup-input-script
require "$SETUP_STARTING_STRUCTURE" --setup-starting-structure
require "$SETUP_POTENTIAL_FILE" --setup-potential-file
require "$ANALYSIS_PARENT_DIR" --analysis-parent-dir
require "$TRAJECTORY_INPUT_SCRIPT" --trajectory-input-script
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

STAGE1_DIR="$(pwd)"
RUN_ID="$(basename "$STAGE1_DIR")"

if [[ -e "$STAGE1_DIR/input_files" || -e "$STAGE1_DIR/run" ]]; then
  echo "Error: $STAGE1_DIR already has input_files/ or run/ — refusing to overwrite." >&2
  exit 1
fi

############################
# Stage 1: setup run
############################

mkdir -p "$STAGE1_DIR/input_files" "$STAGE1_DIR/run"

cp "$SETUP_INPUT_SCRIPT" "$STAGE1_DIR/input_files/in.input"
ln -s "$(realpath "$SETUP_STARTING_STRUCTURE")" "$STAGE1_DIR/input_files/start.data"
ln -s "$(realpath "$SETUP_POTENTIAL_FILE")" "$STAGE1_DIR/input_files/$(basename "$SETUP_POTENTIAL_FILE")"

cp "$LAMMPS_TEMPLATE" "$STAGE1_DIR/lammps_submit.slurm"
cp "$LAMMPS_TEMPLATE" "$STAGE1_DIR/run/lammps_submit.slurm"

setup_sbatch_args=()
[[ -n "$SETUP_NODES" ]]      && setup_sbatch_args+=(--nodes="$SETUP_NODES")
[[ -n "$SETUP_NTASKS" ]]     && setup_sbatch_args+=(--ntasks="$SETUP_NTASKS")
[[ -n "$SETUP_TIME" ]]       && setup_sbatch_args+=(--time="$SETUP_TIME")
[[ -n "$SETUP_JOB_NAME" ]]   && setup_sbatch_args+=(--job-name="$SETUP_JOB_NAME")
[[ -n "$SETUP_CONSTRAINT" ]] && setup_sbatch_args+=(--constraint="$SETUP_CONSTRAINT")

echo "Submitting stage 1 (setup) from $STAGE1_DIR/run ..."
JOBID1="$(cd "$STAGE1_DIR/run" && sbatch --parsable \
  "${setup_sbatch_args[@]+"${setup_sbatch_args[@]}"}" lammps_submit.slurm)"
echo "  stage 1 job id: $JOBID1"

############################
# Stage 2: trajectory run
############################

mkdir -p "$ANALYSIS_PARENT_DIR"
ANALYSIS_PARENT_DIR="$(cd "$ANALYSIS_PARENT_DIR" && pwd)"
STAGE2_DIR="$ANALYSIS_PARENT_DIR/${RUN_ID}_distribution_analysis"

if [[ -e "$STAGE2_DIR" ]]; then
  echo "Error: $STAGE2_DIR already exists — refusing to overwrite." >&2
  exit 1
fi

mkdir -p "$STAGE2_DIR/input_files" "$STAGE2_DIR/run"

cp "$TRAJECTORY_INPUT_SCRIPT" "$STAGE2_DIR/input_files/in.input"
ln -s "$STAGE1_DIR/run/$STAGE1_OUTPUT" "$STAGE2_DIR/input_files/start.data"
ln -s "$(realpath "$SETUP_POTENTIAL_FILE")" "$STAGE2_DIR/input_files/$(basename "$SETUP_POTENTIAL_FILE")"

cp "$LAMMPS_TEMPLATE" "$STAGE2_DIR/lammps_submit.slurm"
cp "$LAMMPS_TEMPLATE" "$STAGE2_DIR/run/lammps_submit.slurm"

trajectory_sbatch_args=()
[[ -n "$TRAJECTORY_NODES" ]]      && trajectory_sbatch_args+=(--nodes="$TRAJECTORY_NODES")
[[ -n "$TRAJECTORY_NTASKS" ]]     && trajectory_sbatch_args+=(--ntasks="$TRAJECTORY_NTASKS")
[[ -n "$TRAJECTORY_TIME" ]]       && trajectory_sbatch_args+=(--time="$TRAJECTORY_TIME")
[[ -n "$TRAJECTORY_JOB_NAME" ]]   && trajectory_sbatch_args+=(--job-name="$TRAJECTORY_JOB_NAME")
[[ -n "$TRAJECTORY_CONSTRAINT" ]] && trajectory_sbatch_args+=(--constraint="$TRAJECTORY_CONSTRAINT")

echo "Submitting stage 2 (trajectory run) from $STAGE2_DIR/run, dependent on job $JOBID1 ..."
JOBID2="$(cd "$STAGE2_DIR/run" && sbatch --parsable --dependency=afterok:"$JOBID1" \
  "${trajectory_sbatch_args[@]+"${trajectory_sbatch_args[@]}"}" lammps_submit.slurm)"
echo "  stage 2 job id: $JOBID2"

############################
# Stage 2b: distribution analysis (reads the trajectory as input; never modifies it)
############################

cp "$ANALYSIS_TEMPLATE_DIR/dsf.py" \
   "$ANALYSIS_TEMPLATE_DIR/rdf_freud.py" \
   "$ANALYSIS_TEMPLATE_DIR/bad_freud.py" \
   "$ANALYSIS_TEMPLATE_DIR/dsf_submit.slurm" \
   "$STAGE2_DIR/"
ln -s "$STAGE2_DIR/run/$DUMP_FILE" "$STAGE2_DIR/$DUMP_FILE"

echo "Submitting distribution analysis from $STAGE2_DIR, dependent on job $JOBID2 ..."
JOBID3="$(cd "$STAGE2_DIR" && sbatch --parsable --dependency=afterok:"$JOBID2" \
  --export=ALL,TRAJ="$DUMP_FILE",RUN_DSF="$RUN_DSF",RUN_RDF="$RUN_RDF",RUN_BAD="$RUN_BAD" \
  dsf_submit.slurm)"
echo "  distribution-analysis job id: $JOBID3"

enabled_scripts=""
[[ "$RUN_DSF" == "1" ]] && enabled_scripts+="dsf.py "
[[ "$RUN_RDF" == "1" ]] && enabled_scripts+="rdf_freud.py "
[[ "$RUN_BAD" == "1" ]] && enabled_scripts+="bad_freud.py "
enabled_scripts="${enabled_scripts% }"

cat <<SUMMARY

Pipeline submitted:
  stage 1 (setup)                 : $STAGE1_DIR  (job $JOBID1)
  stage 2 (trajectory + analysis) : $STAGE2_DIR  (job $JOBID2, analysis job $JOBID3)
  analysis scripts enabled        : ${enabled_scripts:-none}

$(if [[ -n "$enabled_scripts" ]]; then
  echo "Reminder: edit $STAGE2_DIR/{${enabled_scripts// /,}}'s physics config"
  echo "(ELEMENTS, R_CUTOFF, DT, WINDOW_SIZE, ...) before job $JOBID3 actually runs —"
  echo "the dependency only guarantees it won't start early."
fi)
SUMMARY
