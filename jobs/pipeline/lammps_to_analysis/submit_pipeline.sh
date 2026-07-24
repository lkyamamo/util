#!/bin/bash
# submit_pipeline.sh — run a LAMMPS setup+trajectory job, then distribution
# analysis (rdf/bad/dsf) on its dump, via sbatch --dependency.
#
# Run from inside the LAMMPS run directory; its basename is the run id.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: submit_pipeline.sh [options]

Run from inside the LAMMPS run directory, whose name is the run id
(e.g. cd .../runs/water/N5184/0149 && /path/to/submit_pipeline.sh ...).

Required:
  --input-script FILE          LAMMPS .input script (setup + trajectory generation)
  --starting-structure FILE    Starting structure .data file, symlinked as start.data
  --potential-file FILE        Potential file, symlinked by basename
  --analysis-parent-dir DIR    Where <run_id>_distribution_analysis/ is created

Note: --input-script's in.input must write its trajectory dump, flat (no
subdirectory), to dump.lammpstrj (the standardized output name) — this is
what the copied dsf.py/rdf_freud.py/bad_freud.py read. See OH-therm.input/
b-SiO-therm.input in this directory for examples.

Optional:
  --analysis-template-dir DIR   Distribution-analysis scripts to copy into stage 2
                                 (default: <repo>/analysis/distributions/20260608_GrNrBaSqw)
  --run-dsf {0|1}                Run dsf.py (static/dynamic structure factor) (default: 1)
  --run-rdf {0|1}                Run rdf_freud.py (radial distribution function) (default: 1)
  --run-bad {0|1}                Run bad_freud.py (bond angle distribution) (default: 1)

  --nodes N               --analysis-nodes N
  --ntasks N               --analysis-ntasks N
  --time HH:MM:SS           --analysis-time HH:MM:SS
  --job-name NAME           --analysis-job-name NAME
  --constraint NAME         --analysis-constraint NAME
      Overrides for the corresponding #SBATCH directives — the left column
      overrides jobs/slurm/lammps_submit.slurm (the LAMMPS run), the right
      column overrides dsf_submit.slurm (the analysis job). Omit any of
      these to leave the template's own value in effect.

  -h, --help                    Show this help
EOF
}

REPO_ROOT="$HOME/util"
LAMMPS_TEMPLATE="$REPO_ROOT/jobs/slurm/lammps_submit.slurm"
DUMP_FILE="dump.lammpstrj"

INPUT_SCRIPT=""
STARTING_STRUCTURE=""
POTENTIAL_FILE=""
ANALYSIS_PARENT_DIR=""
ANALYSIS_TEMPLATE_DIR="$REPO_ROOT/analysis/distributions/20260608_GrNrBaSqw"
RUN_DSF="1"
RUN_RDF="1"
RUN_BAD="1"
NODES=""
NTASKS=""
TIME=""
JOB_NAME=""
CONSTRAINT=""
ANALYSIS_NODES=""
ANALYSIS_NTASKS=""
ANALYSIS_TIME=""
ANALYSIS_JOB_NAME=""
ANALYSIS_CONSTRAINT=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-script) INPUT_SCRIPT="$2"; shift 2 ;;
    --starting-structure) STARTING_STRUCTURE="$2"; shift 2 ;;
    --potential-file) POTENTIAL_FILE="$2"; shift 2 ;;
    --analysis-parent-dir) ANALYSIS_PARENT_DIR="$2"; shift 2 ;;
    --analysis-template-dir) ANALYSIS_TEMPLATE_DIR="$2"; shift 2 ;;
    --run-dsf) RUN_DSF="$2"; shift 2 ;;
    --run-rdf) RUN_RDF="$2"; shift 2 ;;
    --run-bad) RUN_BAD="$2"; shift 2 ;;
    --nodes) NODES="$2"; shift 2 ;;
    --ntasks) NTASKS="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --job-name) JOB_NAME="$2"; shift 2 ;;
    --constraint) CONSTRAINT="$2"; shift 2 ;;
    --analysis-nodes) ANALYSIS_NODES="$2"; shift 2 ;;
    --analysis-ntasks) ANALYSIS_NTASKS="$2"; shift 2 ;;
    --analysis-time) ANALYSIS_TIME="$2"; shift 2 ;;
    --analysis-job-name) ANALYSIS_JOB_NAME="$2"; shift 2 ;;
    --analysis-constraint) ANALYSIS_CONSTRAINT="$2"; shift 2 ;;
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
require "$INPUT_SCRIPT" --input-script
require "$STARTING_STRUCTURE" --starting-structure
require "$POTENTIAL_FILE" --potential-file
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

STAGE1_DIR="$(pwd)"
RUN_ID="$(basename "$STAGE1_DIR")"

if [[ -e "$STAGE1_DIR/input_files" || -e "$STAGE1_DIR/run" ]]; then
  echo "Error: $STAGE1_DIR already has input_files/ or run/ — refusing to overwrite." >&2
  exit 1
fi

############################
# Stage 1: LAMMPS run (setup + trajectory generation)
############################

mkdir -p "$STAGE1_DIR/input_files" "$STAGE1_DIR/run"

cp "$INPUT_SCRIPT" "$STAGE1_DIR/input_files/in.input"
ln -s "$(realpath "$STARTING_STRUCTURE")" "$STAGE1_DIR/input_files/start.data"
ln -s "$(realpath "$POTENTIAL_FILE")" "$STAGE1_DIR/input_files/$(basename "$POTENTIAL_FILE")"

cp "$LAMMPS_TEMPLATE" "$STAGE1_DIR/lammps_submit.slurm"
cp "$LAMMPS_TEMPLATE" "$STAGE1_DIR/run/lammps_submit.slurm"

sbatch_args=()
[[ -n "$NODES" ]]      && sbatch_args+=(--nodes="$NODES")
[[ -n "$NTASKS" ]]     && sbatch_args+=(--ntasks="$NTASKS")
[[ -n "$TIME" ]]       && sbatch_args+=(--time="$TIME")
[[ -n "$JOB_NAME" ]]   && sbatch_args+=(--job-name="$JOB_NAME")
[[ -n "$CONSTRAINT" ]] && sbatch_args+=(--constraint="$CONSTRAINT")

echo "Submitting LAMMPS run from $STAGE1_DIR/run ..."
JOBID1="$(cd "$STAGE1_DIR/run" && sbatch --parsable \
  "${sbatch_args[@]+"${sbatch_args[@]}"}" lammps_submit.slurm)"
echo "  LAMMPS job id: $JOBID1"

############################
# Stage 2: distribution analysis (reads the trajectory as input; never modifies it)
############################

mkdir -p "$ANALYSIS_PARENT_DIR"
ANALYSIS_PARENT_DIR="$(cd "$ANALYSIS_PARENT_DIR" && pwd)"
STAGE2_DIR="$ANALYSIS_PARENT_DIR/${RUN_ID}_distribution_analysis"

if [[ -e "$STAGE2_DIR" ]]; then
  echo "Error: $STAGE2_DIR already exists — refusing to overwrite." >&2
  exit 1
fi

mkdir -p "$STAGE2_DIR"

cp "$ANALYSIS_TEMPLATE_DIR/dsf.py" \
   "$ANALYSIS_TEMPLATE_DIR/rdf_freud.py" \
   "$ANALYSIS_TEMPLATE_DIR/bad_freud.py" \
   "$ANALYSIS_TEMPLATE_DIR/dsf_submit.slurm" \
   "$STAGE2_DIR/"
ln -s "$STAGE1_DIR/run/$DUMP_FILE" "$STAGE2_DIR/$DUMP_FILE"

analysis_sbatch_args=()
[[ -n "$ANALYSIS_NODES" ]]      && analysis_sbatch_args+=(--nodes="$ANALYSIS_NODES")
[[ -n "$ANALYSIS_NTASKS" ]]     && analysis_sbatch_args+=(--ntasks="$ANALYSIS_NTASKS")
[[ -n "$ANALYSIS_TIME" ]]       && analysis_sbatch_args+=(--time="$ANALYSIS_TIME")
[[ -n "$ANALYSIS_JOB_NAME" ]]   && analysis_sbatch_args+=(--job-name="$ANALYSIS_JOB_NAME")
[[ -n "$ANALYSIS_CONSTRAINT" ]] && analysis_sbatch_args+=(--constraint="$ANALYSIS_CONSTRAINT")

echo "Submitting distribution analysis from $STAGE2_DIR, dependent on job $JOBID1 ..."
JOBID2="$(cd "$STAGE2_DIR" && sbatch --parsable --dependency=afterok:"$JOBID1" \
  "${analysis_sbatch_args[@]+"${analysis_sbatch_args[@]}"}" \
  --export=ALL,TRAJ="$DUMP_FILE",RUN_DSF="$RUN_DSF",RUN_RDF="$RUN_RDF",RUN_BAD="$RUN_BAD" \
  dsf_submit.slurm)"
echo "  distribution-analysis job id: $JOBID2"

enabled_scripts=""
[[ "$RUN_DSF" == "1" ]] && enabled_scripts+="dsf.py "
[[ "$RUN_RDF" == "1" ]] && enabled_scripts+="rdf_freud.py "
[[ "$RUN_BAD" == "1" ]] && enabled_scripts+="bad_freud.py "
enabled_scripts="${enabled_scripts% }"

cat <<SUMMARY

Pipeline submitted:
  LAMMPS run (setup + trajectory) : $STAGE1_DIR  (job $JOBID1)
  distribution analysis           : $STAGE2_DIR  (job $JOBID2)
  analysis scripts enabled        : ${enabled_scripts:-none}

$(if [[ -n "$enabled_scripts" ]]; then
  echo "Reminder: edit $STAGE2_DIR/{${enabled_scripts// /,}}'s physics config"
  echo "(ELEMENTS, R_CUTOFF, DT, WINDOW_SIZE, ...) before job $JOBID2 actually runs —"
  echo "the dependency only guarantees it won't start early."
fi)
SUMMARY
