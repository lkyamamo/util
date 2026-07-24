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
  None of these use commas (sbatch --export is comma-delimited and silently
  truncates any value containing one), so no quoting/encoding is needed
  beyond normal shell quoting of the whole flag value.

  --nodes N               --analysis-nodes N
  --ntasks N               --analysis-ntasks N
  --time HH:MM:SS           --analysis-time HH:MM:SS
  --job-name NAME           --analysis-job-name NAME
  --constraint NAME         --analysis-constraint NAME
                           --analysis-cpus-per-task N
      Overrides for the corresponding #SBATCH directives — the left column
      overrides jobs/slurm/lammps_submit.slurm (the LAMMPS run), the right
      column overrides dsf_submit.slurm (the analysis job). Omit any of
      these to leave the template's own value in effect. --analysis-cpus-per-task
      has no LAMMPS-side counterpart; it also drives OMP_NUM_THREADS for
      dsf.py/rdf_freud.py/bad_freud.py (dsf_submit.slurm sets
      OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK).

  --force REASON                  Overwrite existing input_files/run/ (stage 1) or an
                                  existing <run_id>_distribution_analysis/ (stage 2)
                                  instead of refusing to run. REASON is required (a
                                  short explanation of why you're overwriting) and is
                                  logged, with a timestamp and the exact path removed,
                                  to both stderr and
                                  jobs/pipeline/lammps_to_analysis/overwrite.log.

  -h, --help                    Show this help
EOF
}

REPO_ROOT="$HOME/util"
LAMMPS_TEMPLATE="$REPO_ROOT/jobs/slurm/lammps_submit.slurm"
DUMP_FILE="dump.lammpstrj"
FORCE="0"
FORCE_REASON=""
LOG_FILE="$REPO_ROOT/jobs/pipeline/lammps_to_analysis/overwrite.log"

log_overwrite() {
  local msg
  msg="[$(date "+%Y-%m-%dT%H:%M:%S%z")] $1"
  echo "WARNING: $msg" >&2
  echo "$msg" >> "$LOG_FILE"
}

# general input parameters
INPUT_SCRIPT="/scratch1/lkyamamo/test-runs/potential-verify/runs/0168/OH-therm.input"
STARTING_STRUCTURE="/home1/lkyamamo/util/starting-structures/ICE_CUBIC.data"
POTENTIAL_FILE="/scratch1/lkyamamo/test-runs/potential-verify/potentials/20260723_OH.vashishta"
ANALYSIS_PARENT_DIR="/scratch1/lkyamamo/test-runs/potential-verify/analysis"
ANALYSIS_TEMPLATE_DIR="$REPO_ROOT/analysis/distributions/20260608_GrNrBaSqw"
RUN_DSF="0"
RUN_RDF="1"
RUN_BAD="1"

# trajectory creation slurm parameters
NODES="1"
NTASKS="128"
TIME="1:00:00"
JOB_NAME="water-setup"
CONSTRAINT="epyc-9554"

# analysis slurm parameters
ANALYSIS_NODES="1"
ANALYSIS_NTASKS="1"
ANALYSIS_TIME="1:00:00"
ANALYSIS_JOB_NAME="water-distributions"
ANALYSIS_CONSTRAINT="epyc-9554"
ANALYSIS_CPUS_PER_TASK=""

# distribution code parameters
RDF_R_MAX="8"
RDF_BINS_VAL="800"
BAD_ELEMENTS="O;H"
BAD_R_CUTOFF=""
BAD_R_MINCUT=""
BAD_TRIPLET_CUTOFFS="H-O-H:H-O-H:1.2:0.5:1.2:0.5|H--O--H:H-O-H:2.4:1.49:2.4:1.49|O-H--O:O-H-O:1.2:0.5:2.4:1.49"
BAD_BINS_VAL="360"
DSF_DT=""
DSF_N_FRAMES=""
DSF_STRIDE=""
DSF_WINDOW_SIZE=""
DSF_Q_MAX=""
DSF_N_Q_BINS=""


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
    --force) FORCE="1"; FORCE_REASON="$2"; shift 2 ;;
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

if [[ "$FORCE" == "1" && -z "$FORCE_REASON" ]]; then
  echo 'Error: --force requires a non-empty reason, e.g. --force "re-running after fixing potential file"' >&2
  exit 1
fi

STAGE1_DIR="$(pwd)"
RUN_ID="$(basename "$STAGE1_DIR")"

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
   "$ANALYSIS_TEMPLATE_DIR/dsf_submit.slurm" \
   "$STAGE2_DIR/"
ln -s "$STAGE1_DIR/run/$DUMP_FILE" "$STAGE2_DIR/$DUMP_FILE"

analysis_sbatch_args=()
[[ -n "$ANALYSIS_NODES" ]]      && analysis_sbatch_args+=(--nodes="$ANALYSIS_NODES")
[[ -n "$ANALYSIS_NTASKS" ]]     && analysis_sbatch_args+=(--ntasks="$ANALYSIS_NTASKS")
[[ -n "$ANALYSIS_TIME" ]]       && analysis_sbatch_args+=(--time="$ANALYSIS_TIME")
[[ -n "$ANALYSIS_JOB_NAME" ]]   && analysis_sbatch_args+=(--job-name="$ANALYSIS_JOB_NAME")
[[ -n "$ANALYSIS_CONSTRAINT" ]] && analysis_sbatch_args+=(--constraint="$ANALYSIS_CONSTRAINT")
[[ -n "$ANALYSIS_CPUS_PER_TASK" ]] && analysis_sbatch_args+=(--cpus-per-task="$ANALYSIS_CPUS_PER_TASK")

# None of these values may contain a comma — sbatch --export is
# comma-delimited and silently truncates anything after an embedded one.
export_vars="ALL,TRAJ=$DUMP_FILE,RUN_DSF=$RUN_DSF,RUN_RDF=$RUN_RDF,RUN_BAD=$RUN_BAD"
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

echo "Submitting distribution analysis from $STAGE2_DIR, dependent on job $JOBID1 ..."
JOBID2="$(cd "$STAGE2_DIR" && sbatch --parsable --dependency=afterok:"$JOBID1" \
  "${analysis_sbatch_args[@]+"${analysis_sbatch_args[@]}"}" \
  --export="$export_vars" \
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
  echo "Reminder: any physics config (ELEMENTS, R_CUTOFF, DT, WINDOW_SIZE, ...) not"
  echo "passed via --rdf-*/--bad-*/--dsf-* flags is using each script's own default —"
  echo "check $STAGE2_DIR/{${enabled_scripts// /,}} if that's not what you want."
fi)
SUMMARY
