#!/bin/bash
# submit_dielectric_pipeline.sh — run a LAMMPS dielectric-production job (one
# per temperature), then the dielectric MPI calc on its chunked dumps, via
# sbatch --dependency, then aggregate all temperatures into one CSV.
#
# Run from inside a sweep directory; its basename is the sweep id. A sweep
# of a single temperature works the same way as a multi-temperature sweep.
set -euo pipefail

usage() {
  cat <<'EOF'
Usage: submit_dielectric_pipeline.sh [options]

Run from inside a sweep directory, whose name is the sweep id (e.g. cd
.../dielectric-sweeps/water-OH-0001 && /path/to/submit_dielectric_pipeline.sh ...).
Creates one T<value>/ subdirectory per requested temperature.

Requires submit_dielectric_pipeline.conf next to this script (same
directory), which sets defaults for every parameter below — the script
refuses to run without it. Copy submit_dielectric_pipeline.conf.example to
submit_dielectric_pipeline.conf and edit it for your setup;
submit_dielectric_pipeline.conf itself is gitignored since it holds your own
local paths. Any flag below overrides its config value for that one
invocation.

Required:
  --input-script FILE          LAMMPS .input template referencing ${TARGET_TEMP}
                                (and ${N_TIMES}/${NVT_LENGTH}/${DUMP_EVERY}), e.g.
                                dielectric-therm.input in this directory. Assumes
                                --starting-structure is already equilibrated at
                                each requested temperature (no deform/anneal ramp).
  --starting-structure FILE    Starting structure .data file, symlinked as start.data
  --potential-file FILE        Potential file, symlinked into input_files/
  --analysis-parent-dir DIR    Where <sweep_id>_T<value>_dielectric_calc/ dirs are created
  --temperatures LIST          SEMICOLON-separated list of target temperatures
                                (K), e.g. "303;323;343". A single value (e.g.
                                "303") runs just that one temperature.
  --la FLOAT                   Box length a (Angstrom) — feeds 2.dipole_std.py's
                                dielectric-constant prefactor
  --lb FLOAT                   Box length b (Angstrom)
  --lc FLOAT                   Box length c (Angstrom)

Note: --potential-link-name NAME (optional but often required in practice)
sets the symlink's name in input_files/ (e.g. "OH.usc"). It must match
whatever your in.input's pair_coeff line expects. Defaults to
--potential-file's own basename if omitted.

Optional:
  --analysis-template-dir DIR   Dielectric-calc scripts to copy into stage 2
                                 (default: <repo>/simulation/lammps/20260617_dielectric_multi_traj)

  Dielectric production (stage 1) config — LAMMPS -var overrides:
    --dielectric-n-chunks INT       Number of chunked dielectric.N.custom
                                    files (LAMMPS ${N_TIMES}), e.g. 128
                                    (default: 128)
    --dielectric-chunk-length INT   NVT steps per chunk (LAMMPS ${NVT_LENGTH}),
                                    e.g. 468750 (default: 468750)
    --dielectric-dump-every INT     Steps between dumped frames within a chunk
                                    (LAMMPS ${DUMP_EVERY}), e.g. 10 (default: 10)

  Dielectric calc (stage 2) config — same meaning as 1.calc_mpi.py/2.dipole_std.py's
  own flags/args in simulation/lammps/20260617_dielectric_multi_traj/:
    --dielectric-cutoff FLOAT        O-H bond cutoff (Angstrom) (default: 1.2)
    --dielectric-type-o INT          LAMMPS atom type for oxygen (default: 1)
    --dielectric-type-h INT          LAMMPS atom type for hydrogen (default: 2)
    --dielectric-averaging-method STR  windowed|cumulative|hybrid|binned
                                        (default: cumulative)

  --nodes N               --analysis-nodes N
  --ntasks N                --analysis-time HH:MM:SS
  --time HH:MM:SS           --analysis-job-name NAME
  --job-name NAME           --analysis-constraint NAME
  --constraint NAME         --analysis-nodelist NODES
  --nodelist NODES
      Overrides for the corresponding #SBATCH directives — the left column
      overrides jobs/slurm/lammps_dielectric_submit.slurm (stage 1, per
      temperature), the right column overrides dielectric_submit.slurm
      (stage 2, per temperature). Omit any of these to leave the template's
      own value in effect. Stage 2's rank count is always auto-derived from
      the number of dielectric.N.custom files found (must equal
      --dielectric-n-chunks), so it has no --analysis-ntasks override.

  --skip-trajectory                Skip stage 1 (LAMMPS dielectric production) for
                                  every requested temperature and assume each
                                  T<value>/run/dumps/ already exists and is
                                  populated (i.e. you already ran this pipeline, or
                                  LAMMPS directly, from here). Verified per
                                  temperature before stage 2 runs; the pipeline
                                  aborts with an error if a T<value>/run/dumps/ is
                                  missing or empty. --input-script/
                                  --starting-structure/--potential-file, the
                                  dielectric-production (stage 1) config, and the
                                  trajectory-creation slurm overrides (--nodes/
                                  --ntasks/--time/--job-name/--constraint/
                                  --nodelist) are ignored in this mode, and stage 2
                                  is submitted without a --dependency on stage 1
                                  (nothing to wait on).

  --force REASON                  Overwrite an existing T<value>/ (stage 1) or
                                  an existing <sweep_id>_T<value>_dielectric_calc/
                                  (stage 2) instead of refusing to run. REASON is
                                  required (a short explanation of why you're
                                  overwriting) and is logged, with a timestamp
                                  and the exact path removed, to both stderr and
                                  jobs/pipeline/lammps_to_dielectric/overwrite.log.

  --interactive                   Run all stages in the foreground via plain
                                  bash instead of submitting them with sbatch.
                                  Assumes you already have an interactive
                                  allocation (e.g. via salloc) — this does not
                                  request one. --nodes/--time/--job-name/
                                  --constraint/--nodelist (and their
                                  --analysis-* counterparts) are ignored in
                                  this mode, since no new job is submitted.
                                  --ntasks still applies (default 64 if not
                                  given) since it drives srun -n for stage 1.
                                  Temperatures run sequentially, one after
                                  another.

  -h, --help                    Show this help

All terminal output from this script (not the SLURM jobs themselves) is
also appended to submit_dielectric_pipeline.log in the sweep directory (cwd)
each time it's invoked.
EOF
}

INTERACTIVE="0"
SKIP_TRAJECTORY="0"
REPO_ROOT="$HOME/util"
LAMMPS_TEMPLATE="$REPO_ROOT/jobs/slurm/lammps_dielectric_submit.slurm"
FORCE="0"
FORCE_REASON=""
LOG_FILE="$REPO_ROOT/jobs/pipeline/lammps_to_dielectric/overwrite.log"
AGGREGATE_SCRIPT="$REPO_ROOT/jobs/pipeline/lammps_to_dielectric/aggregate_dielectric_vs_temperature.py"
AGGREGATE_TEMPLATE="$REPO_ROOT/jobs/pipeline/lammps_to_dielectric/aggregate_submit.slurm"

log_overwrite() {
  local msg
  msg="[$(date "+%Y-%m-%dT%H:%M:%S%z")] $1"
  echo "WARNING: $msg" >&2
  echo "$msg" >> "$LOG_FILE"
}

# All parameter defaults (general input, dielectric-production config,
# dielectric-calc config, trajectory-creation slurm, analysis slurm) live in
# submit_dielectric_pipeline.conf next to this script, not in the script
# itself, so personal paths never need to be edited into (or committed from)
# tracked code. See submit_dielectric_pipeline.conf.example.
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG_FILE="$SCRIPT_DIR/submit_dielectric_pipeline.conf"
if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "Error: required config file not found: $CONFIG_FILE" >&2
  echo "Copy $SCRIPT_DIR/submit_dielectric_pipeline.conf.example to $CONFIG_FILE and edit it for your setup." >&2
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
    --temperatures) TEMPERATURES="$2"; shift 2 ;;
    --la) LA="$2"; shift 2 ;;
    --lb) LB="$2"; shift 2 ;;
    --lc) LC="$2"; shift 2 ;;
    --dielectric-n-chunks) DIELECTRIC_N_CHUNKS="$2"; shift 2 ;;
    --dielectric-chunk-length) DIELECTRIC_CHUNK_LENGTH="$2"; shift 2 ;;
    --dielectric-dump-every) DIELECTRIC_DUMP_EVERY="$2"; shift 2 ;;
    --dielectric-cutoff) DIELECTRIC_CUTOFF="$2"; shift 2 ;;
    --dielectric-type-o) DIELECTRIC_TYPE_O="$2"; shift 2 ;;
    --dielectric-type-h) DIELECTRIC_TYPE_H="$2"; shift 2 ;;
    --dielectric-averaging-method) DIELECTRIC_AVERAGING_METHOD="$2"; shift 2 ;;
    --nodes) NODES="$2"; shift 2 ;;
    --ntasks) NTASKS="$2"; shift 2 ;;
    --time) TIME="$2"; shift 2 ;;
    --job-name) JOB_NAME="$2"; shift 2 ;;
    --constraint) CONSTRAINT="$2"; shift 2 ;;
    --nodelist) NODELIST="$2"; shift 2 ;;
    --analysis-nodes) ANALYSIS_NODES="$2"; shift 2 ;;
    --analysis-time) ANALYSIS_TIME="$2"; shift 2 ;;
    --analysis-job-name) ANALYSIS_JOB_NAME="$2"; shift 2 ;;
    --analysis-constraint) ANALYSIS_CONSTRAINT="$2"; shift 2 ;;
    --analysis-nodelist) ANALYSIS_NODELIST="$2"; shift 2 ;;
    --force) FORCE="1"; FORCE_REASON="$2"; shift 2 ;;
    --interactive) INTERACTIVE="1"; shift 1 ;;
    --skip-trajectory) SKIP_TRAJECTORY="1"; shift 1 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage; exit 1 ;;
  esac
done

SWEEP_DIR="$(pwd)"
SWEEP_ID="$(basename "$SWEEP_DIR")"
PIPELINE_LOG="$SWEEP_DIR/submit_dielectric_pipeline.log"
exec > >(tee -a "$PIPELINE_LOG") 2>&1
echo "=== submit_dielectric_pipeline.sh started $(date "+%Y-%m-%dT%H:%M:%S%z") — sweep id: $SWEEP_ID ==="
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
require "$TEMPERATURES" --temperatures
require "$LA" --la
require "$LB" --lb
require "$LC" --lc
if [[ "$missing" -ne 0 ]]; then
  usage
  exit 1
fi

if [[ "$FORCE" == "1" && -z "$FORCE_REASON" ]]; then
  echo 'Error: --force requires a non-empty reason, e.g. --force "re-running after fixing potential file"' >&2
  exit 1
fi

if [[ "$INTERACTIVE" == "1" && -z "${SLURM_JOB_ID:-}" ]]; then
  echo "Error: --interactive requires an active Slurm allocation (run 'salloc ...' first)." >&2
  echo "No \$SLURM_JOB_ID is set in this shell, so srun would fail to allocate resources —" >&2
  echo "and lammps_dielectric_submit.slurm's unconditional 'exit 0' would silently mask that failure" >&2
  echo "instead of stopping the pipeline before stage 2." >&2
  exit 1
fi

mkdir -p "$ANALYSIS_PARENT_DIR"
ANALYSIS_PARENT_DIR="$(cd "$ANALYSIS_PARENT_DIR" && pwd)"

IFS=';' read -ra TEMP_LIST <<< "$TEMPERATURES"
if [[ "${#TEMP_LIST[@]}" -eq 0 ]]; then
  echo "Error: --temperatures produced no entries: '$TEMPERATURES'" >&2
  exit 1
fi

############################
# Per-temperature stage 1 + stage 2
############################

STAGE1_JOBIDS=()
STAGE2_JOBIDS=()
ENTRIES_FILE="$(mktemp)"
trap 'rm -f "$ENTRIES_FILE"' EXIT

for T in "${TEMP_LIST[@]}"; do
  [[ -z "$T" ]] && continue
  echo ""
  echo "--- Temperature $T ---"

  STAGE1_DIR="$SWEEP_DIR/T${T}"
  STAGE2_DIR="$ANALYSIS_PARENT_DIR/${SWEEP_ID}_T${T}_dielectric_calc"

  if [[ "$SKIP_TRAJECTORY" == "1" ]]; then
    if [[ ! -d "$STAGE1_DIR/run/dumps" ]] || [[ -z "$(ls -A "$STAGE1_DIR/run/dumps" 2>/dev/null)" ]]; then
      echo "Error: --skip-trajectory was given but $STAGE1_DIR/run/dumps does not exist or is empty (T=$T)." >&2
      echo "Re-run without --skip-trajectory to generate the trajectory, or fix the path above." >&2
      exit 1
    fi
    echo "--skip-trajectory: found existing dumps in $STAGE1_DIR/run/dumps — skipping stage 1 for T=$T."
  else
    if [[ -e "$STAGE1_DIR/input_files" || -e "$STAGE1_DIR/run" ]]; then
      if [[ "$FORCE" == "1" ]]; then
        log_overwrite "--force ($FORCE_REASON): removing existing $STAGE1_DIR/input_files and/or $STAGE1_DIR/run (sweep id: $SWEEP_ID, T=$T)"
        rm -rf "$STAGE1_DIR/input_files" "$STAGE1_DIR/run"
      else
        echo "Error: $STAGE1_DIR already has input_files/ or run/ — refusing to overwrite (use --force to override)." >&2
        exit 1
      fi
    fi
  fi

  if [[ -e "$STAGE2_DIR" ]]; then
    if [[ "$FORCE" == "1" ]]; then
      log_overwrite "--force ($FORCE_REASON): removing existing $STAGE2_DIR"
      rm -rf "$STAGE2_DIR"
    else
      echo "Error: $STAGE2_DIR already exists — refusing to overwrite (use --force to override)." >&2
      exit 1
    fi
  fi

  if [[ "$SKIP_TRAJECTORY" != "1" ]]; then
    ##########
    # Stage 1: LAMMPS dielectric production
    ##########

    mkdir -p "$STAGE1_DIR/input_files" "$STAGE1_DIR/run"

    cp "$INPUT_SCRIPT" "$STAGE1_DIR/input_files/in.input"
    ln -s "$(realpath "$STARTING_STRUCTURE")" "$STAGE1_DIR/input_files/start.data"
    ln -s "$(realpath "$POTENTIAL_FILE")" "$STAGE1_DIR/input_files/${POTENTIAL_LINK_NAME:-$(basename "$POTENTIAL_FILE")}"

    cp "$LAMMPS_TEMPLATE" "$STAGE1_DIR/lammps_submit.slurm"
    cp "$LAMMPS_TEMPLATE" "$STAGE1_DIR/run/lammps_submit.slurm"
  fi

  ##########
  # Stage 2: dielectric MPI calc (reads stage 1's dumps; never modifies them)
  ##########

  mkdir -p "$STAGE2_DIR/logs"
  cp "$ANALYSIS_TEMPLATE_DIR/1.calc_mpi.py" \
     "$ANALYSIS_TEMPLATE_DIR/2.dipole_std.py" \
     "$ANALYSIS_TEMPLATE_DIR/dielectric_submit.slurm" \
     "$STAGE2_DIR/"
  ln -s "$STAGE1_DIR/run/dumps" "$STAGE2_DIR/dumps"

  stage1_export="TARGET_TEMP=$T,N_TIMES=$DIELECTRIC_N_CHUNKS,NVT_LENGTH=$DIELECTRIC_CHUNK_LENGTH,DUMP_EVERY=$DIELECTRIC_DUMP_EVERY"
  stage2_export="ALL,DUMP_DIR=$STAGE2_DIR/dumps,DUMP_EVERY=$DIELECTRIC_DUMP_EVERY,CUTOFF=$DIELECTRIC_CUTOFF,TYPE_O=$DIELECTRIC_TYPE_O,TYPE_H=$DIELECTRIC_TYPE_H,TEMPERATURE=$T,LA=$LA,LB=$LB,LC=$LC,AVERAGING_METHOD=$DIELECTRIC_AVERAGING_METHOD"

  if [[ "$INTERACTIVE" == "1" ]]; then
    if [[ "$SKIP_TRAJECTORY" == "1" ]]; then
      echo "--skip-trajectory: not running LAMMPS for T=$T, using existing dumps in $STAGE1_DIR/run/dumps."
    else
      echo "Running LAMMPS dielectric production interactively in $STAGE1_DIR/run ..."
      (
        export SLURM_SUBMIT_DIR="$STAGE1_DIR/run"
        export SLURM_NTASKS="${NTASKS:-64}"
        export TARGET_TEMP="$T" N_TIMES="$DIELECTRIC_N_CHUNKS" NVT_LENGTH="$DIELECTRIC_CHUNK_LENGTH" DUMP_EVERY="$DIELECTRIC_DUMP_EVERY"
        cd "$STAGE1_DIR/run" && bash lammps_submit.slurm
      )
      echo "LAMMPS dielectric production finished for T=$T."
    fi

    echo "Running dielectric calc interactively in $STAGE2_DIR ..."
    (
      export SLURM_SUBMIT_DIR="$STAGE2_DIR"
      IFS=',' read -ra _kv_pairs <<< "${stage2_export#ALL,}"
      for pair in "${_kv_pairs[@]}"; do export "$pair"; done
      cd "$STAGE2_DIR" && bash dielectric_submit.slurm
    )
    echo "Dielectric calc finished for T=$T."
  else
    JOBID1=""
    if [[ "$SKIP_TRAJECTORY" == "1" ]]; then
      echo "--skip-trajectory: not submitting a LAMMPS job for T=$T, using existing dumps in $STAGE1_DIR/run/dumps."
    else
      sbatch_args=()
      [[ -n "$NODES" ]]      && sbatch_args+=(--nodes="$NODES")
      [[ -n "$NTASKS" ]]     && sbatch_args+=(--ntasks="$NTASKS")
      [[ -n "$TIME" ]]       && sbatch_args+=(--time="$TIME")
      [[ -n "$JOB_NAME" ]]   && sbatch_args+=(--job-name="${JOB_NAME}-T${T}")
      [[ -n "$CONSTRAINT" ]] && sbatch_args+=(--constraint="$CONSTRAINT")
      [[ -n "$NODELIST" ]]   && sbatch_args+=(--nodelist="$NODELIST")

      echo "Submitting LAMMPS dielectric production from $STAGE1_DIR/run ..."
      JOBID1="$(cd "$STAGE1_DIR/run" && sbatch --parsable \
        "${sbatch_args[@]+"${sbatch_args[@]}"}" \
        --export="ALL,$stage1_export" \
        lammps_submit.slurm)"
      echo "  LAMMPS job id: $JOBID1"
      STAGE1_JOBIDS+=("$JOBID1")
    fi

    analysis_sbatch_args=()
    [[ -n "$ANALYSIS_NODES" ]]      && analysis_sbatch_args+=(--nodes="$ANALYSIS_NODES")
    [[ -n "$ANALYSIS_TIME" ]]       && analysis_sbatch_args+=(--time="$ANALYSIS_TIME")
    [[ -n "$ANALYSIS_JOB_NAME" ]]   && analysis_sbatch_args+=(--job-name="${ANALYSIS_JOB_NAME}-T${T}")
    [[ -n "$ANALYSIS_CONSTRAINT" ]] && analysis_sbatch_args+=(--constraint="$ANALYSIS_CONSTRAINT")
    [[ -n "$ANALYSIS_NODELIST" ]]   && analysis_sbatch_args+=(--nodelist="$ANALYSIS_NODELIST")

    dependency_args=()
    [[ -n "$JOBID1" ]] && dependency_args+=(--dependency=afterok:"$JOBID1")

    echo "Submitting dielectric calc from $STAGE2_DIR${JOBID1:+, dependent on job $JOBID1} ..."
    JOBID2="$(cd "$STAGE2_DIR" && sbatch --parsable \
      "${dependency_args[@]+"${dependency_args[@]}"}" \
      "${analysis_sbatch_args[@]+"${analysis_sbatch_args[@]}"}" \
      --export="$stage2_export" \
      dielectric_submit.slurm)"
    echo "  dielectric-calc job id: $JOBID2"
    STAGE2_JOBIDS+=("$JOBID2")
  fi

  echo "${T}:${STAGE2_DIR}/dipole_output/summary.txt" >> "$ENTRIES_FILE"
done

############################
# Aggregation: temperature vs. dielectric constant
############################

SUMMARY_CSV="$SWEEP_DIR/dielectric_vs_temperature.csv"

if [[ "$INTERACTIVE" == "1" ]]; then
  echo ""
  echo "Aggregating all temperatures into $SUMMARY_CSV ..."
  entry_args=()
  while IFS= read -r line; do
    [[ -n "$line" ]] && entry_args+=(--entry "$line")
  done < "$ENTRIES_FILE"
  python3 "$AGGREGATE_SCRIPT" "${entry_args[@]}" --output "$SUMMARY_CSV"

  cat <<SUMMARY

Pipeline completed interactively:
  Sweep directory      : $SWEEP_DIR
  Temperatures         : ${TEMP_LIST[*]}
  Aggregated summary    : $SUMMARY_CSV
SUMMARY
else
  dependency="afterok"
  for jobid in "${STAGE2_JOBIDS[@]}"; do
    dependency+=":$jobid"
  done

  # Preserve the entries file past this script's exit (its own trap removes
  # it) so the aggregate job — which runs later, asynchronously — can read it.
  PERSISTENT_ENTRIES_FILE="$SWEEP_DIR/.dielectric_entries"
  cp "$ENTRIES_FILE" "$PERSISTENT_ENTRIES_FILE"

  echo ""
  echo "Submitting aggregation job, dependent on: ${STAGE2_JOBIDS[*]} ..."
  JOBID3="$(cd "$SWEEP_DIR" && sbatch --parsable --dependency="$dependency" \
    --export="ALL,ENTRIES_FILE=$PERSISTENT_ENTRIES_FILE,OUTPUT_CSV=$SUMMARY_CSV,AGGREGATE_SCRIPT=$AGGREGATE_SCRIPT" \
    "$AGGREGATE_TEMPLATE")"
  echo "  aggregation job id: $JOBID3"

  cat <<SUMMARY

Pipeline submitted:
  Sweep directory      : $SWEEP_DIR
  Temperatures         : ${TEMP_LIST[*]}
  LAMMPS jobs           : ${STAGE1_JOBIDS[*]}
  Dielectric-calc jobs   : ${STAGE2_JOBIDS[*]}
  Aggregation job        : $JOBID3
  Aggregated summary (once $JOBID3 finishes) : $SUMMARY_CSV
SUMMARY
fi
