#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/voxel_%A_%a.out

# GNU time (supports -v/-o) isn't always on PATH by default on CARC — load
# the module if available, BEFORE activating the venv. module load can
# prepend to PATH; loading it after the venv would risk shadowing the
# venv's python and breaking numpy/h5py imports. Falls through silently
# if module/time is absent.
module load gcc/13.3.0 time/1.9 >/dev/null 2>&1 || true

source "$VENV_PATH/bin/activate"

# VOXEL_MEMORY_PROFILE is set in submit_voxel_pipeline.sh (Configuration —
# diagnostics) and passed through via --export. When set to 1, per-layer
# and per-smoothing-pass memory usage (RSS) is logged to this task's .err log.

mkdir -p "$OUTPUT_DIR"

DUMP=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$PENDING_FILE")

if [ -z "$DUMP" ]; then
    echo "ERROR: No dump file found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

DUMP_PREFIX="${DUMP_GLOB%%\**}"
DUMP_SUFFIX="${DUMP_GLOB##*\*}"
BASENAME=$(basename "$DUMP")
FRAME_ID="${BASENAME#$DUMP_PREFIX}"
FRAME_ID="${FRAME_ID%$DUMP_SUFFIX}"

rm -f "${OUTPUT_DIR}/output_${FRAME_ID}.h5.tmp"

echo "Task ${SLURM_ARRAY_TASK_ID}: processing ${DUMP} (frame ${FRAME_ID})"

SUMMARY_LOG="$SCRIPT_DIR/logs/summary_${SLURM_ARRAY_JOB_ID}.log"
TIME_LOG="$SCRIPT_DIR/logs/time_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.log"

TIME_BIN=$(command -v time 2>/dev/null || true)
if [ -n "$TIME_BIN" ]; then
    TIME_CMD=("$TIME_BIN" -v -o "$TIME_LOG")
else
    TIME_CMD=()
fi

"${TIME_CMD[@]}" python "$SCRIPT_DIR/voxel_analysis.py" "$DUMP" "${OUTPUT_DIR}/output_${FRAME_ID}.h5" "$HYDRONIUM_Y_CENTER" "$HYDRONIUM_Z_CENTER" "$BUBBLE_RADIUS" \
    2>"$SCRIPT_DIR/logs/voxel_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
EXIT_CODE=$?

STDERR_CONTENT=$(cat "$SCRIPT_DIR/logs/voxel_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err")
if [ ${#TIME_CMD[@]} -gt 0 ]; then
    MAX_RSS_KB=$(grep "Maximum resident set size" "$TIME_LOG" | awk '{print $NF}')
else
    MAX_RSS_KB="n/a"
fi

if [ $EXIT_CODE -ne 0 ]; then
    STATUS="ERROR (exit $EXIT_CODE)"
elif echo "$STDERR_CONTENT" | grep -qi "warning"; then
    STATUS="WARNING"
else
    STATUS="OK"
fi

TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
echo "${TIMESTAMP}  task=${SLURM_ARRAY_TASK_ID}  frame=${FRAME_ID}  ${STATUS}  max_rss=${MAX_RSS_KB}KB" >> "$SUMMARY_LOG"

# also print stderr to the per-task stdout log for full context
if [ -n "$STDERR_CONTENT" ]; then
    echo "--- stderr ---"
    echo "$STDERR_CONTENT"
    echo "--------------"
fi

exit $EXIT_CODE
