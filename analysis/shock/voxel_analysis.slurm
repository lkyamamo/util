#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/voxel_%A_%a.out

source "$VENV_PATH/bin/activate"

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

python "$SCRIPT_DIR/voxel_analysis.py" "$DUMP" "${OUTPUT_DIR}/output_${FRAME_ID}.h5" \
    2>"$SCRIPT_DIR/logs/voxel_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err"
EXIT_CODE=$?

STDERR_CONTENT=$(cat "$SCRIPT_DIR/logs/voxel_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err")

if [ $EXIT_CODE -ne 0 ]; then
    STATUS="ERROR (exit $EXIT_CODE)"
elif echo "$STDERR_CONTENT" | grep -qi "warning"; then
    STATUS="WARNING"
else
    STATUS="OK"
fi

TIMESTAMP=$(date '+%Y-%m-%d %H:%M:%S')
echo "${TIMESTAMP}  task=${SLURM_ARRAY_TASK_ID}  frame=${FRAME_ID}  ${STATUS}" >> "$SUMMARY_LOG"

# also print stderr to the per-task stdout log for full context
if [ -n "$STDERR_CONTENT" ]; then
    echo "--- stderr ---"
    echo "$STDERR_CONTENT"
    echo "--------------"
fi

exit $EXIT_CODE
