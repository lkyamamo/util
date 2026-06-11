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

echo "Task ${SLURM_ARRAY_TASK_ID}: processing ${DUMP} (frame ${FRAME_ID})"

python "$SCRIPT_DIR/voxel_analysis.py" "$DUMP" "${OUTPUT_DIR}/output_${FRAME_ID}.h5"
