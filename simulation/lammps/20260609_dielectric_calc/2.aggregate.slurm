#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00

# Single job per chunk — runs after that chunk's calc array tasks complete.
# Required: CHUNK_ID, CHUNK_NUM_FRAMES, OUTPUT_DIR, SCRIPT_DIR, VENV_PATH (optional)

module purge
module load usc

if [ -n "$VENV_PATH" ]; then
    source "$VENV_PATH/bin/activate"
fi

CHUNK_DIR="$OUTPUT_DIR/chunk_${CHUNK_ID}"
CHUNK_OUTPUT="$OUTPUT_DIR/dipole_output_chunk_${CHUNK_ID}.txt"
CHUNK_WARN="$OUTPUT_DIR/dipole_warnings_chunk_${CHUNK_ID}.log"

# 1. Skip if this chunk is already aggregated.
if [ -f "$CHUNK_OUTPUT" ]; then
    echo "Chunk ${CHUNK_ID} already aggregated — skipping"
    exit 0
fi

# 2. Require per-frame dipole files in the chunk directory.
shopt -s nullglob
dipole_files=("$CHUNK_DIR"/dipole_*.txt)
if [ ${#dipole_files[@]} -eq 0 ]; then
    echo "ERROR: No dipole_*.txt files in ${CHUNK_DIR}" >&2
    exit 1
fi

# 3. Run aggregation.
echo "Aggregating chunk ${CHUNK_ID} at $(date)"
echo "  input : ${CHUNK_DIR}"
echo "  frames: ${CHUNK_NUM_FRAMES} (expected)"

agg_args=(
    --input-dir "$CHUNK_DIR"
    --output "$CHUNK_OUTPUT"
    --warn-log "$CHUNK_WARN"
)
if [ -n "$CHUNK_NUM_FRAMES" ]; then
    agg_args+=(--expected-frames "$CHUNK_NUM_FRAMES")
fi

python "$SCRIPT_DIR/2.aggregate.py" "${agg_args[@]}"

# 4. Remove per-frame files after successful aggregation.
rm -f "$CHUNK_DIR"/dipole_*.txt "$CHUNK_DIR"/warn_*.log

echo "Chunk ${CHUNK_ID} aggregate finished at $(date)"
