#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00

# Per-chunk array job — SLURM_ARRAY_TASK_ID is the index within this chunk (0 .. CHUNK_NUM_FRAMES-1).
# CHUNK_ID is still required: it maps the local task ID to the global frame index and output paths.
# Required: CHUNK_ID, CHUNK_SIZE, OUTPUT_DIR, SCRIPT_DIR, DUMP_DIR, DUMP_GLOB, START_TIMESTEP,
#           DUMP_EVERY, CUTOFF, TYPE_O, TYPE_H, VENV_PATH (optional)

module purge
module load usc

if [ -n "$VENV_PATH" ]; then
    source "$VENV_PATH/bin/activate"
fi

: "${CHUNK_ID:?CHUNK_ID must be set via sbatch --export (see submit_pipeline.sh)}"

FRAME_INDEX=$((CHUNK_ID * CHUNK_SIZE + SLURM_ARRAY_TASK_ID))
TIMESTEP=$((START_TIMESTEP + FRAME_INDEX * DUMP_EVERY))
CHUNK_DIR="$OUTPUT_DIR/chunk_${CHUNK_ID}"
CHUNK_OUTPUT="$OUTPUT_DIR/dipole_output_chunk_${CHUNK_ID}.txt"
FRAME_OUTPUT="$CHUNK_DIR/dipole_${TIMESTEP}.txt"
DUMP="$DUMP_DIR/${DUMP_GLOB/\*/$TIMESTEP}"

# 1. If chunk is not yet aggregated, ensure the chunk work directory exists.
if [ -f "$CHUNK_OUTPUT" ]; then
    exit 0
fi
mkdir -p "$CHUNK_DIR"

# 2. Skip if this frame was already calculated.
if [ -f "$FRAME_OUTPUT" ]; then
    exit 0
fi

# 3. Require the dump file.
if [ ! -f "$DUMP" ]; then
    echo "ERROR: Dump file not found for frame ${FRAME_INDEX} (timestep ${TIMESTEP}): ${DUMP}" >&2
    exit 1
fi

# 4. Run the calculation.
echo "Frame ${FRAME_INDEX} (chunk ${CHUNK_ID}, timestep ${TIMESTEP}): ${DUMP}"

python "$SCRIPT_DIR/1.calc_dipole_from_dump.py" "$DUMP" \
    --cutoff "$CUTOFF" \
    --type-o "$TYPE_O" \
    --type-h "$TYPE_H" \
    --output "$FRAME_OUTPUT" \
    --log    "$CHUNK_DIR/warn_${TIMESTEP}.log"
