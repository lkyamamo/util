#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3800M
#SBATCH --time=02:00:00

# Per-group array job — SLURM_ARRAY_TASK_ID is the group index.
# Each task processes GROUP_SIZE consecutive frames sequentially (one frame in
# memory at a time), writing per-frame results to a scratch directory, then
# aggregates them into a single group output file and cleans up.
#
# Required exports: DUMP_GLOB, DUMP_DIR, OUTPUT_DIR, SCRIPT_DIR,
#   START_TIMESTEP, DUMP_EVERY, GROUP_SIZE, NUM_FRAMES,
#   CUTOFF, TYPE_O, TYPE_H, VENV_PATH (optional)

module purge
module load usc

if [ -n "$VENV_PATH" ]; then
    source "$VENV_PATH/bin/activate"
fi

: "${GROUP_SIZE:?}"
: "${NUM_FRAMES:?}"

GROUP_ID=$SLURM_ARRAY_TASK_ID
GROUP_OUTPUT="$OUTPUT_DIR/groups/dipole_group_${GROUP_ID}.txt"

mkdir -p "$OUTPUT_DIR/groups"

python "$SCRIPT_DIR/1.calc_group.py" \
    --group-id       "$GROUP_ID" \
    --group-size     "$GROUP_SIZE" \
    --num-frames     "$NUM_FRAMES" \
    --start-timestep "$START_TIMESTEP" \
    --dump-every     "$DUMP_EVERY" \
    --dump-dir       "$DUMP_DIR" \
    --dump-glob      "$DUMP_GLOB" \
    --cutoff         "$CUTOFF" \
    --type-o         "$TYPE_O" \
    --type-h         "$TYPE_H" \
    --output         "$GROUP_OUTPUT"
