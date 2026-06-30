#!/bin/bash

# =============================================================================
# Configuration — paths and job sizing
# =============================================================================
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DUMP_DIR="/path/to/dump/dir"    # directory containing dump files
DUMP_GLOB="dump.*"              # adjust to match your dump file naming
OUTPUT_DIR="$SCRIPT_DIR/voxel_output"
FINAL_H5="$SCRIPT_DIR/trajectory.h5"
VENV_PATH="$SCRIPT_DIR/venv"    # path to virtual environment
BATCH_SIZE=64

# =============================================================================
# Configuration — hydronium detection
# =============================================================================
# Rod center (y, z coordinates in Å — rotational symmetry axis of the system)
HYDRONIUM_Y_CENTER=0.0
HYDRONIUM_Z_CENTER=0.0

# =============================================================================
# Configuration — crater analysis
# =============================================================================
CRATER_INITIAL_CUTOFF=5.0    # Å below reference plane to classify a bin as crater (first pass)
CRATER_SECONDARY_CUTOFF=3.0  # Å from fitted sphere surface to keep a point (second pass)

# Sphere fit center: fixed y/z coordinates (Å) — axis of rotational symmetry of the jet
SPHERE_Y_CENTER=0.0
SPHERE_Z_CENTER=0.0

# =============================================================================
# Configuration — diagnostics
# =============================================================================
# Set to 1 to log per-layer and per-smoothing-pass memory usage (RSS) from
# voxel_analysis.py to each task's .err log. Off by default — adds print
# overhead per layer, only enable for memory-scaling investigations.
VOXEL_MEMORY_PROFILE=0

N=$(ls "$DUMP_DIR"/$DUMP_GLOB 2>/dev/null | wc -l)

if [ "$N" -eq 0 ]; then
    echo "ERROR: No files matching '$DUMP_GLOB' found in '$DUMP_DIR'" >&2
    exit 1
fi

mkdir -p "$SCRIPT_DIR/logs"
mkdir -p "$OUTPUT_DIR"

# --- Build list of unprocessed dump files ---
DUMP_PREFIX="${DUMP_GLOB%%\**}"
DUMP_SUFFIX="${DUMP_GLOB##*\*}"
PENDING_FILE="$SCRIPT_DIR/pending_dumps.txt"
> "$PENDING_FILE"

skipped=0
for DUMP in $(ls "$DUMP_DIR"/$DUMP_GLOB | sort); do
    BASENAME=$(basename "$DUMP")
    FRAME_ID="${BASENAME#$DUMP_PREFIX}"
    FRAME_ID="${FRAME_ID%$DUMP_SUFFIX}"
    if [ -f "$OUTPUT_DIR/output_${FRAME_ID}.h5" ]; then
        skipped=$((skipped + 1))
    else
        echo "$DUMP" >> "$PENDING_FILE"
    fi
done

N_PENDING=$(wc -l < "$PENDING_FILE")
echo "Found $N dump files total — $skipped already processed, $N_PENDING pending"

# Merge should only run once all voxels are analyzed and trajectory.h5
# doesn't already exist (avoid redundant/overwriting re-merges).
if [ "$N_PENDING" -eq 0 ]; then
    if [ -f "$FINAL_H5" ]; then
        echo "All frames already processed and $FINAL_H5 already exists. Nothing to do."
        exit 0
    fi

    echo "All frames already processed but $FINAL_H5 does not exist — submitting merge."
    MERGE_JOB=$(sbatch --parsable \
        --export=OUTPUT_DIR=$OUTPUT_DIR,FINAL_H5=$FINAL_H5,SCRIPT_DIR=$SCRIPT_DIR,VENV_PATH=$VENV_PATH,CRATER_INITIAL_CUTOFF=$CRATER_INITIAL_CUTOFF,CRATER_SECONDARY_CUTOFF=$CRATER_SECONDARY_CUTOFF,SPHERE_Y_CENTER=$SPHERE_Y_CENTER,SPHERE_Z_CENTER=$SPHERE_Z_CENTER \
        "$SCRIPT_DIR/merge_h5.slurm")

    echo "Merge job ID: $MERGE_JOB"
    exit 0
fi

if [ -f "$FINAL_H5" ]; then
    echo "WARNING: $FINAL_H5 already exists but new dump files are pending."
    echo "         The merge job will overwrite it once the array job completes."
fi

ARRAY_MAX=$((N_PENDING - 1))
echo "Submitting array 0-${ARRAY_MAX} (max ${BATCH_SIZE} simultaneous)"

ARRAY_JOB=$(sbatch --parsable \
    --array=0-${ARRAY_MAX}%${BATCH_SIZE} \
    --export=PENDING_FILE=$PENDING_FILE,DUMP_GLOB=$DUMP_GLOB,OUTPUT_DIR=$OUTPUT_DIR,SCRIPT_DIR=$SCRIPT_DIR,VENV_PATH=$VENV_PATH,HYDRONIUM_Y_CENTER=$HYDRONIUM_Y_CENTER,HYDRONIUM_Z_CENTER=$HYDRONIUM_Z_CENTER,VOXEL_MEMORY_PROFILE=$VOXEL_MEMORY_PROFILE \
    "$SCRIPT_DIR/voxel_analysis.slurm")

if [ -z "$ARRAY_JOB" ]; then
    echo "ERROR: Array job submission failed" >&2
    exit 1
fi

echo "Array job ID: $ARRAY_JOB"

MERGE_JOB=$(sbatch --parsable \
    --dependency=afterok:$ARRAY_JOB \
    --export=OUTPUT_DIR=$OUTPUT_DIR,FINAL_H5=$FINAL_H5,SCRIPT_DIR=$SCRIPT_DIR,VENV_PATH=$VENV_PATH,CRATER_INITIAL_CUTOFF=$CRATER_INITIAL_CUTOFF,CRATER_SECONDARY_CUTOFF=$CRATER_SECONDARY_CUTOFF,SPHERE_Y_CENTER=$SPHERE_Y_CENTER,SPHERE_Z_CENTER=$SPHERE_Z_CENTER \
    "$SCRIPT_DIR/merge_h5.slurm")

echo "Merge job ID: $MERGE_JOB (runs after array completes)"
