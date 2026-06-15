#!/bin/bash

# =============================================================================
# Dielectric constant pipeline — chunked submission
#
# 1. Global calc array  — 1.calc.slurm (one task per frame, %MAX_SIMULTANEOUS)
# 2. Per-chunk aggregate — 2.aggregate.slurm (depends on chunk calc task range)
# 3. Final job          — 3.final.slurm (combine + dipoleStd)
#
# Usage (run on HPC login node):
#   bash submit_pipeline.sh
#
# Edit the Configuration section below before running.
# =============================================================================

# --- Configuration ---
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
DUMP_DIR="../dumps"
DUMP_GLOB="dielectric.*.custom"
OUTPUT_DIR="$SCRIPT_DIR/dipole_output"
START_TIMESTEP=0
END_TIMESTEP=60000000
DUMP_EVERY=10
CHUNK_SIZE=100000

CUTOFF=1.2
TYPE_O=1
TYPE_H=2

TEMPERATURE=303.0
LA=37.2514
LB=37.2514
LC=37.2514
AVERAGING_METHOD=cumulative

VENV_PATH="/home1/lkyamamo/venv/struc_analysis"

# Pipeline-wide cap on concurrent calc array tasks
MAX_SIMULTANEOUS=64

# =============================================================================
# SETUP
# =============================================================================

if [ ! -f "$DUMP_DIR/${DUMP_GLOB/\*/$START_TIMESTEP}" ]; then
    echo "ERROR: No files matching '${DUMP_GLOB/\*/$START_TIMESTEP}' found in '$DUMP_DIR'" >&2
    exit 1
fi

mkdir -p "$SCRIPT_DIR/logs" "$OUTPUT_DIR"

NUM_FRAMES=$(( (END_TIMESTEP - START_TIMESTEP) / DUMP_EVERY + 1 ))
NUM_CHUNKS=$(( (NUM_FRAMES + CHUNK_SIZE - 1) / CHUNK_SIZE ))

CALC_EXPORT="DUMP_GLOB=$DUMP_GLOB,DUMP_DIR=$DUMP_DIR,OUTPUT_DIR=$OUTPUT_DIR"
CALC_EXPORT="${CALC_EXPORT},SCRIPT_DIR=$SCRIPT_DIR,START_TIMESTEP=$START_TIMESTEP,DUMP_EVERY=$DUMP_EVERY"
CALC_EXPORT="${CALC_EXPORT},CHUNK_SIZE=$CHUNK_SIZE,CUTOFF=$CUTOFF,TYPE_O=$TYPE_O,TYPE_H=$TYPE_H"
CALC_EXPORT="${CALC_EXPORT},VENV_PATH=$VENV_PATH"

echo "Frames: $NUM_FRAMES  Chunks: $NUM_CHUNKS  Chunk size: $CHUNK_SIZE"
echo "Max simultaneous calc tasks: $MAX_SIMULTANEOUS"

# =============================================================================
# SUBMIT
# =============================================================================

echo "Submitting global calc array 0-$((NUM_FRAMES - 1)) (max ${MAX_SIMULTANEOUS} simultaneous)"

CALC_JOB=$(sbatch --parsable \
    --array=0-$((NUM_FRAMES - 1))%${MAX_SIMULTANEOUS} \
    --output="$SCRIPT_DIR/logs/calc_%A_%a.out" \
    --export="$CALC_EXPORT" \
    "$SCRIPT_DIR/1.calc.slurm")

if [ -z "$CALC_JOB" ]; then
    echo "ERROR: Calc array submission failed" >&2
    exit 1
fi

echo "Calc array job ID: $CALC_JOB"

AGG_JOBS=()
for ((chunk=0; chunk<NUM_CHUNKS; chunk++)); do
    TASK_START=$((chunk * CHUNK_SIZE))
    CHUNK_NUM_FRAMES=$CHUNK_SIZE
    if (( chunk == NUM_CHUNKS - 1 )); then
        CHUNK_NUM_FRAMES=$((NUM_FRAMES - TASK_START))
    fi
    TASK_END=$((TASK_START + CHUNK_NUM_FRAMES - 1))

    AGG_EXPORT="CHUNK_ID=$chunk,CHUNK_NUM_FRAMES=$CHUNK_NUM_FRAMES"
    AGG_EXPORT="${AGG_EXPORT},OUTPUT_DIR=$OUTPUT_DIR,SCRIPT_DIR=$SCRIPT_DIR,VENV_PATH=$VENV_PATH"

    # Run aggregate only after all calc array tasks for this chunk succeed.
    AGG_DEP="afterok:${CALC_JOB}_${TASK_START}-${TASK_END}"
    echo "  dependency: ${AGG_DEP}"

    AGG_JOB=$(sbatch --parsable \
        --dependency="$AGG_DEP" \
        --output="$SCRIPT_DIR/logs/agg_chunk${chunk}_%j.out" \
        --export="$AGG_EXPORT" \
        "$SCRIPT_DIR/2.aggregate.slurm")

    if [ -z "$AGG_JOB" ]; then
        echo "ERROR: Aggregate submission failed for chunk $chunk" >&2
        exit 1
    fi

    AGG_JOBS+=("$AGG_JOB")
    echo "Chunk ${chunk} aggregate job ID: $AGG_JOB (tasks ${TASK_START}-${TASK_END})"
done

# Run final job only after all chunk aggregate jobs succeed.
FINAL_DEP=""
for agg_id in "${AGG_JOBS[@]}"; do
    FINAL_DEP="${FINAL_DEP}afterok:${agg_id},"
done
FINAL_DEP="${FINAL_DEP%,}"
echo "Final dependency: ${FINAL_DEP}"

FINAL_EXPORT="OUTPUT_DIR=$OUTPUT_DIR,NUM_CHUNKS=$NUM_CHUNKS,SCRIPT_DIR=$SCRIPT_DIR"
FINAL_EXPORT="${FINAL_EXPORT},TEMPERATURE=$TEMPERATURE,LA=$LA,LB=$LB,LC=$LC"
FINAL_EXPORT="${FINAL_EXPORT},AVERAGING_METHOD=$AVERAGING_METHOD,VENV_PATH=$VENV_PATH"

FINAL_JOB=$(sbatch --parsable \
    --dependency="$FINAL_DEP" \
    --output="$SCRIPT_DIR/logs/final_%j.out" \
    --export="$FINAL_EXPORT" \
    "$SCRIPT_DIR/3.final.slurm")

if [ -z "$FINAL_JOB" ]; then
    echo "ERROR: Final job submission failed" >&2
    exit 1
fi

echo "Final job ID: $FINAL_JOB (combine + dipoleStd)"
echo ""

echo "=== scontrol dependency check ==="
echo "--- calc array ($CALC_JOB) ---"
scontrol show job "$CALC_JOB" | grep -E '^JobId=|JobName=|JobState=|Dependency='
for ((chunk=0; chunk<NUM_CHUNKS; chunk++)); do
    echo "--- aggregate chunk ${chunk} (${AGG_JOBS[$chunk]}) ---"
    scontrol show job "${AGG_JOBS[$chunk]}" | grep -E '^JobId=|JobName=|JobState=|Dependency='
done
echo "--- final ($FINAL_JOB) ---"
scontrol show job "$FINAL_JOB" | grep -E '^JobId=|JobName=|JobState=|Dependency='

echo ""
echo "Monitor calc:    squeue -j $CALC_JOB"
echo "Monitor final:   squeue -j $FINAL_JOB"
echo "Chunk outputs:   $OUTPUT_DIR/dipole_output_chunk_*.txt"
echo "Combined output: $OUTPUT_DIR/dipole_output.txt"
echo "Dielectric CSV:  $OUTPUT_DIR/dipole_output_timestep_data_${AVERAGING_METHOD}.csv"
