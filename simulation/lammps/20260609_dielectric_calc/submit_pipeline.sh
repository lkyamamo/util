#!/bin/bash

# =============================================================================
# Dielectric constant pipeline — chunked submission
#
# 1. Per-chunk calc array — 1.calc.slurm (wave-limited: MAX_ACTIVE_CHUNKS at a time)
# 2. Per-chunk aggregate  — 2.aggregate.slurm (afterok on that chunk's calc only)
# 3. Final job            — 3.final.slurm (combine + dipoleStd)
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

# Wave submission: at most MAX_ACTIVE_CHUNKS calc arrays running at once.
# Each calc array runs at most MAX_SIMULTANEOUS tasks (% throttle).
# Aggregate jobs depend only on their chunk's calc array (not on wave limits).
MAX_ACTIVE_CHUNKS=4
MAX_SIMULTANEOUS=16

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

MAX_ARRAY_SIZE=$(scontrol show config 2>/dev/null | awk -F= '/^MaxArraySize/ {gsub(/ /,"",$2); print $2; exit}')
if [ -n "$MAX_ARRAY_SIZE" ] && [ "$MAX_ARRAY_SIZE" -lt "$CHUNK_SIZE" ]; then
    echo "ERROR: CHUNK_SIZE=$CHUNK_SIZE exceeds cluster MaxArraySize=$MAX_ARRAY_SIZE" >&2
    echo "Reduce CHUNK_SIZE or ask admins to raise MaxArraySize." >&2
    exit 1
fi

CALC_EXPORT_BASE="DUMP_GLOB=$DUMP_GLOB,DUMP_DIR=$DUMP_DIR,OUTPUT_DIR=$OUTPUT_DIR"
CALC_EXPORT_BASE="${CALC_EXPORT_BASE},SCRIPT_DIR=$SCRIPT_DIR,START_TIMESTEP=$START_TIMESTEP,DUMP_EVERY=$DUMP_EVERY"
CALC_EXPORT_BASE="${CALC_EXPORT_BASE},CHUNK_SIZE=$CHUNK_SIZE,CUTOFF=$CUTOFF,TYPE_O=$TYPE_O,TYPE_H=$TYPE_H"
CALC_EXPORT_BASE="${CALC_EXPORT_BASE},VENV_PATH=$VENV_PATH"

echo "Frames: $NUM_FRAMES  Chunks: $NUM_CHUNKS  Chunk size: $CHUNK_SIZE"
echo "Wave size: ${MAX_ACTIVE_CHUNKS} calc arrays  Tasks per chunk: ${MAX_SIMULTANEOUS}"
[ -n "$MAX_ARRAY_SIZE" ] && echo "Cluster MaxArraySize: $MAX_ARRAY_SIZE"

# =============================================================================
# SUBMIT — one calc array + one aggregate per chunk; then final
# =============================================================================

CALC_JOBS=()
AGG_JOBS=()

for ((chunk=0; chunk<NUM_CHUNKS; chunk++)); do
    CHUNK_NUM_FRAMES=$CHUNK_SIZE
    if (( chunk == NUM_CHUNKS - 1 )); then
        CHUNK_NUM_FRAMES=$((NUM_FRAMES - chunk * CHUNK_SIZE))
    fi
    CHUNK_ARRAY_MAX=$((CHUNK_NUM_FRAMES - 1))

    if [ -n "$MAX_ARRAY_SIZE" ] && [ "$CHUNK_NUM_FRAMES" -gt "$MAX_ARRAY_SIZE" ]; then
        echo "ERROR: Chunk ${chunk} has ${CHUNK_NUM_FRAMES} frames > MaxArraySize=${MAX_ARRAY_SIZE}" >&2
        exit 1
    fi

    CALC_EXPORT="${CALC_EXPORT_BASE},CHUNK_ID=$chunk"

    # Calc waves: chunk N's calc waits for the previous wave's calc arrays to finish.
    # Aggregate jobs are submitted immediately and only depend on their own calc array.
    CALC_DEP=""
    wave=$((chunk / MAX_ACTIVE_CHUNKS))
    if (( wave > 0 )); then
        prev_wave_start=$(( (wave - 1) * MAX_ACTIVE_CHUNKS ))
        prev_wave_end=$(( wave * MAX_ACTIVE_CHUNKS - 1 ))
        for ((c=prev_wave_start; c<=prev_wave_end && c<NUM_CHUNKS; c++)); do
            CALC_DEP="${CALC_DEP}afterok:${CALC_JOBS[$c]},"
        done
        CALC_DEP="${CALC_DEP%,}"
    fi

    if (( chunk % MAX_ACTIVE_CHUNKS == 0 )); then
        wave_end=$((chunk + MAX_ACTIVE_CHUNKS - 1))
        (( wave_end >= NUM_CHUNKS )) && wave_end=$((NUM_CHUNKS - 1))
        echo "Submitting calc wave $wave (chunks ${chunk}-${wave_end})"
    fi
    echo "  chunk ${chunk} calc array 0-${CHUNK_ARRAY_MAX} (max ${MAX_SIMULTANEOUS} simultaneous)"
    [ -n "$CALC_DEP" ] && echo "  calc dependency: ${CALC_DEP}"

    if [ -n "$CALC_DEP" ]; then
        CALC_JOB=$(sbatch --parsable \
            --dependency="$CALC_DEP" \
            --array=0-${CHUNK_ARRAY_MAX}%${MAX_SIMULTANEOUS} \
            --output="$SCRIPT_DIR/logs/calc_chunk${chunk}_%A_%a.out" \
            --export="$CALC_EXPORT" \
            "$SCRIPT_DIR/1.calc.slurm")
    else
        CALC_JOB=$(sbatch --parsable \
            --array=0-${CHUNK_ARRAY_MAX}%${MAX_SIMULTANEOUS} \
            --output="$SCRIPT_DIR/logs/calc_chunk${chunk}_%A_%a.out" \
            --export="$CALC_EXPORT" \
            "$SCRIPT_DIR/1.calc.slurm")
    fi

    if [ -z "$CALC_JOB" ]; then
        echo "ERROR: Calc array submission failed for chunk $chunk" >&2
        exit 1
    fi

    CALC_JOBS+=("$CALC_JOB")
    echo "  calc job ID: $CALC_JOB"

    AGG_EXPORT="CHUNK_ID=$chunk,CHUNK_NUM_FRAMES=$CHUNK_NUM_FRAMES"
    AGG_EXPORT="${AGG_EXPORT},OUTPUT_DIR=$OUTPUT_DIR,SCRIPT_DIR=$SCRIPT_DIR,VENV_PATH=$VENV_PATH"

    # Run aggregate only after all calc array tasks for this chunk succeed.
    AGG_DEP="afterok:${CALC_JOB}"
    echo "  aggregate dependency: ${AGG_DEP}"

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
    echo "  aggregate job ID: $AGG_JOB"
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
for ((chunk=0; chunk<NUM_CHUNKS; chunk++)); do
    echo "--- calc chunk ${chunk} (${CALC_JOBS[$chunk]}) ---"
    scontrol show job "${CALC_JOBS[$chunk]}" | grep -E '^JobId=|JobName=|JobState=|Dependency='
    echo "--- aggregate chunk ${chunk} (${AGG_JOBS[$chunk]}) ---"
    scontrol show job "${AGG_JOBS[$chunk]}" | grep -E '^JobId=|JobName=|JobState=|Dependency='
done
echo "--- final ($FINAL_JOB) ---"
scontrol show job "$FINAL_JOB" | grep -E '^JobId=|JobName=|JobState=|Dependency='

echo ""
echo "Monitor jobs:    squeue -u \$USER"
echo "Final job:       squeue -j $FINAL_JOB"
echo "Chunk outputs:   $OUTPUT_DIR/dipole_output_chunk_*.txt"
echo "Combined output: $OUTPUT_DIR/dipole_output.txt"
echo "Dielectric CSV:  $OUTPUT_DIR/dipole_output_timestep_data_${AVERAGING_METHOD}.csv"
