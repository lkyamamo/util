#!/bin/bash
# Dielectric constant pipeline — submit script.
#
# Computes the job array size from configuration, submits the calc array,
# then submits the final job with an afterok dependency.
#
# Usage (from this directory on the HPC login node):
#   bash 0.submit.sh
#
# Edit the Configuration section below before running.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# --- Configuration ---
DUMP_DIR="../dumps"
DUMP_GLOB="dielectric.*.custom"
OUTPUT_DIR="$SCRIPT_DIR/dipole_output"
START_TIMESTEP=0
END_TIMESTEP=60000000
DUMP_EVERY=10

CUTOFF=1.2
TYPE_O=1
TYPE_H=2

TEMPERATURE=303.0
LA=37.2514
LB=37.2514
LC=37.2514
AVERAGING_METHOD=cumulative

VENV_PATH="/home1/lkyamamo/venv/struc_analysis"

# ONE_SHOT=true : ignore GROUP_SIZE; distribute frames evenly across
#   MAX_SIMULTANEOUS tasks (one per core) so the entire run fits in a
#   single array submission with no throttle.
# ONE_SHOT=false: use GROUP_SIZE directly; throttle array to MAX_SIMULTANEOUS.
ONE_SHOT=true
MAX_SIMULTANEOUS=128
GROUP_SIZE=5000          # used only when ONE_SHOT=false

# --- Derived ---
NUM_FRAMES=$(( (END_TIMESTEP - START_TIMESTEP) / DUMP_EVERY + 1 ))

if [ "$ONE_SHOT" = true ]; then
    NUM_GROUPS=$MAX_SIMULTANEOUS
    GROUP_SIZE=$(( (NUM_FRAMES + NUM_GROUPS - 1) / NUM_GROUPS ))
    ARRAY_THROTTLE=""
else
    NUM_GROUPS=$(( (NUM_FRAMES + GROUP_SIZE - 1) / GROUP_SIZE ))
    ARRAY_THROTTLE="%${MAX_SIMULTANEOUS}"
fi

ARRAY_MAX=$(( NUM_GROUPS - 1 ))

mkdir -p "$SCRIPT_DIR/logs" "$OUTPUT_DIR/groups"

echo "ONE_SHOT: ${ONE_SHOT}"
echo "Frames: ${NUM_FRAMES}  Groups: ${NUM_GROUPS}  Group size: ${GROUP_SIZE}"
echo "Array: 0-${ARRAY_MAX}${ARRAY_THROTTLE}"

CALC_EXPORT="DUMP_GLOB=$DUMP_GLOB,DUMP_DIR=$DUMP_DIR,OUTPUT_DIR=$OUTPUT_DIR"
CALC_EXPORT="${CALC_EXPORT},SCRIPT_DIR=$SCRIPT_DIR,START_TIMESTEP=$START_TIMESTEP"
CALC_EXPORT="${CALC_EXPORT},DUMP_EVERY=$DUMP_EVERY,GROUP_SIZE=$GROUP_SIZE"
CALC_EXPORT="${CALC_EXPORT},NUM_FRAMES=$NUM_FRAMES,CUTOFF=$CUTOFF"
CALC_EXPORT="${CALC_EXPORT},TYPE_O=$TYPE_O,TYPE_H=$TYPE_H,VENV_PATH=$VENV_PATH"

CALC_JOB=$(sbatch --parsable \
    --array=0-${ARRAY_MAX}${ARRAY_THROTTLE} \
    --output="$SCRIPT_DIR/logs/calc_%A_%a.out" \
    --export="$CALC_EXPORT" \
    "$SCRIPT_DIR/1.calc_group.slurm")
echo "Calc array job: ${CALC_JOB}"

FINAL_EXPORT="OUTPUT_DIR=$OUTPUT_DIR,NUM_GROUPS=$NUM_GROUPS,SCRIPT_DIR=$SCRIPT_DIR"
FINAL_EXPORT="${FINAL_EXPORT},TEMPERATURE=$TEMPERATURE,LA=$LA,LB=$LB,LC=$LC"
FINAL_EXPORT="${FINAL_EXPORT},AVERAGING_METHOD=$AVERAGING_METHOD,VENV_PATH=$VENV_PATH"

FINAL_JOB=$(sbatch --parsable \
    --dependency=afterok:${CALC_JOB} \
    --output="$SCRIPT_DIR/logs/final_%j.out" \
    --export="$FINAL_EXPORT" \
    "$SCRIPT_DIR/2.final.slurm")
echo "Final job: ${FINAL_JOB}"
