#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/merge_%j.out

source "$VENV_PATH/bin/activate"

TIME_LOG="$SCRIPT_DIR/logs/time_merge_${SLURM_JOB_ID}.log"

/usr/bin/time -v -o "$TIME_LOG" \
    python "$SCRIPT_DIR/merge_h5.py" "$OUTPUT_DIR" "$FINAL_H5" \
    --initial-cutoff "$CRATER_INITIAL_CUTOFF" \
    --secondary-cutoff "$CRATER_SECONDARY_CUTOFF" \
    --sphere-y-center "$SPHERE_Y_CENTER" \
    --sphere-z-center "$SPHERE_Z_CENTER"

MAX_RSS_KB=$(grep "Maximum resident set size" "$TIME_LOG" | awk '{print $NF}')
echo "Peak memory usage: ${MAX_RSS_KB} KB"
