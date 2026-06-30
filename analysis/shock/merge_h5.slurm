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

if command -v /usr/bin/time >/dev/null 2>&1; then
    TIME_CMD=(/usr/bin/time -v -o "$TIME_LOG")
else
    TIME_CMD=()
fi

"${TIME_CMD[@]}" python "$SCRIPT_DIR/merge_h5.py" "$OUTPUT_DIR" "$FINAL_H5" \
    --initial-cutoff "$CRATER_INITIAL_CUTOFF" \
    --secondary-cutoff "$CRATER_SECONDARY_CUTOFF" \
    --sphere-y-center "$SPHERE_Y_CENTER" \
    --sphere-z-center "$SPHERE_Z_CENTER"

if [ ${#TIME_CMD[@]} -gt 0 ]; then
    MAX_RSS_KB=$(grep "Maximum resident set size" "$TIME_LOG" | awk '{print $NF}')
    echo "Peak memory usage: ${MAX_RSS_KB} KB"
else
    echo "/usr/bin/time not available — peak memory not logged (use seff/sacct after the job completes)"
fi
