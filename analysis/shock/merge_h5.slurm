#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/merge_%j.out

# GNU time (supports -v/-o) isn't always on PATH by default on CARC — load
# the module if available, BEFORE activating the venv. module load can
# prepend to PATH; loading it after the venv would risk shadowing the
# venv's python and breaking numpy/h5py imports. Falls through silently
# if module/time is absent.
module load gcc/13.3.0 time/1.9 >/dev/null 2>&1 || true

source "$VENV_PATH/bin/activate"

TIME_LOG="$SCRIPT_DIR/logs/time_merge_${SLURM_JOB_ID}.log"

TIME_BIN=$(command -v time 2>/dev/null || true)
if [ -n "$TIME_BIN" ]; then
    TIME_CMD=("$TIME_BIN" -v -o "$TIME_LOG")
else
    TIME_CMD=()
fi

"${TIME_CMD[@]}" python "$SCRIPT_DIR/merge_h5.py" "$OUTPUT_DIR" "$FINAL_H5" \
    --initial-cutoff "$CRATER_INITIAL_CUTOFF" \
    --secondary-cutoff "$CRATER_SECONDARY_CUTOFF" \
    --sphere-y-center "$SPHERE_Y_CENTER" \
    --sphere-z-center "$SPHERE_Z_CENTER"

if [ -n "$TIME_BIN" ]; then
    MAX_RSS_KB=$(grep "Maximum resident set size" "$TIME_LOG" | awk '{print $NF}')
    echo "Peak memory usage: ${MAX_RSS_KB} KB"
else
    echo "GNU time not available — peak memory not logged (use seff/sacct after the job completes)"
fi

############################
# Ledger
############################

/apps/spack/2406/apps/linux-rocky8-x86_64_v3/gcc-13.3.0/python-3.11.9-x74mtjf/bin/python3 /home1/lkyamamo/util/jobs/ledger/bin/simledger hook analysis --dir "$DUMP_DIR" --type voxel_shock_density \
    --output "$FINAL_H5" \
    || echo "WARNING: ledger hook failed" >&2
