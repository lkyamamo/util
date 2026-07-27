#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=logs/merge_%j.out

source "$VENV_PATH/bin/activate"

python "$SCRIPT_DIR/merge_h5.py" "$OUTPUT_DIR" "$FINAL_H5" \
    --initial-cutoff "$CRATER_INITIAL_CUTOFF" \
    --secondary-cutoff "$CRATER_SECONDARY_CUTOFF" \
    --sphere-y-center "$SPHERE_Y_CENTER" \
    --sphere-z-center "$SPHERE_Z_CENTER"

############################
# Ledger
############################

python3 /home1/lkyamamo/util/jobs/ledger/ledger.py add-analysis --dir "$DUMP_DIR" voxel_shock_density \
    --output "$FINAL_H5" \
    || echo "WARNING: ledger add-analysis failed" >&2
