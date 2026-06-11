#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/voxel_%A_%a.out

source "$VENV_PATH/bin/activate"

mkdir -p "$OUTPUT_DIR"
DUMP=$(ls "$DUMP_DIR"/$DUMP_GLOB | sort | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")
python "$SCRIPT_DIR/voxel_analysis.py" "$DUMP" "${OUTPUT_DIR}/output_${SLURM_ARRAY_TASK_ID}.h5"
