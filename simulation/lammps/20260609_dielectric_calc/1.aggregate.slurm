#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:30:00
#SBATCH --output=logs/aggregate_%j.out

# --- Required environment variables (set by submit_pipeline.sh via --export) ---
# OUTPUT_DIR  : directory containing per-frame dipole_*.txt and warn_*.log files
# SCRIPT_DIR  : directory containing 1.aggregate.py

echo "Aggregate job starting at $(date)"
echo "OUTPUT_DIR: $OUTPUT_DIR"

python "$SCRIPT_DIR/1.aggregate.py" \
    --input-dir  "$OUTPUT_DIR" \
    --output     "$OUTPUT_DIR/dipole_output.txt" \
    --warn-log   "$OUTPUT_DIR/dipole_warnings.log"

echo "Aggregate job finished at $(date)"
