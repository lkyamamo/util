#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/dipole_%A_%a.out

# --- Required environment variables (set by submit_pipeline.sh via --export) ---
# DUMP_DIR     : directory containing dump files
# DUMP_GLOB    : glob pattern matching dump files (e.g. "dump.*.lammpstrj")
# OUTPUT_DIR   : directory where per-frame dipole_*.txt and warn_*.log are written
# SCRIPT_DIR   : directory containing this script and 1.calc_dipole_from_dump.py
# CUTOFF       : O-H bond cutoff in Angstroms (e.g. 1.2)
# TYPE_O       : LAMMPS atom type integer for oxygen (e.g. 1)
# TYPE_H       : LAMMPS atom type integer for hydrogen (e.g. 2)

mkdir -p "$OUTPUT_DIR"

# Map SLURM_ARRAY_TASK_ID to a dump file (sort ensures consistent ordering)
DUMP=$(ls "$DUMP_DIR"/$DUMP_GLOB | sort | sed -n "$((SLURM_ARRAY_TASK_ID + 1))p")

if [ -z "$DUMP" ]; then
    echo "ERROR: No dump file found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

echo "Task ${SLURM_ARRAY_TASK_ID}: processing ${DUMP}"

python "$SCRIPT_DIR/1.calc_dipole_from_dump.py" "$DUMP" \
    --cutoff "$CUTOFF" \
    --type-o "$TYPE_O" \
    --type-h "$TYPE_H" \
    --output "$OUTPUT_DIR/dipole_${SLURM_ARRAY_TASK_ID}.txt" \
    --log    "$OUTPUT_DIR/warn_${SLURM_ARRAY_TASK_ID}.log"
