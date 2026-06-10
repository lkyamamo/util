#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
#SBATCH --output=logs/dipole_%A_%a.out

# --- Required environment variables (set by submit_pipeline.sh via --export) ---
# PENDING_FILE : path to pending_dumps.txt listing unprocessed dump files (one per line)
# DUMP_GLOB    : glob pattern (e.g. "dump.*.lammpstrj") — used for timestep extraction
# OUTPUT_DIR   : directory where per-frame dipole_*.txt and warn_*.log are written
# SCRIPT_DIR   : directory containing this script and 1.calc_dipole_from_dump.py
# CUTOFF       : O-H bond cutoff in Angstroms (e.g. 1.2)
# TYPE_O       : LAMMPS atom type integer for oxygen (e.g. 1)
# TYPE_H       : LAMMPS atom type integer for hydrogen (e.g. 2)
# VENV_PATH    : (optional) path to the uv/virtualenv to activate

# --- Environment setup ---
module purge
module load usc

if [ -n "$VENV_PATH" ]; then
    source "$VENV_PATH/bin/activate"
fi

mkdir -p "$OUTPUT_DIR"

# Read dump path for this task from the pre-filtered pending list
DUMP=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" "$PENDING_FILE")

if [ -z "$DUMP" ]; then
    echo "ERROR: No dump file found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
    exit 1
fi

# Extract timestep from filename using the glob prefix/suffix
DUMP_PREFIX="${DUMP_GLOB%%\**}"
DUMP_SUFFIX="${DUMP_GLOB##*\*}"
BASENAME=$(basename "$DUMP")
TIMESTEP="${BASENAME#$DUMP_PREFIX}"
TIMESTEP="${TIMESTEP%$DUMP_SUFFIX}"

echo "Task ${SLURM_ARRAY_TASK_ID}: processing ${DUMP} (timestep ${TIMESTEP})"

python "$SCRIPT_DIR/1.calc_dipole_from_dump.py" "$DUMP" \
    --cutoff "$CUTOFF" \
    --type-o "$TYPE_O" \
    --type-h "$TYPE_H" \
    --output "$OUTPUT_DIR/dipole_${TIMESTEP}.txt" \
    --log    "$OUTPUT_DIR/warn_${TIMESTEP}.log"
