#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=logs/dipoleStd_%j.out

# --- Required environment variables (set by submit_pipeline.sh via --export) ---
# OUTPUT_DIR       : directory containing dipole_output.txt (written by 1.aggregate.py)
# SCRIPT_DIR       : directory containing 2.dipoleStd.py
# TEMPERATURE      : simulation temperature in Kelvin
# LA, LB, LC       : box dimensions in Angstroms
# AVERAGING_METHOD : windowed | cumulative | hybrid | binned (default: binned)
# VENV_PATH        : (optional) path to the virtualenv to activate

# --- Environment setup ---
module purge
module load usc

if [ -n "$VENV_PATH" ]; then
    source "$VENV_PATH/bin/activate"
fi

DIPOLE_FILE="$OUTPUT_DIR/dipole_output.txt"

if [ ! -f "$DIPOLE_FILE" ]; then
    echo "ERROR: dipole_output.txt not found at $DIPOLE_FILE" >&2
    exit 1
fi

echo "dipoleStd job starting at $(date)"
echo "Input : $DIPOLE_FILE"
echo "T=$TEMPERATURE  la=$LA  lb=$LB  lc=$LC  method=${AVERAGING_METHOD:-binned}"

python "$SCRIPT_DIR/2.dipoleStd.py" \
    "$DIPOLE_FILE" \
    "$TEMPERATURE" \
    "$LA" \
    "$LB" \
    "$LC" \
    "${AVERAGING_METHOD:-binned}"

echo "dipoleStd job finished at $(date)"
