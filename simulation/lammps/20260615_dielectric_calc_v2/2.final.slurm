#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3800M
#SBATCH --time=01:00:00

# Runs after all calc_group array tasks complete.
# Concatenates group dipole files in order, then runs dipoleStd.
#
# Required: OUTPUT_DIR, NUM_GROUPS, SCRIPT_DIR, TEMPERATURE, LA, LB, LC,
#           AVERAGING_METHOD, VENV_PATH (optional)

set -euo pipefail

module purge
module load usc

if [ -n "$VENV_PATH" ]; then
    source "$VENV_PATH/bin/activate"
fi

echo "Final job starting at $(date)"
echo "OUTPUT_DIR: $OUTPUT_DIR"
echo "NUM_GROUPS: $NUM_GROUPS"

# 1. Require all group dipole outputs.
missing=()
for ((i=0; i<NUM_GROUPS; i++)); do
    if [ ! -f "$OUTPUT_DIR/groups/dipole_group_${i}.txt" ]; then
        missing+=("$i")
    fi
done
if [ ${#missing[@]} -gt 0 ]; then
    echo "ERROR: Missing group dipole output(s) for group index: ${missing[*]}" >&2
    exit 1
fi

# 2. Concatenate group files in order.
DIPOLE_FILE="$OUTPUT_DIR/dipole_output.txt"
echo "Concatenating ${NUM_GROUPS} group files..."
: > "$DIPOLE_FILE"
for ((i=0; i<NUM_GROUPS; i++)); do
    cat "$OUTPUT_DIR/groups/dipole_group_${i}.txt" >> "$DIPOLE_FILE"
done
echo "  Wrote: $DIPOLE_FILE ($(wc -l < "$DIPOLE_FILE") lines)"

# 3. Run dipoleStd.
echo ""
echo "Running dipoleStd on ${DIPOLE_FILE}"
echo "T=$TEMPERATURE  la=$LA  lb=$LB  lc=$LC  method=${AVERAGING_METHOD:-binned}"

python "$SCRIPT_DIR/2.dipole_std.py" \
    "$DIPOLE_FILE" \
    "$TEMPERATURE" \
    "$LA" \
    "$LB" \
    "$LC" \
    "${AVERAGING_METHOD:-binned}"

echo "Final job finished at $(date)"
echo "Dielectric CSV: $OUTPUT_DIR/dipole_output_timestep_data_${AVERAGING_METHOD:-binned}.csv"
