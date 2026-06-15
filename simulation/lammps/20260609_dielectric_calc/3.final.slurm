#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00

# Runs after all chunk aggregate jobs complete.
# Required: OUTPUT_DIR, NUM_CHUNKS, SCRIPT_DIR, TEMPERATURE, LA, LB, LC,
#           AVERAGING_METHOD, VENV_PATH (optional)

set -euo pipefail

module purge
module load usc

if [ -n "$VENV_PATH" ]; then
    source "$VENV_PATH/bin/activate"
fi

echo "Final job starting at $(date)"
echo "OUTPUT_DIR: $OUTPUT_DIR"
echo "NUM_CHUNKS: $NUM_CHUNKS"

# 1. Require all chunk dipole outputs; report any that are missing.
missing=()
for ((i=0; i<NUM_CHUNKS; i++)); do
    if [ ! -f "$OUTPUT_DIR/dipole_output_chunk_${i}.txt" ]; then
        missing+=("$i")
    fi
done
if [ ${#missing[@]} -gt 0 ]; then
    echo "ERROR: Missing chunk dipole output(s) for chunk index: ${missing[*]}" >&2
    exit 1
fi

# 2. Combine chunk dipole files.
echo "Combining ${NUM_CHUNKS} chunk dipole files..."
: > "$OUTPUT_DIR/dipole_output.txt"
for ((i=0; i<NUM_CHUNKS; i++)); do
    cat "$OUTPUT_DIR/dipole_output_chunk_${i}.txt" >> "$OUTPUT_DIR/dipole_output.txt"
done
echo "  Wrote: $OUTPUT_DIR/dipole_output.txt ($(wc -l < "$OUTPUT_DIR/dipole_output.txt") lines)"

# 3. Combine chunk warning logs.
echo "Combining ${NUM_CHUNKS} chunk warning logs..."
: > "$OUTPUT_DIR/dipole_warnings.log"
for ((i=0; i<NUM_CHUNKS; i++)); do
    warn_chunk="$OUTPUT_DIR/dipole_warnings_chunk_${i}.log"
    if [ -f "$warn_chunk" ]; then
        cat "$warn_chunk" >> "$OUTPUT_DIR/dipole_warnings.log"
    fi
done
echo "  Wrote: $OUTPUT_DIR/dipole_warnings.log"

# 4. Run dipoleStd.
DIPOLE_FILE="$OUTPUT_DIR/dipole_output.txt"
echo ""
echo "Running dipoleStd on ${DIPOLE_FILE}"
echo "T=$TEMPERATURE  la=$LA  lb=$LB  lc=$LC  method=${AVERAGING_METHOD:-binned}"

python "$SCRIPT_DIR/3.dipoleStd.py" \
    "$DIPOLE_FILE" \
    "$TEMPERATURE" \
    "$LA" \
    "$LB" \
    "$LC" \
    "${AVERAGING_METHOD:-binned}"

echo "Final job finished at $(date)"
echo "Dielectric CSV: $OUTPUT_DIR/dipole_output_timestep_data_${AVERAGING_METHOD:-binned}.csv"
