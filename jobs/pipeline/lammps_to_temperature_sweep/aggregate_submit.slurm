#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:05:00
#SBATCH --output=STREAM_OUTPUT
#SBATCH --job-name=dielectric-aggregate

# Combines one route's per-temperature results into a single CSV. Submitted
# once per route by submit_temperature_sweep.sh, with --dependency=afterok on
# every terminal job of that route: the dielectric calcs, or the msd jobs.
# AGGREGATE_SCRIPT is what makes it route-agnostic.
#
# ENTRIES_FILE: one "CELSIUS:KELVIN:PATH" per line (avoids the comma
# restriction of sbatch --export for an arbitrarily long list).
# OUTPUT_CSV / AGGREGATE_SCRIPT: set via sbatch --export.
# VENV_PATH: virtualenv to activate, set via --export from the driver, so the
# aggregators run under the same python as everything else in the pipeline.

set -euo pipefail

if [ -z "${ENTRIES_FILE:-}" ] || [ -z "${OUTPUT_CSV:-}" ] || [ -z "${AGGREGATE_SCRIPT:-}" ]; then
    echo "ERROR: ENTRIES_FILE, OUTPUT_CSV, and AGGREGATE_SCRIPT must be set via sbatch --export." >&2
    exit 1
fi

if [ -n "${VENV_PATH:-}" ]; then
    source "$VENV_PATH/bin/activate"
fi

entry_args=()
while IFS= read -r line; do
    [ -n "$line" ] && entry_args+=(--entry "$line")
done < "$ENTRIES_FILE"

python "$AGGREGATE_SCRIPT" "${entry_args[@]}" --output "$OUTPUT_CSV"
