#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --partition=priya
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:05:00
#SBATCH --output=STREAM_OUTPUT
#SBATCH --job-name=dielectric-aggregate

# Combines all per-temperature dielectric_submit.slurm summaries into one
# CSV. Submitted by submit_dielectric_pipeline.sh with --dependency=afterok
# on every per-temperature stage-2 job.
#
# ENTRIES_FILE: one "TEMPERATURE:SUMMARY_FILE" per line (avoids the comma
# restriction of sbatch --export for an arbitrarily long list).
# OUTPUT_CSV / AGGREGATE_SCRIPT: set via sbatch --export.

set -euo pipefail

if [ -z "${ENTRIES_FILE:-}" ] || [ -z "${OUTPUT_CSV:-}" ] || [ -z "${AGGREGATE_SCRIPT:-}" ]; then
    echo "ERROR: ENTRIES_FILE, OUTPUT_CSV, and AGGREGATE_SCRIPT must be set via sbatch --export." >&2
    exit 1
fi

entry_args=()
while IFS= read -r line; do
    [ -n "$line" ] && entry_args+=(--entry "$line")
done < "$ENTRIES_FILE"

python3 "$AGGREGATE_SCRIPT" "${entry_args[@]}" --output "$OUTPUT_CSV"
