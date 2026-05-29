#!/bin/bash
# =============================================================================
# bbox_subset.slurm — LAMMPS trajectory bounding-box subset (parallel)
#
# Usage:
#   sbatch bbox_subset.slurm
#
# What this script does:
#   1. Counts total frames in INPUT.
#   2. Launches N_WORKERS (= $SLURM_NTASKS) parallel Python subset workers,
#      each handling a disjoint frame range and writing a part file.
#   3. Waits for all workers; checks each exit code individually.
#   4. Runs merge with V1-V7 validation in the same allocation.
#   5. (Optional) Removes part files after a successful merge.
#
# Worker parallelism note:
#   All workers share one node and each re-read the source file from the start,
#   skipping to their frame range with fast readline() calls (no coordinate
#   parsing). Acceptable on shared filesystems for typical trajectory sizes.
#   For multi-node parallelism, a Slurm array variant can be added later
#   without changing the Python tool.
#
# Edit the two sections marked ── EDIT ── before submitting.
# =============================================================================

# ── EDIT: Slurm resource directives ─────────────────────────────────────────
#SBATCH --job-name=bbox_subset
#SBATCH --nodes=1
#SBATCH --ntasks=8              # <── number of parallel workers; also sets N_WORKERS below
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=priya
#SBATCH --account=priyav_216
#SBATCH --output=bbox_subset_%j.log
#SBATCH --error=bbox_subset_%j.err
# ── end Slurm directives ────────────────────────────────────────────────────

# ── EDIT: job configuration ──────────────────────────────────────────────────
INPUT="/scratch1/lkyamamo/finalized-bubble-collapse/analysis/small_interface/0108_analysis/SiOH.custom"     # source trajectory (read-only)
OUT_DIR="/scratch1/lkyamamo/finalized-bubble-collapse/analysis/small_interface/0108_analysis"         # part files + merged output written here
FORMAT="lammps_custom"                               # auto | lammps_custom | lammps_xyz

# Bounding box — absolute Cartesian coordinates in Å
BBOX_XMIN=200.0
BBOX_XMAX=250.0
BBOX_YMIN=0.0
BBOX_YMAX=428.1
BBOX_ZMIN=0.0
BBOX_ZMAX=428.1

# Set to 1 to delete part files after a successful merge; 0 to keep them.
CLEANUP_PARTS=1
# ── end job configuration ────────────────────────────────────────────────────

# ── derived (do not edit below this line) ───────────────────────────────────
set -euo pipefail

N_WORKERS=$SLURM_NTASKS
EXT="${INPUT##*.}"
STEM="$(basename "${INPUT%.*}")"
OUTPUT="${OUT_DIR}/${STEM}_bbox.${EXT}"

# SLURM_SUBMIT_DIR is the directory from which sbatch was invoked — i.e. where
# the .slurm file and all sibling Python files (lammps_bbox_subset.py,
# frame_filter.py, formats/) live.  Slurm stages the script into its own spool
# directory, so BASH_SOURCE[0] would resolve there instead; do not use it.
# Submit from visualization/lammps_bbox_subset/ so SCRIPT_DIR is correct.
SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
SCRIPT="${SCRIPT_DIR}/lammps_bbox_subset.py"

module purge
module load python

echo "============================================================"
echo "bbox_subset — job ${SLURM_JOB_ID}"
echo "Input:    ${INPUT}"
echo "Output:   ${OUTPUT}"
echo "Out dir:  ${OUT_DIR}"
echo "Format:   ${FORMAT}"
echo "Workers:  ${N_WORKERS}"
echo "Bbox:     x[${BBOX_XMIN}, ${BBOX_XMAX}]  y[${BBOX_YMIN}, ${BBOX_YMAX}]  z[${BBOX_ZMIN}, ${BBOX_ZMAX}]"
echo "Start:    $(date)"
echo "============================================================"

mkdir -p "${OUT_DIR}"

# Count total frames (uses line-count arithmetic via wc -l — no parsing)
TOTAL=$(python3 "${SCRIPT}" count-frames --format "${FORMAT}" "${INPUT}")
echo "Total frames in source: ${TOTAL}"

CHUNK=$(( (TOTAL + N_WORKERS - 1) / N_WORKERS ))
echo "Chunk size per worker:  ${CHUNK}"

# ── parallel subset workers ──────────────────────────────────────────────────
pids=()
worker_info=()

for (( i=0; i<N_WORKERS; i++ )); do
    START=$(( i * CHUNK ))
    END=$(( START + CHUNK ))
    (( END > TOTAL )) && END=$TOTAL
    (( START >= END )) && continue      # last chunk may be empty if TOTAL % N_WORKERS == 0

    PART="${OUT_DIR}/part$(printf '%04d' ${i}).${EXT}"

    python3 "${SCRIPT}" subset "${INPUT}" \
        --format "${FORMAT}" \
        --frames "${START}" "${END}" \
        --xmin "${BBOX_XMIN}" --xmax "${BBOX_XMAX}" \
        --ymin "${BBOX_YMIN}" --ymax "${BBOX_YMAX}" \
        --zmin "${BBOX_ZMIN}" --zmax "${BBOX_ZMAX}" \
        --output "${PART}" &

    pids+=($!)
    worker_info+=("worker ${i} frames ${START}-${END} → $(basename "${PART}")")
    echo "Launched: ${worker_info[-1]} (PID ${pids[-1]})"
done

echo "Waiting for ${#pids[@]} workers..."

# Check each worker's exit code individually — do not use 'wait -n'
failed=0
for j in "${!pids[@]}"; do
    if wait "${pids[$j]}"; then
        echo "OK:     ${worker_info[$j]}"
    else
        echo "FAILED: ${worker_info[$j]} (PID ${pids[$j]})" >&2
        failed=1
    fi
done

if [ "${failed}" -ne 0 ]; then
    echo "One or more workers failed — aborting merge. Part files kept in ${OUT_DIR}." >&2
    exit 1
fi

echo "All workers completed successfully."

# ── merge with V1-V7 validation ──────────────────────────────────────────────
echo "Running merge..."

python3 "${SCRIPT}" merge \
    --format "${FORMAT}" \
    --parts "${OUT_DIR}"/part*.${EXT} \
    --output "${OUTPUT}" \
    --expected-parts "${N_WORKERS}"

echo "Merge complete: ${OUTPUT}"

# ── optional part cleanup ─────────────────────────────────────────────────────
if [ "${CLEANUP_PARTS}" -eq 1 ]; then
    echo "Removing part files..."
    rm -f "${OUT_DIR}"/part*.${EXT}
    echo "Part files removed."
else
    echo "Part files retained in ${OUT_DIR} (set CLEANUP_PARTS=1 to remove)."
fi

echo "============================================================"
echo "Done: $(date)"
echo "Output: ${OUTPUT}"
echo "============================================================"
