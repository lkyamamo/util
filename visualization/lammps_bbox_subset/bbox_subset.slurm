#!/bin/bash
# =============================================================================
# bbox_subset.slurm — LAMMPS trajectory bounding-box subset (parallel)
#
# Usage:
#   cd visualization/lammps_bbox_subset/
#   sbatch bbox_subset.slurm
#
# What this script does:
#   1. Builds a binary byte-offset frame index (source.idx) — O(1) seeking.
#   2. Counts total frames (TOTAL) in the source trajectory.
#   3. Launches N_WORKERS parallel subset workers (round-robin assignment):
#        worker W processes frames W, W+N, W+2N, …
#        each frame written to FRAMES_DIR/frame_{i:07d}.{ext}
#   4. Waits for all workers; checks exit codes individually.
#   5. Merges: sorts frame files by index, validates, concatenates → OUTPUT.
#   6. (Optional) Removes per-frame files after a successful merge.
#
# Parallelism model:
#   Workers are assigned frames mod N_WORKERS (round-robin).  Because the
#   frame index already gives every worker O(1) seek to any frame, there is
#   no sequential skip overhead regardless of frame number.
#
# Submit from visualization/lammps_bbox_subset/ so $SLURM_SUBMIT_DIR resolves
# to the directory containing lammps_bbox_subset.py and its siblings.
# =============================================================================

# ── EDIT: Slurm resource directives ─────────────────────────────────────────
#SBATCH --job-name=bbox_subset
#SBATCH --nodes=1
#SBATCH --ntasks=8              # <── number of parallel workers
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --partition=priya
#SBATCH --account=priyav_216
#SBATCH --output=bbox_subset_%j.log
#SBATCH --error=bbox_subset_%j.err
# ── end Slurm directives ────────────────────────────────────────────────────

# ── EDIT: job configuration ──────────────────────────────────────────────────
INPUT="/scratch1/lkyamamo/finalized-bubble-collapse/analysis/small_interface/0108_analysis/SiOH.custom"
OUT_DIR="/scratch1/lkyamamo/finalized-bubble-collapse/analysis/small_interface/0108_analysis"
FORMAT="lammps_custom"                # auto | lammps_custom | lammps_xyz

# Bounding box — absolute Cartesian coordinates in Å
BBOX_XMIN=200.0
BBOX_XMAX=250.0
BBOX_YMIN=0.0
BBOX_YMAX=428.1
BBOX_ZMIN=0.0
BBOX_ZMAX=428.1

# Set to 1 to delete per-frame files after a successful merge; 0 to keep them.
CLEANUP_PARTS=1
# ── end job configuration ────────────────────────────────────────────────────

set -euo pipefail

N_WORKERS=$SLURM_NTASKS
EXT="${INPUT##*.}"
STEM="$(basename "${INPUT%.*}")"
OUTPUT="${OUT_DIR}/${STEM}_bbox.${EXT}"

# SLURM_SUBMIT_DIR is the directory from which sbatch was invoked.  Slurm
# stages the script into a spool directory, so BASH_SOURCE[0] would resolve
# there.  Submit from visualization/lammps_bbox_subset/ so SCRIPT_DIR is set
# correctly and all Python siblings (lammps_bbox_subset.py, formats/, …) are
# reachable.
SCRIPT_DIR="${SLURM_SUBMIT_DIR}"
SCRIPT="${SCRIPT_DIR}/lammps_bbox_subset.py"

# Per-frame output files live in a subdirectory to keep OUT_DIR tidy.
FRAMES_DIR="${OUT_DIR}/frames_${SLURM_JOB_ID}"

module purge
module load python

echo "============================================================"
echo "bbox_subset — job ${SLURM_JOB_ID}"
echo "Input:    ${INPUT}"
echo "Output:   ${OUTPUT}"
echo "Frames:   ${FRAMES_DIR}"
echo "Format:   ${FORMAT}"
echo "Workers:  ${N_WORKERS}"
echo "Bbox:     x[${BBOX_XMIN}, ${BBOX_XMAX}]  y[${BBOX_YMIN}, ${BBOX_YMAX}]  z[${BBOX_ZMIN}, ${BBOX_ZMAX}]"
echo "Start:    $(date)"
echo "============================================================"

mkdir -p "${OUT_DIR}" "${FRAMES_DIR}"

# ── Step 1: build byte-offset frame index ────────────────────────────────────
# The index is a compact binary file of uint64 byte offsets — one per frame.
# Building it requires a full sequential pass over the source file (O(file_size)),
# but the source trajectory never changes, so we store the index next to the
# source file and reuse it across all subsequent bbox runs.  The bottleneck is
# paid exactly once.
#
# Naming: ${INPUT}.idx  — stable across jobs, keyed to the source file.
# To force a rebuild (e.g. after the source changes), delete this file first.
INDEX="${INPUT}.idx"

if [ -f "${INDEX}" ]; then
    echo "Reusing existing frame index: ${INDEX}"
else
    echo "Building frame index (one-time cost for this source file)..."
    python3 "${SCRIPT}" build-index --format "${FORMAT}" "${INPUT}" --output "${INDEX}"
fi

# ── Step 2: count total frames ────────────────────────────────────────────────
# count-frames reads the pre-built index (O(1) — just byte-counts the .idx file)
# rather than re-scanning the source.
TOTAL=$(python3 "${SCRIPT}" count-frames --format "${FORMAT}" "${INPUT}")
echo "Total frames: ${TOTAL}"

# ── Step 3: launch parallel round-robin workers ──────────────────────────────
# Worker W writes frames W, W+N_WORKERS, W+2*N_WORKERS, … each to its own
# file: FRAMES_DIR/frame_{i:07d}.{ext}.  No chunk math required.
pids=()
worker_info=()

for (( W=0; W<N_WORKERS; W++ )); do
    python3 "${SCRIPT}" subset "${INPUT}" \
        --format  "${FORMAT}" \
        --frame-index "${INDEX}" \
        --worker  "${W}" \
        --num-workers "${N_WORKERS}" \
        --output-dir  "${FRAMES_DIR}" \
        --xmin "${BBOX_XMIN}" --xmax "${BBOX_XMAX}" \
        --ymin "${BBOX_YMIN}" --ymax "${BBOX_YMAX}" \
        --zmin "${BBOX_ZMIN}" --zmax "${BBOX_ZMAX}" &

    pids+=($!)
    worker_info+=("worker ${W}/${N_WORKERS}")
    echo "Launched: ${worker_info[-1]}  (PID ${pids[-1]})"
done

echo "Waiting for ${#pids[@]} workers..."

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
    echo "One or more workers failed — aborting merge. Frame files kept in ${FRAMES_DIR}." >&2
    exit 1
fi

echo "All workers completed successfully."

# ── Step 4: merge ────────────────────────────────────────────────────────────
# Frame files are named frame_{i:07d}.{ext}; sorting by filename gives the
# correct temporal order.  Merge validates file count (V1), non-empty files
# (V2), and strict timestep monotonicity (V7) of the merged output.
echo "Running merge..."

python3 "${SCRIPT}" merge \
    --format "${FORMAT}" \
    --parts-dir "${FRAMES_DIR}" \
    --output    "${OUTPUT}" \
    --expected-frames "${TOTAL}"

echo "Merge complete: ${OUTPUT}"

# ── Step 5: optional cleanup ─────────────────────────────────────────────────
# Per-frame files are transient — always safe to remove after a successful merge.
# The index (${INPUT}.idx) is intentionally kept: it is reused by future runs
# against the same source trajectory, so the O(file_size) build cost is paid once.
if [ "${CLEANUP_PARTS}" -eq 1 ]; then
    echo "Removing per-frame files..."
    rm -rf "${FRAMES_DIR}"
    echo "Done. Index retained at: ${INDEX}"
else
    echo "Per-frame files retained in ${FRAMES_DIR} (set CLEANUP_PARTS=1 to remove)."
fi

echo "============================================================"
echo "Done: $(date)"
echo "Output: ${OUTPUT}"
echo "============================================================"
