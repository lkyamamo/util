#!/bin/bash
# distribution_run.sh — local runner for the analysis pipeline (all cores,
# no Slurm); mirrors distribution_submit.slurm's run flags/thread wiring.
#
# Usage:
#   ./distribution_run.sh              — use all available cores
#   ./distribution_run.sh 8            — use 8 threads

############################
# Run flags — set 1 to run, 0 to skip
############################

RUN_DSF=1
RUN_RDF=1
RUN_BAD=1
RUN_VDOS=1

############################
# Thread count
############################

N_THREADS=${1:-$(nproc)}
export OMP_NUM_THREADS=$N_THREADS
# vdos.py's fft_periodogram method reads VDOS_THREADS directly (scipy.fft
# workers); its default vacf_cosine_transform method instead benefits from
# OMP_NUM_THREADS above via numpy's underlying BLAS matmul.
export VDOS_THREADS=$N_THREADS
echo "Threads: $OMP_NUM_THREADS"

############################
# Load environment
############################

source /home1/lkyamamo/venv/struc_analysis/bin/activate

############################
# Run scripts
############################

if [ "$RUN_DSF" -eq 1 ]; then
    echo "--- dsf.py ---"
    python dsf.py
fi

if [ "$RUN_RDF" -eq 1 ]; then
    echo "--- rdf_freud.py ---"
    python rdf_freud.py
fi

if [ "$RUN_BAD" -eq 1 ]; then
    echo "--- bad_freud.py ---"
    python bad_freud.py
fi

if [ "$RUN_VDOS" -eq 1 ]; then
    echo "--- vdos.py ---"
    python vdos.py
fi
