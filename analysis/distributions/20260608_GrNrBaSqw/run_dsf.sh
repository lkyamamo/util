#!/bin/bash
# run_dsf.sh — local runner for the analysis pipeline
#
# Usage:
#   ./run_dsf.sh              — use all available cores
#   ./run_dsf.sh 8            — use 8 threads

############################
# Run flags — set 1 to run, 0 to skip
############################

RUN_DSF=1
RUN_RDF=1
RUN_BAD=1

############################
# Thread count
############################

N_THREADS=${1:-$(nproc)}
export OMP_NUM_THREADS=$N_THREADS
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
