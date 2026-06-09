#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --partition=priya
#SBATCH --time=04:00:00
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --output=STREAM_OUTPUT
#SBATCH --job-name=analysis
#SBATCH --constraint=epyc-7513
#SBATCH --mail-type=all
#SBATCH --mail-user=lkyamamo@usc.edu
ulimit -s unlimited
echo "starting analysis **************************************"
date

############################
# Run flags — set 1 to run, 0 to skip
############################

RUN_DSF=1
RUN_RDF=1
RUN_BAD=1

############################
# Paths
############################

RUN_DIR=$SLURM_SUBMIT_DIR
cd "$RUN_DIR"
echo "Running in $(pwd)"

############################
# Load modules + environment
############################

module purge
module load usc
conda activate dynasor_env        # env with dynasor, icc_rt, freud, numpy, matplotlib

############################
# Thread count
############################

# dsf.py (dynasor/numba) uses OMP_NUM_THREADS.
# rdf_freud.py and bad_freud.py (freud/TBB) use the same variable.
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
echo "Threads: $OMP_NUM_THREADS"

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

date
echo "analysis finished **************************************"
echo

exit 0
