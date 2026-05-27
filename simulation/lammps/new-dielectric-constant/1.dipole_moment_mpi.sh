#!/bin/bash
#
# MPI Bond Processing Script for Interactive Use
#
# Usage:
#   ./1.dipole_moment_mpi.sh
#
#   Modify the NUM_PROCS variable below to set the number of MPI processes
#

# Configuration - MODIFY THESE FOR YOUR SYSTEM
INPUT_FILE="../all_lammps.xyz"
START=0
END=10000000
INCREMENT=10
TYPE_A=1
TYPE_B=2
CUTOFF=1.2
LA=37.2514
LB=37.2514
LC=37.2514

# Output file
OUTPUT_FILE="dipole_output.txt"

# Number of MPI processes - MODIFY THIS TO CHANGE THE NUMBER OF PROCESSES
NUM_PROCS=64

# Number of CPU cores per task - MODIFY THIS TO CHANGE CORES PER TASK
CORES_PER_TASK=1

# Log job start
echo "=========================================="
echo "MPI Job starting at $(date)"
echo "Number of MPI processes: ${NUM_PROCS}"
echo "Cores per task: ${CORES_PER_TASK}"
echo "Output: ${OUTPUT_FILE}"
echo "=========================================="

# Load MPI module - IMPORTANT: Only load ONE MPI implementation!
# Option 1: Use OpenMPI
module purge
module load usc
module load openmpi/5.0.5

source ~/venvs/mpi-ompi5/bin/activate

# Run the MPI script
echo "Starting MPI execution..."

# Check if running under SLURM or directly
if [ -n "$SLURM_JOB_ID" ]; then
    # Running under SLURM - use srun
    echo "Running under SLURM (job ID: $SLURM_JOB_ID)"
    echo "Using $NUM_PROCS MPI tasks with $CORES_PER_TASK cores per task"
    srun -n ${NUM_PROCS} -c ${CORES_PER_TASK} --mpi=pmix python process_bonds_inline_3atom_mpi.py \
        ${INPUT_FILE} \
        ${START} \
        ${END} \
        ${INCREMENT} \
        ${TYPE_A} \
        ${TYPE_B} \
        ${CUTOFF} \
        ${LA} \
        ${LB} \
        ${LC} \
        ${OUTPUT_FILE}
else
    # Running directly - use mpirun
    echo "Running directly (interactive mode)"
    echo "Using $NUM_PROCS MPI processes with $CORES_PER_TASK cores per task"
    mpirun -np ${NUM_PROCS} --bind-to core --map-by core:PE=${CORES_PER_TASK} python process_bonds_inline_3atom_mpi_optimized.py \
        ${INPUT_FILE} \
        ${START} \
        ${END} \
        ${INCREMENT} \
        ${TYPE_A} \
        ${TYPE_B} \
        ${CUTOFF} \
        ${LA} \
        ${LB} \
        ${LC} \
        ${OUTPUT_FILE}
fi

# Check exit status
EXIT_CODE=$?
if [ $EXIT_CODE -eq 0 ]; then
    echo "=========================================="
    echo "MPI Job completed successfully at $(date)"
    echo "=========================================="
else
    echo "=========================================="
    echo "MPI Job failed with exit code ${EXIT_CODE} at $(date)"
    echo "=========================================="
    exit $EXIT_CODE
fi
