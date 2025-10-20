#!/bin/bash

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

# Log job start
echo "=========================================="
echo "MPI Job starting at $(date)"
echo "Nodes: ${SLURM_NNODES}"
echo "Total tasks: ${SLURM_NTASKS}"
echo "Tasks per node: ${SLURM_TASKS_PER_NODE}"
echo "CPUs per task: ${SLURM_CPUS_PER_TASK}"
echo "Output: ${OUTPUT_FILE}"
echo "=========================================="

source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate /home1/lkyamamo/.conda/envs/analysis

# Load MPI module - IMPORTANT: Only load ONE MPI implementation!
# Option 1: Use OpenMPI (recommended for SLURM)
module load ver/2506  gcc/14.3.0
module load openmpi/5.0.8

# Option 2: Use MVAPICH2 (uncomment if you prefer)
# module load legacy/CentOS7
# module load mvapich2/2.3.7

# Run the MPI script
echo "Starting MPI execution..."

# Check if running under SLURM or directly
if [ -n "$SLURM_JOB_ID" ]; then
    # Running under SLURM - use srun
    echo "Running under SLURM (job ID: $SLURM_JOB_ID)"
    srun python process_bonds_inline_3atom_mpi.py \
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
    echo "Running directly (not under SLURM)"
    NUM_PROCS=${NUM_PROCS:-1}  # Default to 1 process if not set
    echo "Using $NUM_PROCS MPI processes"
    mpirun -np ${NUM_PROCS} python process_bonds_inline_3atom_mpi.py \
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