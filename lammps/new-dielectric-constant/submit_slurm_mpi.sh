#!/bin/bash
#SBATCH --job-name=dipole_calc_mpi
#SBATCH --time=08:00:00                # Max runtime
#SBATCH --nodes=1                     # Number of nodes
#SBATCH --ntasks=64            # MPI tasks per node
#SBATCH --cpus-per-task=1              # CPUs per task
#SBATCH --mem=0                      # Memory per task
#SBATCH --partition=priya             # Partition name (adjust for your cluster)
#SBATCH --constraint=epyc-7513

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

# Load MPI module - Use compatible OpenMPI version
module load legacy/CentOS7  gcc/12.3.0
module load openmpi/4.1.5 

source /apps/conda/miniforge3/24.11.3/etc/profile.d/conda.sh
conda activate /home1/lkyamamo/.conda/envs/analysis

echo "Testing MPI..."
srun -n 2 python -c "from mpi4py import MPI; print(f'MPI works! Rank {MPI.COMM_WORLD.Get_rank()} of {MPI.COMM_WORLD.Get_size()}')"

# Run the MPI script
echo "Starting Full MPI execution..."
echo "Using ${SLURM_NTASKS} MPI processes"

# Use srun with OpenMPI
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