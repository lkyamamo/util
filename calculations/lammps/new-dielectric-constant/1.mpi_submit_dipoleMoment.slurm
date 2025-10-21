#!/bin/bash
#SBATCH --job-name=dipole_calc_harmonic_pot
#SBATCH --time=6:00:00          	# Max runtime
#SBATCH --nodes=2	                # Number of nodes
#SBATCH --ntasks=128    	       	# MPI tasks per node
#SBATCH --cpus-per-task=1              	# CPUs per task
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0                      # Memory per task
#SBATCH --partition=priya             # Partition name (adjust for your cluster)
#SBATCH --constraint=epyc-7513
#SBATCH --exclusive

# Configuration - MODIFY THESE FOR YOUR SYSTEM
START=350000
END=10350000
INCREMENT=100
TEMPLATE="../dumps/303"

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

module purge
module load usc
module load openmpi/5.0.5
source ~/venvs/mpi-ompi5/bin/activate

which python

srun --mpi=pmix -n $SLURM_NTASKS python -c "from mpi4py import MPI; print(f'MPI works! Rank {MPI.COMM_WORLD.Get_rank()} of {MPI.COMM_WORLD.Get_size()}')"

# Run the MPI script
echo "Starting MPI execution..."
srun --mpi=pmix -n $SLURM_NTASKS python 1.dipoleMoment_mpi.py \
    ${START} \
    ${END} \
    ${INCREMENT} \
    ${OUTPUT_FILE} \
    ${TEMPLATE}

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

