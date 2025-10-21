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
TEMPLATE="303"

# Output file (use absolute path to avoid working directory issues)
OUTPUT_FILE="all_lammps.xyz"

# Note: Let the Python script handle its own path resolution
echo "SLURM working directory: $(pwd)"

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
echo "Template prefix: ${TEMPLATE}"
echo "Output file: ${OUTPUT_FILE}"
echo "Timestep range: ${START} to ${END} (step ${INCREMENT})"
echo "Dumps directory should be: $(dirname "$(pwd)")/dumps"
echo "Checking if dumps directory exists:"
if [ -d "../dumps" ]; then
    echo "✓ Dumps directory exists"
    echo "Checking for sample 303.* files:"
    # Use find with -maxdepth 1 and limit results to avoid performance issues
    find ../dumps -maxdepth 1 -name "303.*" -type f | head -3 | while read file; do
        echo "  Found: $(basename "$file")"
    done
    # Count total 303.* files efficiently
    file_count=$(find ../dumps -maxdepth 1 -name "303.*" -type f | wc -l)
    echo "  Total 303.* files found: $file_count"
else
    echo "✗ Dumps directory not found"
fi
echo "=========================================="

srun --mpi=pmix -n $SLURM_NTASKS python old_to_new_mpi.py \
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

