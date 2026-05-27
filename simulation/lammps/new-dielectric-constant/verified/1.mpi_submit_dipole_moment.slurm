#!/bin/bash
#SBATCH --job-name=dipole_calc_harmonic_pot
#SBATCH --time=12:00:00          	# Max runtime
#SBATCH --nodes=4	                # Number of nodes
#SBATCH --ntasks=256   	       	# MPI tasks per node
#SBATCH --cpus-per-task=1              	# CPUs per task
#SBATCH --ntasks-per-node=64
#SBATCH --mem=0                      # Memory per task
#SBATCH --partition=priya             # Partition name (adjust for your cluster)
#SBATCH --constraint=epyc-7513
#SBATCH --exclusive

# Configuration - MODIFY THESE FOR YOUR SYSTEM
INPUT_FILE="../all_lammps.xyz"
START=270000
END=10270000
INCREMENT=10
TYPE_A=1
TYPE_B=2
CUTOFF=1.2
LA=37.2514
LB=37.2514
LC=37.2514

# Output file
OUTPUT_FILE="dipole_output_mpi.txt"

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
srun --mpi=pmix -n $SLURM_NTASKS python process_bonds_inline_3atom_mpi.py \
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

