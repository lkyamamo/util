#!/bin/bash
#SBATCH --job-name=paraview_render
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem=240G
#SBATCH --time=04:00:00
#SBATCH --partition=compute
#SBATCH --output=render_%j.log
#SBATCH --error=render_%j.err

module purge
module load paraview/5.13.3-osmesa

# Directory where your dump files live
DUMP_DIR=/scratch/$USER/simulation/
SCRIPT_DIR=$SLURM_SUBMIT_DIR

echo "Starting render at $(date)"
echo "Using $SLURM_NTASKS MPI ranks"

mpirun -np $SLURM_NTASKS pvpython \
    $SCRIPT_DIR/render_movie.py \
    --dump-dir $DUMP_DIR \
    --output $SCRIPT_DIR/output/

echo "Finished at $(date)"