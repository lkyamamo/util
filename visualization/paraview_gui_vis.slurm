#!/bin/bash
#SBATCH --job-name=paraview_server
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=64
#SBATCH --cpus-per-task=1
#SBATCH --mem=0
#SBATCH --time=02:00:00
#SBATCH --partition=priya
#SBATCH --output=paraview_%j.log
#SBATCH --constraint=epyc-7513

module purge
module load paraview/5.11.1-server

# Print connection info to log file so you know where to tunnel
echo "============================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Node list: $SLURM_NODELIST"
echo "First node: $(echo $SLURM_NODELIST | cut -d',' -f1)"
echo "Tunnel command:"
echo "  ssh -L 11111:$(echo $SLURM_NODELIST | cut -d',' -f1):11111 $USER@<cluster_login_node>"
echo "============================================"

mpirun -np $SLURM_NTASKS pvserver \
    --server-port=11111 \
    --force-offscreen-rendering