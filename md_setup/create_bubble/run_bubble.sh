#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --partition=priya
#SBATCH --time=01:00:00
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --output=bubble_%j.out
#SBATCH --job-name=create_bubble
#SBATCH --constraint=epyc-7513
#SBATCH --mail-type=all
#SBATCH --mail-user=lkyamamo@usc.edu

set -e
ulimit -s unlimited

# ── config ────────────────────────────────────────────────────────────────────
INPUT_FILE="system.data"
OUTPUT_FILE="system_bubble.data"
SCRIPT_FILE="bubble.in"

CX=400.0; CY=400.0; CZ=881.0
RADIUS=150.0
SHELL_THICKNESS=2.0

O_TYPE=2
H_TYPE=3
SI_TYPE=1

export lmp="/home1/lkyamamo/executables/lammps/lmp_mpi_2019"
VENV_DIR="/home1/lkyamamo/util/md_setup/.venv"
# ─────────────────────────────────────────────────────────────────────────────
# Submit this script from the directory containing INPUT_FILE.
# All relative paths (INPUT_FILE, OUTPUT_FILE, SCRIPT_FILE) resolve against
# SLURM_SUBMIT_DIR, which is the cwd for all three steps.
# ─────────────────────────────────────────────────────────────────────────────

SCRIPT_DIR="$SLURM_SUBMIT_DIR"

echo "starting bubble creation **************************************"
date

module purge
module load usc
module load fftw

source "${VENV_DIR}/bin/activate"

# ── Step 1: generate the LAMMPS input script (serial) ────────────────────────
echo "Generating LAMMPS input script..."
python3 "${SCRIPT_DIR}/generate_lammps_script.py" \
    --input  "$INPUT_FILE"  \
    --output "$OUTPUT_FILE" \
    --script "$SCRIPT_FILE" \
    --cx "$CX" --cy "$CY" --cz "$CZ" \
    --radius "$RADIUS"

# ── Step 2: run LAMMPS in parallel ───────────────────────────────────────────
echo "Running LAMMPS..."
srun --mpi=pmix_v5 -n $SLURM_NTASKS $lmp \
    -log log.lammps \
    -in "$SCRIPT_FILE"

# ── Step 3: Python orphan atom cleanup (serial) ──────────────────────────────
echo "Running orphan atom cleanup..."
python3 "${SCRIPT_DIR}/cleanup_orphan_atoms.py" \
    --original "$INPUT_FILE"  \
    --data     "$OUTPUT_FILE" \
    --cx "$CX" --cy "$CY" --cz "$CZ" \
    --radius "$RADIUS" \
    --shell-thickness "$SHELL_THICKNESS" \
    --o-type "$O_TYPE" \
    --h-type "$H_TYPE" \
    --si-type "$SI_TYPE"

date
echo "bubble creation finished **************************************"

exit 0
