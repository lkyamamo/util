#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition=priya
#SBATCH --time=24:00:00
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --output=msd_%a.log
#SBATCH --job-name=MSD_batch
#SBATCH --constraint=epyc-7513
#SBATCH --mail-type=all
#SBATCH --mail-user=lkyamamo@usc.edu
#SBATCH --array=0-7

TEMPERATURES=(15 30 45 60 75 90 120 150)
T="${TEMPERATURES[$SLURM_ARRAY_TASK_ID]}"

INPUT_DIR="../input_files"
DYNAMICS_DIR="/path/to/util/analysis/dynamics"

RUN_DIR="${SLURM_SUBMIT_DIR}"
INPUT="${INPUT_DIR}/${T}C.custom"
OUT_DIR="${RUN_DIR}/msd_${T}C"
XYZ="${OUT_DIR}/${T}C_calculation.xyz"
SRC_DIR="${DYNAMICS_DIR}/src"

mkdir -p "${OUT_DIR}"

module purge
module load legacy/CentOS7
module load usc gcc/9.2.0 boost
module load python

cd "${SRC_DIR}" || exit 1
make msd

export PYTHONPATH="${DYNAMICS_DIR}:${PYTHONPATH:-}"

python3 -c "
import lammps_custom_to_input_xyz as m
m.input_filename = '${INPUT}'
m.output_filename = '${XYZ}'
m.main()
"

cd "${OUT_DIR}" || exit 1
"${SRC_DIR}/qmd_msd.exe" "${XYZ}" \
  "${CORR_LENGTH}" "${CORR_INTERVAL}" "${TIME_UNIT}" "${MAX_FREQ}"