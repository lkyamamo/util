#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --partition=priya
#SBATCH --time=01:00:00
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --output=msd_%a.log
#SBATCH --job-name=MSD_batch
#SBATCH --constraint=epyc-7513
#SBATCH --mail-type=all
#SBATCH --mail-user=lkyamamo@usc.edu
#SBATCH --array=0-7

############################
# Edit for each submission
############################

TEMPERATURES=(15 30 45 60 75 90 120 150)

INPUT_DIR="../input_files"
DYNAMICS_DIR="/path/to/util/analysis/dynamics"

CORR_LENGTH=5000
CORR_INTERVAL=10
TIME_UNIT=2.5
MAX_FREQ=0.5

############################
# Resolve this task's temperature
############################

T="${TEMPERATURES[$SLURM_ARRAY_TASK_ID]}"

############################
# Resolve paths
############################

RUN_DIR="${SLURM_SUBMIT_DIR}"

if [[ "${INPUT_DIR}" != /* ]]; then
  INPUT_DIR="${RUN_DIR}/${INPUT_DIR}"
fi
INPUT_DIR="$(cd "${INPUT_DIR}" && pwd)"

if [[ "${DYNAMICS_DIR}" != /* ]]; then
  DYNAMICS_DIR="${RUN_DIR}/${DYNAMICS_DIR}"
fi
DYNAMICS_DIR="$(cd "${DYNAMICS_DIR}" && pwd)"
SRC_DIR="${DYNAMICS_DIR}/src"

INPUT="${INPUT_DIR}/${T}C.custom"
OUT_DIR="${RUN_DIR}/msd_${T}C"
XYZ="${OUT_DIR}/${T}C_calculation.xyz"

############################
# Validate
############################

if [[ ! -d "${SRC_DIR}" ]] || [[ ! -f "${SRC_DIR}/Makefile" ]]; then
  echo "ERROR: missing src/ build directory under ${DYNAMICS_DIR}" >&2
  exit 1
fi

if [[ ! -f "${INPUT}" ]]; then
  echo "ERROR: missing input file ${INPUT}" >&2
  exit 1
fi

############################
# Load modules and build
############################

module purge
module load legacy/CentOS7
module load usc gcc/9.2.0 boost
module load python

echo "Temperature:       ${T}C"
echo "Input file:        ${INPUT}"
echo "Output directory:  ${OUT_DIR}"
echo "CORR_LENGTH:       ${CORR_LENGTH}"
echo "CORR_INTERVAL:     ${CORR_INTERVAL}"
echo "TIME_UNIT:         ${TIME_UNIT}"
echo "MAX_FREQ:          ${MAX_FREQ}"
date

cd "${SRC_DIR}" || exit 1
make msd

if [[ ! -x "${SRC_DIR}/qmd_msd.exe" ]]; then
  echo "ERROR: qmd_msd.exe build failed" >&2
  exit 1
fi

export PYTHONPATH="${DYNAMICS_DIR}:${PYTHONPATH:-}"

############################
# Convert dump to XYZ
############################

mkdir -p "${OUT_DIR}"

python3 -c "
import lammps_custom_to_input_xyz as m
m.input_filename = '${INPUT}'
m.output_filename = '${XYZ}'
m.main()
"

if [[ ! -f "${XYZ}" ]]; then
  echo "ERROR: XYZ conversion failed, ${XYZ} not found" >&2
  exit 1
fi

############################
# Run MSD calculation
############################

cd "${OUT_DIR}" || exit 1

"${SRC_DIR}/qmd_msd.exe" "${XYZ}" \
  "${CORR_LENGTH}" "${CORR_INTERVAL}" "${TIME_UNIT}" "${MAX_FREQ}"

echo "Wrote ${OUT_DIR}/msd.dat and ${OUT_DIR}/dos.dat"
date