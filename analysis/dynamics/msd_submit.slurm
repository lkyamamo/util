#!/bin/bash
#SBATCH --account=priyav_216
#SBATCH --nodes=1
#SBATCH --ntasks=64
#SBATCH --partition=priya
#SBATCH --time=01:00:00
#SBATCH --mem=0
#SBATCH --exclusive
#SBATCH --output=STREAM_OUTPUT
#SBATCH --job-name=MSD_batch
#SBATCH --constraint=epyc-7513
#SBATCH --mail-type=all
#SBATCH --mail-user=lkyamamo@usc.edu

ulimit -s unlimited
echo "starting MSD batch **************************************"
date

############################
# Edit for each submission
# Input dumps: {T}C.custom in INPUT_DIR (e.g. 30C.custom)
# Relative paths are resolved from the submit directory ($SLURM_SUBMIT_DIR)
############################

INPUT_DIR="../input_files"
DYNAMICS_DIR="/path/to/util/analysis/dynamics"

TEMPERATURES=(
  15
  30
  45
  60
  75
  90
  120
  150
)

CORR_LENGTH=5000    # correlation length [fs]
CORR_INTERVAL=10    # interval between reference frames [fs]
TIME_UNIT=2.5         # time between MD frames [fs]
MAX_FREQ=0.5          # max frequency for DOS [eV]

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

if [[ ! -d "${SRC_DIR}" ]] || [[ ! -f "${SRC_DIR}/Makefile" ]]; then
  echo "ERROR: missing src/ build directory under ${DYNAMICS_DIR}" >&2
  exit 1
fi

############################
# Load modules and build
############################

module purge
module load legacy/CentOS7
module load usc gcc/9.2.0 boost
module load python

echo "Run directory (outputs):       ${RUN_DIR}"
echo "Input directory ({T}C.custom): ${INPUT_DIR}"
echo "Dynamics directory:            ${DYNAMICS_DIR}"

cd "${SRC_DIR}" || exit 1
make msd

if [[ ! -x "${SRC_DIR}/qmd_msd.exe" ]]; then
  echo "ERROR: qmd_msd.exe build failed" >&2
  exit 1
fi

export PYTHONPATH="${DYNAMICS_DIR}:${PYTHONPATH:-}"

############################
# Process each temperature
############################

for T in "${TEMPERATURES[@]}"; do
  INPUT="${INPUT_DIR}/${T}C.custom"
  OUT_DIR="${RUN_DIR}/msd_${T}C"
  XYZ="${OUT_DIR}/${T}C_calculation.xyz"
  MSD_OUT="${OUT_DIR}/msd.dat"
  DOS_OUT="${OUT_DIR}/dos.dat"

  echo "----------------------------------------"
  echo "Temperature: ${T}C"
  echo "Input: ${INPUT}"
  echo "Output directory: ${OUT_DIR}"

  if [[ ! -f "${INPUT}" ]]; then
    echo "ERROR: missing input file ${INPUT}" >&2
    exit 1
  fi

  mkdir -p "${OUT_DIR}"

  python3 -c "
import lammps_custom_to_input_xyz as m
m.input_filename = '${INPUT}'
m.output_filename = '${XYZ}'
m.main()
"

  # qmd_msd.exe always writes msd.dat and dos.dat to cwd — run in per-T dir
  (
    cd "${OUT_DIR}" || exit 1
    "${SRC_DIR}/qmd_msd.exe" "${XYZ}" \
      "${CORR_LENGTH}" "${CORR_INTERVAL}" "${TIME_UNIT}" "${MAX_FREQ}"
  )

  echo "Wrote ${MSD_OUT} and ${DOS_OUT}"
done

date
echo "MSD batch finished **************************************"

exit 0
