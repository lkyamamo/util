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
TIME_UNIT=2.5       # time between MD frames [fs]
MAX_FREQ=0.5        # max frequency for DOS [eV]

# Max temperatures to run at once (default: SLURM_NTASKS)
MAX_PARALLEL=""

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

MAX_PARALLEL="${MAX_PARALLEL:-${SLURM_NTASKS:-8}}"

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
echo "Max parallel temperatures:     ${MAX_PARALLEL}"

cd "${SRC_DIR}" || exit 1
make msd

if [[ ! -x "${SRC_DIR}/qmd_msd.exe" ]]; then
  echo "ERROR: qmd_msd.exe build failed" >&2
  exit 1
fi

export PYTHONPATH="${DYNAMICS_DIR}:${PYTHONPATH:-}"
export RUN_DIR INPUT_DIR SRC_DIR
export CORR_LENGTH CORR_INTERVAL TIME_UNIT MAX_FREQ

############################
# Per-temperature worker (runs in parallel)
############################

run_temperature() {
  local T="$1"
  local INPUT="${INPUT_DIR}/${T}C.custom"
  local OUT_DIR="${RUN_DIR}/msd_${T}C"
  local XYZ="${OUT_DIR}/${T}C_calculation.xyz"
  local LOG="${OUT_DIR}/run.log"

  mkdir -p "${OUT_DIR}"

  {
    echo "Temperature: ${T}C"
    echo "Input: ${INPUT}"
    echo "Output directory: ${OUT_DIR}"
    date

    python3 -c "
import lammps_custom_to_input_xyz as m
m.input_filename = '${INPUT}'
m.output_filename = '${XYZ}'
m.main()
"

    cd "${OUT_DIR}" || exit 1
    "${SRC_DIR}/qmd_msd.exe" "${XYZ}" \
      "${CORR_LENGTH}" "${CORR_INTERVAL}" "${TIME_UNIT}" "${MAX_FREQ}"

    echo "Wrote ${OUT_DIR}/msd.dat and ${OUT_DIR}/dos.dat"
    date
  } > "${LOG}" 2>&1
}

export -f run_temperature

############################
# Validate inputs, then run in parallel
############################

for T in "${TEMPERATURES[@]}"; do
  if [[ ! -f "${INPUT_DIR}/${T}C.custom" ]]; then
    echo "ERROR: missing input file ${INPUT_DIR}/${T}C.custom" >&2
    exit 1
  fi
done

PIDS=()
FAILED=0

for T in "${TEMPERATURES[@]}"; do
  while (( $(jobs -r | wc -l) >= MAX_PARALLEL )); do
    sleep 2
  done
  echo "Launching ${T}C (running: $(jobs -r | wc -l))"
  run_temperature "${T}" &
  PIDS+=($!)
done

for PID in "${PIDS[@]}"; do
  if ! wait "${PID}"; then
    FAILED=1
  fi
done

if (( FAILED )); then
  echo "ERROR: one or more temperatures failed (see msd_*/run.log)" >&2
  exit 1
fi

date
echo "MSD batch finished **************************************"

exit 0
