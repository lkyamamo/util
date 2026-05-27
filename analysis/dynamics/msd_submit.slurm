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
# Temperature list (edit)
# Input files: {T}C.custom  (e.g. 30C.custom, 40C.custom)
############################

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

############################
# MSD parameters (same units as run.sh)
############################

CORR_LENGTH=5000
CORR_INTERVAL=10
TIME_UNIT=2.5
MAX_FREQ=0.5

############################
# Load modules and build
############################

module purge
module load legacy/CentOS7
module load usc gcc/9.2.0 boost
module load python

cd "${SLURM_SUBMIT_DIR}" || exit 1
echo "Working directory: $(pwd)"

make msd

############################
# Process each temperature
############################

for T in "${TEMPERATURES[@]}"; do
  INPUT="${T}C.custom"
  XYZ="${T}C_calculation.xyz"
  MSD_OUT="msd_${T}C.dat"
  DOS_OUT="dos_${T}C.dat"

  echo "----------------------------------------"
  echo "Temperature: ${T}C"
  echo "Input: ${INPUT}"

  if [[ ! -f "${INPUT}" ]]; then
    echo "ERROR: missing input file ${INPUT}" >&2
    exit 1
  fi

  python3 -c "
import lammps_custom_to_input_xyz as m
m.input_filename = '${INPUT}'
m.output_filename = '${XYZ}'
m.main()
"

  ./qmd_msd.exe "${XYZ}" "${CORR_LENGTH}" "${CORR_INTERVAL}" "${TIME_UNIT}" "${MAX_FREQ}"

  mv -f msd.dat "${MSD_OUT}"
  mv -f dos.dat "${DOS_OUT}"

  echo "Wrote ${MSD_OUT} and ${DOS_OUT}"
done

date
echo "MSD batch finished **************************************"

exit 0
