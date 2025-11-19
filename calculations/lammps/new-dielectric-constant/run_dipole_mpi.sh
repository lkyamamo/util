#!/bin/bash
# Shell script to run the MPI version of 1.dipoleMoment.py

# Check if help is requested first
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    # Show help without checking dependencies
    echo "Usage: $0 [OPTIONS]"
    echo ""
    echo "Options:"
    echo "  -np, --processes NUM    Number of MPI processes (default: 4)"
    echo "  -s, --start NUM         Start timestep (default: 0)"
    echo "  -e, --end NUM           End timestep (default: 1000)"
    echo "  -i, --increment NUM     Timestep increment (default: 10)"
    echo "  -o, --output FILE       Output file (default: dipole_output.txt)"
    echo "  -p, --prefix STRING     File prefix (default: 298)"
    echo "  -h, --help              Show this help message"
    echo ""
    echo "Examples:"
    echo "  $0 -np 64 -s 0 -e 10000 -i 10 -o my_dipole_output.txt"
    echo "  $0 -np 64 -s 0 -e 10000 -i 10 -o my_dipole_output.txt -p 303"
    echo "  $0 -np 64 -s 0 -e 10000 -i 10 -o my_dipole_output.txt -p 273.15"
    echo ""
    echo "File naming:"
    echo "  Input files: {prefix}.{timestep} (e.g., 303.0, 303.10, 303.20, ...)"
    echo ""
    echo "Requirements:"
    echo "  - mpi4py: pip install mpi4py"
    echo "  - MPI: OpenMPI or MPICH"
    exit 0
fi

# Check if mpi4py is available
python -c "import mpi4py" 2>/dev/null
if [ $? -ne 0 ]; then
    echo "Error: mpi4py is not installed!"
    echo "Please install it with: pip install mpi4py"
    exit 1
fi

# Check if mpirun is available
if ! command -v mpirun &> /dev/null; then
    echo "Error: mpirun is not available!"
    echo "Please install MPI (e.g., OpenMPI or MPICH)"
    exit 1
fi

# Default values
NP=4
START=0
END=1000
INCREMENT=10
OUTPUT="dipole_output.txt"
PREFIX="298"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -np|--processes)
            NP="$2"
            shift 2
            ;;
        -s|--start)
            START="$2"
            shift 2
            ;;
        -e|--end)
            END="$2"
            shift 2
            ;;
        -i|--increment)
            INCREMENT="$2"
            shift 2
            ;;
        -o|--output)
            OUTPUT="$2"
            shift 2
            ;;
        -p|--prefix)
            PREFIX="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  -np, --processes NUM    Number of MPI processes (default: 4)"
            echo "  -s, --start NUM         Start timestep (default: 0)"
            echo "  -e, --end NUM           End timestep (default: 1000)"
            echo "  -i, --increment NUM     Timestep increment (default: 10)"
            echo "  -o, --output FILE       Output file (default: dipole_output.txt)"
            echo "  -p, --prefix STRING     File prefix (default: 298)"
            echo "  -h, --help              Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0 -np 64 -s 0 -e 10000 -i 10 -o my_dipole_output.txt"
            echo "  $0 -np 64 -s 0 -e 10000 -i 10 -o my_dipole_output.txt -p 303"
            echo "  $0 -np 64 -s 0 -e 10000 -i 10 -o my_dipole_output.txt -p 273.15"
            echo ""
            echo "File naming:"
            echo "  Input files: {prefix}.{timestep} (e.g., 303.0, 303.10, 303.20, ...)"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            echo "Use -h or --help for usage information"
            exit 1
            ;;
    esac
done

echo "Running MPI dipole moment calculation..."
echo "  Processes: $NP"
echo "  Timesteps: $START to $END (step $INCREMENT)"
echo "  File prefix: $PREFIX"
echo "  Input files: $PREFIX.$START, $PREFIX.$((START+INCREMENT)), ..."
echo "  Output: $OUTPUT"
echo ""

# Run the MPI version
mpirun -np $NP python 1.dipoleMoment_mpi.py $START $END $INCREMENT $OUTPUT $PREFIX

echo ""
echo "MPI calculation complete!"
echo "Check output file: $OUTPUT"
