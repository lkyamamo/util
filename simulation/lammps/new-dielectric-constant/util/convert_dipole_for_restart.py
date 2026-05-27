#!/usr/bin/env python3
"""
Convert old-format dipole output files to new format with timesteps for restart capability.

Old format:
    dipole_x  dipole_y  dipole_z

New format:
    timestep  dipole_x  dipole_y  dipole_z

Usage:
    python convert_dipole_for_restart.py <input_file> <output_file> [start] [increment]
    
Examples:
    # Convert with default start=0, increment=10
    python convert_dipole_for_restart.py dipole_output.txt dipole_output_restart.txt
    
    # Convert with custom start and increment
    python convert_dipole_for_restart.py dipole_output.txt dipole_output_restart.txt 0 10
"""

import sys
import os


def convert_dipole_file(input_file, output_file, start=0, increment=10):
    """
    Convert old-format dipole file to new format with timesteps.
    
    Parameters:
    -----------
    input_file : str
        Path to old-format dipole file
    output_file : str
        Path to new-format dipole file
    start : int
        Starting timestep (default: 0)
    increment : int
        Timestep increment (default: 10)
    """
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist")
        return False
    
    print(f"Converting: {input_file} -> {output_file}")
    print(f"Start timestep: {start}, Increment: {increment}")
    
    try:
        with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
            timestep = start
            line_count = 0
            
            for line in f_in:
                line = line.strip()
                
                # Skip empty lines and comments
                if not line or line.startswith('#'):
                    continue
                
                # Parse dipole values
                parts = line.split()
                if len(parts) >= 3:
                    try:
                        dipole_x = float(parts[0])
                        dipole_y = float(parts[1])
                        dipole_z = float(parts[2])
                        
                        # Write new format with timestep
                        f_out.write(f"{timestep}  {dipole_x:14.6f}  {dipole_y:14.6f}  {dipole_z:14.6f}\n")
                        
                        timestep += increment
                        line_count += 1
                        
                    except ValueError as e:
                        print(f"Warning: Could not parse line: {line}")
                        print(f"  Error: {e}")
                        continue
                else:
                    print(f"Warning: Skipping line with insufficient data: {line}")
            
            print(f"\nConverted {line_count} dipole values")
            print(f"Timestep range: {start} to {start + (line_count - 1) * increment}")
            print(f"Output written to: {output_file}")
            return True
            
    except Exception as e:
        print(f"Error converting file: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    if len(sys.argv) < 3:
        print("Usage:")
        print("  python convert_dipole_for_restart.py <input_file> <output_file> [start] [increment]")
        print("\nExamples:")
        print("  # Convert with default start=0, increment=10")
        print("  python convert_dipole_for_restart.py dipole_output.txt dipole_output_restart.txt")
        print("\n  # Convert with custom start and increment")
        print("  python convert_dipole_for_restart.py dipole_output.txt dipole_output_restart.txt 0 10")
        print("\n  # Convert test output")
        print("  python convert_dipole_for_restart.py test_outputs/new_script_output.txt test_outputs/new_script_output_restart.txt")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    # Parse optional start and increment
    start = int(sys.argv[3]) if len(sys.argv) > 3 else 0
    increment = int(sys.argv[4]) if len(sys.argv) > 4 else 10
    
    # Convert file
    success = convert_dipole_file(input_file, output_file, start, increment)
    
    if success:
        print("\n✓ Conversion successful!")
        print(f"\nYou can now use '{output_file}' for restart with:")
        print(f"  python process_bonds_inline_3atom.py ... {output_file}")
        print(f"  sbatch submit_slurm_mpi.sh")
    else:
        print("\n✗ Conversion failed!")
        sys.exit(1)


if __name__ == "__main__":
    main()

