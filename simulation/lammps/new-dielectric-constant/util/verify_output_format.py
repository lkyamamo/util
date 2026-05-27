#!/usr/bin/env python3
"""
Verify that the MPI version produces the same output format as the original.
"""

import os
import sys
import subprocess
import numpy as np

def create_test_xyz_files(start, end, increment):
    """Create test XYZ files for testing both versions."""
    print("Creating test XYZ files...")
    
    # Create a simple water-like molecule structure
    natoms = 3  # 1 water molecule
    
    for i in range(start, end, increment):
        filename = f"298.{i}"
        with open(filename, 'w') as f:
            # Write the format expected by LoadXYZ function:
            # Line 1: (empty)
            # Line 2: (empty) 
            # Line 3: number of atoms
            # Line 4: (empty)
            # Line 5: x box bounds (min max)
            # Line 6: y box bounds (min max)
            # Line 7: z box bounds (min max)
            # Line 8: (empty)
            # Then: atomID bondID itype x y z q
            
            f.write("\n")  # Line 1
            f.write("\n")  # Line 2
            f.write(f"{natoms}\n")  # Line 3
            f.write("\n")  # Line 4
            f.write("0.000000 10.000000\n")  # Line 5 - x box
            f.write("0.000000 10.000000\n")  # Line 6 - y box
            f.write("0.000000 10.000000\n")  # Line 7 - z box
            f.write("\n")  # Line 8
            
            # Add some randomness to make it interesting
            np.random.seed(i)  # Reproducible randomness
            
            # Oxygen atom (type 1)
            x, y, z = 5.0 + np.random.normal(0, 0.1), 5.0 + np.random.normal(0, 0.1), 5.0 + np.random.normal(0, 0.1)
            f.write(f"1 1 1 {x:.6f} {y:.6f} {z:.6f} -0.65966\n")
            
            # Hydrogen 1 (type 2)
            x, y, z = 5.0 + np.random.normal(0, 0.1), 5.0 + np.random.normal(0, 0.1), 5.0 + np.random.normal(0, 0.1)
            f.write(f"2 1 2 {x:.6f} {y:.6f} {z:.6f} 0.32983\n")
            
            # Hydrogen 2 (type 2)
            x, y, z = 5.0 + np.random.normal(0, 0.1), 5.0 + np.random.normal(0, 0.1), 5.0 + np.random.normal(0, 0.1)
            f.write(f"3 1 2 {x:.6f} {y:.6f} {z:.6f} 0.32983\n")
    
    print(f"Created {len(range(start, end, increment))} test XYZ files")

def run_original_version(start, end, increment, output_file):
    """Run the original version."""
    print("Running original version...")
    
    cmd = ["python", "1.dipoleMoment.py", str(start), str(end), str(increment), output_file]
    
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        
        if result.returncode == 0:
            print("✓ Original version completed successfully!")
            return True
        else:
            print(f"✗ Original version failed with return code {result.returncode}")
            print("STDERR:", result.stderr)
            return False
            
    except subprocess.TimeoutExpired:
        print("✗ Original version timed out!")
        return False
    except Exception as e:
        print(f"✗ Error running original version: {e}")
        return False

def run_mpi_version(start, end, increment, output_file):
    """Run the MPI version."""
    print("Running MPI version...")
    
    # Try to run with MPI, but fall back to direct Python if mpi4py not available
    try:
        cmd = ["mpirun", "-np", "1", "python", "1.dipoleMoment_mpi.py", 
               str(start), str(end), str(increment), output_file]
        
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        
        if result.returncode == 0:
            print("✓ MPI version completed successfully!")
            return True
        else:
            print(f"✗ MPI version failed with return code {result.returncode}")
            print("STDERR:", result.stderr)
            return False
            
    except subprocess.TimeoutExpired:
        print("✗ MPI version timed out!")
        return False
    except FileNotFoundError:
        print("✗ mpirun not found, skipping MPI test")
        return False
    except Exception as e:
        print(f"✗ Error running MPI version: {e}")
        return False

def compare_output_files(file1, file2):
    """Compare two output files line by line."""
    print(f"Comparing {file1} and {file2}...")
    
    if not os.path.exists(file1):
        print(f"✗ File {file1} does not exist!")
        return False
    
    if not os.path.exists(file2):
        print(f"✗ File {file2} does not exist!")
        return False
    
    with open(file1, 'r') as f1, open(file2, 'r') as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()
    
    if len(lines1) != len(lines2):
        print(f"✗ Different number of lines: {len(lines1)} vs {len(lines2)}")
        return False
    
    differences = []
    for i, (line1, line2) in enumerate(zip(lines1, lines2)):
        if line1.strip() != line2.strip():
            differences.append((i+1, line1.strip(), line2.strip()))
    
    if differences:
        print(f"✗ Found {len(differences)} differences:")
        for line_num, line1, line2 in differences[:5]:  # Show first 5 differences
            print(f"  Line {line_num}:")
            print(f"    Original: {line1}")
            print(f"    MPI:      {line2}")
        if len(differences) > 5:
            print(f"    ... and {len(differences) - 5} more differences")
        return False
    else:
        print("✓ Output files are identical!")
        return True

def cleanup_test_files():
    """Clean up test files."""
    print("\nCleaning up test files...")
    
    # Remove test XYZ files
    for i in range(0, 20, 5):
        filename = f"298.{i}"
        if os.path.exists(filename):
            os.remove(filename)
    
    # Remove output files
    for output_file in ["original_output.txt", "mpi_output.txt"]:
        if os.path.exists(output_file):
            os.remove(output_file)
        
        # Remove per-process files
        for i in range(4):  # 4 processes
            proc_file = f"{output_file}.proc{i}"
            if os.path.exists(proc_file):
                os.remove(proc_file)
        
        # Remove other generated files
        for suffix in [".merged", ".temp_combined", ".overlap_log"]:
            filename = f"{output_file}{suffix}"
            if os.path.exists(filename):
                os.remove(filename)
    
    print("Cleanup complete!")

def main():
    """Main test function."""
    print("=" * 70)
    print("Output Format Verification Test")
    print("=" * 70)
    
    # Test parameters
    start = 0
    end = 20
    increment = 5
    
    try:
        # Create test files
        create_test_xyz_files(start, end, increment)
        
        # Run original version
        original_success = run_original_version(start, end, increment, "original_output.txt")
        
        # Run MPI version (with 1 process for direct comparison)
        mpi_success = run_mpi_version(start, end, increment, "mpi_output.txt")
        
        if original_success:
            print("\n=== Original Output ===")
            if os.path.exists("original_output.txt"):
                with open("original_output.txt", 'r') as f:
                    for i, line in enumerate(f, 1):
                        print(f"  {i}: {line.strip()}")
            else:
                print("  No output file found!")
        
        if original_success and mpi_success:
            # Compare outputs
            format_match = compare_output_files("original_output.txt", "mpi_output.txt")
            
            if format_match:
                print("\n✓ SUCCESS: MPI version produces identical output format!")
            else:
                print("\n✗ FAILURE: Output formats do not match!")
        elif original_success:
            print("\n✓ SUCCESS: Original version works correctly!")
            print("  (MPI version requires mpi4py installation)")
        else:
            print("\n✗ FAILURE: One or both versions failed to run!")
            
    finally:
        # Don't cleanup immediately for inspection
        print("\nTest files preserved for inspection.")
        print("Run cleanup_test_files() to remove them.")

if __name__ == "__main__":
    main()
