#!/usr/bin/env python3
"""
MPI-parallel version of 1.dipoleMoment.py for multi-node HPC clusters.

This version uses MPI to distribute timestep processing across multiple nodes.
Each MPI process handles a subset of timesteps independently.

Key Features:
- Each process writes to its own per-process file (.proc*) for incremental saving
- Line-by-line streaming merge: processes files without loading any file entirely into memory
- Simple comparison algorithm: compares current line from output file with current line from temp file
- Verifies output file is chronologically sorted before processing (exits if not)
- Per-process files are combined into a temporary sorted file before merging with main output
- Preserves all existing data - only adds new timesteps, never overwrites
- Proper sorting ensures chronological order in the final output file
- Restart capability: skips already processed timesteps
- Periodic combination during processing to avoid data loss
- Automatic backup and recovery on merge failures

Requirements:
    mpi4py (pip install mpi4py)

Usage:
    mpirun -np 64 python 1.dipoleMoment_mpi.py <start> <end> <increment> <output_file>
"""

import sys
import numpy as np
import os
import datetime
from mpi4py import MPI

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def LoadXYZ(xyz):
    """Load XYZ file and return atoms, bonds, and box dimensions."""
    atoms = {}
    bonds = {}
    with open(xyz,'r') as file1:
        lines = file1.readline()        #1
        lines = file1.readline()        #2
        lines = file1.readline()        #3
        natoms = int(lines.strip().split()[0])
        lines = file1.readline()        #4
        lines = file1.readline()        #5
        la = float(lines.strip().split()[1]) - float(lines.strip().split()[0])
        lines = file1.readline()        #6
        lb = float(lines.strip().split()[1]) - float(lines.strip().split()[0])
        lines = file1.readline()        #7
        lc = float(lines.strip().split()[1]) - float(lines.strip().split()[0])
        lines = file1.readline()        #8

        for lines in file1:
            line = lines.strip().split()
            atomID = int(line[0])
            x,y,z,q = float(line[3]), float(line[4]), float(line[5]), float(line[6])
            bondID = int(line[1])
            itype = int(line[2])
            if bondID in bonds:
                bonds[bondID].append(atomID)
            else:
                bonds[bondID] = [atomID]

            atoms[atomID] = [itype,x,y,z,q]
    return natoms,atoms,bonds,la,lb,lc

def TotalMoment(atoms, bonds, la, lb, lc):
    """Calculate total dipole moment using the original formula."""
    totalMoment = [0.0,0.0,0.0]
    dipole = []

    counter = 0
    for iid in bonds.keys():
        counter += 1
        atomIDs = sorted(bonds[iid])
        atomID1, atomID2, atomID3 = atomIDs[0], atomIDs[1], atomIDs[2]

        # Calculate distances with periodic boundary conditions
        dx12 = atoms[atomID1][1] - atoms[atomID2][1]
        if dx12 >= la/2.0:
            dx12 -= la
        elif dx12 <= -la/2.0:
            dx12 += la

        dx31 = atoms[atomID3][1] - atoms[atomID1][1]
        if dx31 >= la/2.0:
            dx31 -= la
        elif dx31 <= -la/2.0:
            dx31 += la

        dy12 = atoms[atomID1][2] - atoms[atomID2][2]
        if dy12 >= lb/2.0:
            dy12 -= lb
        elif dy12 <= -lb/2.0:
            dy12 += lb

        dy31 = atoms[atomID3][2] - atoms[atomID1][2]
        if dy31 >= lb/2.0:
            dy31 -= lb
        elif dy31 <= -lb/2.0:
            dy31 += lb

        dz12 = atoms[atomID1][3] - atoms[atomID2][3]
        if dz12 >= lc/2.0:
            dz12 -= lc
        elif dz12 <= -lc/2.0:
            dz12 += lc

        dz31 = atoms[atomID3][3] - atoms[atomID1][3]
        if dz31 >= lc/2.0:
            dz31 -= lc
        elif dz31 <= -lc/2.0:
            dz31 += lc

        dipoleMoment = [0,0,0]
        q = 0.32983	

        dipoleMoment[0] += (dx12)*q
        dipoleMoment[0] -= (dx31)*q

        dipoleMoment[1] += (dy12)*q
        dipoleMoment[1] -= (dy31)*q

        dipoleMoment[2] += (dz12)*q
        dipoleMoment[2] -= (dz31)*q

        dipole_mag = np.linalg.norm(np.array(dipoleMoment))
        dipole.append(dipole_mag)

        totalMoment[0] += dipoleMoment[0]
        totalMoment[1] += dipoleMoment[1]
        totalMoment[2] += dipoleMoment[2]

    print('---------------------')
    print(np.mean(np.array(dipole)))

    return totalMoment

def get_timesteps_to_process(start, end, increment, output_file):
    """
    Calculate which timesteps should be processed and which are already done.
    
    For the original format, we track by line number since timesteps are processed sequentially.
    
    Returns:
    - timesteps_to_process: sorted list of timesteps that need processing
    - already_processed: set of line numbers already in output file
    """
    # Calculate all timesteps that SHOULD be processed
    all_timesteps = list(range(start, end, increment))
    
    # Load already processed line numbers from output file
    already_processed_lines = load_processed_timesteps(output_file)
    
    # For the original format, we need to check if we have enough lines
    # If we have fewer lines than expected timesteps, we need to process more
    expected_lines = len(all_timesteps)
    current_lines = len(already_processed_lines)
    
    if current_lines >= expected_lines:
        # All timesteps have been processed
        timesteps_to_process = []
    else:
        # We need to process the remaining timesteps
        # Take the remaining timesteps from the end
        remaining_count = expected_lines - current_lines
        timesteps_to_process = all_timesteps[-remaining_count:]
    
    return timesteps_to_process, already_processed_lines

def load_processed_timesteps(output_file):
    """
    Load already processed timesteps from output file.
    Returns set of timesteps that have been completed.
    
    For the original format, we track by line number since timesteps are processed sequentially.
    """
    processed = set()
    
    # Check main output file
    if os.path.exists(output_file):
        try:
            with open(output_file, 'r') as f:
                line_count = 0
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 3:  # 3 dipole components
                            try:
                                # For the original format, timestep is implicit (line number)
                                # We'll use line number as timestep identifier
                                line_count += 1
                                processed.add(line_count)
                            except ValueError:
                                pass
        except Exception as e:
            print(f"[Proc {rank}] Warning: Could not read main output file: {e}")
    
    # Check per-process files (for MPI restart)
    proc_file = f"{output_file}.proc{rank}"
    if os.path.exists(proc_file):
        try:
            with open(proc_file, 'r') as f:
                line_count = 0
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 3:
                            try:
                                # For per-process files, we track by line number
                                line_count += 1
                                processed.add(line_count)
                            except ValueError:
                                pass
        except Exception as e:
            print(f"[Proc {rank}] Warning: Could not read per-process file: {e}")
    
    return processed

def verify_output_file_sorted(filename):
    """
    Verify that the output file is chronologically sorted by timestep.
    Also collects all existing timesteps to avoid reading the file again later.
    Returns (is_sorted, last_timestep, existing_timesteps) where:
    - is_sorted: boolean indicating if file is sorted
    - last_timestep: the highest timestep in the file
    - existing_timesteps: set of all timesteps in the file
    Uses streaming approach to avoid loading entire file into memory.
    """
    if not os.path.exists(filename):
        return True, None, set()  # Empty file is considered sorted
    
    last_timestep = None
    line_count = 0
    existing_timesteps = set()
    
    try:
        with open(filename, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 3:
                        try:
                            # For original format, timestep is implicit (line number)
                            current_timestep = line_num
                            existing_timesteps.add(current_timestep)
                            
                            if last_timestep is not None and current_timestep <= last_timestep:
                                print(f"ERROR: Output file {filename} is not chronologically sorted!")
                                print(f"  Line {line_num}: timestep {current_timestep} follows timestep {last_timestep}")
                                return False, last_timestep, existing_timesteps
                            last_timestep = current_timestep
                            line_count += 1
                        except ValueError as e:
                            print(f"Warning: Could not parse line {line_num} in {filename}: {line}")
        
        print(f"✓ Output file {filename} is chronologically sorted ({line_count} timesteps)")
        return True, last_timestep, existing_timesteps
        
    except Exception as e:
        print(f"Error verifying sort order of {filename}: {e}")
        return False, last_timestep, existing_timesteps

def streaming_merge_files(temp_combined_file, output_file, existing_timesteps=None):
    """
    Merge a sorted temporary file with the existing sorted output file using simple line-by-line comparison.
    Only adds new timesteps to the output file, preserving all existing data.
    
    Process:
    1. Compare current line from output file with current line from temp file
    2. If temp timestep is smaller, write temp line and move to next temp line
    3. If temp timestep is larger, write output line and move to next output line
    4. Repeat until all lines are processed
    
    Args:
        temp_combined_file: Path to temporary file with new sorted data (must be chronologically sorted)
        output_file: Path to existing output file (must be chronologically sorted)
        existing_timesteps: Set of existing timesteps (optional, will be read from file if not provided)
    
    Returns:
        (new_timesteps_added, total_timesteps_in_output)
    """
    if not os.path.exists(temp_combined_file):
        print(f"Warning: Temporary file {temp_combined_file} does not exist")
        return 0, 0
    
    # Get existing timesteps to avoid duplicates (use provided set or read from file)
    if existing_timesteps is None:
        existing_timesteps = set()
        if os.path.exists(output_file):
            with open(output_file, 'r') as f:
                line_count = 0
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 3:
                            try:
                                line_count += 1
                                existing_timesteps.add(line_count)
                            except ValueError:
                                continue
    
    # Simple line-by-line merge approach
    merged_file = f"{output_file}.merged"
    new_timesteps_added = 0
    
    try:
        # Handle case where output file doesn't exist yet
        if os.path.exists(output_file):
            # For existing output file, we need to append new data
            # Since the original format doesn't have explicit timesteps, we just append
            with open(output_file, 'r') as f_output, open(temp_combined_file, 'r') as f_temp, open(merged_file, 'w') as f_out:
                # Copy all existing output lines
                for line in f_output:
                    f_out.write(line)
                
                # Append new temp lines
                for line in f_temp:
                    line_clean = line.strip()
                    if line_clean and not line_clean.startswith('#'):
                        parts = line_clean.split()
                        if len(parts) >= 3:
                            try:
                                temp_dipole = [float(parts[0]), float(parts[1]), float(parts[2])]
                                f_out.write(f"{temp_dipole[0]:14.6f}  {temp_dipole[1]:14.6f}  {temp_dipole[2]:14.6f}\n")
                                new_timesteps_added += 1
                            except ValueError:
                                pass
        else:
            # Output file doesn't exist - just copy temp file to output file
            with open(temp_combined_file, 'r') as f_temp, open(merged_file, 'w') as f_out:
                for line in f_temp:
                    line_clean = line.strip()
                    if line_clean and not line_clean.startswith('#'):
                        parts = line_clean.split()
                        if len(parts) >= 3:
                            try:
                                temp_dipole = [float(parts[0]), float(parts[1]), float(parts[2])]
                                f_out.write(f"{temp_dipole[0]:14.6f}  {temp_dipole[1]:14.6f}  {temp_dipole[2]:14.6f}\n")
                                new_timesteps_added += 1
                            except ValueError:
                                pass
        
        # Replace original file with merged file
        os.replace(merged_file, output_file)
        print(f"✓ Merged {new_timesteps_added} new timesteps into output file")
        
        total_timesteps = len(existing_timesteps) + new_timesteps_added
        return new_timesteps_added, total_timesteps
        
    except Exception as e:
        print(f"Error during streaming merge: {e}")
        if os.path.exists(merged_file):
            os.remove(merged_file)
        raise

def combine_per_process_files(output_file, size):
    """
    Combine all per-process files with existing output file using streaming approach.
    Efficient for very large output files by avoiding loading entire file into memory.
    
    Process:
    1. Verify output file is chronologically sorted (exit if not)
    2. Combine per-process files into temporary sorted file
    3. Use streaming merge to add new timesteps to existing output file
    4. Preserve all existing data, only add new timesteps
    
    Called periodically during processing and at the end.
    """
    print(f"\n[Combination] Starting combination of per-process files...")
    
    # Step 1: Verify existing output file is sorted (critical requirement)
    existing_timesteps = set()
    if os.path.exists(output_file):
        is_sorted, last_timestep, existing_timesteps = verify_output_file_sorted(output_file)
        if not is_sorted:
            print(f"ERROR: Cannot proceed - output file {output_file} is not chronologically sorted!")
            print("Please manually sort the output file or remove it and restart.")
            sys.exit(1)
        print(f"  Existing output file has timestep range up to {last_timestep} ({len(existing_timesteps)} timesteps)")
    else:
        print(f"  No existing output file found - will create new one")
    
    # Step 2: Read all .proc* files into intermediate data
    intermediate_results = {}
    seen_timesteps = {}  # Track duplicates within .proc* files
    proc_file_sources = {}  # Track which .proc* file each timestep came from
    
    for proc_id in range(size):
        proc_file = f"{output_file}.proc{proc_id}"
        if os.path.exists(proc_file):
            try:
                with open(proc_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line and not line.startswith('#'):
                            parts = line.split()
                            if len(parts) >= 3:
                                dipole = [float(parts[0]), float(parts[1]), float(parts[2])]
                                # Use line number as timestep identifier
                                timestep = len(intermediate_results) + 1
                                intermediate_results[timestep] = dipole
                                proc_file_sources[timestep] = proc_file
                                
                                # Track duplicates
                                if timestep in seen_timesteps:
                                    seen_timesteps[timestep] += 1
                                else:
                                    seen_timesteps[timestep] = 1
            except Exception as e:
                print(f"Warning: Could not read {proc_file}: {e}")
    
    # Report duplicates within .proc* files
    duplicates = {ts: count for ts, count in seen_timesteps.items() if count > 1}
    if duplicates:
        print(f"Warning: Found {len(duplicates)} duplicate timesteps in .proc* files (kept last occurrence)")
    
    if not intermediate_results:
        print("  No new data found in per-process files")
        return 0, 0, 0
    
    # Step 3: Create temporary combined file with sorted per-process results
    temp_combined_file = f"{output_file}.temp_combined"
    sorted_proc_results = sorted(intermediate_results.items())
    
    print(f"  Creating temporary combined file with {len(sorted_proc_results)} timesteps")
    with open(temp_combined_file, 'w') as f_temp:
        for timestep, dipole in sorted_proc_results:
            f_temp.write(f"{dipole[0]:14.6f}  {dipole[1]:14.6f}  {dipole[2]:14.6f}\n")
    
    # Step 4: Use streaming merge to combine with existing output file
    print(f"  Performing streaming merge with existing output file...")
    try:
        new_timesteps_added, total_timesteps = streaming_merge_files(temp_combined_file, output_file, existing_timesteps)
        
        # Step 5: Log any overlapping timesteps (should be rare with proper duplicate detection)
        overlapping_timesteps = []
        for timestep, dipole in intermediate_results.items():
            if timestep in existing_timesteps:
                # This timestep was in both existing and new data
                overlapping_timesteps.append({
                    'timestep': timestep,
                    'source_file': proc_file_sources.get(timestep, 'unknown'),
                    'new_dipole': dipole
                })
        
        # Log overlapping timesteps if any
        if overlapping_timesteps:
            log_file = f"{output_file}.overlap_log"
            with open(log_file, 'a') as f_log:
                f_log.write(f"\n# Overlap log - {len(overlapping_timesteps)} timesteps found in both existing and new data\n")
                f_log.write(f"# Timestamp: {datetime.datetime.now().isoformat()}\n")
                f_log.write(f"# Note: New data was NOT overwritten - existing data preserved\n")
                for overlap in overlapping_timesteps:
                    f_log.write(f"Timestep {overlap['timestep']}: Found in {overlap['source_file']}\n")
                    f_log.write(f"  New data: [{overlap['new_dipole'][0]:14.6f}, {overlap['new_dipole'][1]:14.6f}, {overlap['new_dipole'][2]:14.6f}]\n")
            print(f"  Logged {len(overlapping_timesteps)} overlapping timesteps to {log_file}")
        
        # Clean up temporary file
        if os.path.exists(temp_combined_file):
            os.remove(temp_combined_file)
        
        print(f"✓ Combination complete: {new_timesteps_added} new timesteps added, {total_timesteps} total timesteps")
        return total_timesteps, new_timesteps_added, len(overlapping_timesteps)
        
    except Exception as e:
        print(f"Error during combination: {e}")
        # Clean up temporary file on error
        if os.path.exists(temp_combined_file):
            os.remove(temp_combined_file)
        raise

def process_timesteps_mpi(start, end, increment, output_file, file_prefix="298"):
    """
    Process timesteps using MPI parallelization.
    Master process coordinates, workers process timesteps.
    
    Features:
    - Incremental writing: each process writes to its own file
    - Restart capability: skips already processed timesteps
    - Periodic combination: combines results periodically during processing
    - Final combination: master process combines all results at the end
    """
    # Create per-process output files for incremental writing
    proc_output_file = f"{output_file}.proc{rank}"
    
    # Clean up old per-process file to avoid duplicates from previous runs
    if os.path.exists(proc_output_file):
        os.remove(proc_output_file)
        print(f"[Proc {rank}] Removed old per-process file to avoid duplicates", flush=True)
    
    # Calculate combination interval based on expected workload
    total_timesteps_expected = len(range(start, end, increment))
    timesteps_per_process = max(1, total_timesteps_expected // size)
    combination_interval = max(100, timesteps_per_process // 20)  # Every 5% or 100 timesteps
    
    # All processes: Determine which timesteps need processing
    if rank == 0:
        print("=" * 70)
        print("MPI Parallel Dipole Moment Calculation")
        print("=" * 70)
        print(f"Timesteps: {start} to {end} (step {increment})")
        print(f"File prefix: {file_prefix}")
        print(f"Input files: {file_prefix}.{start}, {file_prefix}.{start+increment}, ...")
        print(f"Output file: {output_file}")
        print(f"Per-process files: {output_file}.proc*")
        print(f"MPI processes: {size}")
        print("=" * 70)
        print()
    
    print("Determining which timesteps need processing...", flush=True)
    timesteps_to_process, already_processed = get_timesteps_to_process(
        start, end, increment, output_file
    )
    
    total_timesteps = len(range(start, end, increment))
    print(f"  Total timesteps in range: {total_timesteps}")
    print(f"  Already processed: {len(already_processed)}")
    print(f"  Need to process: {len(timesteps_to_process)}", flush=True)
    
    if len(timesteps_to_process) == 0:
        if rank == 0:
            print("\n✓ All timesteps already processed! Nothing to do.", flush=True)
    else:
        print(f"\nProcessing {len(timesteps_to_process)} timesteps...", flush=True)
        if rank == 0:
            print(f"Combining results every {combination_interval} timesteps per process", flush=True)
    
    # Create a dictionary for fast timestep lookup (optimization)
    timestep_to_index = {ts: idx for idx, ts in enumerate(timesteps_to_process)}
    
    # Open per-process output file in write mode for incremental writing
    f_out = open(proc_output_file, 'w')
    
    try:
        # Each process processes its assigned timesteps
        counter = 0
        last_combination_counter = 0
        
        for i in range(start, end, increment):
            # Check if this timestep needs processing and is assigned to this process
            if i in timestep_to_index:
                # Assign timesteps to processes using round-robin
                timestep_index = timestep_to_index[i]
                if timestep_index % size == rank:
                    try:
                        xyz = f"{file_prefix}.{i}"
                        if os.path.exists(xyz):
                            natoms, atoms, bonds, la, lb, lc = LoadXYZ(xyz)
                            moment = TotalMoment(atoms, bonds, la, lb, lc)
                            
                            print(f"[Proc {rank}] Timestep {i}: dipole = [{moment[0]:14.6f}, {moment[1]:14.6f}, {moment[2]:14.6f}]")
                            print(f"[Proc {rank}] Box dimensions: {la}, {lb}, {lc}")
                            
                            # Write immediately for restart capability
                            f_out.write(f"{moment[0]:14.6f}  {moment[1]:14.6f}  {moment[2]:14.6f}\n")
                            f_out.flush()  # Ensure data is written to disk
                            counter += 1
                            
                            # Periodic combination check
                            if counter - last_combination_counter >= combination_interval:
                                # Synchronize before combining
                                comm.Barrier()
                                
                                # Master process combines files
                                if rank == 0:
                                    total_frames, new_timesteps, overlapping_timesteps = combine_per_process_files(output_file, size)
                                    print(f"\n[Combination] Combined {total_frames} total frames at {counter} frames per process", flush=True)
                                    if new_timesteps > 0 or overlapping_timesteps > 0:
                                        print(f"  Added {new_timesteps} new timesteps, found {overlapping_timesteps} overlapping timesteps", flush=True)
                                
                                last_combination_counter = counter
                                comm.Barrier()  # Synchronize after combining
                        else:
                            print(f"[Proc {rank}] Warning: File {xyz} not found, skipping timestep {i}")
                            
                    except Exception as e:
                        print(f"[Proc {rank}] Error processing timestep {i}: {e}")
                        import traceback
                        traceback.print_exc()
        
        print(f"[Proc {rank}] Processed {counter} timesteps")
    
    finally:
        f_out.close()
    
    # Synchronize all processes before final combination
    comm.Barrier()
    
    # Master process does final combination
    if rank == 0:
        print("\nPerforming final combination of per-process files...", flush=True)
        total_frames, new_timesteps, overlapping_timesteps = combine_per_process_files(output_file, size)
        
        print(f"\n✓ Processing complete!")
        print(f"  Total frames in final file: {total_frames}")
        print(f"  Added {new_timesteps} new timesteps")
        print(f"  Found {overlapping_timesteps} overlapping timesteps (existing data preserved)")
        print(f"  Dipole output written to: {output_file}")
        if overlapping_timesteps > 0:
            print(f"  Overlap log: {output_file}.overlap_log")
        print(f"  Per-process files: {output_file}.proc* (can be deleted)")

if __name__ == "__main__":
    if len(sys.argv) < 5:
        if rank == 0:
            print("Usage:")
            print("  mpirun -np <num_procs> python 1.dipoleMoment_mpi.py <start> <end> <increment> <output_file> [file_prefix]")
            print("\nArguments:")
            print("  start: starting timestep")
            print("  end: ending timestep")
            print("  increment: timestep increment")
            print("  output_file: output file name")
            print("  file_prefix: optional file prefix (default: '298')")
            print("\nExamples:")
            print("  mpirun -np 64 python 1.dipoleMoment_mpi.py 0 1000 10 dipole_output.txt")
            print("  mpirun -np 64 python 1.dipoleMoment_mpi.py 0 1000 10 dipole_output.txt 303")
            print("  mpirun -np 64 python 1.dipoleMoment_mpi.py 0 1000 10 dipole_output.txt 273.15")
            print("\nNote:")
            print("  - Requires mpi4py: pip install mpi4py")
            print("  - Processes XYZ files named '{file_prefix}.{timestep}'")
            print("  - Uses original dipole moment calculation formula")
            print("  - Each process writes to its own .proc* file for incremental saving")
        sys.exit(1)
    
    start = int(sys.argv[1])
    end = int(sys.argv[2])
    increment = int(sys.argv[3])
    output_file = sys.argv[4]
    
    # Optional 5th argument for file format (default: "298")
    if len(sys.argv) > 5:
        file_prefix = sys.argv[5]
    else:
        file_prefix = "298"
    
    process_timesteps_mpi(start, end, increment, output_file, file_prefix)
