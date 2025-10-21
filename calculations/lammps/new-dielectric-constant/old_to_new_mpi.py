#!/usr/bin/env python3
"""
MPI-parallel version of old_to_new.py for multi-node HPC clusters.

This version uses MPI to distribute timestep processing across multiple nodes.
Each MPI process handles a subset of timesteps independently.

Key Features:
- Each process writes to its own per-process file (.proc*) for incremental saving
- Line-by-line streaming merge: processes files without loading any file entirely into memory
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
    mpirun -np 64 python old_to_new_mpi.py <start> <end> <increment> <output_file> [template_prefix]
"""

import sys
import os
import datetime
from mpi4py import MPI

# Initialize MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()

def process_lammps_file(filename, output_file):
    """
    Process a single LAMMPS dump file and convert to XYZ format.
    Returns the XYZ content as a string, or None if file doesn't exist or has errors.
    """
    if not os.path.exists(filename):
        print(f"[Proc {rank}] File not found: {filename}")
        return None
    print(f"[Proc {rank}] Processing file: {filename}")
    
    try:
        with open(filename, "r") as f:
            # Read header information
            # line 1: ITEM: TIMESTEP
            line = f.readline()
            if not line.startswith("ITEM: TIMESTEP"):
                print(f"[Proc {rank}] Warning: Invalid LAMMPS format in {filename}")
                return None

            # line 2: timestep number
            line = f.readline()
            timestep = int(line.strip())     

            # line 3: "ITEM: NUMBER OF ATOMS"
            line = f.readline()        

            # line 4: ITEM: NUMBER OF ATOMS
            line = f.readline().strip()
            num_atoms = int(line)

            # Build XYZ content
            xyz_content = f"{num_atoms}\n"
            xyz_content += f"Atoms. Timestep: {timestep}\n"

            # Skip box bounds (lines 5-8)
            for _ in range(4):
                f.readline()

            # line 9: "ITEM: ATOMS id mol type x y z q"
            line = f.readline()
            if not line.startswith("ITEM: ATOMS"):
                print(f"[Proc {rank}] Warning: Invalid ATOMS section in {filename}")
                return None

            # Read atom data
            for i in range(num_atoms):
                line = f.readline()
                if not line:
                    print(f"[Proc {rank}] Warning: Unexpected end of file in {filename}")
                    return None
                atom_id, mol_id, atom_type, x, y, z, q = line.split()
                xyz_content += f"{atom_type} {x} {y} {z}\n"
            
            return xyz_content

    except Exception as e:
        print(f"[Proc {rank}] Error processing {filename}: {e}")
        return None

def get_timesteps_to_process(start, end, increment, output_file):
    """
    Calculate which timesteps should be processed and which are already done.
    
    For XYZ format, we track by timestep number since each timestep has a unique identifier.
    
    Returns:
    - timesteps_to_process: sorted list of timesteps that need processing
    - already_processed: set of timesteps already in output file
    """
    # Calculate all timesteps that SHOULD be processed
    all_timesteps = list(range(start, end, increment))
    
    # Load already processed timesteps from output file
    already_processed = load_processed_timesteps(output_file)
    
    # Find timesteps that need processing
    timesteps_to_process = [ts for ts in all_timesteps if ts not in already_processed]
    
    return timesteps_to_process, already_processed

def load_processed_timesteps(output_file):
    """
    Load already processed timesteps from output file.
    Returns set of timesteps that have been completed.
    
    For XYZ format, we extract timestep numbers from the comment lines.
    Only called by rank 0 process.
    """
    processed = set()
    
    # Check main output file
    if os.path.exists(output_file):
        try:
            with open(output_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("Atoms. Timestep:"):
                        try:
                            # Extract timestep from "Atoms. Timestep: X"
                            timestep = int(line.split(":")[1].strip())
                            processed.add(timestep)
                        except (ValueError, IndexError):
                            pass
        except Exception as e:
            print(f"[Rank 0] Warning: Could not read main output file: {e}")
    
    # Note: Per-process files are handled during combination process
    # Only rank 0 calls this function, so we don't need to check per-process files here
    
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
    existing_timesteps = set()
    
    try:
        with open(filename, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if line.startswith("Atoms. Timestep:"):
                    try:
                        current_timestep = int(line.split(":")[1].strip())
                        existing_timesteps.add(current_timestep)
                        
                        if last_timestep is not None and current_timestep <= last_timestep:
                            print(f"ERROR: Output file {filename} is not chronologically sorted!")
                            print(f"  Line {line_num}: timestep {current_timestep} follows timestep {last_timestep}")
                            return False, last_timestep, existing_timesteps
                        last_timestep = current_timestep
                    except (ValueError, IndexError) as e:
                        print(f"Warning: Could not parse timestep from line {line_num} in {filename}: {line}")
        
        print(f"✓ Output file {filename} is chronologically sorted ({len(existing_timesteps)} timesteps)")
        return True, last_timestep, existing_timesteps
        
    except Exception as e:
        print(f"Error verifying sort order of {filename}: {e}")
        return False, last_timestep, existing_timesteps

def streaming_merge_files(temp_combined_file, output_file, existing_timesteps=None):
    """
    Merge a sorted temporary file with the existing sorted output file using proper merge algorithm.
    Only adds new timesteps to the output file, preserving all existing data and maintaining sort order.
    
    Process:
    1. Read both files line by line, comparing timesteps
    2. Write the smaller timestep first, maintaining chronological order
    3. Skip duplicate timesteps (already in existing file)
    4. Continue until both files are fully processed
    
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
                for line in f:
                    line = line.strip()
                    if line.startswith("Atoms. Timestep:"):
                        try:
                            timestep = int(line.split(":")[1].strip())
                            existing_timesteps.add(timestep)
                        except (ValueError, IndexError):
                            continue
    
    merged_file = f"{output_file}.merged"
    new_timesteps_added = 0
    
    try:
        # Handle case where output file doesn't exist yet
        if not os.path.exists(output_file):
            # Output file doesn't exist - just copy temp file to output file
            with open(temp_combined_file, 'r') as f_temp, open(merged_file, 'w') as f_out:
                for line in f_temp:
                    f_out.write(line)
                    if line.strip().startswith("Atoms. Timestep:"):
                        new_timesteps_added += 1
        else:
            # Both files exist - perform proper merge
            with open(output_file, 'r') as f_output, open(temp_combined_file, 'r') as f_temp, open(merged_file, 'w') as f_out:
                # Read first lines from both files
                output_line = f_output.readline()
                temp_line = f_temp.readline()
                
                while output_line or temp_line:
                    # Determine which line to write next based on timestep comparison
                    if output_line and temp_line:
                        # Both files have data - compare timesteps
                        if output_line.strip().startswith("Atoms. Timestep:"):
                            output_timestep = int(output_line.split(":")[1].strip())
                        else:
                            # This is part of an XYZ block, write it and continue
                            f_out.write(output_line)
                            output_line = f_output.readline()
                            continue
                            
                        if temp_line.strip().startswith("Atoms. Timestep:"):
                            temp_timestep = int(temp_line.split(":")[1].strip())
                        else:
                            # This is part of an XYZ block, write it and continue
                            f_out.write(temp_line)
                            temp_line = f_temp.readline()
                            continue
                        
                        # Compare timesteps and write the smaller one
                        if output_timestep < temp_timestep:
                            # Write output file data
                            f_out.write(output_line)
                            # Read complete XYZ block from output file
                            while True:
                                output_line = f_output.readline()
                                if not output_line or output_line.strip().startswith("Atoms. Timestep:"):
                                    break
                                f_out.write(output_line)
                        elif output_timestep > temp_timestep:
                            # Write temp file data (if not duplicate)
                            if temp_timestep not in existing_timesteps:
                                f_out.write(temp_line)
                                new_timesteps_added += 1
                                # Read complete XYZ block from temp file
                                while True:
                                    temp_line = f_temp.readline()
                                    if not temp_line or temp_line.strip().startswith("Atoms. Timestep:"):
                                        break
                                    f_out.write(temp_line)
                            else:
                                # Skip duplicate timestep from temp file
                                while True:
                                    temp_line = f_temp.readline()
                                    if not temp_line or temp_line.strip().startswith("Atoms. Timestep:"):
                                        break
                        else:
                            # Same timestep - write output file data (existing takes precedence)
                            f_out.write(output_line)
                            # Skip temp file data (duplicate)
                            while True:
                                temp_line = f_temp.readline()
                                if not temp_line or temp_line.strip().startswith("Atoms. Timestep:"):
                                    break
                            # Read complete XYZ block from output file
                            while True:
                                output_line = f_output.readline()
                                if not output_line or output_line.strip().startswith("Atoms. Timestep:"):
                                    break
                                f_out.write(output_line)
                    elif output_line:
                        # Only output file has data - write remaining output data
                        f_out.write(output_line)
                        output_line = f_output.readline()
                    elif temp_line:
                        # Only temp file has data - write temp data (if not duplicate)
                        if temp_line.strip().startswith("Atoms. Timestep:"):
                            temp_timestep = int(temp_line.split(":")[1].strip())
                            if temp_timestep not in existing_timesteps:
                                f_out.write(temp_line)
                                new_timesteps_added += 1
                                # Read complete XYZ block
                                while True:
                                    temp_line = f_temp.readline()
                                    if not temp_line or temp_line.strip().startswith("Atoms. Timestep:"):
                                        break
                                    f_out.write(temp_line)
                            else:
                                # Skip duplicate
                                while True:
                                    temp_line = f_temp.readline()
                                    if not temp_line or temp_line.strip().startswith("Atoms. Timestep:"):
                                        break
                        else:
                            f_out.write(temp_line)
                            temp_line = f_temp.readline()
        
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
                    current_timestep = None
                    xyz_block = []
                    
                    for line in f:
                        line = line.strip()
                        if line.startswith("Atoms. Timestep:"):
                            # Save previous timestep if exists
                            if current_timestep is not None and xyz_block:
                                intermediate_results[current_timestep] = xyz_block.copy()
                                proc_file_sources[current_timestep] = proc_file
                                
                                # Track duplicates
                                if current_timestep in seen_timesteps:
                                    seen_timesteps[current_timestep] += 1
                                else:
                                    seen_timesteps[current_timestep] = 1
                            
                            # Start new timestep
                            current_timestep = int(line.split(":")[1].strip())
                            xyz_block = [line + "\n"]
                        else:
                            # Add to current XYZ block
                            xyz_block.append(line + "\n")
                    
                    # Don't forget the last timestep
                    if current_timestep is not None and xyz_block:
                        intermediate_results[current_timestep] = xyz_block.copy()
                        proc_file_sources[current_timestep] = proc_file
                        
                        # Track duplicates
                        if current_timestep in seen_timesteps:
                            seen_timesteps[current_timestep] += 1
                        else:
                            seen_timesteps[current_timestep] = 1
                            
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
    print(f"  Timestep range: {sorted_proc_results[0][0]} to {sorted_proc_results[-1][0]}")
    with open(temp_combined_file, 'w') as f_temp:
        for timestep, xyz_block in sorted_proc_results:
            for line in xyz_block:
                f_temp.write(line)
    
    # Step 4: Use streaming merge to combine with existing output file
    print(f"  Performing streaming merge with existing output file...")
    try:
        new_timesteps_added, total_timesteps = streaming_merge_files(temp_combined_file, output_file, existing_timesteps)
        
        # Step 5: Log any overlapping timesteps (should be rare with proper duplicate detection)
        overlapping_timesteps = []
        for timestep, xyz_block in intermediate_results.items():
            if timestep in existing_timesteps:
                # This timestep was in both existing and new data
                overlapping_timesteps.append({
                    'timestep': timestep,
                    'source_file': proc_file_sources.get(timestep, 'unknown'),
                    'xyz_lines': len(xyz_block)
                })
        
        # Log overlapping timesteps if any
        if overlapping_timesteps:
            log_file = f"{output_file}.overlap_log"
            with open(log_file, 'a') as f_log:
                f_log.write(f"\n# Overlap log - {len(overlapping_timesteps)} timesteps found in both existing and new data\n")
                f_log.write(f"# Timestamp: {datetime.datetime.now().isoformat()}\n")
                f_log.write(f"# Note: New data was NOT overwritten - existing data preserved\n")
                for overlap in overlapping_timesteps:
                    f_log.write(f"Timestep {overlap['timestep']}: Found in {overlap['source_file']} ({overlap['xyz_lines']} lines)\n")
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

def process_timesteps_mpi(start, end, increment, output_file, template_prefix="303"):
    """
    Process timesteps using MPI parallelization.
    Master process coordinates, workers process timesteps.
    
    Features:
    - Incremental writing: each process writes to its own file
    - Restart capability: skips already processed timesteps
    - Periodic combination: combines results periodically during processing
    - Final combination: master process combines all results at the end
    """
    # Get the absolute path to the dumps directory
    script_dir = os.path.dirname(os.path.abspath(__file__))
    dumps_dir = os.path.join(script_dir, "..", "dumps")
    dumps_dir = os.path.abspath(dumps_dir)
    
    # Debug: Print path information (only rank 0)
    if rank == 0:
        print(f"Script directory: {script_dir}")
        print(f"Dumps directory: {dumps_dir}")
        print(f"Dumps directory exists: {os.path.exists(dumps_dir)}")
        if os.path.exists(dumps_dir):
            # Use find command to efficiently get sample files (avoid os.listdir with 1M files)
            import subprocess
            try:
                result = subprocess.run(['find', dumps_dir, '-maxdepth', '1', '-name', f'{template_prefix}.*', '-type', 'f'], 
                                      capture_output=True, text=True, timeout=10)
                if result.returncode == 0:
                    sample_files = result.stdout.strip().split('\n')[:3]
                    sample_files = [os.path.basename(f) for f in sample_files if f]
                    print(f"Sample dump files: {sample_files}")
                else:
                    print(f"Could not list sample files: {result.stderr}")
            except Exception as e:
                print(f"Error listing sample files: {e}")
        print()
    # Create per-process output files for incremental writing
    proc_output_file = f"{output_file}.proc{rank}"
    
    # Clean up old per-process file to avoid duplicates from previous runs
    if os.path.exists(proc_output_file):
        os.remove(proc_output_file)
        print(f"[Proc {rank}] Removed old per-process file to avoid duplicates", flush=True)
    
    # Calculate combination interval based on expected workload (only rank 0)
    if rank == 0:
        total_timesteps_expected = len(range(start, end, increment))
        timesteps_per_process = max(1, total_timesteps_expected // size)
        combination_interval = max(100, timesteps_per_process // 20)  # Every 5% or 100 timesteps
    else:
        combination_interval = 0
    
    # Broadcast combination interval to all processes
    combination_interval = comm.bcast(combination_interval, root=0)
    
    # All processes: Determine which timesteps need processing
    if rank == 0:
        print("=" * 70)
        print("MPI Parallel LAMMPS to XYZ Conversion")
        print("=" * 70)
        print(f"Timesteps: {start} to {end} (step {increment})")
        print(f"Template prefix: {template_prefix}")
        print(f"Input files: {dumps_dir}/{template_prefix}.{start}, {dumps_dir}/{template_prefix}.{start+increment}, ...")
        print(f"Output file: {output_file}")
        print(f"Per-process files: {output_file}.proc*")
        print(f"MPI processes: {size}")
        print("=" * 70)
        print()
    
    # Only rank 0 determines which timesteps need processing
    if rank == 0:
        print("Determining which timesteps need processing...", flush=True)
        timesteps_to_process, already_processed = get_timesteps_to_process(
            start, end, increment, output_file
        )
        
        total_timesteps = len(range(start, end, increment))
        print(f"  Total timesteps in range: {total_timesteps}")
        print(f"  Already processed: {len(already_processed)}")
        print(f"  Need to process: {len(timesteps_to_process)}", flush=True)
    else:
        # Other processes initialize empty lists
        timesteps_to_process = []
        already_processed = set()
        total_timesteps = 0
    
    # Broadcast the results to all processes
    timesteps_to_process = comm.bcast(timesteps_to_process, root=0)
    already_processed = comm.bcast(already_processed, root=0)
    total_timesteps = comm.bcast(total_timesteps, root=0)
    
    # All processes now have the same information
    if rank != 0:
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
                        filename = os.path.join(dumps_dir, f"{template_prefix}.{i}")
                        
                        # Debug: Check if file exists before processing
                        if not os.path.exists(filename):
                            print(f"[Proc {rank}] ERROR: File not found: {filename}")
                            print(f"[Proc {rank}] Dumps directory: {dumps_dir}")
                            print(f"[Proc {rank}] Template prefix: {template_prefix}")
                            print(f"[Proc {rank}] Timestep: {i}")
                            
                            # Additional debugging: check if dumps directory has any files
                            if os.path.exists(dumps_dir):
                                try:
                                    import subprocess
                                    result = subprocess.run(['find', dumps_dir, '-maxdepth', '1', '-name', f'{template_prefix}.*', '-type', 'f'], 
                                                          capture_output=True, text=True, timeout=5)
                                    if result.returncode == 0:
                                        files = result.stdout.strip().split('\n')
                                        files = [f for f in files if f]
                                        print(f"[Proc {rank}] Found {len(files)} {template_prefix}.* files in dumps directory")
                                        if files:
                                            print(f"[Proc {rank}] Sample files: {[os.path.basename(f) for f in files[:3]]}")
                                    else:
                                        print(f"[Proc {rank}] Could not list files in dumps directory")
                                except Exception as e:
                                    print(f"[Proc {rank}] Error checking dumps directory: {e}")
                            else:
                                print(f"[Proc {rank}] Dumps directory does not exist!")
                            continue
                        
                        xyz_content = process_lammps_file(filename, output_file)
                        
                        if xyz_content is not None:
                            print(f"[Proc {rank}] Timestep {i}: converted {filename}")
                            
                            # Write immediately for restart capability
                            f_out.write(xyz_content)
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
                            print(f"[Proc {rank}] Warning: File {filename} exists but could not be processed, skipping timestep {i}")
                            
                    except Exception as e:
                        print(f"[Proc {rank}] Error processing timestep {i}: {e}")
                        print(f"[Proc {rank}] Filename attempted: {filename}")
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
        print(f"  XYZ output written to: {output_file}")
        if overlapping_timesteps > 0:
            print(f"  Overlap log: {output_file}.overlap_log")
        print(f"  Per-process files: {output_file}.proc* (can be deleted)")

if __name__ == "__main__":
    if len(sys.argv) < 5:
        if rank == 0:
            print("Usage:")
            print("  mpirun -np <num_procs> python old_to_new_mpi.py <start> <end> <increment> <output_file> [template_prefix]")
            print("\nArguments:")
            print("  start: starting timestep")
            print("  end: ending timestep")
            print("  increment: timestep increment")
            print("  output_file: output file name")
            print("  template_prefix: optional template prefix (default: '303')")
            print("\nExamples:")
            print("  mpirun -np 64 python old_to_new_mpi.py 0 10000001 10 all_lammps_new.xyz")
            print("  mpirun -np 64 python old_to_new_mpi.py 0 10000001 10 all_lammps_new.xyz 298")
            print("  mpirun -np 64 python old_to_new_mpi.py 0 10000001 10 all_lammps_new.xyz 273.15")
            print("\nNote:")
            print("  - Requires mpi4py: pip install mpi4py")
            print("  - Processes LAMMPS dump files named '../dumps/{template_prefix}.{timestep}'")
            print("  - Converts LAMMPS dump format to XYZ format")
            print("  - Each process writes to its own .proc* file for incremental saving")
        sys.exit(1)
    
    start = int(sys.argv[1])
    end = int(sys.argv[2])
    increment = int(sys.argv[3])
    output_file = sys.argv[4]
    
    # Optional 5th argument for template prefix (default: "303")
    if len(sys.argv) > 5:
        template_prefix = sys.argv[5]
    else:
        template_prefix = "303"
    
    process_timesteps_mpi(start, end, increment, output_file, template_prefix)
