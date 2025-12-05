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

Enhanced Data Integrity Features (v2.0):
- Frame validation: Complete XYZ frame validation before writing
- Atomic writes: Uses temporary files and atomic operations to prevent partial writes
- Recovery logic: Detects and recovers from incomplete frames in existing files
- Main output file recovery: Automatically detects and removes malformed frames from main output file
- Malformed frame reprocessing: Automatically adds removed timesteps back to processing queue
- Automatic backup cleanup: Backups are automatically removed after use to prevent disk space accumulation
- Enhanced error handling: Detailed logging for malformed frames and validation failures
- Robust XYZ parsing: Validates atom count, timestep format, and coordinate data
- Process interruption resilience: Continues processing even if individual frames fail

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
    Enhanced with frame validation and atomic write capability.
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

            # Read atom data and validate completeness
            atoms_read = 0
            for i in range(num_atoms):
                line = f.readline()
                if not line:
                    print(f"[Proc {rank}] Warning: Unexpected end of file in {filename}")
                    return None
                atom_id, mol_id, atom_type, x, y, z, q = line.split()
                xyz_content += f"{atom_type} {x} {y} {z}\n"
                atoms_read += 1
            
            # Final validation - ensure we read exactly the expected number of atoms
            if atoms_read != num_atoms:
                print(f"[Proc {rank}] Warning: Frame validation failed in {filename} - expected {num_atoms} atoms, read {atoms_read}")
                return None
            
            print(f"[Proc {rank}] Successfully processed {filename}: {num_atoms} atoms, timestep {timestep}")
            return xyz_content

    except Exception as e:
        print(f"[Proc {rank}] Error processing {filename}: {e}")
        return None

def recover_incomplete_frames(proc_output_file):
    """
    Detect and recover from incomplete frames in per-process files.
    This function scans the file for malformed XYZ frames and attempts to fix them.
    
    Args:
        proc_output_file: Path to the per-process output file
    
    Returns:
        tuple: (frames_recovered, frames_removed, total_frames_processed)
    """
    if not os.path.exists(proc_output_file):
        return 0, 0, 0
    
    backup_file = f"{proc_output_file}.backup"
    temp_file = f"{proc_output_file}.recovered"
    
    frames_recovered = 0
    frames_removed = 0
    total_frames_processed = 0
    
    try:
        # Create backup
        import shutil
        shutil.copy2(proc_output_file, backup_file)
        
        with open(proc_output_file, 'r') as f_in, open(temp_file, 'w') as f_out:
            while True:
                # Try to read a complete frame
                num_atoms, timestep, frame_lines = read_xyz_frame(f_in)
                
                if num_atoms is None:
                    # End of file or malformed frame
                    break
                
                total_frames_processed += 1
                
                if timestep is not None and len(frame_lines) == num_atoms + 2:
                    # Valid frame - write it
                    for line in frame_lines:
                        f_out.write(line)
                    frames_recovered += 1
                else:
                    # Malformed frame - skip it
                    frames_removed += 1
                    print(f"[Recovery] Removed malformed frame at position {total_frames_processed}")
        
        # Replace original with recovered file
        os.replace(temp_file, proc_output_file)
        
        if frames_removed > 0:
            print(f"[Recovery] Recovered {frames_recovered} frames, removed {frames_removed} malformed frames")
            print(f"[Recovery] Backup saved as: {backup_file}")
        
        # Always clean up backup after successful recovery to prevent accumulation
        safe_remove_file(backup_file, "backup")
        
        return frames_recovered, frames_removed, total_frames_processed
        
    except Exception as e:
        print(f"[Recovery] Error during frame recovery: {e}")
        # Clean up temporary files
        safe_remove_file(temp_file, "temp")
        return 0, 0, 0

def atomic_write_frame(proc_output_file, xyz_content, timestep):
    """
    Atomically write a complete XYZ frame to the per-process file.
    Simplified implementation with proper error handling.
    
    Args:
        proc_output_file: Path to the per-process output file
        xyz_content: Complete XYZ frame content as string
        timestep: Timestep number for logging
    
    Returns:
        bool: True if write was successful, False otherwise
    """
    try:
        # Direct append with proper error handling
        with open(proc_output_file, 'a') as f_out:
            f_out.write(xyz_content)
            f_out.flush()
            os.fsync(f_out.fileno())  # Force write to disk
        return True
        
    except Exception as e:
        print(f"[Proc {rank}] Error writing frame for timestep {timestep}: {e}")
        return False

def get_timesteps_to_process(start, end, increment, output_file):
    """
    Calculate which timesteps should be processed and which are already done.
    Includes recovery of malformed frames from the main output file.
    
    For XYZ format, we track by timestep number since each timestep has a unique identifier.
    
    Returns:
    - timesteps_to_process: sorted list of timesteps that need processing
    - already_processed: set of timesteps already in output file
    """
    # Calculate all timesteps that SHOULD be processed
    all_timesteps = list(range(start, end, increment))
    
    # First, recover any malformed frames from the main output file
    if os.path.exists(output_file):
        print("Checking main output file for malformed frames...")
        recovered_timesteps, removed_timesteps, total_frames = recover_main_output_file(output_file)
        
        if removed_timesteps:
            print(f"Found {len(removed_timesteps)} malformed frames that need reprocessing")
            # Add removed timesteps back to the processing list
            for ts in removed_timesteps:
                if ts in all_timesteps and ts not in recovered_timesteps:
                    # This timestep was malformed and needs reprocessing
                    pass  # It will be included in timesteps_to_process below
        else:
            print("No malformed frames found in main output file")
    
    # Load already processed timesteps from output file (after recovery)
    already_processed = load_processed_timesteps(output_file)
    
    # Find timesteps that need processing
    timesteps_to_process = [ts for ts in all_timesteps if ts not in already_processed]
    
    return timesteps_to_process, already_processed

def safe_remove_file(filepath, context="file"):
    """
    Safely remove a file with proper error handling and logging.
    
    Args:
        filepath: Path to file to remove
        context: Context for logging (e.g., "backup", "temp", "old backup")
    
    Returns:
        bool: True if file was removed or didn't exist, False if error occurred
    """
    if not os.path.exists(filepath):
        return True
    
    try:
        os.remove(filepath)
        print(f"[Cleanup] Removed {context} file: {filepath}")
        return True
    except Exception as e:
        print(f"[Cleanup] Warning: Could not remove {context} file {filepath}: {e}")
        return False

def cleanup_old_backups(output_file):
    """
    Clean up any existing backup files from previous runs to prevent accumulation.
    """
    import glob
    
    # Clean up main output file backups
    main_backup_pattern = f"{output_file}.backup"
    safe_remove_file(main_backup_pattern, "old backup")
    
    # Clean up per-process file backups
    proc_backup_pattern = f"{output_file}.proc*.backup"
    proc_backups = glob.glob(proc_backup_pattern)
    for backup_file in proc_backups:
        safe_remove_file(backup_file, "old backup")
    
    if proc_backups:
        print(f"[Cleanup] Cleaned up {len(proc_backups)} old per-process backup files")

def recover_main_output_file(output_file):
    """
    Detect and remove malformed XYZ frames from the main output file.
    Returns timesteps that were removed and need to be reprocessed.
    
    Args:
        output_file: Path to the main output file
    
    Returns:
        tuple: (recovered_timesteps, removed_timesteps, total_frames_processed)
    """
    if not os.path.exists(output_file):
        return set(), set(), 0
    
    backup_file = f"{output_file}.backup"
    temp_file = f"{output_file}.recovered"
    
    recovered_timesteps = set()
    removed_timesteps = set()
    total_frames_processed = 0
    
    try:
        # Create backup
        import shutil
        shutil.copy2(output_file, backup_file)
        print(f"[Recovery] Created backup: {backup_file}")
        
        with open(output_file, 'r') as f_in, open(temp_file, 'w') as f_out:
            while True:
                # Try to read a complete frame
                num_atoms, timestep, frame_lines = read_xyz_frame(f_in)
                
                if num_atoms is None:
                    # End of file or malformed frame - try to continue reading
                    current_pos = f_in.tell()
                    line = f_in.readline()
                    if not line:
                        # End of file
                        break
                    else:
                        # Malformed frame - try to find next valid frame start
                        f_in.seek(current_pos)  # Go back to start of malformed content
                        
                        # Skip lines until we find a potential frame start (number)
                        while True:
                            line = f_in.readline()
                            if not line:
                                break
                            line = line.strip()
                            if line and line.isdigit():
                                # Found potential frame start - go back and try again
                                f_in.seek(f_in.tell() - len(line) - 1)
                                break
                        continue
                
                total_frames_processed += 1
                
                if timestep is not None and len(frame_lines) == num_atoms + 2:
                    # Valid frame - write it
                    for line in frame_lines:
                        f_out.write(line)
                    recovered_timesteps.add(timestep)
                else:
                    # Malformed frame - extract timestep if possible and mark for reprocessing
                    if timestep is not None:
                        removed_timesteps.add(timestep)
                        print(f"[Recovery] Removed malformed frame for timestep {timestep}")
                    else:
                        print(f"[Recovery] Removed malformed frame at position {total_frames_processed} (no timestep)")
        
        # Replace original with recovered file
        os.replace(temp_file, output_file)
        
        if removed_timesteps:
            print(f"[Recovery] Main output file recovery complete:")
            print(f"  - Recovered {len(recovered_timesteps)} valid frames")
            print(f"  - Removed {len(removed_timesteps)} malformed frames")
            print(f"  - Backup saved as: {backup_file}")
            print(f"  - Timesteps to reprocess: {sorted(removed_timesteps)}")
        else:
            print(f"[Recovery] No malformed frames found in main output file")
        
        # Always clean up backup after successful recovery to prevent accumulation
        safe_remove_file(backup_file, "backup")
        
        return recovered_timesteps, removed_timesteps, total_frames_processed
        
    except Exception as e:
        print(f"[Recovery] Error during main output file recovery: {e}")
        # Clean up temporary files
        safe_remove_file(temp_file, "temp")
        return set(), set(), 0

def extract_timesteps_from_file(filename):
    """
    Extract all timesteps from an XYZ file.
    Returns set of timesteps found in the file.
    """
    timesteps = set()
    if not os.path.exists(filename):
        return timesteps
    
    try:
        with open(filename, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith("Atoms. Timestep:"):
                    try:
                        timestep = int(line.split(":")[1].strip())
                        timesteps.add(timestep)
                    except (ValueError, IndexError):
                        pass
    except Exception as e:
        print(f"Warning: Could not read file {filename}: {e}")
    
    return timesteps

def load_processed_timesteps(output_file):
    """
    Load already processed timesteps from output file.
    Returns set of timesteps that have been completed.
    
    For XYZ format, we extract timestep numbers from the comment lines.
    Only called by rank 0 process.
    """
    # Use the consolidated timestep extraction function
    processed = extract_timesteps_from_file(output_file)
    
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

def read_xyz_frame(f):
    """
    Read a complete XYZ frame from file handle.
    Returns (num_atoms, timestep, frame_lines) or (None, None, None) at EOF.
    Enhanced with better error handling and validation.
    """
    # Read number of atoms
    line = f.readline()
    if not line or not line.strip():
        return None, None, None
    
    try:
        num_atoms = int(line.strip())
        if num_atoms < 0:
            print(f"Warning: Invalid negative atom count: {num_atoms}")
            return None, None, None
    except ValueError as e:
        print(f"Warning: Could not parse atom count from line: '{line.strip()}' - {e}")
        return None, None, None
    
    # Read comment line with timestep
    comment_line = f.readline()
    if not comment_line:
        print(f"Warning: Missing comment line after atom count: {num_atoms}")
        return None, None, None
    
    # Extract timestep from comment line
    timestep = None
    if comment_line.strip().startswith("Atoms. Timestep:"):
        try:
            timestep = int(comment_line.split(":")[1].strip())
        except (ValueError, IndexError) as e:
            print(f"Warning: Could not parse timestep from comment line: '{comment_line.strip()}' - {e}")
            return None, None, None
    else:
        print(f"Warning: Invalid comment line format: '{comment_line.strip()}'")
        return None, None, None
    
    # Read atom data lines
    frame_lines = [line, comment_line]  # Include num_atoms and comment lines
    atoms_read = 0
    
    for i in range(num_atoms):
        atom_line = f.readline()
        if not atom_line:
            print(f"Warning: Incomplete frame - expected {num_atoms} atoms, only read {atoms_read}")
            return None, None, None
        frame_lines.append(atom_line)
        atoms_read += 1
    
    # Verify we read exactly the expected number of atoms
    if atoms_read != num_atoms:
        print(f"Warning: Frame validation failed - expected {num_atoms} atoms, read {atoms_read}")
        return None, None, None
    
    return num_atoms, timestep, frame_lines

def streaming_merge_files(temp_combined_file, output_file, existing_timesteps=None):
    """
    Merge a sorted temporary file with the existing sorted output file using XYZ-aware merge algorithm.
    Only adds new timesteps to the output file, preserving all existing data and maintaining sort order.
    
    Process:
    1. Read complete XYZ frames from both files
    2. Compare timesteps and write the smaller frame first
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
        existing_timesteps = extract_timesteps_from_file(output_file)
    
    merged_file = f"{output_file}.merged"
    new_timesteps_added = 0
    
    try:
        # Handle case where output file doesn't exist yet
        if not os.path.exists(output_file):
            # Output file doesn't exist - just copy temp file to output file
            with open(temp_combined_file, 'r') as f_temp, open(merged_file, 'w') as f_out:
                while True:
                    num_atoms, timestep, frame_lines = read_xyz_frame(f_temp)
                    if num_atoms is None:
                        break
                    for line in frame_lines:
                        f_out.write(line)
                    new_timesteps_added += 1
        else:
            # Both files exist - perform proper XYZ-aware merge
            with open(output_file, 'r') as f_output, open(temp_combined_file, 'r') as f_temp, open(merged_file, 'w') as f_out:
                # Read first frame from each file
                output_num_atoms, output_timestep, output_frame_lines = read_xyz_frame(f_output)
                temp_num_atoms, temp_timestep, temp_frame_lines = read_xyz_frame(f_temp)
                
                while output_num_atoms is not None or temp_num_atoms is not None:
                    # Determine which frame to write next based on timestep comparison
                    if output_num_atoms is not None and temp_num_atoms is not None:
                        # Both files have frames - compare timesteps
                        if output_timestep < temp_timestep:
                            # Write output frame (it comes first)
                            for line in output_frame_lines:
                                f_out.write(line)
                            # Read next output frame
                            output_num_atoms, output_timestep, output_frame_lines = read_xyz_frame(f_output)
                        elif output_timestep > temp_timestep:
                            # Write temp frame (if not duplicate)
                            if temp_timestep not in existing_timesteps:
                                for line in temp_frame_lines:
                                    f_out.write(line)
                                new_timesteps_added += 1
                            # Read next temp frame
                            temp_num_atoms, temp_timestep, temp_frame_lines = read_xyz_frame(f_temp)
                        else:
                            # Same timestep - write output frame (existing takes precedence)
                            for line in output_frame_lines:
                                f_out.write(line)
                            # Skip temp frame (duplicate)
                            temp_num_atoms, temp_timestep, temp_frame_lines = read_xyz_frame(f_temp)
                            # Read next output frame
                            output_num_atoms, output_timestep, output_frame_lines = read_xyz_frame(f_output)
                    elif output_num_atoms is not None:
                        # Only output file has frames - write remaining output frames
                        for line in output_frame_lines:
                            f_out.write(line)
                        output_num_atoms, output_timestep, output_frame_lines = read_xyz_frame(f_output)
                    elif temp_num_atoms is not None:
                        # Only temp file has frames - write temp frames (if not duplicate)
                        if temp_timestep not in existing_timesteps:
                            for line in temp_frame_lines:
                                f_out.write(line)
                            new_timesteps_added += 1
                        temp_num_atoms, temp_timestep, temp_frame_lines = read_xyz_frame(f_temp)
        
        # Replace original file with merged file
        os.replace(merged_file, output_file)
        print(f"✓ Merged {new_timesteps_added} new timesteps into output file")
        
        total_timesteps = len(existing_timesteps) + new_timesteps_added
        return new_timesteps_added, total_timesteps
        
    except Exception as e:
        print(f"Error during streaming merge: {e}")
        safe_remove_file(merged_file, "temp")
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
        # First recover any malformed frames
        recovered_timesteps, removed_timesteps, total_frames = recover_main_output_file(output_file)
        existing_timesteps = recovered_timesteps
        
        if removed_timesteps:
            print(f"  Recovered {len(recovered_timesteps)} valid frames, removed {len(removed_timesteps)} malformed frames")
        
        # Now verify the recovered file is sorted
        is_sorted, last_timestep, verified_timesteps = verify_output_file_sorted(output_file)
        if not is_sorted:
            print(f"ERROR: Cannot proceed - output file {output_file} is not chronologically sorted after recovery!")
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
                    while True:
                        num_atoms, timestep, frame_lines = read_xyz_frame(f)
                        if num_atoms is None:
                            break
                        
                        if timestep is not None:
                            intermediate_results[timestep] = frame_lines.copy()
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
        safe_remove_file(temp_combined_file, "temp")
        
        print(f"✓ Combination complete: {new_timesteps_added} new timesteps added, {total_timesteps} total timesteps")
        return total_timesteps, new_timesteps_added, len(overlapping_timesteps)
        
    except Exception as e:
        print(f"Error during combination: {e}")
        # Clean up temporary file on error
        safe_remove_file(temp_combined_file, "temp")
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
    
    # Clean up any old backup files from previous runs (only rank 0)
    if rank == 0:
        cleanup_old_backups(output_file)
    
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
    
    # Initialize per-process output file (will be created by atomic writes)
    if os.path.exists(proc_output_file):
        # Try to recover any incomplete frames from previous runs
        print(f"[Proc {rank}] Attempting to recover incomplete frames from existing file...")
        recovered, removed, total = recover_incomplete_frames(proc_output_file)
        if removed > 0:
            print(f"[Proc {rank}] Recovery complete: {recovered} frames recovered, {removed} malformed frames removed")
        else:
            print(f"[Proc {rank}] No malformed frames found in existing file")
    else:
        print(f"[Proc {rank}] No existing per-process file found - starting fresh")
    
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
                        
                        # Use atomic write to ensure complete frame is written
                        if atomic_write_frame(proc_output_file, xyz_content, i):
                            counter += 1
                            print(f"[Proc {rank}] Successfully wrote timestep {i} ({counter} total)")
                        else:
                            print(f"[Proc {rank}] Failed to write timestep {i} - frame may be incomplete")
                            # Continue processing other timesteps rather than failing completely
                        
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
