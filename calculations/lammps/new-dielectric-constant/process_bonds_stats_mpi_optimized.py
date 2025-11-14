#!/usr/bin/env python3
"""
OPTIMIZED MPI-parallel version for collecting bond statistics from molecular dynamics trajectories.

This version uses MPI to distribute frame processing across multiple nodes.
Each MPI process handles a subset of frames independently and collects bond statistics.

Key Optimizations for Large-Scale Processing (1M+ frames):
- Incremental statistics: Only stores running sums, not individual values
- Memory-efficient: O(1) memory usage regardless of frame count
- Optimized MPI communication: Only sends summary statistics, not raw data
- Streaming processing: No large arrays in memory
- Periodic statistics reporting: Shows progress during long runs

Requirements:
    mpi4py (pip install mpi4py)

Usage:
    mpirun -np 64 python process_bonds_stats_mpi_optimized.py <input_file> <start> <end> <increment> 
           <type_A> <type_B> <cutoff> <la> <lb> <lc> <output_file>
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

def read_frame(f):
    """Read a single frame from file handle. Returns (atoms, num_atoms, timestep) or (None, None, None) at EOF."""
    atoms = {}
    
    line = f.readline()
    if not line or not line.strip():
        return None, None, None
    
    num_atoms = int(line.strip())
    timestep_line = f.readline().strip()
    timestep = int(timestep_line.split(':')[1].strip())
    
    for i in range(1, num_atoms + 1):
        line = f.readline()
        if not line:
            raise ValueError(f"Unexpected end of file while reading atom {i}")
        
        parts = line.strip().split()
        atoms[i] = [int(parts[0]), float(parts[1]), float(parts[2]), float(parts[3])]
    
    return atoms, num_atoms, timestep


def create_3atom_bonds(atoms, type_A, type_B, cutoff, box_dims):
    """
    Create 3-atom molecular bonds (like water molecules).
    
    Algorithm:
    1. Find all type_A atoms (e.g., oxygen)
    2. For each type_A, find closest 2 type_B atoms (e.g., hydrogens)
    3. Create a bond with these 3 atoms
    
    Returns:
    --------
    bonds : dict
        {bondID: [atomID1, atomID2, atomID3]}
        where atomID1 is type_A, atomID2 and atomID3 are type_B
    bond_stats : dict
        Statistics about bond formation including:
        - 'bonds_created': number of bonds successfully created
        - 'bonds_missing': number of type_A atoms that couldn't form bonds
        - 'atoms_within_cutoff': total count of type_B atoms within cutoff of any type_A
        - 'type_A_count': total number of type_A atoms
        - 'type_B_count': total number of type_B atoms
    """
    bonds = {}
    bond_id = 1
    cutoff_sq = cutoff * cutoff
    
    # Statistics tracking
    bonds_created = 0
    bonds_missing = 0
    atoms_within_cutoff = 0
    type_A_count = 0
    type_B_count = 0
    
    atom_ids = list(atoms.keys())
    la, lb, lc = box_dims
    
    # Group atoms by type
    type_A_atoms = []
    type_B_atoms = []
    
    for atom_id in atom_ids:
        atom_type, _, _, _ = atoms[atom_id]
        if atom_type == type_A:
            type_A_atoms.append(atom_id)
            type_A_count += 1
        elif atom_type == type_B:
            type_B_atoms.append(atom_id)
            type_B_count += 1
    
    # For each type_A atom, find closest 2 type_B atoms
    for atom_id1 in type_A_atoms:
        _, x1, y1, z1 = atoms[atom_id1]
        
        # Find distances to all type_B atoms
        distances = []
        type_B_within_cutoff_count = 0  # Count for this specific type_A atom
        
        for atom_id2 in type_B_atoms:
            _, x2, y2, z2 = atoms[atom_id2]
            
            # Calculate distance with periodic boundary conditions
            dx = x2 - x1
            dy = y2 - y1
            dz = z2 - z1
            
            if dx >= la/2.0: dx -= la
            elif dx <= -la/2.0: dx += la
            if dy >= lb/2.0: dy -= lb
            elif dy <= -lb/2.0: dy += lb
            if dz >= lc/2.0: dz -= lc
            elif dz <= -lc/2.0: dz += lc
            
            dist_sq = dx*dx + dy*dy + dz*dz
            distances.append((dist_sq, atom_id2))
            
            # Count type B atoms within cutoff distance for this type_A atom
            if dist_sq <= cutoff_sq:
                type_B_within_cutoff_count += 1
        
        # Add to total count (for average calculation)
        atoms_within_cutoff += type_B_within_cutoff_count
        
        # Sort by distance and take closest 2
        distances.sort()
        
        # Check if closest 2 are within cutoff
        if len(distances) >= 2:
            dist1, atom_id2 = distances[0]
            dist2, atom_id3 = distances[1]
            
            if dist1 <= cutoff_sq and dist2 <= cutoff_sq:
                # Create 3-atom bond: [type_A, type_B1, type_B2]
                bonds[bond_id] = [atom_id1, atom_id2, atom_id3]
                bond_id += 1
                bonds_created += 1
            else:
                bonds_missing += 1
                print(f"Warning: Bond {bond_id} not created: dist_sq1 = {dist1}, dist_sq2 = {dist2}")
            
            # Check for third closest atom if it exists
            if len(distances) >= 3:
                dist3, atom_id4 = distances[2]
                if dist3 <= cutoff_sq:
                    print(f"Warning: Third closest atom {atom_id4} is within cutoff for atom {atom_id1}")

        else:
            bonds_missing += 1
            print(f"No 2 closest atoms found for atom {atom_id1}")
    
    # Calculate average type B atoms within cutoff per type A atom (molecule)
    avg_type_B_within_cutoff = atoms_within_cutoff / type_A_count if type_A_count > 0 else 0.0
    
    # Create statistics dictionary
    bond_stats = {
        'bonds_created': bonds_created,
        'bonds_missing': bonds_missing,
        'atoms_within_cutoff': avg_type_B_within_cutoff,  # Now average per molecule
        'type_A_count': type_A_count,
        'type_B_count': type_B_count
    }
    
    return bonds, bond_stats


def process_frame(atoms, num_atoms, timestep, type_A, type_B, cutoff, box_dims):
    """Process a single frame and return bond statistics."""
    bonds, bond_stats = create_3atom_bonds(atoms, type_A, type_B, cutoff, box_dims)
    
    # Total bonds = created + missing (all potential bonds)
    total_bonds = bond_stats['bonds_created'] + bond_stats['bonds_missing']
    
    return timestep, total_bonds, bond_stats


class SimpleAverage:
    """
    Simple average calculator for memory efficiency.
    Only tracks sum and count for mean calculation.
    Memory usage: O(1) regardless of number of samples.
    """
    
    def __init__(self):
        self.count = 0
        self.sum = 0.0
    
    def add_value(self, value):
        """Add a new value to the average."""
        self.count += 1
        self.sum += value
    
    def get_mean(self):
        """Get the mean value."""
        return self.sum / self.count if self.count > 0 else 0.0
    
    def get_count(self):
        """Get the number of values added."""
        return self.count


def get_timesteps_to_process(start, end, increment, output_file):
    """
    Calculate which timesteps should be processed and which are already done.
    
    Returns:
    - timesteps_to_process: sorted list of timesteps that need processing
    - already_processed: set of timesteps already in output file
    """
    # Calculate all timesteps that SHOULD be processed
    all_timesteps = set(range(start, end+increment, increment))
    
    # Load already processed timesteps from output file
    already_processed = load_processed_timesteps(output_file)
    
    # Calculate which timesteps NEED processing
    timesteps_to_process = sorted(all_timesteps - already_processed)
    
    return timesteps_to_process, already_processed


def load_processed_timesteps(output_file):
    """
    Load already processed timesteps from output file.
    Returns set of timesteps that have been completed.
    
    Also checks per-process files (.proc*) if they exist.
    """
    processed = set()
    
    # Check main output file
    if os.path.exists(output_file):
        try:
            with open(output_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 7:  # timestep + 6 bond statistics
                            try:
                                timestep = int(parts[0])
                                processed.add(timestep)
                            except ValueError:
                                pass
        except Exception as e:
            print(f"[Proc {rank}] Warning: Could not read main output file: {e}")
    
    # Check per-process files (for MPI restart)
    proc_file = f"{output_file}.proc{rank}"
    if os.path.exists(proc_file):
        try:
            with open(proc_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 7:
                            try:
                                timestep = int(parts[0])
                                processed.add(timestep)
                            except ValueError:
                                pass
        except Exception as e:
            print(f"[Proc {rank}] Warning: Could not read per-process file: {e}")
    
    return processed


def gather_simple_statistics(local_stats, comm, rank, size):
    """
    Gather simple statistics from all MPI processes and calculate global averages.
    Only calculates the requested statistics: average type B neighbors, average failed bonds,
    average created bonds, and percentage of failed bonds.
    
    Returns:
        dict: Global statistics with only the requested metrics
    """
    # Prepare summary data for each statistic
    summary_data = {}
    for key in ['bonds_created', 'bonds_missing', 'atoms_within_cutoff']:
        stats = local_stats[key]
        summary_data[key] = {
            'count': stats.get_count(),
            'sum': stats.sum
        }
    
    # Gather summary data from all processes
    all_summaries = comm.gather(summary_data, root=0)
    
    if rank == 0 and all_summaries:
        # Combine statistics from all processes
        combined_stats = {}
        for key in ['bonds_created', 'bonds_missing', 'atoms_within_cutoff']:
            total_count = 0
            total_sum = 0.0
            
            for proc_summary in all_summaries:
                if proc_summary and key in proc_summary:
                    data = proc_summary[key]
                    total_count += data['count']
                    total_sum += data['sum']
            
            # Calculate combined average
            mean_val = total_sum / total_count if total_count > 0 else 0.0
            combined_stats[key] = {
                'count': total_count,
                'mean': mean_val,
                'total_sum': total_sum
            }
        
        # Calculate percentage of failed bonds over all possible bonds
        total_created = combined_stats['bonds_created']['total_sum']
        total_missing = combined_stats['bonds_missing']['total_sum']
        total_possible = total_created + total_missing
        failed_percentage = (total_missing / total_possible * 100) if total_possible > 0 else 0.0
        
        return {
            'avg_type_b_neighbors': combined_stats['atoms_within_cutoff']['mean'],
            'avg_failed_bonds': combined_stats['bonds_missing']['mean'],
            'avg_created_bonds': combined_stats['bonds_created']['mean'],
            'failed_bond_percentage': failed_percentage,
            'total_frames': combined_stats['bonds_created']['count']
        }
    else:
        return None


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
                    if len(parts) >= 7:
                        try:
                            current_timestep = int(parts[0])
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
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 7:
                            try:
                                timestep = int(parts[0])
                                existing_timesteps.add(timestep)
                            except ValueError:
                                continue
    
    # Simple line-by-line merge approach
    merged_file = f"{output_file}.merged"
    new_timesteps_added = 0
    
    try:
        # Handle case where output file doesn't exist yet
        if os.path.exists(output_file):
            with open(output_file, 'r') as f_output, open(temp_combined_file, 'r') as f_temp, open(merged_file, 'w') as f_out:
                # Read first line from each file
                output_line = f_output.readline()
                temp_line = f_temp.readline()
                
                # Parse timestep from output line
                output_timestep = None
                if output_line:
                    line = output_line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 7:
                            try:
                                output_timestep = int(parts[0])
                            except ValueError:
                                pass
                
                # Main merge loop
                while output_line and temp_line:
                    # Parse temp line
                    temp_line_clean = temp_line.strip()
                    temp_timestep = None
                    temp_stats = None
                    
                    if temp_line_clean and not temp_line_clean.startswith('#'):
                        parts = temp_line_clean.split()
                        if len(parts) >= 7:
                            try:
                                temp_timestep = int(parts[0])
                                temp_stats = parts[1:7]  # All 6 statistics
                            except ValueError:
                                pass
                    
                    # Skip temp line if invalid or already exists
                    if temp_timestep is None or temp_timestep in existing_timesteps:
                        temp_line = f_temp.readline()
                        continue
                    
                    # Compare timesteps and decide what to write
                    if output_timestep is None or temp_timestep < output_timestep:
                        # Write temp line (it comes before current output line)
                        f_out.write(f"{temp_timestep}  {'  '.join(temp_stats)}\n")
                        new_timesteps_added += 1
                        temp_line = f_temp.readline()  # Move to next temp line
                        
                    elif temp_timestep == output_timestep:
                        # Timestep already exists, skip temp line
                        temp_line = f_temp.readline()
                        
                    else:
                        # Write output line (it comes before current temp line)
                        f_out.write(output_line)
                        output_line = f_output.readline()  # Move to next output line
                        
                        # Parse timestep from new output line
                        output_timestep = None
                        if output_line:
                            line = output_line.strip()
                            if line and not line.startswith('#'):
                                parts = line.split()
                                if len(parts) >= 7:
                                    try:
                                        output_timestep = int(parts[0])
                                    except ValueError:
                                        pass
                
                # Write remaining output lines
                while output_line:
                    f_out.write(output_line)
                    output_line = f_output.readline()
                
                # Write remaining temp lines (skip duplicates)
                while temp_line:
                    temp_line_clean = temp_line.strip()
                    if temp_line_clean and not temp_line_clean.startswith('#'):
                        parts = temp_line_clean.split()
                        if len(parts) >= 7:
                            try:
                                temp_timestep = int(parts[0])
                                temp_stats = parts[1:7]
                                if temp_timestep not in existing_timesteps:
                                    f_out.write(f"{temp_timestep}  {'  '.join(temp_stats)}\n")
                                    new_timesteps_added += 1
                            except ValueError:
                                pass
                    temp_line = f_temp.readline()
        else:
            # Output file doesn't exist - just copy temp file to output file
            with open(temp_combined_file, 'r') as f_temp, open(merged_file, 'w') as f_out:
                for line in f_temp:
                    line_clean = line.strip()
                    if line_clean and not line_clean.startswith('#'):
                        parts = line_clean.split()
                        if len(parts) >= 7:
                            try:
                                temp_timestep = int(parts[0])
                                temp_stats = parts[1:7]
                                f_out.write(f"{temp_timestep}  {'  '.join(temp_stats)}\n")
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
                            if len(parts) >= 7:
                                timestep = int(parts[0])
                                stats = parts[1:7]  # All 6 statistics
                                intermediate_results[timestep] = stats  # Last write wins
                                proc_file_sources[timestep] = proc_file  # Track source
                                
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
        for timestep, stats in sorted_proc_results:
            f_temp.write(f"{timestep}  {'  '.join(stats)}\n")
    
    # Step 4: Use streaming merge to combine with existing output file
    print(f"  Performing streaming merge with existing output file...")
    try:
        new_timesteps_added, total_timesteps = streaming_merge_files(temp_combined_file, output_file, existing_timesteps)
        
        # Step 5: Log any overlapping timesteps (should be rare with proper duplicate detection)
        overlapping_timesteps = []
        for timestep, stats in intermediate_results.items():
            if timestep in existing_timesteps:
                # This timestep was in both existing and new data
                overlapping_timesteps.append({
                    'timestep': timestep,
                    'source_file': proc_file_sources.get(timestep, 'unknown'),
                    'new_stats': stats
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
                    f_log.write(f"  New data: {'  '.join(overlap['new_stats'])}\n")
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


def process_concatenated_file_mpi(filename, start, end, increment, type_A, type_B, 
                                   cutoff, box_dims, output_file):
    """
    Process concatenated frames using MPI parallelization.
    Master process coordinates, workers process frames.
    
    Features:
    - Incremental writing: each process writes to its own file
    - Restart capability: skips already processed timesteps
    - Periodic combination: combines results periodically during processing
    - Final combination: master process combines all results at the end
    - Memory-efficient running statistics: O(1) memory usage regardless of frame count
    - Periodic statistics reporting: shows progress during long runs
    """
    # Create per-process output files for incremental writing
    proc_output_file = f"{output_file}.proc{rank}"
    
    # Clean up old per-process file to avoid duplicates from previous runs
    if os.path.exists(proc_output_file):
        os.remove(proc_output_file)
        print(f"[Proc {rank}] Removed old per-process file to avoid duplicates", flush=True)
    
    # Calculate combination interval based on expected workload
    # For 1M frames with 64 processes: ~15,625 frames per process
    # Combine every 5% of expected work or every 1000 frames, whichever is larger
    total_frames_expected = len(range(start, end, increment))
    frames_per_process = max(1, total_frames_expected // size)
    combination_interval = max(1000, frames_per_process // 20)  # Every 5% or 1000 frames
    
    # Statistics reporting interval (show progress during long runs)
    stats_report_interval = max(100, frames_per_process // 100)  # Every 1% or 100 frames
    
    # All processes: Determine which timesteps need processing
    if rank == 0:
        print("=" * 70)
        print("MPI Parallel Bond Statistics Collection - 3-Atom Molecular Bonds (OPTIMIZED)")
        print("=" * 70)
        print(f"Input file: {filename}")
        print(f"Timesteps: {start} to {end} (step {increment})")
        print(f"Bond types: {type_A} (central) + 2*{type_B} (terminal)")
        print(f"Cutoff: {cutoff}")
        print(f"Box dimensions: {box_dims[0]} x {box_dims[1]} x {box_dims[2]}")
        print(f"Output file: {output_file}")
        print(f"Per-process files: {output_file}.proc*")
        print(f"MPI processes: {size}")
        print(f"Memory optimization: Incremental statistics (O(1) memory)")
        print(f"Statistics reporting: Every {stats_report_interval} frames per process")
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
            print(f"Combining results every {combination_interval} frames per process", flush=True)
    
    # Create a dictionary for fast timestep lookup (optimization)
    timestep_to_index = {ts: idx for idx, ts in enumerate(timesteps_to_process)}
    
    # Initialize simple average tracking for only the required statistics
    running_stats = {
        'bonds_created': SimpleAverage(),
        'bonds_missing': SimpleAverage(),
        'atoms_within_cutoff': SimpleAverage()
    }
    
    # Open per-process output file in write mode for incremental writing
    f_out = open(proc_output_file, 'w')
    
    # Initialize counters before try block to ensure they're always defined
    counter = 0
    frame_count = 0
    last_combination_counter = 0
    last_stats_report_counter = 0
    
    try:
        # Each process reads through the file and processes its assigned timesteps
        
        with open(filename, 'r') as f_in:
            while True:
                atoms, num_atoms, timestep = read_frame(f_in)
                if atoms is None:
                    break
                
                frame_count += 1
                
                # Check if this timestep needs processing and is assigned to this process
                if timestep in timestep_to_index:
                    # Assign timesteps to processes using round-robin
                    timestep_index = timestep_to_index[timestep]
                    if timestep_index % size == rank:
                        try:
                            result = process_frame(atoms, num_atoms, timestep, type_A, type_B,
                                                  cutoff, box_dims)
                            if result is not None:
                                timestep, num_bonds, bond_stats = result
                                
                                print(f"[Proc {rank}] Timestep {timestep}: {num_bonds} bonds, "
                                      f"created: {bond_stats['bonds_created']}, "
                                      f"missing: {bond_stats['bonds_missing']}, "
                                      f"within_cutoff: {bond_stats['atoms_within_cutoff']}, "
                                      f"type_A: {bond_stats['type_A_count']}, "
                                      f"type_B: {bond_stats['type_B_count']}")
                                
                                # Accumulate only the required statistics
                                running_stats['bonds_created'].add_value(bond_stats['bonds_created'])
                                running_stats['bonds_missing'].add_value(bond_stats['bonds_missing'])
                                running_stats['atoms_within_cutoff'].add_value(bond_stats['atoms_within_cutoff'])
                                
                                # Write immediately with timestep for restart capability
                                # Format: timestep num_bonds bonds_created bonds_missing atoms_within_cutoff type_A_count type_B_count
                                f_out.write(f"{timestep}  {num_bonds}  {bond_stats['bonds_created']}  "
                                          f"{bond_stats['bonds_missing']}  {bond_stats['atoms_within_cutoff']}  "
                                          f"{bond_stats['type_A_count']}  {bond_stats['type_B_count']}\n")
                                f_out.flush()  # Ensure data is written to disk
                                counter += 1
                                
                                # Periodic statistics reporting (show progress)
                                if counter - last_stats_report_counter >= stats_report_interval:
                                    print(f"[Proc {rank}] Progress: {counter} frames processed, "
                                          f"avg created: {running_stats['bonds_created'].get_mean():.2f}, "
                                          f"avg failed: {running_stats['bonds_missing'].get_mean():.2f}")
                                    last_stats_report_counter = counter
                                
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
                                
                        except Exception as e:
                            print(f"[Proc {rank}] Error processing timestep {timestep}: {e}")
                            import traceback
                            traceback.print_exc()
        
        print(f"[Proc {rank}] Processed {counter} frames")
    
    finally:
        f_out.close()
    
    # Synchronize all processes before final combination
    comm.Barrier()
    
    # Gather and print running statistics from all processes
    print(f"\n" + "="*70)
    print("BOND STATISTICS SUMMARY")
    print("="*70)
    
    # Gather statistics from all processes
    global_stats = gather_simple_statistics(running_stats, comm, rank, size)
    
    if rank == 0 and global_stats:
        print(f"Total frames processed: {global_stats['total_frames']}")
        print(f"")
        print(f"BOND STATISTICS SUMMARY:")
        print(f"  Average type B atoms within cutoff per molecule: {global_stats['avg_type_b_neighbors']:8.2f}")
        print(f"  Average failed bonds per frame:         {global_stats['avg_failed_bonds']:8.2f}")
        print(f"  Average created bonds per frame:        {global_stats['avg_created_bonds']:8.2f}")
        print(f"  Percentage of failed bonds:             {global_stats['failed_bond_percentage']:8.2f}%")
        print("="*70)
    elif rank == 0:
        print("No bond statistics data available.")
        print("="*70)
    
    # Synchronize after statistics printing
    comm.Barrier()
    
    # Master process does final combination
    if rank == 0:
        print("\nPerforming final combination of per-process files...", flush=True)
        total_frames, new_timesteps, overlapping_timesteps = combine_per_process_files(output_file, size)
        
        print(f"\n✓ Processing complete!")
        print(f"  Total frames in final file: {total_frames}")
        print(f"  Added {new_timesteps} new timesteps")
        print(f"  Found {overlapping_timesteps} overlapping timesteps (existing data preserved)")
        print(f"  Output written to: {output_file}")
        if overlapping_timesteps > 0:
            print(f"  Overlap log: {output_file}.overlap_log")
        print(f"  Per-process files: {output_file}.proc* (can be deleted)")
        print(f"\nOutput format:")
        print(f"  timestep  num_bonds  bonds_created  bonds_missing  atoms_within_cutoff  type_A_count  type_B_count")


if __name__ == "__main__":
    if len(sys.argv) < 12:
        if rank == 0:
            print("Usage:")
            print("  mpirun -np <num_procs> python process_bonds_stats_mpi_optimized.py <input_file> <start> <end> <increment> "
                  "<type_A> <type_B> <cutoff> <la> <lb> <lc> <output_file>")
            print("\nExample:")
            print("  mpirun -np 64 python process_bonds_stats_mpi_optimized.py trajectory.dat 0 1000 10 1 2 1.2 "
                  "35.0 35.0 35.0 bond_stats_output")
            print("\nNote:")
            print("  - Requires mpi4py: pip install mpi4py")
            print("  - Creates 3-atom molecular bonds (type_A + 2*type_B)")
            print("  - Collects bond statistics instead of dipole moments")
            print("  - Memory optimized for large-scale processing (1M+ frames)")
            print("  - Output file automatically gets .statistics extension")
            print("  - Output format: timestep num_bonds bonds_created bonds_missing atoms_within_cutoff type_A_count type_B_count")
        sys.exit(1)
    
    filename = sys.argv[1]
    start = int(sys.argv[2])
    end = int(sys.argv[3])
    increment = int(sys.argv[4])
    type_A = int(sys.argv[5])
    type_B = int(sys.argv[6])
    cutoff = float(sys.argv[7])
    la = float(sys.argv[8])
    lb = float(sys.argv[9])
    lc = float(sys.argv[10])
    output_file = sys.argv[11]
    
    # Add .statistics extension to make it clear this is a statistics file
    if not output_file.endswith('.statistics'):
        output_file = f"{output_file}.statistics"
    
    process_concatenated_file_mpi(filename, start, end, increment, type_A, type_B,
                                   cutoff, [la, lb, lc], output_file)
