#!/usr/bin/env python3
"""
Combine MPI per-process output files into a single output file.

This script reads all .proc* files and combines them into the main output file,
sorted by timestep. Also combines statistics files (.statistics.proc*) into a
single statistics file with summary. Useful for:
- Manual combination if needed
- Recovery after job crash
- Recombining after manual edits

Usage:
    python combine_mpi_output.py <output_file>

Example:
    python combine_mpi_output.py dipole_output.txt

This will combine:
- dipole_output.txt.proc* → dipole_output.txt
- dipole_output.txt.statistics.proc* → dipole_output.txt.statistics
"""

import sys
import os
import glob
import datetime
from pathlib import Path


def read_dipole_file(filename):
    """Read a dipole output file and return a dictionary of {timestep: dipole}."""
    results = {}
    if not os.path.exists(filename):
        return results
    
    try:
        with open(filename, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if line and not line.startswith('#'):
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            timestep = int(parts[0])
                            dipole = [float(parts[1]), float(parts[2]), float(parts[3])]
                            results[timestep] = dipole
                        except ValueError as e:
                            print(f"Warning: Could not parse line {line_num} in {filename}: {line}")
        return results
    except Exception as e:
        print(f"Error reading {filename}: {e}")
        return {}


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
                        if len(parts) >= 4:
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
                        if len(parts) >= 4:
                            try:
                                output_timestep = int(parts[0])
                            except ValueError:
                                pass
                
                # Main merge loop
                while output_line and temp_line:
                    # Parse temp line
                    temp_line_clean = temp_line.strip()
                    temp_timestep = None
                    temp_dipole = None
                    
                    if temp_line_clean and not temp_line_clean.startswith('#'):
                        parts = temp_line_clean.split()
                        if len(parts) >= 4:
                            try:
                                temp_timestep = int(parts[0])
                                temp_dipole = [float(parts[1]), float(parts[2]), float(parts[3])]
                            except ValueError:
                                pass
                    
                    # Skip temp line if invalid or already exists
                    if temp_timestep is None or temp_timestep in existing_timesteps:
                        temp_line = f_temp.readline()
                        continue
                    
                    # Compare timesteps and decide what to write
                    if output_timestep is None or temp_timestep < output_timestep:
                        # Write temp line (it comes before current output line)
                        f_out.write(f"{temp_timestep}  {temp_dipole[0]:14.6f}  {temp_dipole[1]:14.6f}  {temp_dipole[2]:14.6f}\n")
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
                                if len(parts) >= 4:
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
                        if len(parts) >= 4:
                            try:
                                temp_timestep = int(parts[0])
                                temp_dipole = [float(parts[1]), float(parts[2]), float(parts[3])]
                                if temp_timestep not in existing_timesteps:
                                    f_out.write(f"{temp_timestep}  {temp_dipole[0]:14.6f}  {temp_dipole[1]:14.6f}  {temp_dipole[2]:14.6f}\n")
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
                        if len(parts) >= 4:
                            try:
                                temp_timestep = int(parts[0])
                                temp_dipole = [float(parts[1]), float(parts[2]), float(parts[3])]
                                f_out.write(f"{temp_timestep}  {temp_dipole[0]:14.6f}  {temp_dipole[1]:14.6f}  {temp_dipole[2]:14.6f}\n")
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


def combine_mpi_files(output_file):
    """
    Combine all .proc* files with existing output file, avoiding duplicates.
    Uses streaming merge approach for memory efficiency with large files.
    
    Args:
        output_file: Base name of output file (e.g., 'dipole_output.txt')
    
    Returns:
        Number of frames in final combined file
    """
    # Find all per-process files
    proc_files = sorted(glob.glob(f"{output_file}.proc*"))
    
    if not proc_files:
        print(f"No per-process files found for {output_file}")
        print(f"Looking for files matching: {output_file}.proc*")
        
        # Check if there's an existing output file to report on
        if os.path.exists(output_file):
            print(f"\nReading existing output file...")
            existing_results = read_dipole_file(output_file)
            if existing_results:
                print(f"  ✓ Found {len(existing_results)} timesteps in existing file")
                print(f"  Timestep range: {min(existing_results.keys())} to {max(existing_results.keys())}")
                return len(existing_results)
            else:
                print(f"  ✓ Existing file is empty")
        else:
            print(f"  ✓ No existing output file found")
        
        return 0
    
    print(f"Found {len(proc_files)} per-process files:")
    for proc_file in proc_files:
        print(f"  - {proc_file}")
    
    # Step 1: Read all .proc* files into intermediate data
    print("\nStep 1: Reading per-process files...")
    intermediate_results = {}
    seen_timesteps = {}  # Track duplicates within .proc* files
    
    for proc_file in proc_files:
        try:
            with open(proc_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 4:
                            try:
                                timestep = int(parts[0])
                                dipole = [float(parts[1]), float(parts[2]), float(parts[3])]
                                intermediate_results[timestep] = dipole  # Last write wins
                                
                                # Track duplicates
                                if timestep in seen_timesteps:
                                    seen_timesteps[timestep] += 1
                                else:
                                    seen_timesteps[timestep] = 1
                            except ValueError as e:
                                print(f"Warning: Could not parse line {line_num} in {proc_file}: {line}")
            print(f"  ✓ {proc_file}: read successfully")
        except Exception as e:
            print(f"Error reading {proc_file}: {e}")
    
    if not intermediate_results:
        print("\nNo valid data found in per-process files!")
        return 0
    
    # Report duplicates within .proc* files
    duplicates = {ts: count for ts, count in seen_timesteps.items() if count > 1}
    if duplicates:
        print(f"\n⚠ Found {len(duplicates)} duplicate timesteps in .proc* files (kept last occurrence)")
    
    # Step 2: Read existing output file (if it exists)
    print("\nStep 2: Reading existing output file...")
    existing_results = read_dipole_file(output_file)
    existing_timesteps = set(existing_results.keys()) if existing_results else set()
    
    if existing_results:
        print(f"  ✓ Found {len(existing_results)} timesteps in existing file")
    else:
        print(f"  ✓ No existing file (will create new one)")
    
    # Step 3: Create temporary combined file with sorted per-process results
    temp_combined_file = f"{output_file}.temp_combined"
    sorted_proc_results = sorted(intermediate_results.items())
    
    print(f"\nStep 3: Creating temporary combined file with {len(sorted_proc_results)} timesteps")
    with open(temp_combined_file, 'w') as f_temp:
        for timestep, dipole in sorted_proc_results:
            f_temp.write(f"{timestep}  {dipole[0]:14.6f}  {dipole[1]:14.6f}  {dipole[2]:14.6f}\n")
    
    # Step 4: Use streaming merge to combine with existing output file
    print(f"\nStep 4: Performing streaming merge with existing output file...")
    try:
        new_timesteps_added, total_timesteps = streaming_merge_files(temp_combined_file, output_file, existing_timesteps)
        
        # Clean up temporary file
        if os.path.exists(temp_combined_file):
            os.remove(temp_combined_file)
        
        print(f"\n✓ Combination complete!")
        print(f"  Added {new_timesteps_added} new timesteps")
        print(f"  Total frames in final file: {total_timesteps}")
        print(f"  Output file: {output_file}")
        
        return total_timesteps
        
    except Exception as e:
        print(f"Error during combination: {e}")
        # Clean up temporary file on error
        if os.path.exists(temp_combined_file):
            os.remove(temp_combined_file)
        raise


def combine_statistics_files(stats_file):
    """
    Combine all .statistics.proc* files into a single statistics file.
    Similar to combine_mpi_files but for statistics.
    
    Format of statistics lines: timestep  frames_processed  avg_created  avg_failed  avg_type_B_per_mol  failed_percentage
    """
    # Find all per-process statistics files
    proc_stats_files = sorted(glob.glob(f"{stats_file}.proc*"))
    
    if not proc_stats_files:
        print(f"No per-process statistics files found for {stats_file}")
        print(f"Looking for files matching: {stats_file}.proc*")
        
        # Check if there's an existing statistics file
        if os.path.exists(stats_file):
            print(f"\nReading existing statistics file...")
            with open(stats_file, 'r') as f:
                lines = [l for l in f if l.strip() and not l.strip().startswith('#')]
                if lines:
                    print(f"  ✓ Found {len(lines)} statistics entries in existing file")
                else:
                    print(f"  ✓ Existing file is empty")
        else:
            print(f"  ✓ No existing statistics file found")
        
        return 0
    
    print(f"\nFound {len(proc_stats_files)} per-process statistics files:")
    for proc_file in proc_stats_files:
        print(f"  - {proc_file}")
    
    # Read all per-process statistics files
    print("\nReading per-process statistics files...")
    all_stats_entries = []
    
    for proc_file in proc_stats_files:
        try:
            with open(proc_file, 'r', buffering=8192) as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 6:
                            try:
                                timestep = int(parts[0])
                                frames_processed = int(parts[1])
                                avg_created = float(parts[2])
                                avg_failed = float(parts[3])
                                avg_type_B_per_mol = float(parts[4])
                                failed_percentage = float(parts[5])
                                all_stats_entries.append({
                                    'timestep': timestep,
                                    'frames_processed': frames_processed,
                                    'avg_created': avg_created,
                                    'avg_failed': avg_failed,
                                    'avg_type_B_per_mol': avg_type_B_per_mol,
                                    'failed_percentage': failed_percentage
                                })
                            except ValueError as e:
                                print(f"Warning: Could not parse line {line_num} in {proc_file}: {line}")
            print(f"  ✓ {proc_file}: read successfully")
        except Exception as e:
            print(f"Error reading {proc_file}: {e}")
    
    if not all_stats_entries:
        print("\nNo valid statistics data found in per-process files!")
        return 0
    
    # Read existing statistics file to merge
    print("\nReading existing statistics file...")
    existing_entries = []
    if os.path.exists(stats_file):
        try:
            with open(stats_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('# SUMMARY'):
                        break  # Stop at summary section
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if len(parts) >= 6:
                            try:
                                timestep = int(parts[0])
                                frames_processed = int(parts[1])
                                avg_created = float(parts[2])
                                avg_failed = float(parts[3])
                                avg_type_B_per_mol = float(parts[4])
                                failed_percentage = float(parts[5])
                                existing_entries.append({
                                    'timestep': timestep,
                                    'frames_processed': frames_processed,
                                    'avg_created': avg_created,
                                    'avg_failed': avg_failed,
                                    'avg_type_B_per_mol': avg_type_B_per_mol,
                                    'failed_percentage': failed_percentage
                                })
                            except ValueError:
                                continue
            if existing_entries:
                print(f"  ✓ Found {len(existing_entries)} statistics entries in existing file")
            else:
                print(f"  ✓ No existing entries (will create new file)")
        except Exception as e:
            print(f"Warning: Could not read existing statistics file: {e}")
    else:
        print(f"  ✓ No existing file (will create new one)")
    
    # Merge existing and new entries (keep unique timesteps, new data wins)
    existing_timesteps = {entry['timestep'] for entry in existing_entries}
    merged_entries = existing_entries.copy()
    
    for entry in all_stats_entries:
        if entry['timestep'] not in existing_timesteps:
            merged_entries.append(entry)
    
    # Sort merged entries by timestep
    merged_entries.sort(key=lambda x: x['timestep'])
    
    # Calculate aggregated statistics for summary
    if merged_entries:
        # Group by timestep and aggregate (in case multiple processes reported at same timestep)
        timestep_groups = {}
        for entry in merged_entries:
            ts = entry['timestep']
            if ts not in timestep_groups:
                timestep_groups[ts] = []
            timestep_groups[ts].append(entry)
        
        # For summary, use the latest timestep's aggregated values
        latest_timestep = max(timestep_groups.keys())
        latest_entries = timestep_groups[latest_timestep]
        
        # Aggregate across all processes for the latest timestep
        total_frames_latest = sum(e['frames_processed'] for e in latest_entries)
        total_created_latest = sum(e['avg_created'] * e['frames_processed'] for e in latest_entries)
        total_failed_latest = sum(e['avg_failed'] * e['frames_processed'] for e in latest_entries)
        total_type_B_latest = sum(e['avg_type_B_per_mol'] * e['frames_processed'] for e in latest_entries)
        
        avg_created_overall = total_created_latest / total_frames_latest if total_frames_latest > 0 else 0.0
        avg_failed_overall = total_failed_latest / total_frames_latest if total_frames_latest > 0 else 0.0
        avg_type_B_per_mol_overall = total_type_B_latest / total_frames_latest if total_frames_latest > 0 else 0.0
        failed_percentage_overall = (avg_failed_overall / (avg_created_overall + avg_failed_overall) * 100) if (avg_created_overall + avg_failed_overall) > 0 else 0.0
        total_frames = total_frames_latest
    else:
        avg_created_overall = avg_failed_overall = avg_type_B_per_mol_overall = failed_percentage_overall = 0.0
        total_frames = 0
    
    # Write merged statistics file with summary
    temp_stats_file = f"{stats_file}.temp"
    try:
        with open(temp_stats_file, 'w') as f:
            # Write header
            f.write("# Bond Statistics - Periodic Reports\n")
            f.write(f"# Generated: {datetime.datetime.now().isoformat()}\n")
            f.write("# Format: timestep  frames_processed  avg_created  avg_failed  avg_type_B_per_mol  failed_percentage\n")
            f.write("#\n")
            
            # Write data entries
            for entry in merged_entries:
                f.write(f"{entry['timestep']}  {entry['frames_processed']}  "
                       f"{entry['avg_created']:.6f}  {entry['avg_failed']:.6f}  "
                       f"{entry['avg_type_B_per_mol']:.6f}  {entry['failed_percentage']:.6f}\n")
            
            # Write summary section
            f.write("\n" + "="*70 + "\n")
            f.write("# SUMMARY - Aggregated Statistics Across All Timesteps\n")
            f.write("="*70 + "\n")
            f.write(f"# Total frames processed: {total_frames}\n")
            f.write(f"# Total statistics entries: {len(merged_entries)}\n")
            f.write("#\n")
            f.write("AGGREGATED STATISTICS:\n")
            f.write(f"  Average type B atoms within cutoff per molecule: {avg_type_B_per_mol_overall:8.2f}\n")
            f.write(f"  Average failed bonds per frame:         {avg_failed_overall:8.2f}\n")
            f.write(f"  Average created bonds per frame:        {avg_created_overall:8.2f}\n")
            f.write(f"  Percentage of failed bonds:             {failed_percentage_overall:8.2f}%\n")
            f.write("="*70 + "\n")
        
        # Replace original file
        os.replace(temp_stats_file, stats_file)
        print(f"\n✓ Combined {len(merged_entries)} statistics entries (added {len(all_stats_entries)} new entries)")
        print(f"  Summary: {total_frames} total frames, avg created: {avg_created_overall:.2f}, avg failed: {avg_failed_overall:.2f}")
        
        return len(merged_entries)
        
    except Exception as e:
        print(f"Error writing statistics file: {e}")
        if os.path.exists(temp_stats_file):
            os.remove(temp_stats_file)
        raise


def main():
    """Main function."""
    if len(sys.argv) != 2:
        print("Usage: python combine_mpi_output.py <output_file>")
        print("\nExample:")
        print("  python combine_mpi_output.py dipole_output.txt")
        print("\nThis will combine:")
        print("  - dipole_output.txt.proc* → dipole_output.txt")
        print("  - dipole_output.txt.statistics.proc* → dipole_output.txt.statistics")
        sys.exit(1)
    
    output_file = sys.argv[1]
    stats_file = f"{output_file}.statistics"
    
    print("=" * 70)
    print("MPI Output File Combiner")
    print("=" * 70)
    print(f"Output file: {output_file}")
    print(f"Statistics file: {stats_file}")
    print()
    
    # Note: Script now merges with existing file instead of overwriting
    # No need to ask for confirmation - it will merge intelligently
    
    # Combine dipole files (will merge with existing if it exists)
    num_frames = combine_mpi_files(output_file)
    
    # Combine statistics files (will merge with existing if it exists)
    num_stats_entries = combine_statistics_files(stats_file)
    
    if num_frames == 0 and num_stats_entries == 0:
        print("\nNo data to combine.")
        sys.exit(1)
    
    print("\n" + "=" * 70)
    print("✓ Success!")
    print("=" * 70)
    print(f"\nYou can now delete the per-process files:")
    print(f"  rm {output_file}.proc*")
    print(f"  rm {stats_file}.proc*")
    print()


if __name__ == "__main__":
    main()

