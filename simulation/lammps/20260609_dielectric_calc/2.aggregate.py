#!/usr/bin/env python3
"""
Aggregate per-frame dipole output files into a single sorted file
compatible with 2.dipoleStd.py.

Reads all dipole_*.txt files from the per-frame output directory, sorts
entries by timestep, and writes dipole_output.txt with the format:
    timestep  Mx  My  Mz

Also consolidates all warn_*.log files into a single dipole_warnings.log,
sorted by timestep, and prints a summary of how many frames had unassigned
hydrogen atoms.

Usage:
    python 1.aggregate.py
        --input-dir <dir>
        --output    <dipole_output.txt>
        --warn-log  <dipole_warnings.log>
"""

import argparse
import glob
import os
import re
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description="Aggregate per-frame dipole files into a single sorted output."
    )
    parser.add_argument("--input-dir", required=True,
                        help="Directory containing dipole_*.txt and warn_*.log files")
    parser.add_argument("--output", required=True,
                        help="Path for the combined dipole output file")
    parser.add_argument("--warn-log", required=True,
                        help="Path for the consolidated warnings log")
    parser.add_argument("--expected-frames", type=int, default=None,
                        help="If set, exit with error unless exactly this many frames are collected")
    return parser.parse_args()


def collect_dipole_files(input_dir):
    """Return sorted list of (timestep, line) tuples from all dipole_*.txt files."""
    pattern = os.path.join(input_dir, "dipole_*.txt")
    files = glob.glob(pattern)

    if not files:
        print(f"ERROR: No dipole_*.txt files found in '{input_dir}'", file=sys.stderr)
        sys.exit(1)

    rows = []
    missing = 0

    for fpath in files:
        try:
            with open(fpath) as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split()
                    if len(parts) < 4:
                        print(f"WARNING: Skipping malformed line in {fpath}: '{line}'",
                              file=sys.stderr)
                        continue
                    timestep = int(parts[0])
                    rows.append((timestep, line))
        except Exception as e:
            print(f"WARNING: Could not read '{fpath}': {e}", file=sys.stderr)
            missing += 1

    if missing:
        print(f"WARNING: {missing} file(s) could not be read.", file=sys.stderr)

    rows.sort(key=lambda r: r[0])
    return rows


def collect_warn_files(input_dir):
    """
    Return list of (timestep, line) tuples from all warn_*.log files,
    sorted by timestep. Lines without a parseable timestep are appended last.
    """
    pattern = os.path.join(input_dir, "warn_*.log")
    files = glob.glob(pattern)

    entries = []

    for fpath in files:
        try:
            with open(fpath) as f:
                for line in f:
                    line = line.rstrip("\n")
                    if not line:
                        continue
                    # Extract timestep from "WARN [timestep=NNNN]" format
                    m = re.search(r"timestep=(\d+)", line)
                    ts = int(m.group(1)) if m else -1
                    entries.append((ts, line))
        except Exception as e:
            print(f"WARNING: Could not read '{fpath}': {e}", file=sys.stderr)

    entries.sort(key=lambda e: e[0])
    return entries


def main():
    args = parse_args()

    print(f"Aggregating dipole files from: {args.input_dir}")

    rows = collect_dipole_files(args.input_dir)
    print(f"  Collected {len(rows)} frames")

    if args.expected_frames is not None and len(rows) != args.expected_frames:
        print(f"ERROR: Expected {args.expected_frames} frames, got {len(rows)}",
              file=sys.stderr)
        sys.exit(1)

    with open(args.output, "w") as out:
        for _ts, line in rows:
            out.write(line + "\n")

    print(f"  Wrote: {args.output}")

    # --- Warnings ---
    warn_entries = collect_warn_files(args.input_dir)

    if warn_entries:
        frames_with_warnings = sum(1 for ts, _ in warn_entries if ts >= 0)
        print(f"\n  WARNING SUMMARY: {frames_with_warnings} frame(s) had unassigned H atoms.")
        print(f"  See: {args.warn_log}")
        with open(args.warn_log, "w") as wf:
            wf.write(f"# Consolidated warnings from {args.input_dir}\n")
            wf.write(f"# {frames_with_warnings} frame(s) with unassigned H atoms\n")
            for _ts, line in warn_entries:
                wf.write(line + "\n")
    else:
        print("  No unassigned H atom warnings — all frames fully bonded.")
        # Write an empty log so downstream scripts can always expect the file
        with open(args.warn_log, "w") as wf:
            wf.write("# No unassigned H atom warnings\n")

    print("\nDone.")
    print(f"  Dipole output : {args.output}")
    print(f"  Warnings log  : {args.warn_log}")
    print(f"\nNext step: python 2.dipoleStd.py {args.output} <T> <la> <lb> <lc>")


if __name__ == "__main__":
    main()
