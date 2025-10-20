#!/usr/bin/env python3
import argparse
import glob
import os
import re
import csv
import sys

def append_numbered_csvs(input_dir: str, output_file: str):
    # Find all files matching the pattern
    pattern = os.path.join(input_dir, "*.thermodynamics.csv")
    files = glob.glob(pattern)
    if not files:
        print(f"No files found matching pattern: {pattern}", file=sys.stderr)
        return

    # Extract leading number and sort
    def extract_number(path):
        name = os.path.basename(path)
        m = re.match(r"(\d+)\.thermodynamics\.csv$", name)
        if not m:
            raise ValueError(f"Filename does not match expected format: {name}")
        return int(m.group(1))

    files_sorted = sorted(files, key=extract_number)

    with open(output_file, 'w', newline='') as fout:
        writer = None

        for csv_path in files_sorted:
            with open(csv_path, 'r', newline='') as fin:
                reader = csv.reader(fin)
                try:
                    header = next(reader)
                except StopIteration:
                    # skip empty files
                    continue

                # Write header once
                if writer is None:
                    writer = csv.writer(fout)
                    writer.writerow(header)

                # Write all data rows
                for row in reader:
                    writer.writerow(row)

    print(f"Appended {len(files_sorted)} files into {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Append thermodynamics CSVs named '<number>.thermodynamics.csv' in numeric order.")
    parser.add_argument(
        "input_dir",
        help="Directory containing '*.thermodynamics.csv' files.")
    parser.add_argument(
        "output",
        help="Path to the combined output CSV.")
    args = parser.parse_args()

    append_numbered_csvs(args.input_dir, args.output)
