#!/usr/bin/env python3
"""
Combine per-temperature dielectric_submit.slurm summaries into one CSV.

Each --entry TEMPERATURE:SUMMARY_FILE points at a dipole_output/summary.txt
written by dielectric_submit.slurm (the teed stdout of 2.dipole_std.py),
which ends with a line of the form:

    eps_x =    12.345678, eps_y =    12.345678, eps_z =    12.345678, eps_total =    12.345678

Usage: called by submit_dielectric_pipeline.sh — not run standalone.
"""

import argparse
import csv
import re
import sys

EPS_LINE_RE = re.compile(
    r"eps_x\s*=\s*([-\d.]+),\s*eps_y\s*=\s*([-\d.]+),\s*"
    r"eps_z\s*=\s*([-\d.]+),\s*eps_total\s*=\s*([-\d.]+)"
)


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--entry", action="append", required=True,
        metavar="TEMPERATURE:SUMMARY_FILE",
        help="Repeatable. One per temperature.",
    )
    p.add_argument("--output", required=True)
    return p.parse_args()


def parse_summary(path):
    with open(path) as f:
        lines = f.readlines()

    for line in reversed(lines):
        m = EPS_LINE_RE.search(line)
        if m:
            eps_x, eps_y, eps_z, eps_total = (float(g) for g in m.groups())
            return eps_x, eps_y, eps_z, eps_total

    raise ValueError(f"No eps_x/eps_y/eps_z/eps_total line found in {path}")


def main():
    args = parse_args()

    rows = []
    errors = []
    for entry in args.entry:
        temperature_str, _, summary_path = entry.partition(":")
        if not summary_path:
            errors.append(f"Malformed --entry (expected TEMPERATURE:PATH): {entry}")
            continue
        try:
            eps_x, eps_y, eps_z, eps_total = parse_summary(summary_path)
        except (OSError, ValueError) as e:
            errors.append(f"{entry}: {e}")
            continue
        rows.append((float(temperature_str), eps_x, eps_y, eps_z, eps_total))

    if errors:
        for e in errors:
            print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    rows.sort(key=lambda r: r[0])

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["temperature", "eps_x", "eps_y", "eps_z", "eps_total"])
        writer.writerows(rows)

    print(f"Wrote {len(rows)} row(s) to {args.output}")


if __name__ == "__main__":
    main()
