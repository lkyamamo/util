#!/usr/bin/env python3
"""
Combine per-temperature msd.py diffusion results into one CSV.

Each --entry CELSIUS:KELVIN:ANALYSIS_DIR points at a per-temperature MSD
analysis directory, in which msd.py wrote a dated diffusion CSV:

    <YYYYMMDD>_diffusion.csv
    label,D_1e-5_cm2_s
    O,0.123456
    H,0.234567
    total,0.180000

The newest dated file wins, matching how the rest of the pipeline treats
re-run analyses (outputs are added to a directory, never replaced).

Reading msd.py's CSV rather than scraping its stdout is deliberate: the
D values are a structured output, and a regex over log text would break the
first time someone reformats a print statement.

Usage: called by submit_temperature_sweep.sh — not run standalone.
"""

import argparse
import csv
import sys
from pathlib import Path

DIFFUSION_GLOB = "[0-9]*_diffusion.csv"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument(
        "--entry", action="append", required=True,
        metavar="CELSIUS:KELVIN:ANALYSIS_DIR",
        help="Repeatable. One per temperature.",
    )
    p.add_argument("--output", required=True)
    return p.parse_args()


def parse_diffusion_dir(path):
    """{label: D} from the newest <date>_diffusion.csv in an analysis dir."""
    directory = Path(path)
    if not directory.is_dir():
        raise ValueError(f"not a directory: {directory}")

    candidates = sorted(directory.glob(DIFFUSION_GLOB))
    if not candidates:
        raise ValueError(
            f"no {DIFFUSION_GLOB} in {directory} — did msd.py run there?"
        )
    newest = candidates[-1]

    values = {}
    with open(newest, newline="") as f:
        reader = csv.DictReader(f)
        if reader.fieldnames is None or "label" not in reader.fieldnames:
            raise ValueError(f"{newest}: missing a 'label' column")
        value_field = next(
            (c for c in reader.fieldnames if c.startswith("D_")), None
        )
        if value_field is None:
            raise ValueError(f"{newest}: no D_* column among {reader.fieldnames}")
        for row in reader:
            values[row["label"]] = float(row[value_field])

    if not values:
        raise ValueError(f"{newest}: no rows")
    return values


def label_sort_key(label):
    """Element labels alphabetically, 'total' last."""
    return (1, "") if label == "total" else (0, label)


def main():
    args = parse_args()

    rows = []
    errors = []
    for entry in args.entry:
        parts = entry.split(":")
        if len(parts) != 3:
            errors.append(
                f"Malformed --entry (expected CELSIUS:KELVIN:DIR): {entry}"
            )
            continue
        celsius_str, kelvin_str, analysis_dir = parts
        try:
            values = parse_diffusion_dir(analysis_dir)
        except (OSError, ValueError) as e:
            errors.append(f"{entry}: {e}")
            continue
        rows.append((float(celsius_str), float(kelvin_str), values))

    if errors:
        for e in errors:
            print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    # Columns come from the union across every temperature, so one run that
    # saw an extra species does not silently drop it from the table. A
    # temperature missing a label gets a blank cell, not a zero — zero is a
    # physically meaningful diffusion coefficient and would be a lie here.
    labels = sorted({lbl for _, _, v in rows for lbl in v}, key=label_sort_key)

    rows.sort(key=lambda r: r[0])

    with open(args.output, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(
            ["temperature_C", "temperature_K"] + [f"D_{lbl}" for lbl in labels]
        )
        for celsius, kelvin, values in rows:
            writer.writerow(
                [celsius, kelvin]
                + [values.get(lbl, "") for lbl in labels]
            )

    print(f"Wrote {len(rows)} row(s) to {args.output}")
    print("D values are in 1e-5 cm^2/s (msd.py's Einstein-relation fit).")


if __name__ == "__main__":
    main()
