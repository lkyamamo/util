#!/usr/bin/env python3
"""
bubble_setup.py

Generate and run a LAMMPS input script to create spherical bubble(s) in a
silica/water LAMMPS data file, with optional filler liquid insertion.

LAMMPS handles all atom deletion, broken-molecule cleanup, and topology
filtering natively and in parallel.  Python only builds the script and
reads a few header lines to detect atom-type counts.

Usage: edit the CONFIG section and run as a script,
       or import and call create_bubble() directly.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

# ─── CONFIG ──────────────────────────────────────────────────────────────────
INPUT_FILE = "/Users/loganyamamoto/Desktop/Research/grants/geo_sciences/bubble_collapse/data/systems/sioh_cov/a-SiO/2_surface/N113.4M/0005.data"

# None → <stem>_bubble<suffix> written beside INPUT_FILE
OUTPUT_FILE: str | None = None

# (x, y, z) bubble center positions in Angstroms — must be in the water medium
BUBBLE_CENTERS: list[tuple[float, float, float]] = [
    (400.0, 400.0, 881.0),
]

BUBBLE_RADIUS = 150.0  # Angstroms

# Optional filler LAMMPS data file.
# Must share the same coordinate system as INPUT_FILE.
FILLER_FILE: str | None = None  # e.g. "ethanol.data"

# atom_style used in both data files — must match what the files were written with
ATOM_STYLE = "full"

# LAMMPS executable.  For MPI runs use e.g. "mpirun -np 16 lammps"
LAMMPS_EXEC = "lammps"

# If True, save the generated .in script beside the output file for inspection.
SAVE_SCRIPT = True
# ─────────────────────────────────────────────────────────────────────────────


def default_output_path(input_file: str) -> str:
    """<stem>_bubble<suffix> in the same directory as input_file."""
    p = Path(input_file)
    return str(p.with_name(f"{p.stem}_bubble{p.suffix}"))


def read_header(data_file: str) -> dict:
    """
    Stream just the header block of a LAMMPS data file.

    Returns a dict with integer values for any of:
        n_atoms, n_bonds, n_angles, n_dihedrals, n_impropers,
        n_atom_types, n_bond_types, n_angle_types,
        n_dihedral_types, n_improper_types
    """
    info: dict[str, int] = {}
    with open(data_file) as fh:
        next(fh)  # skip comment line
        for line in fh:
            stripped = line.split("#")[0].strip()
            if not stripped:
                continue
            # first alphabetic section keyword signals end of header
            if stripped[0].isalpha() and "types" not in stripped:
                break
            parts = stripped.split()
            if len(parts) >= 2:
                if parts[1] == "atoms":          info["n_atoms"]          = int(parts[0])
                elif parts[1] == "bonds":        info["n_bonds"]          = int(parts[0])
                elif parts[1] == "angles":       info["n_angles"]         = int(parts[0])
                elif parts[1] == "dihedrals":    info["n_dihedrals"]      = int(parts[0])
                elif parts[1] == "impropers":    info["n_impropers"]      = int(parts[0])
            if "atom types"     in stripped: info["n_atom_types"]     = int(parts[0])
            elif "bond types"   in stripped: info["n_bond_types"]     = int(parts[0])
            elif "angle types"  in stripped: info["n_angle_types"]    = int(parts[0])
            elif "dihedral types" in stripped: info["n_dihedral_types"] = int(parts[0])
            elif "improper types" in stripped: info["n_improper_types"] = int(parts[0])
    return info


def _bubble_regions(centers: list[tuple], radius: float) -> str:
    """
    Return LAMMPS region commands defining:
        'bubble'         — union of all bubble spheres (side in)
        'outside_bubble' — atoms outside ALL bubble spheres (side out)

    For a single center these are simply a pair of sphere regions.
    For N centers, 'bubble' is a union and 'outside_bubble' is the
    intersection of N 'side out' spheres (= outside every sphere).
    """
    r = float(radius)
    lines: list[str] = []

    if len(centers) == 1:
        x, y, z = centers[0]
        lines.append(f"region bubble         sphere {x} {y} {z} {r} side in")
        lines.append(f"region outside_bubble sphere {x} {y} {z} {r} side out")
    else:
        for i, (x, y, z) in enumerate(centers):
            lines.append(f"region bubble_{i} sphere {x} {y} {z} {r} side in")
        in_ids  = " ".join(f"bubble_{i}"  for i in range(len(centers)))
        lines.append(f"region bubble union {len(centers)} {in_ids}")

        for i, (x, y, z) in enumerate(centers):
            lines.append(f"region outside_{i} sphere {x} {y} {z} {r} side out")
        out_ids = " ".join(f"outside_{i}" for i in range(len(centers)))
        lines.append(f"region outside_bubble intersect {len(centers)} {out_ids}")

    return "\n".join(lines)


def generate_script(
    input_file:     str,
    output_file:    str,
    centers:        list[tuple],
    radius:         float,
    atom_style:     str,
    filler_file:    str | None,
) -> str:
    """Build and return the full LAMMPS input script as a string."""

    lines: list[str] = [
        "# Generated by bubble_setup.py — do not edit by hand",
        "",
        f"atom_style {atom_style}",
        "",
        f"read_data {input_file}",
        "",
        "# ── bubble region(s) ──────────────────────────────────────────",
        _bubble_regions(centers, radius),
        "",
        "# ── hollow out water (mol yes cleans up broken molecules) ─────",
        "delete_atoms region bubble mol yes",
        "",
    ]

    if filler_file:
        base_n_types   = read_header(input_file).get("n_atom_types")
        filler_n_types = read_header(filler_file).get("n_atom_types")

        if base_n_types is None:
            raise ValueError(f"Could not read n_atom_types from {input_file}")
        if filler_n_types is None:
            raise ValueError(f"Could not read n_atom_types from {filler_file}")

        filler_lo = base_n_types + 1
        filler_hi = base_n_types + filler_n_types

        lines += [
            "# ── filler insertion ──────────────────────────────────────────",
            f"# Filler atom types will be {filler_lo}–{filler_hi} after append",
            f"read_data {filler_file} add append",
            "",
            f"group filler        type {filler_lo} {filler_hi}",
            "group outside_atoms region outside_bubble",
            "",
            "# Filler atoms that landed outside all bubble(s)",
            "group filler_outside intersect filler outside_atoms",
            "",
            "# Expand to whole molecules that straddle the boundary, then remove",
            "group filler_broken  molecule filler_outside",
            "delete_atoms group filler_broken",
            "",
            "group filler         delete",
            "group outside_atoms  delete",
            "group filler_outside delete",
            "group filler_broken  delete",
            "",
        ]

    lines += [
        "# ── write result ──────────────────────────────────────────────",
        f"write_data {output_file}",
        "",
    ]

    return "\n".join(lines)


def run_lammps(script: str, script_path: Path, lammps_exec: str) -> None:
    """Write the script to disk and execute LAMMPS."""
    script_path.write_text(script)
    print(f"Script written to: {script_path}")

    # Support "mpirun -np 8 lammps" — split on whitespace
    cmd = lammps_exec.split() + ["-in", str(script_path)]
    print(f"Running: {' '.join(cmd)}\n")

    result = subprocess.run(cmd)
    if result.returncode != 0:
        raise RuntimeError(
            f"LAMMPS exited with code {result.returncode}. "
            "Check output above for errors."
        )


def create_bubble(
    input_file:  str,
    centers:     list[tuple[float, float, float]],
    radius:      float,
    output_file: str | None = None,
    atom_style:  str        = "full",
    filler_file: str | None = None,
    lammps_exec: str        = "lammps",
    save_script: bool       = True,
) -> None:
    """
    Create spherical bubble(s) in a LAMMPS silica/water data file.

    Generates a LAMMPS input script and runs it.  LAMMPS handles all
    heavy lifting: atom deletion, broken-molecule removal (mol yes),
    topology filtering, and parallel I/O.

    Args:
        input_file:  LAMMPS data file (silica slab + water).
        centers:     (x, y, z) bubble center(s) in Angstroms.
                     Each center must lie within the water medium —
                     not inside the silica slab, and at least `radius`
                     away from the silica surface.
        radius:      Bubble radius in Angstroms.
        output_file: Output path.  Defaults to <stem>_bubble<ext>.
        atom_style:  Must match the atom_style used in both data files.
        filler_file: Optional LAMMPS data file to carve and insert.
                     Must share the same coordinate system as input_file.
        lammps_exec: LAMMPS binary.  Use "mpirun -np N lammps" for MPI.
        save_script: If True, save the generated .in file beside output.
    """
    if output_file is None:
        output_file = default_output_path(input_file)

    script = generate_script(
        input_file, output_file, centers, radius, atom_style, filler_file
    )

    print("─" * 64)
    print(script)
    print("─" * 64)

    out = Path(output_file)
    script_path = out.with_suffix(".in") if save_script \
                  else Path(output_file).parent / "_bubble_setup_tmp.in"

    try:
        run_lammps(script, script_path, lammps_exec)
    finally:
        if not save_script:
            script_path.unlink(missing_ok=True)

    print(f"\nDone.  Output: {output_file}")


if __name__ == "__main__":
    create_bubble(
        input_file=INPUT_FILE,
        centers=BUBBLE_CENTERS,
        radius=BUBBLE_RADIUS,
        output_file=OUTPUT_FILE,
        atom_style=ATOM_STYLE,
        filler_file=FILLER_FILE,
        lammps_exec=LAMMPS_EXEC,
        save_script=SAVE_SCRIPT,
    )
