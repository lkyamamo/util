#!/usr/bin/env python3
"""Compute the mass density of a user-defined region in a LAMMPS dump frame.

Usage:
    python density.py --local-path PATH [options]

Primary options:
    --local-path PATH   Local LAMMPS dump file (single frame is read; if the
                         file has multiple frames, use --frame to pick one).
    --frame N            Frame index to read (default 0).

Region (box mode, default):
    --xlo/--xhi/--ylo/--yhi/--zlo/--zhi
                         Region bounds in Angstrom. Any bound left unset
                         defaults to that frame's box extent on that axis,
                         so e.g. only --zlo/--zhi gives a full-width slab.

Region (expression mode, overrides box mode):
    --expression EXPR    Arbitrary OVITO selection expression (e.g. a
                         sphere or cylinder). Since volume cannot be
                         derived from an arbitrary expression, --volume
                         must also be given.
    --volume V           Region volume in Angstrom^3 (required with
                         --expression).

Assumes atom positions in the dump are already wrapped into the
simulation box (no unwrap/PBC-rewrap step is performed).
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Optional

import numpy as np
from ovito.io import import_file
from ovito.modifiers import ExpressionSelectionModifier


# ---------------------------------------------------------------------------
# Element -> mass table (amu). Extend as needed.
# ---------------------------------------------------------------------------

ELEMENT_TO_MASS: Dict[str, float] = {
    "H": 1.008,
    "O": 15.999,
    "Si": 28.085,
    "Na": 22.990,
    "Cl": 35.453,
    "C": 12.011,
    "N": 14.007,
    "Ca": 40.078,
    "Al": 26.982,
    "Mg": 24.305,
    "K": 39.098,
    "S": 32.06,
}

AMU_A3_TO_G_CM3 = 1.6605  # amu/Angstrom^3 -> g/cm^3 (LAMMPS metal units)


# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------

def _resolve_type_mass_lookup(data) -> np.ndarray:
    """Return an array mapping particle type id -> mass (amu), indexed by id.

    Particle type names come from the dump's 'element' column (OVITO maps
    that column to particle type names directly).
    """
    ptypes = data.particles.particle_types.types
    mass_by_id: Dict[int, float] = {}
    unresolved = []

    for t in ptypes:
        name = t.name.strip()
        mass = ELEMENT_TO_MASS.get(name)
        if mass is None:
            for symbol, m in ELEMENT_TO_MASS.items():
                if symbol.lower() == name.lower():
                    mass = m
                    break
        if mass is None:
            unresolved.append((t.id, name))
        else:
            mass_by_id[t.id] = mass

    if unresolved:
        raise ValueError(
            f"Unresolved element(s) not in ELEMENT_TO_MASS: {unresolved}. "
            "Add them to ELEMENT_TO_MASS at the top of this file."
        )

    lookup = np.zeros(max(mass_by_id) + 1)
    for tid, mass in mass_by_id.items():
        lookup[tid] = mass
    return lookup


def compute_region_density(
    pipeline,
    *,
    frame: int = 0,
    xlo: Optional[float] = None,
    xhi: Optional[float] = None,
    ylo: Optional[float] = None,
    yhi: Optional[float] = None,
    zlo: Optional[float] = None,
    zhi: Optional[float] = None,
    expression: Optional[str] = None,
    volume: Optional[float] = None,
) -> Dict[str, float]:
    """Return density_g_cm3, n_atoms, volume_ang3, mass_sum_amu for the region."""
    if expression is not None:
        if volume is None:
            raise ValueError("--volume is required when --expression is used.")
        if any(b is not None for b in (xlo, xhi, ylo, yhi, zlo, zhi)):
            raise ValueError("--expression cannot be combined with box bound flags.")
        pipeline.modifiers.append(ExpressionSelectionModifier(expression=expression))
        data = pipeline.compute(frame)
        selected = np.asarray(data.particles["Selection"].array, dtype=bool)
        region_volume = float(volume)
    else:
        data = pipeline.compute(frame)
        cell = data.cell.matrix
        box_lo = np.array([cell[0, 3], cell[1, 3], cell[2, 3]])
        box_hi = box_lo + np.array([cell[0, 0], cell[1, 1], cell[2, 2]])

        bounds_lo = np.array([
            xlo if xlo is not None else box_lo[0],
            ylo if ylo is not None else box_lo[1],
            zlo if zlo is not None else box_lo[2],
        ])
        bounds_hi = np.array([
            xhi if xhi is not None else box_hi[0],
            yhi if yhi is not None else box_hi[1],
            zhi if zhi is not None else box_hi[2],
        ])

        pos = np.asarray(data.particles.positions)
        selected = np.all((pos >= bounds_lo) & (pos <= bounds_hi), axis=1)
        region_volume = float(np.prod(bounds_hi - bounds_lo))

    mass_lookup = _resolve_type_mass_lookup(data)
    ptype = np.asarray(data.particles["Particle Type"].array)
    atom_masses = mass_lookup[ptype]

    mass_sum = float(atom_masses[selected].sum())
    n_atoms = int(selected.sum())
    density = (mass_sum / region_volume) * AMU_A3_TO_G_CM3 if region_volume > 0 else float("nan")

    return {
        "density_g_cm3": density,
        "n_atoms": n_atoms,
        "volume_ang3": region_volume,
        "mass_sum_amu": mass_sum,
    }


# ---------------------------------------------------------------------------
# Argument parsing / entry point
# ---------------------------------------------------------------------------

def _parse_args(argv: Optional[list] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute region density from a LAMMPS dump (OVITO).")
    p.add_argument("--local-path", type=Path, required=True, help="Local LAMMPS dump file.")
    p.add_argument("--frame", type=int, default=None,
                    help="Frame index to read. If omitted and the file has "
                         "multiple frames, you will be prompted for one.")

    p.add_argument("--xlo", type=float, default=None)
    p.add_argument("--xhi", type=float, default=None)
    p.add_argument("--ylo", type=float, default=None)
    p.add_argument("--yhi", type=float, default=None)
    p.add_argument("--zlo", type=float, default=None)
    p.add_argument("--zhi", type=float, default=None)

    p.add_argument("--expression", type=str, default=None,
                    help="OVITO selection expression; overrides box bounds.")
    p.add_argument("--volume", type=float, default=None,
                    help="Region volume in Angstrom^3 (required with --expression).")
    return p.parse_args(argv)


def _resolve_frame(pipeline, requested_frame: Optional[int]) -> int:
    if requested_frame is not None:
        return requested_frame

    num_frames = pipeline.source.num_frames
    if num_frames <= 1:
        return 0

    while True:
        raw = input(f"File has {num_frames} frames. Enter frame index to compute density on (0-{num_frames - 1}): ")
        try:
            frame = int(raw)
        except ValueError:
            print("Please enter an integer.")
            continue
        if not (0 <= frame < num_frames):
            print(f"Frame index must be between 0 and {num_frames - 1}.")
            continue
        return frame


def main(argv: Optional[list] = None) -> None:
    args = _parse_args(argv)

    local_path = args.local_path.expanduser().resolve()
    pipeline = import_file(str(local_path))
    frame = _resolve_frame(pipeline, args.frame)

    result = compute_region_density(
        pipeline,
        frame=frame,
        xlo=args.xlo, xhi=args.xhi,
        ylo=args.ylo, yhi=args.yhi,
        zlo=args.zlo, zhi=args.zhi,
        expression=args.expression,
        volume=args.volume,
    )

    print(f"density_g_cm3={result['density_g_cm3']:.6f}")
    print(f"n_atoms={result['n_atoms']}")
    print(f"volume_ang3={result['volume_ang3']:.6f}")
    print(f"mass_sum_amu={result['mass_sum_amu']:.6f}")


if __name__ == "__main__":
    main()
