#!/usr/bin/env python3
"""Report box dimensions, current density, and target-density dimensions for a
LAMMPS data file.

Same reporting as box_size.py, but the composition and box come from an actual
LAMMPS data file instead of a unit cell definition plus repeats: atom counts
per type, per-type masses (Masses section), and the box are all read from the
file with ASE.

Usage
-----
    python box_size_from_data.py system.data
    python box_size_from_data.py system.data --target-density 1.0
    python box_size_from_data.py system.data --target-density 1.0 2.5 --fix xy

Notes
-----
- The atom style is auto-detected from the "Atoms # style" comment; override
  with --atom-style if the comment is missing or wrong.
- Triclinic boxes are handled: volume comes from the full cell matrix, and
  tilt factors are reported but left untouched by the target-density solve.
- Element symbols are ASE's inference from the Masses section and are shown
  as labels only; every number below comes from the masses themselves.
"""

import argparse
import sys
from pathlib import Path

import numpy as np
from ase.io import read

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from box_size import (  # noqa: E402  (path setup must precede the import)
    AMU_TO_G,
    cubic_dimension_A,
    density_g_cm3,
    third_dimension_A,
)

AXES = ("x", "y", "z")


def read_lammps_data(path, atom_style=None):
    """Read a LAMMPS data file and pull out box, per-type masses, and counts.

    Returns a dict with: n_atoms, lengths (lx, ly, lz), tilt (xy, xz, yz),
    volume_A3, and composition as a list of (type_id, symbol, mass_amu, count)
    sorted by type id.
    """
    kwargs = {"atom_style": atom_style} if atom_style else {}
    atoms = read(str(path), format="lammps-data", index=0, **kwargs)
    if isinstance(atoms, list):
        atoms = atoms[0]

    # ASE stores a LAMMPS cell in lower-triangular form:
    #   row 0 = (lx, 0, 0), row 1 = (xy, ly, 0), row 2 = (xz, yz, lz)
    cell = atoms.cell.array
    lengths = tuple(float(cell[i, i]) for i in range(3))
    tilt = (float(cell[1, 0]), float(cell[2, 0]), float(cell[2, 1]))

    masses = atoms.get_masses()
    symbols = atoms.get_chemical_symbols()
    type_ids = atoms.arrays.get("type")
    if type_ids is None:  # no Atoms type column recovered; group by mass instead
        type_ids = np.unique(masses, return_inverse=True)[1] + 1

    composition = []
    for type_id in np.unique(type_ids):
        sel = type_ids == type_id
        type_masses = np.unique(masses[sel])
        if type_masses.size != 1:
            raise ValueError(
                f"{path}: atom type {type_id} has multiple masses {type_masses}"
            )
        mass = float(type_masses[0])
        if mass <= 0:
            raise ValueError(
                f"{path}: atom type {type_id} has no mass; the file needs a "
                f"Masses section"
            )
        symbol = symbols[int(np.flatnonzero(sel)[0])]
        composition.append((int(type_id), symbol, mass, int(sel.sum())))

    return {
        "n_atoms": len(atoms),
        "lengths": lengths,
        "tilt": tilt,
        "volume_A3": float(atoms.get_volume()),
        "composition": composition,
    }


def total_mass_amu(composition):
    return sum(mass * count for _, _, mass, count in composition)


def isotropic_scale(current_density, target_density):
    """Uniform scale factor on all three edges that yields target_density."""
    return (current_density / target_density) ** (1 / 3)


def report(info, target_densities, fixed_pairs):
    lengths = info["lengths"]
    volume = info["volume_A3"]
    composition = info["composition"]

    mass_amu = total_mass_amu(composition)
    mass_g = mass_amu * AMU_TO_G
    current_density = density_g_cm3(mass_g, volume)

    print("box:")
    for axis, length in zip(AXES, lengths):
        print(f"  l{axis} = {length:.4f} A")
    if any(info["tilt"]):
        xy, xz, yz = info["tilt"]
        print(f"  tilt: xy = {xy:.4f}, xz = {xz:.4f}, yz = {yz:.4f}")
    print(f"  volume: {volume:.4f} A^3")

    print("composition:")
    for type_id, symbol, mass, count in composition:
        frac = count * mass / mass_amu
        print(
            f"  type {type_id} ({symbol}): {count} atoms of mass {mass:g} amu "
            f"({100 * frac:.2f}% of mass)"
        )
    print(f"total atom count: {info['n_atoms']}")
    print(f"total mass: {mass_g:.6e} g ({mass_amu:.4f} amu)")
    print(f"current density: {current_density:.6f} g/cm^3")
    print(f"number density: {info['n_atoms'] / volume:.6f} atoms/A^3")

    results = {}
    for target in target_densities:
        scale = isotropic_scale(current_density, target)
        scaled = tuple(scale * length for length in lengths)
        cubic = cubic_dimension_A(mass_g, target)

        print(f"\ntarget density: {target:.6f} g/cm^3")
        print(
            f"  isotropic scale: {scale:.6f}  ->  "
            f"{scaled[0]:.4f} x {scaled[1]:.4f} x {scaled[2]:.4f} A"
        )
        print(f"  cubic supercell dimension: {cubic:.4f} A")
        print("  holding two edges fixed, solving the third:")
        solved = {}
        for pair in fixed_pairs:
            free = next(a for a in AXES if a not in pair)
            i, j, k = (AXES.index(a) for a in (*pair, free))
            third = third_dimension_A(mass_g, target, lengths[i], lengths[j])
            dims = list(lengths)
            dims[k] = third
            solved[free] = third
            print(
                f"    fix l{pair[0]}, l{pair[1]} -> l{free} = {third:.4f} A  "
                f"({dims[0]:.4f} x {dims[1]:.4f} x {dims[2]:.4f} A)"
            )
        results[target] = {
            "scale": scale,
            "scaled_lengths": scaled,
            "cubic": cubic,
            "solved": solved,
        }

    return {
        "mass_g": mass_g,
        "volume_A3": volume,
        "current_density_g_cm3": current_density,
        "targets": results,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Report box dimensions, current density, and target-density "
            "dimensions for a LAMMPS data file."
        )
    )
    parser.add_argument("data_file", help="LAMMPS data file")
    parser.add_argument(
        "--target-density",
        type=float,
        nargs="+",
        default=[],
        metavar="G_CM3",
        help="One or more target mass densities in g/cm^3.",
    )
    parser.add_argument(
        "--fix",
        action="append",
        choices=["xy", "xz", "yz"],
        metavar="{xy,xz,yz}",
        help=(
            "Axis pair to hold fixed while solving for the third edge. "
            "Repeatable; default is all three pairs."
        ),
    )
    parser.add_argument(
        "--atom-style",
        metavar="STYLE",
        help="Override the atom style auto-detected from the file.",
    )
    args = parser.parse_args(argv)

    path = Path(args.data_file)
    info = read_lammps_data(path, atom_style=args.atom_style)
    fixed_pairs = [tuple(p) for p in (args.fix or ["xy", "xz", "yz"])]

    print(f"file: {path}")
    return report(info, args.target_density, fixed_pairs)


if __name__ == "__main__":
    main()
