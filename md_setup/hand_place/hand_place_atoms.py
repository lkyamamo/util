#!/usr/bin/env python3
"""
hand_place_atoms.py

Adds individual atoms to a LAMMPS data file at positions defined relative to
atoms that are already there — e.g. hanging an H off a specific surface O to
make a silanol. Placements are hard-coded in the PLACEMENTS list below; each
entry is (target_atom_id, new_atom_type, (dx, dy, dz)) and creates one atom at
position(target) + offset.

The data-file reading, wrapping, distance, and writing helpers live in
../lammps_data.py, shared with the other md_setup scripts that add atoms.

Algorithm
---------
0.  Detect the atom style from the file's "Atoms # <style>" comment (pymatgen
    needs to be told the style; it does not read the comment itself).
1.  Read the data file with pymatgen's LammpsData (Atoms/Masses/Velocities as
    pandas tables indexed by atom id / type id).
2.  Validate the whole placement list up front: targets resolve, new types have
    masses, no topology sections, no duplicate ids.
3.  Apply placements in list order, one at a time. Each new atom is inserted
    immediately, so a later placement may target an atom an earlier one created.
    New ids continue from the file's highest id; existing ids are never touched.
    A placement that lands outside the box is wrapped back in and warned about,
    with the image indices it moved through folded into its image flags.
4.  Append Masses rows for any new types and zero Velocities rows for the new
    atoms, then write the output file.

Assumptions
-----------
- Intended for hand placement on small systems; the whole file is held in
  memory and distance checks are O(N x placements).
- The atom style comes from the "Atoms # <style>" comment. Files without it
  must set ATOM_STYLE, since guessing from the column count cannot tell
  molecular / bond / angle / charge apart — they are all six columns wide.
- New atoms inherit the target atom's molecule ID and image flags; charge (for
  atom_style charge/full) defaults to 0.0. Edit the written file if you need
  something else.
- Wrapping is per-axis and assumes an orthogonal box; on a triclinic cell the
  warning says so, and the image flags should be checked by hand. Files with no
  image flag columns are still wrapped, but the indices cannot be recorded.
- New atom types must extend the existing range contiguously, because the
  "N atom types" header is written as the number of Masses rows.
- Files with Bonds/Angles/Dihedrals are rejected: adding an atom would leave
  those sections stale and nothing here fixes them.
- pymatgen reformats the whole file — its own column spacing and fixed-decimal
  coordinates — and drops non-standard header lines such as
  "extra bond per atom". The output is a valid read_data input, not a minimal
  diff of the input.

Usage:
    edit the Configuration block below, then
    python3 hand_place_atoms.py
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import numpy as np
from pymatgen.io.lammps.data import LammpsData

# The shared data-file helpers live at the md_setup root, one level up.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import lammps_data as ld  # noqa: E402

# ── Configuration ─────────────────────────────────────────────────────────────
INPUT_FILE  = "system.data"
OUTPUT_FILE = "system_placed.data"

# Atom style. None reads it from the file's "Atoms # <style>" comment, which is
# the normal case. Set it only for files that lack the comment, or to override
# one that is wrong; a set value that disagrees with the comment is an error.
# One of: atomic, charge, molecular, full, bond, angle.
ATOM_STYLE = None

# (target_atom_id, new_atom_type, (dx, dy, dz)) — applied in order.
# A target id may be one created by an earlier entry.
PLACEMENTS = [
    (1, 2, (0.0, 0.0, 0.9572)),
]

# Masses in amu for any type id beyond the file's existing "N atom types".
NEW_TYPE_MASSES: dict[int, float] = {}

# Warn if a placed atom lands closer than this (Å) to an existing atom.
# None disables the check.
MIN_SEPARATION = 0.5
# ──────────────────────────────────────────────────────────────────────────────

_COORDS = ld.COORDS


# ── validation ────────────────────────────────────────────────────────────────

def _validate(data: LammpsData, placements: list, new_type_masses: dict[int, float]) -> None:
    """
    Check everything that can be checked before any atom is placed, and raise on
    the first problem found. Positions are not known yet, so the box-bounds and
    separation checks happen during placement instead.
    """
    if data.topology:
        sections = ", ".join(sorted(data.topology))
        raise ValueError(
            f"{INPUT_FILE} has topology sections ({sections}); adding atoms would "
            f"leave them inconsistent. Remove them or extend this script."
        )

    if data.atoms.index.has_duplicates:
        dupes = sorted(set(data.atoms.index[data.atoms.index.duplicated()]))
        raise ValueError(f"{INPUT_FILE} has duplicate atom ids: {dupes}")

    known_ids   = set(int(i) for i in data.atoms.index)
    known_types = set(int(t) for t in data.masses.index)
    next_new_id = int(data.atoms.index.max()) + 1
    next_type   = max(known_types) + 1

    for position, (target_id, new_type, offset) in enumerate(placements, start=1):
        where = f"placement {position} {(target_id, new_type, offset)}"

        if len(offset) != 3:
            raise ValueError(f"{where}: offset must have three components")

        if target_id not in known_ids:
            raise ValueError(
                f"{where}: target atom {target_id} does not exist in {INPUT_FILE} "
                f"and is not created by an earlier placement"
            )
        known_ids.add(next_new_id)
        next_new_id += 1

        # Only a type seen for the first time extends the range; later placements
        # may reuse a type an earlier one introduced.
        if new_type not in known_types:
            if new_type != next_type:
                raise ValueError(
                    f"{where}: type {new_type} would leave a gap in the type range "
                    f"(next available is {next_type}); atom types must be contiguous"
                )
            if new_type not in new_type_masses:
                raise ValueError(
                    f"{where}: type {new_type} is new; add its mass to "
                    f"NEW_TYPE_MASSES, e.g. NEW_TYPE_MASSES = {{{new_type}: 1.008}}"
                )
            known_types.add(new_type)
            next_type += 1


# ── placement ─────────────────────────────────────────────────────────────────

def place_atoms(data: LammpsData, placements: list,
                new_type_masses: dict[int, float]) -> list[tuple[int, int, np.ndarray]]:
    """
    Apply the placements in order, mutating `data` in place. Returns the
    (atom_id, type, position) of each atom added.
    """
    columns = list(data.atoms.columns)
    next_id = int(data.atoms.index.max()) + 1
    placed: list[tuple[int, int, np.ndarray]] = []

    for position_index, (target_id, new_type, offset) in enumerate(placements, start=1):
        where  = f"placement {position_index} {(target_id, new_type, offset)}"
        target = data.atoms.loc[target_id]
        raw_position = target[_COORDS].to_numpy(dtype=float) + np.asarray(offset, dtype=float)

        new_position, image_shift = ld.wrap_into_box(data, raw_position)
        if image_shift.any():
            warnings.warn(
                ld.format_wrap_warning(where, next_id, raw_position, new_position,
                                       image_shift, data, columns),
                stacklevel=2,
            )

        if MIN_SEPARATION is not None:
            neighbor_id, _, distance = ld.nearest_existing(data, new_position)[0]
            if distance < MIN_SEPARATION:
                warnings.warn(
                    f"{where}: new atom {next_id} is {distance:.4f} A from atom "
                    f"{neighbor_id} (below MIN_SEPARATION = {MIN_SEPARATION})",
                    stacklevel=2,
                )

        if new_type not in data.masses.index:
            ld.ensure_mass(data, new_type, new_type_masses[new_type])

        ld.add_atom(data, new_type, new_position, image_shift=image_shift,
                    template_row=target, atom_id=next_id)

        placed.append((next_id, int(new_type), new_position))
        next_id += 1

    return placed


# ── output ────────────────────────────────────────────────────────────────────

def write_data_file(data: LammpsData, output_file: str, input_file: str,
                    atom_style: str) -> None:
    """Write the data file with a first line that records where it came from."""
    ld.write_data_file(
        data, output_file,
        f"LAMMPS data file via hand_place_atoms.py, from {input_file}",
        atom_style,
    )


def _report(data: LammpsData, placed: list[tuple[int, int, np.ndarray]],
            atoms_before: int, types_before: int) -> None:
    print(f"placed {len(placed)} atom(s):")
    for atom_id, atom_type, position in placed:
        x, y, z = position
        print(f"  id {atom_id:>8}  type {atom_type}  at {x:.4f} {y:.4f} {z:.4f}")
    print(f"atoms:      {atoms_before} -> {len(data.atoms)}")
    print(f"atom types: {types_before} -> {len(data.masses)}")
    if len(data.masses) > types_before:
        for atom_type in sorted(data.masses.index)[types_before:]:
            print(f"  new type {atom_type}: mass {float(data.masses.loc[atom_type, 'mass'])} amu")
    print(f"velocities: {'padded with zeros' if data.velocities is not None else 'no section'}")


# ── entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    if not PLACEMENTS:
        print("PLACEMENTS is empty; nothing to do.")
        return

    data, atom_style = ld.load(INPUT_FILE, ATOM_STYLE, setting_name="ATOM_STYLE")
    atoms_before = len(data.atoms)
    types_before = len(data.masses)

    _validate(data, PLACEMENTS, NEW_TYPE_MASSES)
    placed = place_atoms(data, PLACEMENTS, NEW_TYPE_MASSES)

    write_data_file(data, OUTPUT_FILE, INPUT_FILE, atom_style)
    _report(data, placed, atoms_before, types_before)
    print(f"wrote {OUTPUT_FILE}")


if __name__ == "__main__":
    try:
        main()
    except (ValueError, KeyError, OSError) as error:
        sys.exit(f"hand_place_atoms.py: {error}")
