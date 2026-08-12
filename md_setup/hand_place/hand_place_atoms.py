#!/usr/bin/env python3
"""
hand_place_atoms.py

Adds individual atoms to a LAMMPS data file at positions defined relative to
atoms that are already there — e.g. hanging an H off a specific surface O to
make a silanol. Placements are hard-coded in the PLACEMENTS list below; each
entry is (target_atom_id, new_atom_type, (dx, dy, dz)) and creates one atom at
position(target) + offset.

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
import pandas as pd
from pymatgen.io.lammps.data import ATOMS_HEADERS, LammpsData

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

_COORDS = ["x", "y", "z"]


# ── atom style ────────────────────────────────────────────────────────────────

def detect_atom_style(input_file: str) -> str | None:
    """
    Read the atom style from the "Atoms # <style>" comment that LAMMPS writes.
    Returns None if the Atoms section carries no comment; raises if the comment
    names a style pymatgen cannot parse.
    """
    with open(input_file) as handle:
        for line in handle:
            keyword, _, comment = line.partition("#")
            if keyword.strip() != "Atoms":
                continue
            words = comment.split()
            if not words:
                return None
            style = words[0]
            if style not in ATOMS_HEADERS:
                supported = ", ".join(sorted(ATOMS_HEADERS))
                raise ValueError(
                    f"{input_file}: Atoms section declares atom style '{style}', "
                    f"which is not supported (supported: {supported})"
                )
            return style
    raise ValueError(f"{input_file}: no Atoms section found")


def resolve_atom_style(input_file: str, configured: str | None) -> str:
    """Reconcile the style in the file with the one set in the config block."""
    detected = detect_atom_style(input_file)

    if configured is None:
        if detected is None:
            raise ValueError(
                f"{input_file}: the Atoms section has no '# <style>' comment, so "
                f"the atom style cannot be detected; set ATOM_STYLE in the "
                f"configuration block"
            )
        print(f"atom style: {detected} (from the Atoms section comment)")
        return detected

    if detected is not None and detected != configured:
        raise ValueError(
            f"{input_file}: ATOM_STYLE is '{configured}' but the Atoms section "
            f"says '{detected}'; fix whichever is wrong (leave ATOM_STYLE = None "
            f"to trust the file)"
        )
    print(f"atom style: {configured} (from ATOM_STYLE)")
    return configured


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

def _new_atom_row(target: pd.Series, new_type: int, new_position: np.ndarray,
                  image_shift: np.ndarray, columns: list[str]) -> dict:
    """
    Build an Atoms row for a placed atom: the given type and position, with
    molecule ID inherited from the target, charge zeroed, and image flags taken
    from the target plus whatever wrapping the placement needed.
    """
    row: dict = {"type": int(new_type)}
    row.update(dict(zip(_COORDS, (float(c) for c in new_position))))
    if "molecule-ID" in columns:
        row["molecule-ID"] = int(target["molecule-ID"])
    if "q" in columns:
        row["q"] = 0.0
    for axis, flag in enumerate(("nx", "ny", "nz")):
        if flag in columns:
            row[flag] = int(target[flag]) + int(image_shift[axis])
    return row


def _wrap_into_box(data: LammpsData, position: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Wrap a position into the box under PBC. Returns the wrapped position and the
    image indices it moved through, such that wrapped + image x L = the original
    position — the same convention as LAMMPS' image flags.
    """
    bounds  = np.asarray(data.box.bounds, dtype=float)
    lengths = bounds[:, 1] - bounds[:, 0]
    image_shift = np.floor((position - bounds[:, 0]) / lengths).astype(int)
    return position - image_shift * lengths, image_shift


def _warn_wrapped(where: str, atom_id: int, raw: np.ndarray, wrapped: np.ndarray,
                  image_shift: np.ndarray, data: LammpsData, columns: list[str]) -> None:
    """Report a placement that landed outside the box and had to be wrapped."""
    moved = ", ".join(
        f"{axis} {raw[i]:.4f} -> {wrapped[i]:.4f}"
        for i, axis in enumerate(_COORDS) if image_shift[i]
    )
    images = " ".join(str(int(i)) for i in image_shift)
    message = (
        f"{where}: new atom {atom_id} landed outside the box and was wrapped "
        f"({moved}); image indices {images}"
    )
    if not any(flag in columns for flag in ("nx", "ny", "nz")):
        message += (
            " — the file has no image flag columns, so the indices cannot be "
            "recorded and the unwrapped position is lost"
        )
    if data.box.tilt is not None and any(data.box.tilt):
        message += (
            " — the box is triclinic and this wrap ignores tilt, so check the "
            "result by hand"
        )
    warnings.warn(message, stacklevel=3)


def _closest_existing(data: LammpsData, position: np.ndarray) -> tuple[int, float]:
    """Nearest existing atom to a position under the minimum-image convention."""
    bounds  = np.asarray(data.box.bounds, dtype=float)
    lengths = bounds[:, 1] - bounds[:, 0]
    deltas  = data.atoms[_COORDS].to_numpy(dtype=float) - position
    deltas -= lengths * np.round(deltas / lengths)  # minimum image; ignores tilt
    distances = np.linalg.norm(deltas, axis=1)
    nearest   = int(np.argmin(distances))
    return int(data.atoms.index[nearest]), float(distances[nearest])


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

        new_position, image_shift = _wrap_into_box(data, raw_position)
        if image_shift.any():
            _warn_wrapped(where, next_id, raw_position, new_position,
                          image_shift, data, columns)

        if MIN_SEPARATION is not None:
            neighbor_id, distance = _closest_existing(data, new_position)
            if distance < MIN_SEPARATION:
                warnings.warn(
                    f"{where}: new atom {next_id} is {distance:.4f} A from atom "
                    f"{neighbor_id} (below MIN_SEPARATION = {MIN_SEPARATION})",
                    stacklevel=2,
                )

        if new_type not in data.masses.index:
            data.masses.loc[new_type] = {"mass": float(new_type_masses[new_type])}
            data.masses.sort_index(inplace=True)

        data.atoms.loc[next_id] = _new_atom_row(target, new_type, new_position,
                                                image_shift, columns)
        if data.velocities is not None:
            data.velocities.loc[next_id] = {"vx": 0.0, "vy": 0.0, "vz": 0.0}

        placed.append((next_id, int(new_type), new_position))
        next_id += 1

    return placed


# ── output ────────────────────────────────────────────────────────────────────

def _restore_int_columns(data: LammpsData) -> None:
    """
    Keep integer columns integral. pymatgen formats only x/y/z/q/v* explicitly
    and prints everything else with pandas' default repr, so a `type` column
    that pandas upcast to float would be written as "1.0" and break read_data.
    """
    integral = [c for c in ("type", "molecule-ID", "nx", "ny", "nz") if c in data.atoms.columns]
    data.atoms[integral] = data.atoms[integral].astype(int)


def write_data_file(data: LammpsData, output_file: str, input_file: str,
                    atom_style: str) -> None:
    """
    Write the data file, then restore the two things pymatgen's writer drops:
    an informative first line, and the "# <style>" comment on the Atoms section
    that LAMMPS, ASE, and box_size_from_data.py use to detect the atom style.
    """
    _restore_int_columns(data)
    data.write_file(output_file, distance=10, charge=8)  # wider than the 6/4 defaults

    path = Path(output_file)
    lines = path.read_text().splitlines(keepends=True)
    lines[0] = f"LAMMPS data file via hand_place_atoms.py, from {input_file}\n"
    for index, line in enumerate(lines):
        if line.strip() == "Atoms":
            lines[index] = f"Atoms  # {atom_style}\n"
            break
    path.write_text("".join(lines))


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

    atom_style = resolve_atom_style(INPUT_FILE, ATOM_STYLE)
    data = LammpsData.from_file(INPUT_FILE, atom_style=atom_style)
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
