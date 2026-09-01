#!/usr/bin/env python3
"""
lammps_data.py

Shared helpers for reading, editing, and writing LAMMPS data files on top of
pymatgen's LammpsData. Extracted from hand_place_atoms.py so that every script
in md_setup that adds atoms to an existing structure agrees on the fiddly
parts: detecting the atom style, keeping the Atoms columns consistent, wrapping
new positions under PBC while recording image flags, measuring distances under
the minimum-image convention, and restoring the two things pymatgen's writer
drops (an informative first line and the "Atoms # <style>" comment).

Nothing here mutates a LammpsData except `ensure_mass` and `add_atom`; the
policy decisions (where an atom goes, what counts as a clash, what to do about
one) stay in the calling script.

This module lives at the md_setup root so scripts in the per-tool
subdirectories can import it:

    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from lammps_data import ...

Assumptions inherited from pymatgen and shared by every caller:
- The atom style comes from the "Atoms # <style>" comment. Files without it
  must be told the style, since the column count cannot tell molecular / bond /
  angle / charge apart - they are all six columns wide.
- Wrapping and minimum-image distances are per-axis and assume an orthogonal
  box. On a triclinic cell they ignore tilt; `is_triclinic` is provided so
  callers can warn.
- pymatgen reformats the whole file. Output is a valid read_data input, not a
  minimal diff of the input.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from pymatgen.io.lammps.data import ATOMS_HEADERS, LammpsData

COORDS = ["x", "y", "z"]
IMAGE_FLAG_COLUMNS = ("nx", "ny", "nz")
INTEGER_COLUMNS = ("type", "molecule-ID", "nx", "ny", "nz")

# Masses (amu) for the handful of elements md_setup deals with, for scripts that
# want to sanity-check a type against the element it is supposed to be. Compare
# with a tolerance; the exact isotope-averaged value in a file varies.
ELEMENT_MASSES = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
    "Na": 22.990,
    "Al": 26.982,
    "Si": 28.086,
    "P": 30.974,
    "S": 32.06,
    "Ca": 40.078,
}


# -- atom style ---------------------------------------------------------------

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


def resolve_atom_style(input_file: str, configured: str | None,
                       setting_name: str = "the atom style setting",
                       quiet: bool = False) -> str:
    """
    Reconcile the style declared in the file with the one the caller was given.
    `setting_name` names the caller's knob (a config constant, a CLI flag) so
    the error messages point at something the user can actually edit.
    """
    detected = detect_atom_style(input_file)

    if configured is None:
        if detected is None:
            raise ValueError(
                f"{input_file}: the Atoms section has no '# <style>' comment, so "
                f"the atom style cannot be detected; set {setting_name}"
            )
        if not quiet:
            print(f"atom style: {detected} (from the Atoms section comment)")
        return detected

    if detected is not None and detected != configured:
        raise ValueError(
            f"{input_file}: {setting_name} is '{configured}' but the Atoms "
            f"section says '{detected}'; fix whichever is wrong (leave it unset "
            f"to trust the file)"
        )
    if not quiet:
        print(f"atom style: {configured} (from {setting_name})")
    return configured


def load(input_file: str, atom_style: str | None = None,
         setting_name: str = "the atom style setting",
         quiet: bool = False) -> tuple[LammpsData, str]:
    """Resolve the atom style and read the file. Returns (data, atom_style)."""
    style = resolve_atom_style(input_file, atom_style, setting_name, quiet)
    return LammpsData.from_file(input_file, atom_style=style), style


# -- geometry -----------------------------------------------------------------

def box_lengths(data: LammpsData) -> np.ndarray:
    bounds = np.asarray(data.box.bounds, dtype=float)
    return bounds[:, 1] - bounds[:, 0]


def is_triclinic(data: LammpsData) -> bool:
    return data.box.tilt is not None and any(data.box.tilt)


def minimum_image(delta: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """Wrap a displacement into [-L/2, L/2) per axis. Ignores tilt."""
    return delta - lengths * np.round(delta / lengths)


def displacement(data: LammpsData, start, end) -> np.ndarray:
    """The shortest vector from `start` to `end` under PBC. Ignores tilt."""
    delta = np.asarray(end, dtype=float) - np.asarray(start, dtype=float)
    return minimum_image(delta, box_lengths(data))


def distance(data: LammpsData, a, b) -> float:
    """Minimum-image distance between two positions. Ignores tilt."""
    return float(np.linalg.norm(displacement(data, a, b)))


def wrap_into_box(data: LammpsData, position: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """
    Wrap a position into the box under PBC. Returns the wrapped position and the
    image indices it moved through, such that wrapped + image x L = the original
    position - the same convention as LAMMPS' image flags.
    """
    bounds = np.asarray(data.box.bounds, dtype=float)
    lengths = bounds[:, 1] - bounds[:, 0]
    image_shift = np.floor((position - bounds[:, 0]) / lengths).astype(int)
    return position - image_shift * lengths, image_shift


def position_of(data: LammpsData, atom_id: int) -> np.ndarray:
    """Coordinates of one atom, by id."""
    if atom_id not in data.atoms.index:
        raise ValueError(f"atom id {atom_id} is not in the data file")
    return data.atoms.loc[atom_id, COORDS].to_numpy(dtype=float)


def nearest_existing(data: LammpsData, position: np.ndarray, count: int = 1,
                     exclude: set[int] | None = None) -> list[tuple[int, int, float]]:
    """
    The `count` nearest existing atoms to a position under the minimum-image
    convention, as a list of (atom_id, atom_type, distance), closest first.
    `exclude` drops atom ids from the search - for an atom being deliberately
    walked toward a target, the target is not a clash.
    """
    atoms = data.atoms
    if exclude:
        atoms = atoms.drop(index=[i for i in exclude if i in atoms.index])
    if atoms.empty:
        return []
    lengths = box_lengths(data)
    deltas = atoms[COORDS].to_numpy(dtype=float) - np.asarray(position, dtype=float)
    deltas -= lengths * np.round(deltas / lengths)  # minimum image; ignores tilt
    distances = np.linalg.norm(deltas, axis=1)
    order = np.argsort(distances)[:count]
    return [
        (int(atoms.index[i]), int(atoms.iloc[i]["type"]), float(distances[i]))
        for i in order
    ]


def format_wrap_warning(where: str, atom_id: int, raw: np.ndarray, wrapped: np.ndarray,
                        image_shift: np.ndarray, data: LammpsData,
                        columns: list[str]) -> str:
    """Message for a new atom that landed outside the box and had to be wrapped."""
    moved = ", ".join(
        f"{axis} {raw[i]:.4f} -> {wrapped[i]:.4f}"
        for i, axis in enumerate(COORDS) if image_shift[i]
    )
    images = " ".join(str(int(i)) for i in image_shift)
    message = (
        f"{where}: new atom {atom_id} landed outside the box and was wrapped "
        f"({moved}); image indices {images}"
    )
    if not any(flag in columns for flag in IMAGE_FLAG_COLUMNS):
        message += (
            " - the file has no image flag columns, so the indices cannot be "
            "recorded and the unwrapped position is lost"
        )
    if is_triclinic(data):
        message += (
            " - the box is triclinic and this wrap ignores tilt, so check the "
            "result by hand"
        )
    return message


# -- types and masses ---------------------------------------------------------

def ensure_mass(data: LammpsData, atom_type: int, mass: float) -> bool:
    """
    Add a Masses row for `atom_type` if it has none. Returns True if a row was
    added. New types must extend the range contiguously, because the
    "N atom types" header is written as the number of Masses rows.
    """
    if atom_type in data.masses.index:
        return False
    next_type = int(max(data.masses.index)) + 1
    if atom_type != next_type:
        raise ValueError(
            f"atom type {atom_type} would leave a gap in the type range "
            f"(next available is {next_type}); atom types must be contiguous"
        )
    data.masses.loc[atom_type] = {"mass": float(mass)}
    data.masses.sort_index(inplace=True)
    return True


# -- adding atoms -------------------------------------------------------------

def next_atom_id(data: LammpsData) -> int:
    return int(data.atoms.index.max()) + 1


def next_molecule_id(data: LammpsData) -> int:
    if "molecule-ID" not in data.atoms.columns:
        return 0
    return int(data.atoms["molecule-ID"].max()) + 1


def new_atom_row(data: LammpsData, atom_type: int, position: np.ndarray,
                 image_shift: np.ndarray | None = None,
                 template_row: pd.Series | None = None,
                 molecule_id: int | None = None,
                 charge: float = 0.0) -> dict:
    """
    Build an Atoms row for a new atom. Columns the style needs but the caller
    did not specify are filled from `template_row` when one is given (this is
    how hand-placed atoms inherit a molecule ID and image flags from the atom
    they hang off), otherwise from `molecule_id` / `charge` / `image_shift`.

    Raises if the row would not cover every column of the Atoms table, rather
    than letting pandas insert NaN and write an unreadable data file.
    """
    columns = list(data.atoms.columns)
    if image_shift is None:
        image_shift = np.zeros(3, dtype=int)

    row: dict = {"type": int(atom_type)}
    row.update(dict(zip(COORDS, (float(c) for c in position))))

    if "molecule-ID" in columns:
        if template_row is not None:
            row["molecule-ID"] = int(template_row["molecule-ID"])
        elif molecule_id is not None:
            row["molecule-ID"] = int(molecule_id)
        else:
            raise ValueError(
                "the Atoms section has a molecule-ID column; pass molecule_id "
                "or a template_row"
            )
    if "q" in columns:
        row["q"] = float(template_row["q"]) if template_row is not None else float(charge)
    for axis, flag in enumerate(IMAGE_FLAG_COLUMNS):
        if flag in columns:
            base = int(template_row[flag]) if template_row is not None else 0
            row[flag] = base + int(image_shift[axis])

    missing = [c for c in columns if c not in row]
    if missing:
        raise ValueError(
            f"cannot build an Atoms row: no value for column(s) {missing}; "
            f"the atom style may need support adding here"
        )
    return row


def add_atom(data: LammpsData, atom_type: int, position: np.ndarray,
             image_shift: np.ndarray | None = None,
             template_row: pd.Series | None = None,
             molecule_id: int | None = None,
             charge: float = 0.0,
             atom_id: int | None = None) -> int:
    """
    Append one atom and pad the Velocities section with zeros if there is one.
    Returns the new atom id. The position is used as given - wrap it with
    `wrap_into_box` first if it might sit outside the box.
    """
    new_id = next_atom_id(data) if atom_id is None else int(atom_id)
    if new_id in data.atoms.index:
        raise ValueError(f"atom id {new_id} already exists")
    data.atoms.loc[new_id] = new_atom_row(
        data, atom_type, position, image_shift, template_row, molecule_id, charge
    )
    if data.velocities is not None:
        data.velocities.loc[new_id] = {"vx": 0.0, "vy": 0.0, "vz": 0.0}
    return new_id


# -- output -------------------------------------------------------------------

def restore_int_columns(data: LammpsData) -> None:
    """
    Keep integer columns integral. pymatgen formats only x/y/z/q/v* explicitly
    and prints everything else with pandas' default repr, so a `type` column
    that pandas upcast to float would be written as "1.0" and break read_data.
    """
    integral = [c for c in INTEGER_COLUMNS if c in data.atoms.columns]
    data.atoms[integral] = data.atoms[integral].astype(int)


def write_data_file(data: LammpsData, output_file: str, header_comment: str,
                    atom_style: str) -> None:
    """
    Write the data file, then restore the two things pymatgen's writer drops:
    an informative first line, and the "# <style>" comment on the Atoms section
    that LAMMPS, ASE, and box_size_from_data.py use to detect the atom style.
    """
    restore_int_columns(data)
    data.write_file(output_file, distance=10, charge=8)  # wider than the 6/4 defaults

    path = Path(output_file)
    lines = path.read_text().splitlines(keepends=True)
    lines[0] = f"{header_comment}\n"
    for index, line in enumerate(lines):
        if line.strip() == "Atoms":
            lines[index] = f"Atoms  # {atom_style}\n"
            break
    path.write_text("".join(lines))
