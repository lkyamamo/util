#!/usr/bin/env python3
"""
water_placement.py

Putting a water molecule into a LAMMPS data file: resolving its atom types,
copying a structure and adding the three atoms, and reporting what each of them
nearly ran into.

This is the format-dependent counterpart to water_path/water_geometry.py, which
builds and aims the molecule but knows nothing about files. It is built on
lammps_data, so it is LAMMPS-only by construction.

It used to live in generate_water_path.py. That script has since moved onto
structure_io, an ASE-backed layer that writes POSCARs as well as data files, and
resolves types through structure_io.add_water instead. The ring tools stayed on
lammps_data, so these kept versions live here rather than being deleted.
"""

from __future__ import annotations

import copy
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "water_path"))
import lammps_data as ld  # noqa: E402
from water_geometry import WATER_LABELS  # noqa: E402

# ── Configuration ─────────────────────────────────────────────────────────────
# The atom types the water is written as. Every frame adds atoms of both: one
# oxygen and two hydrogens.
#
# These are not guessed from the input file. The type ids have to match the
# ordering of the element list on your pair_coeff line, which lives in the
# LAMMPS input script rather than in the data file, so set them to whatever that
# ordering says and change them per system.
#
# A type the input file already defines is used as it stands, and its mass here
# is only cross-checked against the file's. A type the file does not define is
# added to Masses with the mass here, which requires it to extend the file's
# type range contiguously - so with a 2-type silica file (Si, O), a new hydrogen
# type must be 3.
WATER_TYPES = {
    "O": {"type": 2, "mass": 15.9994},
    "H": {"type": 3, "mass": 1.00784},
}

# Tolerance in amu when sanity-checking a mass against the element it should be.
MASS_TOLERANCE = 0.2
# ──────────────────────────────────────────────────────────────────────────────

# Two atom types have to be settled before any frame is written - not just the
# hydrogen type that a silica file is most likely to be missing.
WATER_ELEMENTS = ("O", "H")
ELEMENT_WORDS = {"O": "oxygen", "H": "hydrogen"}

# The nearest-neighbour columns every manifest carries. Callers prepend their own
# coordinate columns and append "clash" and "data_file".
NEAREST_COLUMNS = [
    "o_nearest_id", "o_nearest_type", "o_nearest_distance",
    "h1_nearest_id", "h1_nearest_type", "h1_nearest_distance",
    "h2_nearest_id", "h2_nearest_type", "h2_nearest_distance",
]



def template_offsets(template_file: str, atom_style: str | None) -> tuple[np.ndarray, np.ndarray]:
    """
    Read the two O->H offsets out of a data file holding a single water
    molecule, so the seeded water has whatever geometry that file was built or
    equilibrated with. The heaviest of the three atoms is taken as the oxygen.
    """
    water, _ = ld.load(template_file, atom_style, setting_name="--atom-style", quiet=True)
    if len(water.atoms) != 3:
        raise ValueError(
            f"{template_file}: expected a single water molecule (3 atoms), "
            f"found {len(water.atoms)}"
        )
    masses = {int(t): float(m) for t, m in water.masses["mass"].items()}
    by_mass = sorted(water.atoms.index, key=lambda i: masses[int(water.atoms.loc[i, "type"])])
    hydrogen_ids, oxygen_id = by_mass[:2], by_mass[2]

    oxygen = ld.position_of(water, oxygen_id)
    offsets = [ld.displacement(water, oxygen, ld.position_of(water, h)) for h in hydrogen_ids]
    lengths = [np.linalg.norm(o) for o in offsets]
    if not all(0.5 < length < 1.5 for length in lengths):
        raise ValueError(
            f"{template_file}: O-H distances are {lengths[0]:.3f} and {lengths[1]:.3f} A, "
            f"which do not look like a water molecule"
        )
    return offsets[0], offsets[1]


def apply_water_types(data) -> tuple[int, int, list[str]]:
    """
    Register the WATER_TYPES entries against the input file, and report what
    each element ended up as. Both types are settled before any frame is
    written, because every frame adds atoms of both.
    """
    types = {element: int(WATER_TYPES[element]["type"]) for element in WATER_ELEMENTS}
    masses = {element: float(WATER_TYPES[element]["mass"]) for element in WATER_ELEMENTS}
    notes: list[str] = []

    if types["O"] == types["H"]:
        raise ValueError(
            f"WATER_TYPES gives the oxygen and the hydrogens the same atom type "
            f"({types['O']}); they must be different types"
        )

    for element in WATER_ELEMENTS:
        if types[element] in data.masses.index:
            word = ELEMENT_WORDS[element]
            existing = float(data.masses.loc[types[element], "mass"])
            notes.append(
                f"water {word} uses existing atom type {types[element]} (mass {existing})")
            if abs(existing - masses[element]) > MASS_TOLERANCE:
                notes.append(
                    f"  WARNING: WATER_TYPES gives {word} a mass of {masses[element]}, but "
                    f"type {types[element]} in the input file has mass {existing} - check "
                    f"that this is the type you meant"
                )

    # Added in ascending type order so a file missing both types still ends up
    # with a contiguous range, whatever order WATER_TYPES lists them in.
    missing = [e for e in WATER_ELEMENTS if types[e] not in data.masses.index]
    for element in sorted(missing, key=lambda e: types[e]):
        word = ELEMENT_WORDS[element]
        next_type = int(max(data.masses.index)) + 1
        if types[element] != next_type:
            raise ValueError(
                f"WATER_TYPES puts {word} at atom type {types[element]}, which the input "
                f"file does not define and which would leave a gap in the type range "
                f"(next available is {next_type}); atom types must be contiguous"
            )
        ld.ensure_mass(data, types[element], masses[element])
        notes.append(
            f"added atom type {types[element]} for {word} (mass {masses[element]}) - "
            f"the potential's element list and pair_coeff line must be extended to match"
        )

    notes.append(
        f"each frame adds 3 atoms of 2 types: 1 oxygen of type {types['O']} and "
        f"2 hydrogens of type {types['H']}"
    )
    return types["O"], types["H"], notes


def build_frame(data, positions: dict[str, np.ndarray], oxygen_type: int,
                hydrogen_type: int, charges: dict[str, float], molecule_id: int,
                where: str) -> tuple[object, dict[str, int], int, list[str]]:
    """
    Copy the structure and add the three water atoms.
    Returns (frame, {label: atom id}, atoms wrapped, warnings).

    A path that crosses a periodic boundary wraps as a matter of course and the
    wrapped structure is physically identical, so wrapping is only counted, not
    warned about - except on a triclinic box, where the per-axis wrap ignores
    tilt and the result is worth checking by hand.
    """
    frame = copy.deepcopy(data)
    columns = list(frame.atoms.columns)
    types = {"O": oxygen_type, "H1": hydrogen_type, "H2": hydrogen_type}
    ids: dict[str, int] = {}
    messages: list[str] = []
    wrapped_count = 0

    for label in WATER_LABELS:
        wrapped, image_shift = ld.wrap_into_box(frame, positions[label])
        atom_id = ld.add_atom(frame, types[label], wrapped, image_shift=image_shift,
                              molecule_id=molecule_id, charge=charges[label])
        ids[label] = atom_id
        if image_shift.any():
            wrapped_count += 1
            if ld.is_triclinic(frame):
                messages.append(ld.format_wrap_warning(
                    f"{where} {label}", atom_id, positions[label], wrapped,
                    image_shift, frame, columns))

    return frame, ids, wrapped_count, messages


def nearest_report(data, positions: dict[str, np.ndarray],
                   exclude: dict[str, set[int]] | None = None
                   ) -> dict[str, tuple[int, int, float]]:
    """
    The nearest existing atom to each of the three water atoms separately, as
    {label: (atom id, atom type, distance)}.

    Separately, because a single overall minimum hides a hydrogen buried in a
    wall whenever the oxygen happens to be closer to something else. Measured
    against the original structure, so the water is not compared with itself.
    """
    exclude = exclude or {}
    return {label: ld.nearest_existing(data, positions[label],
                                       exclude=exclude.get(label))[0]
            for label in WATER_LABELS}


def nearest_columns(nearest: dict[str, tuple[int, int, float]]) -> dict:
    """Flatten a nearest_report into the NEAREST_COLUMNS manifest fields."""
    row = {}
    for label in WATER_LABELS:
        near_id, near_type, near_distance = nearest[label]
        prefix = label.lower()
        row[f"{prefix}_nearest_id"] = near_id
        row[f"{prefix}_nearest_type"] = near_type
        row[f"{prefix}_nearest_distance"] = round(near_distance, 4)
    return row
