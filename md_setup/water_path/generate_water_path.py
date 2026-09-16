#!/usr/bin/env python3
"""
generate_water_path.py

Seed a series of structures for a stepwise hydrolysis study: a water molecule is
placed in the pore/interstice of a silica structure and its oxygen is walked
along a straight line toward a targeted silicon, one file per step. Each file is
an independent starting point for its own constrained minimization / MD run,
where the O-Si distance is restrained and everything else - the hydrogens
included - is free to relax.

LAMMPS data files and VASP POSCARs are both supported, and the format is set by
the required --format flag. It applies to the input and the output together: a
data file in gives data files out, a POSCAR in gives POSCARs out.

What this script is and is not
------------------------------
The frames are *seeds*, not a reaction path. The two O-H offsets are held rigid
along the line, which is only a starting guess: the real hydrolysis reaction is
the protons rearranging, so by the late frames the rigid hydrogens are pointing
somewhere unphysical. That is fine as long as every frame gets a thorough
relaxation with the hydrogens unconstrained, and it is *not* fine if you plan to
read energies straight off these geometries. If you want the hydrogens to follow
the reaction coordinate, the frames have to be chained - frame i+1 seeded from
the relaxed output of frame i - which is a different script (a driver that
shells out to LAMMPS between steps).

Along the same lines: a straight line to the silicon nucleus is a crude
coordinate. Water attacks silicon along the back side of an existing Si-O bond,
so prefer --backside-of <bridging O id>, which starts the water on the far side
of the silicon from that oxygen and walks it straight in. The per-frame nearest
-neighbour report is there to catch the other failure mode, a line that drives
the oxygen through an atom that happens to be in the way.

The roll of the water about that line is a free parameter - the approach
direction fixes where the oxygen points but not which way the hydrogens are
turned around it. --h-rotation sets that roll in degrees, and the nearest
-neighbour column is how you tell whether a given one puts a hydrogen into the
pore wall.

Output
------
Frames are numbered from 1 to the frame count, so --n-steps N writes frame 1 (the
water at --start-distance) through frame N+1 (the water at --final-distance).

  --format lammps   <i>.data, ready for read_data.
  --format poscar   <i>/POSCAR, one directory per frame. VASP reads a file named
                    literally POSCAR, so each frame is ready to run once an
                    INCAR, KPOINTS, and POTCAR are dropped in beside it.
  manifest.csv      one row per frame, at the top of --outdir: the O-Si distance
                    the run should restrain to, the ids of the three new atoms,
                    and the closest approach of each to an existing atom (the
                    target silicon is skipped for the oxygen, whose distance to
                    it is the O-Si column, but kept for the hydrogens).

Atom ids
--------
Every id in the manifest is an index into the frame that was written, not into
the input file. Under --format lammps the two agree, because atom ids survive a
write. Under --format poscar they do not: a POSCAR groups its sites by species
and states the counts in a header, so the water has to join the end of the O and
H runs to keep that grouping, and every atom after an insertion point shifts
down. --silicon-id and --backside-of are read as 1-based line numbers in the
*input* coordinate block; the manifest's si_id_out column gives the silicon's
index in the output, which is the one to name in an ICONST constraint.

Atom types and species
----------------------
Every frame adds one oxygen and two hydrogens.

Under --format lammps both come from the WATER_TYPES dictionary in the
configuration block below - set it to match your potential's element list before
running. The elements behind the input's types are known, from its Masses table,
so a type WATER_TYPES names is checked against what the file says that type
actually is. A type the input does not define is added to Masses, and needs a
matching entry on the pair_coeff line, which this script cannot write for you.

Under --format poscar there are no atom types: the water is written as species O
and H, and WATER_TYPES is ignored. A species the input does not already contain
becomes a new group on the species line, and the POTCAR has to be extended to
match - the same caveat as the pair_coeff line above.

What is not carried over
------------------------
Velocities and LAMMPS image flags are dropped. Every frame is an independent
starting point for its own relaxation, so a velocity or an image count inherited
from whatever run produced the input has no meaning in it. Positions are used
exactly as the input states them; only the three new atoms are wrapped into the
cell. Everything else - selective dynamics, charges, molecule ids, the cell, the
species order - is carried through unchanged.

Usage examples
--------------
  # back-side attack on the Si-O bond between Si 5 and bridging O 91
  python generate_water_path.py --format lammps --input silica.data \\
      --outdir frames --silicon-id 5 --backside-of 91 --start-distance 4.5 \\
      --final-distance 1.9 --n-steps 10

  # the same path in VASP, writing frames/1/POSCAR ... frames/11/POSCAR
  python generate_water_path.py --format poscar --input POSCAR \\
      --outdir frames --silicon-id 5 --backside-of 91 --start-distance 4.5 \\
      --final-distance 1.9 --n-steps 10

  # the same, set up for a CG minimization of the water against a rigid
  # framework: every atom from the input is frozen, only the water relaxes
  python generate_water_path.py --format poscar --input POSCAR \\
      --outdir frames --silicon-id 5 --backside-of 91 --freeze-substrate

  # explicit starting point, water geometry taken from an existing data file
  python generate_water_path.py --format lammps --input silica.data \\
      --outdir frames --silicon-id 5 --initial-position 7.0 7.0 7.0 \\
      --water-template ../isolated_water.data --dry-run
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

# The shared structure helpers live at the md_setup root, one level up. One
# module reads and writes both formats, on top of ASE, so everything below is
# format-blind apart from the few places that have to know (`water_types`, and
# the flag guards in parse_args).
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import structure_io as sio  # noqa: E402


# ── Configuration ─────────────────────────────────────────────────────────────
# The atom types the water is written as, under --format lammps only. Every frame
# adds atoms of both: one oxygen and two hydrogens. A POSCAR has no atom types,
# so under --format poscar this is ignored and the water is written as species O
# and H.
#
# These are not guessed from the input file. The type ids have to match the
# ordering of the element list on your pair_coeff line, which lives in the
# LAMMPS input script rather than in the data file, so set them to whatever that
# ordering says and change them per system.
#
# A type the input file already defines has to agree with what the file says that
# type is, and the run stops if it does not. A type the file does not define is
# added to Masses, which requires it to extend the file's type range contiguously
# - so with a 2-type silica file (Si, O), a new hydrogen type must be 3. Masses
# are not set here; they come from the element.
WATER_TYPES = {"O": 2, "H": 3}
# ──────────────────────────────────────────────────────────────────────────────


# ── vector helpers ────────────────────────────────────────────────────────────

def normalize(vector: np.ndarray) -> np.ndarray:
    vector = np.asarray(vector, dtype=float)
    length = np.linalg.norm(vector)
    if length < 1e-9:
        raise ValueError("cannot normalize a zero-length vector")
    return vector / length


def perpendicular_to(direction: np.ndarray, hint: np.ndarray | None = None) -> np.ndarray:
    """
    A unit vector perpendicular to `direction`, deterministically chosen so that
    two runs with the same inputs produce the same water orientation. `hint`
    picks which perpendicular, when the caller cares; a hint parallel to
    `direction` is ignored.
    """
    direction = normalize(direction)
    if hint is not None:
        residual = np.asarray(hint, dtype=float)
        residual = residual - np.dot(residual, direction) * direction
        if np.linalg.norm(residual) > 1e-6:
            return normalize(residual)
    axis = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(direction, axis)) > 0.9:
        axis = np.array([0.0, 1.0, 0.0])
    return normalize(np.cross(direction, axis))


# ── water geometry ────────────────────────────────────────────────────────────

def canonical_offsets(oh_length: float, hoh_angle_deg: float) -> tuple[np.ndarray, np.ndarray]:
    """Two O->H offsets in a local frame, bisector along +z and both H in the xz plane."""
    half = np.radians(hoh_angle_deg) / 2.0
    bisector = np.array([0.0, 0.0, 1.0])
    spread = np.array([1.0, 0.0, 0.0])
    h1 = oh_length * (np.cos(half) * bisector + np.sin(half) * spread)
    h2 = oh_length * (np.cos(half) * bisector - np.sin(half) * spread)
    return h1, h2


def template_offsets(fmt: str, template_file: str,
                     atom_style: str | None) -> tuple[np.ndarray, np.ndarray]:
    """
    Read the two O->H offsets out of a file holding a single water molecule, so
    the seeded water has whatever geometry that file was built or equilibrated
    with. The template is read in the same --format as everything else.
    """
    water = sio.load(template_file, fmt, atom_style, quiet=True)
    if sio.atom_count(water) != 3:
        raise ValueError(
            f"{template_file}: expected a single water molecule (3 atoms), "
            f"found {sio.atom_count(water)}"
        )
    oxygen_id, hydrogen_ids = sio.identify_water(water)

    oxygen = sio.position_of(water, oxygen_id)
    offsets = [sio.displacement(water, oxygen, sio.position_of(water, h))
               for h in hydrogen_ids]
    lengths = [np.linalg.norm(o) for o in offsets]
    if not all(0.5 < length < 1.5 for length in lengths):
        raise ValueError(
            f"{template_file}: O-H distances are {lengths[0]:.3f} and {lengths[1]:.3f} A, "
            f"which do not look like a water molecule"
        )
    return offsets[0], offsets[1]


def orient_offsets(h1: np.ndarray, h2: np.ndarray, bisector: np.ndarray,
                   in_plane_hint: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """
    Rigidly rotate a pair of O->H offsets so their bisector points along
    `bisector`, preserving both bond lengths and the H-O-H angle exactly. The
    molecular plane is fixed by `in_plane_hint`, whose component perpendicular
    to the bisector becomes the direction the two hydrogens spread along.
    """
    h1 = np.asarray(h1, dtype=float)
    h2 = np.asarray(h2, dtype=float)

    source_bisector = normalize(h1) + normalize(h2)
    if np.linalg.norm(source_bisector) < 1e-6:
        raise ValueError("the two O-H offsets are anti-parallel; there is no bisector")
    source_bisector = normalize(source_bisector)
    source_spread = normalize(h1 - np.dot(h1, source_bisector) * source_bisector)

    target_bisector = normalize(bisector)
    target_spread = perpendicular_to(target_bisector, in_plane_hint)

    return tuple(
        np.dot(offset, source_bisector) * target_bisector
        + np.dot(offset, source_spread) * target_spread
        for offset in (h1, h2)
    )


def rotate_about_bisector(h1: np.ndarray, h2: np.ndarray,
                          degrees: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Spin the two O->H offsets about their own H-O-H bisector by `degrees`.

    The bisector is the rotation axis, so it does not move: the oxygen keeps
    facing wherever the orientation aimed it, and what turns is the plane the
    two hydrogens lie in. Under --orientation away/toward the bisector is the
    O-Si axis, so
    this is the free parameter that nothing else fixes - which way the water is
    rolled around its line of approach.

    The rotation is counterclockwise about the bisector, seen looking down the
    bisector toward the oxygen (equivalently, the right-hand rule with the thumb
    along the bisector). It is always counterclockwise: the angle is reduced
    into [0, 360) first, so a negative value becomes the counterclockwise turn
    that lands in the same place rather than a clockwise one.

    Bond lengths and the H-O-H angle are untouched - this is a rigid rotation.
    """
    angle = np.radians(float(degrees) % 360.0)
    axis = normalize(normalize(h1) + normalize(h2))

    def rotated(vector):
        # Rodrigues' rotation formula about the unit `axis`.
        return (vector * np.cos(angle)
                + np.cross(axis, vector) * np.sin(angle)
                + axis * np.dot(axis, vector) * (1.0 - np.cos(angle)))

    return rotated(np.asarray(h1, dtype=float)), rotated(np.asarray(h2, dtype=float))


def resolve_offsets(fmt: str, args, atom_style: str | None,
                    approach: np.ndarray) -> tuple[np.ndarray, np.ndarray, str]:
    """
    Work out the two rigid O->H offsets for the whole path. `approach` is the
    unit vector pointing from the starting position toward the silicon.
    """
    if args.h1_offset is not None:
        h1, h2 = np.array(args.h1_offset), np.array(args.h2_offset)
        source = "given explicitly"
    else:
        if args.water_template is not None:
            h1, h2 = template_offsets(fmt, args.water_template, atom_style)
            source = f"from {args.water_template}"
        else:
            h1, h2 = canonical_offsets(args.oh_length, args.hoh_angle)
            source = (f"built from --oh-length {args.oh_length} and "
                      f"--hoh-angle {args.hoh_angle}")

        if args.orientation == "away":
            # Bisector anti-parallel to the approach: the hydrogens trail behind
            # the oxygen and its lone pairs face the silicon, which is the
            # geometry a nucleophilic attack starts from.
            h1, h2 = orient_offsets(h1, h2, -approach)
            source += ", hydrogens pointing away from the target"
        elif args.orientation == "toward":
            # Bisector along the approach: the hydrogens lead and reach the
            # silicon before the oxygen does. This is not the nucleophilic
            # attack geometry - it is the arrangement to use when the proton,
            # not the oxygen, is what should arrive first.
            h1, h2 = orient_offsets(h1, h2, approach)
            source += ", hydrogens pointing toward the target"
        else:
            # Molecular plane perpendicular to the approach: both hydrogens stay
            # off the line of travel, which keeps them clear of a tight channel.
            bisector = perpendicular_to(approach)
            h1, h2 = orient_offsets(h1, h2, bisector,
                                    in_plane_hint=np.cross(approach, bisector))
            source += ", molecular plane perpendicular to the path"

    # Applied last, and to every source of offsets, so --h-rotation always means
    # the same thing: a roll about the bisector the orientation just set. The
    # bisector is the axis, so this cannot undo the orientation above.
    rotation = float(args.h_rotation) % 360.0
    if rotation:
        h1, h2 = rotate_about_bisector(h1, h2, rotation)
        source += f", rolled {rotation:g} deg counterclockwise about the bisector"

    return h1, h2, source


# ── the path ──────────────────────────────────────────────────────────────────

def step_distances(start_distance: float, final_distance: float, n_steps: int,
                   exponent: float) -> list[float]:
    """
    O-Si distances for the n_steps + 1 frames, from `start_distance` down to
    `final_distance`. An exponent above 1 bunches the frames toward the end of
    the path, where the energy changes fastest.
    """
    distances = []
    for step in range(n_steps + 1):
        fraction = step / n_steps
        if exponent != 1.0:
            fraction = 1.0 - (1.0 - fraction) ** exponent
        distances.append(start_distance + fraction * (final_distance - start_distance))
    return distances


def resolve_start(frame, args) -> tuple[np.ndarray, list[str]]:
    """The starting oxygen position, plus any notes to print about how it was chosen."""
    silicon = sio.position_of(frame, args.silicon_id)
    notes = []

    if args.initial_position is not None:
        return np.array(args.initial_position, dtype=float), notes

    neighbour = sio.position_of(frame, args.backside_of)
    bond_length = sio.distance(frame, neighbour, silicon)
    if bond_length > 2.5:
        notes.append(
            f"atom {args.backside_of} is {bond_length:.3f} A from silicon "
            f"{args.silicon_id}, which is too far to be a bonded neighbour - "
            f"check that this is the bond you meant to attack"
        )
    # Continue through the silicon, away from the named neighbour.
    outward = normalize(sio.displacement(frame, neighbour, silicon))
    notes.append(
        f"starting {args.start_distance:.3f} A from silicon {args.silicon_id}, "
        f"opposite atom {args.backside_of} (Si-neighbour distance {bond_length:.3f} A)"
    )
    return silicon + args.start_distance * outward, notes


# ── atom types ────────────────────────────────────────────────────────────────

# Every frame adds one oxygen and two hydrogens, so two atom types have to be
# settled before any frame is written - not just the hydrogen type that a silica
# file is most likely to be missing.
WATER_ELEMENTS = ("O", "H")
ELEMENT_WORDS = {"O": "oxygen", "H": "hydrogen"}

# The three atoms a frame adds, in the order they are written and reported.
WATER_LABELS = ("O", "H1", "H2")
# Which element each of them is. Under --format poscar this is the whole story;
# under --format lammps it is the key into WATER_TYPES.
WATER_SPECIES = {"O": "O", "H1": "H", "H2": "H"}


def resolve_water_types(frame) -> tuple[list[str], list[str]]:
    """
    Reconcile WATER_TYPES with the atom types the input file actually uses, and
    return the type->element list every frame is written with (position i is
    type i+1) plus notes to print.

    LAMMPS-only: this is entirely about numeric atom types and the pair_coeff
    line, neither of which a POSCAR has.

    The elements behind the input's types are known - they come from its Masses
    table - so this checks the stronger thing the old mass comparison was a proxy
    for: that the type WATER_TYPES names really is the element it claims.
    """
    types = {element: int(WATER_TYPES[element]) for element in WATER_ELEMENTS}
    notes: list[str] = []

    if types["O"] == types["H"]:
        raise ValueError(
            f"WATER_TYPES gives the oxygen and the hydrogens the same atom type "
            f"({types['O']}); they must be different types"
        )

    # New elements are appended in ascending configured type order, so a file
    # missing both still ends up with a contiguous range whatever order
    # WATER_TYPES lists them in.
    specorder = list(frame.specorder)
    present = set(specorder)
    for element in sorted(WATER_ELEMENTS, key=lambda e: types[e]):
        if element not in present:
            specorder.append(element)

    for element in WATER_ELEMENTS:
        word = ELEMENT_WORDS[element]
        actual = specorder.index(element) + 1
        if actual != types[element]:
            if element in present:
                raise ValueError(
                    f"WATER_TYPES puts {word} at atom type {types[element]}, but "
                    f"{frame.source} already has {element} as type {actual}; set "
                    f"WATER_TYPES to match the file, or attack a different file"
                )
            raise ValueError(
                f"WATER_TYPES puts {word} at atom type {types[element]}, which "
                f"{frame.source} does not define and which would leave a gap in the "
                f"type range (next available is {actual}); atom types must be "
                f"contiguous"
            )
        if element in present:
            notes.append(f"water {word} uses existing atom type {actual}")
        else:
            notes.append(
                f"added atom type {actual} for {word} - the potential's element list "
                f"and pair_coeff line must be extended to match"
            )

    notes.append(
        f"each frame adds 3 atoms of 2 types: 1 oxygen of type {types['O']} and "
        f"2 hydrogens of type {types['H']}"
    )
    return specorder, notes


def water_spec(frame, args) -> tuple[dict, list[str]]:
    """
    Describe the three atoms every frame adds, in whatever terms the format
    works in, and report what was settled. This is the one place that has to know
    which format is in play: a LAMMPS file needs atom types, charges, and a
    molecule id, and a POSCAR needs none of them because it names elements
    outright.
    """
    spec: dict = {"species": dict(WATER_SPECIES)}

    if frame.fmt != sio.LAMMPS:
        notes = ["each frame adds 3 atoms: 1 oxygen and 2 hydrogens, written as "
                 "species O and H"]
        added = sio.new_species(frame, spec, WATER_LABELS)
        if added:
            plural = len(added) > 1
            notes.append(
                f"{' and '.join(added)} {'are' if plural else 'is'} not in "
                f"{args.input}, so {'they become new groups' if plural else 'it becomes a new group'} "
                f"on the species line - the POTCAR must be extended to match"
            )
        return spec, notes

    specorder, notes = resolve_water_types(frame)
    spec["specorder"] = specorder
    spec["molecule_id"] = sio.next_molecule_id(frame)
    spec["charges"] = {
        "O": 0.0 if args.o_charge is None else args.o_charge,
        "H1": 0.0 if args.h_charge is None else args.h_charge,
        "H2": 0.0 if args.h_charge is None else args.h_charge,
    }

    has_charges = "initial_charges" in frame.atoms.arrays
    if has_charges and args.o_charge is None and args.h_charge is None:
        print(f"NOTE: atom style {frame.atom_style} carries charges and the water is "
              f"being written with q = 0; correct for a potential that computes charges "
              f"itself, otherwise pass --o-charge / --h-charge.", file=sys.stderr)
    elif not has_charges and (args.o_charge is not None or args.h_charge is not None):
        raise ValueError(
            f"--o-charge/--h-charge were given but atom style {frame.atom_style} has no "
            f"charge column, so the charges cannot be written"
        )
    return spec, notes


# ── writing ───────────────────────────────────────────────────────────────────


def type_cell(label) -> str:
    """
    How an atom's kind reads in the report table: "t2" for a LAMMPS atom type,
    "O" for a POSCAR species, which needs no prefix to be recognisable.
    """
    return f"t{label}" if isinstance(label, int) else str(label)


MANIFEST_COLUMNS = [
    "frame", "o_si_distance", "fraction", "o_x", "o_y", "o_z",
    # Every id here indexes the frame that was written, not the input file. Under
    # --format poscar those differ, so si_id_out carries the target silicon's
    # index in the output - the one to name in an ICONST constraint.
    "o_id", "h1_id", "h2_id", "si_id_out", "o_type", "h_type",
    # Nearest existing atom to each of the three; the target silicon is skipped
    # for the oxygen only (see the search in the frame loop).
    "o_nearest_id", "o_nearest_type", "o_nearest_distance",
    "h1_nearest_id", "h1_nearest_type", "h1_nearest_distance",
    "h2_nearest_id", "h2_nearest_type", "h2_nearest_distance",
    "clash", "data_file",
]


# ── main ──────────────────────────────────────────────────────────────────────

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--format", required=True, choices=sorted(sio.FORMATS),
                        help="the file format to read and write. 'lammps' reads a data "
                             "file and writes <i>.data; 'poscar' reads a VASP POSCAR and "
                             "writes <i>/POSCAR, one directory per frame. Input and output "
                             "are always the same format")
    parser.add_argument("--input", required=True,
                        help="input structure: a LAMMPS data file or a POSCAR, per --format")
    parser.add_argument("--outdir", required=True, help="directory to write the frames into")
    parser.add_argument("--atom-style", default=None,
                        help="override the style in the file's 'Atoms # <style>' comment "
                             "(--format lammps only)")

    parser.add_argument("--silicon-id", type=int, required=True,
                        help="atom id of the silicon being attacked. Under --format poscar "
                             "there are no atom ids, so this is the 1-based line number in "
                             "the input file's coordinate block")

    start = parser.add_mutually_exclusive_group(required=True)
    start.add_argument("--initial-position", nargs=3, type=float, metavar=("X", "Y", "Z"),
                       help="explicit starting position for the water oxygen")
    start.add_argument("--backside-of", type=int, metavar="ATOM_ID",
                       help="start on the far side of the silicon from this atom - normally "
                            "the bridging oxygen whose Si-O bond is being broken. Numbered "
                            "like --silicon-id")
    parser.add_argument("--start-distance", type=float, default=4.0,
                        help="O-Si distance at the first frame, with --backside-of "
                             "(default: 4.0 A)")
    parser.add_argument("--final-distance", type=float, default=1.9,
                        help="O-Si distance at the last frame (default: 1.9 A, near the "
                             "pentacoordinate intermediate)")
    parser.add_argument("--n-steps", type=int, default=8,
                        help="number of steps along the path; writes n+1 frames, 1.data at "
                             "--start-distance through <n+1>.data at --final-distance "
                             "(default: 8)")
    parser.add_argument("--spacing-exponent", type=float, default=1.0,
                        help="1.0 spaces the frames evenly in distance; above 1.0 bunches "
                             "them toward the silicon (default: 1.0)")

    parser.add_argument("--freeze-substrate", action="store_true",
                        help="write selective dynamics freezing every atom that came from "
                             "the input file, leaving only the three water atoms free - the "
                             "setup for a CG minimization of the water against a rigid "
                             "framework (--format poscar only). This replaces any selective "
                             "dynamics the input already had; without it those flags are "
                             "carried through unchanged")

    parser.add_argument("--o-charge", type=float, default=None,
                        help="charge on the water oxygen, for atom styles with a q column "
                             "(default: 0.0, which is right for a potential that computes "
                             "charges itself, e.g. ReaxFF/QEq, and wrong for a fixed-charge "
                             "model)")
    parser.add_argument("--h-charge", type=float, default=None,
                        help="charge on each water hydrogen; see --o-charge")

    parser.add_argument("--water-template", default=None, metavar="FILE",
                        help="take the O-H lengths and H-O-H angle from a data file holding "
                             "a single water molecule")
    parser.add_argument("--oh-length", type=float, default=0.9572,
                        help="O-H bond length in A (default: 0.9572)")
    parser.add_argument("--hoh-angle", type=float, default=104.52,
                        help="H-O-H angle in degrees (default: 104.52)")
    parser.add_argument("--orientation", choices=("away", "toward", "perpendicular"),
                        default="away",
                        help="which way the hydrogens face along the path. 'away' (default) "
                             "points them away from the target silicon, so the oxygen's lone "
                             "pairs lead - the geometry a nucleophilic attack starts from. "
                             "'toward' points them at the silicon, so a proton arrives first "
                             "instead of the oxygen. 'perpendicular' puts the molecular plane "
                             "across the path, keeping both hydrogens off the line of travel")
    parser.add_argument("--h-rotation", type=float, default=0.0, metavar="DEGREES",
                        help="roll the hydrogens about the H-O-H bisector by this angle in "
                             "degrees, counterclockwise looking down the bisector toward the "
                             "oxygen (default: 0). The bisector is the rotation axis, so the "
                             "orientation above is unchanged and only the plane the hydrogens "
                             "lie in turns; under away/toward that is the roll about the O-Si "
                             "axis, which nothing else fixes. Always counterclockwise - a "
                             "negative angle is taken as the counterclockwise turn to the same "
                             "place")
    parser.add_argument("--h1-offset", nargs=3, type=float, default=None, metavar=("DX", "DY", "DZ"),
                        help="explicit O->H1 offset, overriding every other geometry option")
    parser.add_argument("--h2-offset", nargs=3, type=float, default=None, metavar=("DX", "DY", "DZ"),
                        help="explicit O->H2 offset; required with --h1-offset")

    parser.add_argument("--min-separation", type=float, default=1.5,
                        help="warn when a water atom lands closer than this to an existing "
                             "atom (default: 1.5 A; anything under this is a hard clash for "
                             "a reactive potential)")
    parser.add_argument("--fail-on-clash", action="store_true",
                        help="exit non-zero if any frame is below --min-separation")
    parser.add_argument("--dry-run", action="store_true",
                        help="report the path and the clash check without writing any files")

    args = parser.parse_args(argv)
    described = sio.DESCRIPTION[args.format]

    # Flags that describe something only one of the formats has. Rejecting them
    # outright beats writing a file that quietly ignored them.
    if sio.STYLE_OPTION[args.format] is None and args.atom_style is not None:
        parser.error(f"--atom-style does not apply to --format {args.format}: a "
                     f"{described} has no atom style")
    if not sio.SUPPORTS_CHARGES[args.format] and (args.o_charge is not None
                                                  or args.h_charge is not None):
        parser.error(f"--o-charge/--h-charge do not apply to --format {args.format}: a "
                     f"{described} has no per-atom charge column")
    if not sio.SUPPORTS_FREEZING[args.format] and args.freeze_substrate:
        parser.error(f"--freeze-substrate does not apply to --format {args.format}: a "
                     f"{described} has no selective dynamics. Freeze atoms in the LAMMPS "
                     f"input script instead, with a group and fix setforce 0 0 0")

    # --format is what picks the reader, so a mismatch would otherwise surface as
    # a confusing traceback from deep inside ASE.
    for path, what in ((args.input, "--input"),
                       (args.water_template, "--water-template")):
        if path is not None and not sio.looks_like(path, args.format):
            parser.error(f"{what} {path} does not look like a {described}, which is "
                         f"what --format {args.format} expects")

    if (args.h1_offset is None) != (args.h2_offset is None):
        parser.error("--h1-offset and --h2-offset must be given together")
    given = list(argv) if argv is not None else sys.argv[1:]
    if args.initial_position is not None and any(
            arg == "--start-distance" or arg.startswith("--start-distance=") for arg in given):
        parser.error("--start-distance applies to --backside-of; with --initial-position "
                     "the starting distance comes from the position itself")
    if args.n_steps < 1:
        parser.error("--n-steps must be at least 1")
    if args.final_distance <= 0:
        parser.error("--final-distance must be positive")
    if args.spacing_exponent <= 0:
        parser.error("--spacing-exponent must be positive")
    return args


def main(argv=None) -> int:
    args = parse_args(argv)

    topology = sio.find_topology(args.input, args.format)
    if topology:
        raise ValueError(
            f"{args.input} has topology sections ({', '.join(topology)}); adding atoms "
            f"would leave them inconsistent. This script is for reactive/pairwise "
            f"potentials, where water needs no Bonds or Angles."
        )

    data = sio.load(args.input, args.format, args.atom_style)

    silicon = sio.position_of(data, args.silicon_id)
    target_note = sio.check_target(data, args.silicon_id, "Si")
    if target_note:
        print(f"NOTE: {target_note}.", file=sys.stderr)

    start_position, start_notes = resolve_start(data, args)
    to_silicon = sio.displacement(data, start_position, silicon)
    start_distance = float(np.linalg.norm(to_silicon))
    if start_distance <= args.final_distance:
        raise ValueError(
            f"the starting position is {start_distance:.3f} A from silicon "
            f"{args.silicon_id}, which is already inside --final-distance "
            f"({args.final_distance:.3f} A); there is no path to walk"
        )
    approach = to_silicon / start_distance

    spec, type_notes = water_spec(data, args)
    h1_offset, h2_offset, geometry_note = resolve_offsets(
        args.format, args, args.atom_style, approach)

    # Input ids to ids in the frames about to be written. Only the layout of the
    # input decides this, so it is the same for every frame.
    to_output = sio.index_map(data, spec, WATER_LABELS)
    silicon_out = to_output[args.silicon_id]

    if args.freeze_substrate:
        already = sio.frozen_count(data)
        type_notes.append(
            f"freezing all {data.n_input} atoms from {args.input}; only the 3 water "
            f"atoms are free to move"
            + (f" (this replaces the {already} already frozen in the input)"
               if already else "")
        )

    for note in start_notes + type_notes:
        print(f"  {note}")
    print(f"  water geometry {geometry_note}")
    print(f"  O-H {np.linalg.norm(h1_offset):.4f} / {np.linalg.norm(h2_offset):.4f} A, "
          f"H-O-H {np.degrees(np.arccos(np.dot(normalize(h1_offset), normalize(h2_offset)))):.2f} deg")
    # Named relative to --outdir, so "1.data" or "1/POSCAR" depending on format.
    first = sio.frame_path("", 1, args.format)
    last = sio.frame_path("", args.n_steps + 1, args.format)
    print(f"  walking the oxygen from {start_distance:.3f} A to "
          f"{args.final_distance:.3f} A of silicon {args.silicon_id} "
          f"in {args.n_steps} steps, writing {first} through {last}")
    if silicon_out != args.silicon_id:
        print(f"  the water shifts the numbering: silicon {args.silicon_id} is atom "
              f"{silicon_out} in every frame written, and the manifest's ids all index "
              f"the frames, not {args.input}")

    outdir = Path(args.outdir)
    if not args.dry_run:
        outdir.mkdir(parents=True, exist_ok=True)

    distances = step_distances(start_distance, args.final_distance,
                               args.n_steps, args.spacing_exponent)
    span = start_distance - args.final_distance
    rows = []
    messages: list[str] = []
    clashes = 0
    wrapped_atoms = 0

    print()
    print(f"{'':>5}  {'':>7}  {'':>5}  "
          f"{'nearest existing atom to each, in A (target silicon skipped for O only)':^70}")
    print(f"{'frame':>5}  {'O-Si':>7}  {'frac':>5}  "
          f"{'O':>22}  {'H1':>22}  {'H2':>22}  file")

    # Frames are numbered from 1, so the file names line up with the frame count
    # rather than with the step count.
    for index, target_distance in enumerate(distances, start=1):
        # Measured back from the silicon so the last frame sits at exactly
        # --final-distance, whatever rounding the start position carries.
        oxygen = silicon - target_distance * approach
        positions = {"O": oxygen, "H1": oxygen + h1_offset, "H2": oxygen + h2_offset}

        relative_path = sio.frame_path("", index, args.format)
        frame_file = outdir / relative_path
        where = f"frame {index}"

        frame, ids, _, wrapped = sio.add_water(
            data, positions, WATER_LABELS, spec, args.freeze_substrate)
        wrapped_atoms += wrapped

        # The nearest existing atom to each of the three separately - a single
        # overall minimum hides a hydrogen buried in a wall whenever the oxygen
        # happens to be closer to something else. Searched in the frame with the
        # water already in it, excluding the water, so the ids reported are the
        # ones the written file uses and no atom is compared with itself.
        #
        # The target silicon is skipped for the oxygen only. That distance is
        # the coordinate being driven, it is already the O-Si column, and
        # leaving it in would mask every other contact on the last few frames.
        # For the hydrogens it is an ordinary neighbour and has to stay in:
        # under --orientation toward they reach the silicon before the oxygen
        # does, and excluding it would hide exactly the collision that mode
        # risks.
        water_ids = set(ids.values())
        nearest = {
            label: sio.nearest_existing(
                frame, positions[label],
                exclude=water_ids | {silicon_out} if label == "O" else water_ids,
            )[0]
            for label in WATER_LABELS
        }

        clashing = [label for label in WATER_LABELS
                    if nearest[label][2] < args.min_separation]
        if clashing:
            clashes += 1
            for label in clashing:
                near_id, near_type, near_distance = nearest[label]
                messages.append(
                    f"{where}: {label} is {near_distance:.3f} A from atom {near_id} "
                    f"(type {near_type}), below --min-separation {args.min_separation}"
                )

        wrapped_oxygen = sio.position_of(frame, ids["O"])
        fraction = (start_distance - target_distance) / span if span else 0.0

        if not args.dry_run:
            sio.write_frame(
                frame, frame_file,
                f"{sio.DESCRIPTION[args.format]} via generate_water_path.py, from "
                f"{args.input} - frame {index}/{len(distances)}, "
                f"O-Si {target_distance:.4f} A")

        row = {
            "frame": index,
            "o_si_distance": round(target_distance, 4),
            "fraction": round(fraction, 4),
            "o_x": round(float(wrapped_oxygen[0]), 4),
            "o_y": round(float(wrapped_oxygen[1]), 4),
            "o_z": round(float(wrapped_oxygen[2]), 4),
            "o_id": ids["O"], "h1_id": ids["H1"], "h2_id": ids["H2"],
            "si_id_out": silicon_out,
            "o_type": sio.species_label(frame, ids["O"]),
            "h_type": sio.species_label(frame, ids["H1"]),
            "clash": int(bool(clashing)),
            "data_file": str(relative_path),
        }
        for label in WATER_LABELS:
            near_id, near_type, near_distance = nearest[label]
            prefix = label.lower()
            row[f"{prefix}_nearest_id"] = near_id
            row[f"{prefix}_nearest_type"] = near_type
            row[f"{prefix}_nearest_distance"] = round(near_distance, 4)
        rows.append(row)

        cells = "  ".join(
            f"{f'{nearest[label][2]:.3f} id {nearest[label][0]} ({type_cell(nearest[label][1])})':>22}"
            for label in WATER_LABELS
        )
        marker = "  CLASH " + ",".join(clashing) if clashing else ""
        print(f"{index:>5}  {target_distance:>7.3f}  {fraction:>5.3f}  {cells}  "
              f"{relative_path}{marker}")

    if not args.dry_run:
        manifest = outdir / "manifest.csv"
        with open(manifest, "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS)
            writer.writeheader()
            writer.writerows(rows)

    print()
    if wrapped_atoms:
        print(f"{wrapped_atoms} of {3 * len(rows)} water atoms crossed a periodic boundary "
              f"and were wrapped back into the box; the structures are unchanged.")

    sys.stdout.flush()  # keep the warnings below the table when output is piped
    if messages:
        print(f"{len(messages)} warning(s):", file=sys.stderr)
        for message in messages:
            print(f"  {message}", file=sys.stderr)
        print(file=sys.stderr)

    if args.dry_run:
        print(f"dry run: nothing written to {outdir}/")
    else:
        print(f"wrote {sio.frame_path('', 1, args.format)} through "
              f"{sio.frame_path('', len(rows), args.format)} and manifest.csv to {outdir}/")
    print(f"Every frame is a seed, not a relaxed geometry: "
          f"{sio.restraint_advice(args.format)}.")

    if clashes and args.fail_on_clash:
        print(f"{clashes} frame(s) below --min-separation", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, OSError) as error:
        sys.exit(f"generate_water_path.py: {error}")
