#!/usr/bin/env python3
"""
generate_water_path.py

Seed a series of LAMMPS data files for a stepwise hydrolysis study: a water
molecule is placed in the pore/interstice of a silica structure and its oxygen
is walked along a straight line toward a targeted silicon, one data file per
step. Each file is an independent starting point for its own constrained
minimization / MD run, where the O-Si distance is restrained and everything
else - the hydrogens included - is free to relax.

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

Output
------
  <i>.data       the structure, ready for read_data. Frames are numbered from 1
                 to the frame count, so --n-steps N writes 1.data (the water at
                 --start-distance) through <N+1>.data (the water at
                 --final-distance).
  manifest.csv   one row per frame: the O-Si distance the run should restrain
                 to, the ids and types of the three new atoms, and the closest
                 approach to an existing atom other than the target silicon.

Every frame adds atoms of two types, one oxygen and two hydrogens. Both come
from the WATER_TYPES dictionary in the configuration block below - set it to
match your potential's element list before running. A type the input file does
not already define is added to Masses, and needs a matching entry on the
pair_coeff line, which this script cannot write for you.

Usage examples
--------------
  # back-side attack on the Si-O bond between Si 5 and bridging O 91
  python generate_water_path.py --input silica.data --outdir frames \\
      --silicon-id 5 --backside-of 91 --start-distance 4.5 \\
      --final-distance 1.9 --n-steps 10

  # explicit starting point, water geometry taken from an existing data file
  python generate_water_path.py --input silica.data --outdir frames \\
      --silicon-id 5 --initial-position 7.0 7.0 7.0 \\
      --water-template ../isolated_water.data --dry-run
"""

from __future__ import annotations

import argparse
import copy
import csv
import sys
from pathlib import Path

import numpy as np

# The shared data-file helpers live at the md_setup root, one level up.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import lammps_data as ld  # noqa: E402


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


def resolve_offsets(args, atom_style: str | None,
                    approach: np.ndarray) -> tuple[np.ndarray, np.ndarray, str]:
    """
    Work out the two rigid O->H offsets for the whole path. `approach` is the
    unit vector pointing from the starting position toward the silicon.
    """
    if args.h1_offset is not None:
        return np.array(args.h1_offset), np.array(args.h2_offset), "given explicitly"

    if args.water_template is not None:
        h1, h2 = template_offsets(args.water_template, atom_style)
        source = f"from {args.water_template}"
    else:
        h1, h2 = canonical_offsets(args.oh_length, args.hoh_angle)
        source = f"built from --oh-length {args.oh_length} and --hoh-angle {args.hoh_angle}"

    if args.orientation == "lone-pair":
        # Bisector anti-parallel to the approach: the hydrogens trail behind the
        # oxygen and its lone pairs face the silicon, which is the geometry a
        # nucleophilic attack starts from.
        h1, h2 = orient_offsets(h1, h2, -approach)
        source += ", hydrogens pointing away from the silicon"
    else:
        # Molecular plane perpendicular to the approach: both hydrogens stay off
        # the line of travel, which keeps them clear of a tight channel.
        bisector = perpendicular_to(approach)
        h1, h2 = orient_offsets(h1, h2, bisector, in_plane_hint=np.cross(approach, bisector))
        source += ", molecular plane perpendicular to the path"

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


def resolve_start(data, args) -> tuple[np.ndarray, list[str]]:
    """The starting oxygen position, plus any notes to print about how it was chosen."""
    silicon = ld.position_of(data, args.silicon_id)
    notes = []

    if args.initial_position is not None:
        return np.array(args.initial_position, dtype=float), notes

    neighbour = ld.position_of(data, args.backside_of)
    bond_length = ld.distance(data, neighbour, silicon)
    if bond_length > 2.5:
        notes.append(
            f"atom {args.backside_of} is {bond_length:.3f} A from silicon "
            f"{args.silicon_id}, which is too far to be a bonded neighbour - "
            f"check that this is the bond you meant to attack"
        )
    # Continue through the silicon, away from the named neighbour.
    outward = normalize(ld.displacement(data, neighbour, silicon))
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


# ── writing ───────────────────────────────────────────────────────────────────

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

    for label in ("O", "H1", "H2"):
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


MANIFEST_COLUMNS = [
    "frame", "o_si_distance", "fraction", "o_x", "o_y", "o_z",
    "o_id", "h1_id", "h2_id", "o_type", "h_type",
    "closest_atom_id", "closest_atom_type", "closest_distance", "closest_to",
    "clash", "data_file",
]


# ── main ──────────────────────────────────────────────────────────────────────

def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--input", required=True, help="input LAMMPS data file")
    parser.add_argument("--outdir", required=True, help="directory to write the frames into")
    parser.add_argument("--atom-style", default=None,
                        help="override the style in the file's 'Atoms # <style>' comment")

    parser.add_argument("--silicon-id", type=int, required=True,
                        help="atom id of the silicon being attacked")

    start = parser.add_mutually_exclusive_group(required=True)
    start.add_argument("--initial-position", nargs=3, type=float, metavar=("X", "Y", "Z"),
                       help="explicit starting position for the water oxygen")
    start.add_argument("--backside-of", type=int, metavar="ATOM_ID",
                       help="start on the far side of the silicon from this atom - normally "
                            "the bridging oxygen whose Si-O bond is being broken")
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
    parser.add_argument("--orientation", choices=("lone-pair", "perpendicular"),
                        default="lone-pair",
                        help="lone-pair points the hydrogens away from the silicon; "
                             "perpendicular puts the molecular plane across the path "
                             "(default: lone-pair)")
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

    data, atom_style = ld.load(args.input, args.atom_style, setting_name="--atom-style")
    if data.topology:
        sections = ", ".join(sorted(data.topology))
        raise ValueError(
            f"{args.input} has topology sections ({sections}); adding atoms would leave "
            f"them inconsistent. This script is for reactive/pairwise potentials, where "
            f"water needs no Bonds or Angles."
        )
    if ld.is_triclinic(data):
        print("NOTE: the box is triclinic; wrapping and distances here ignore tilt.",
              file=sys.stderr)

    silicon = ld.position_of(data, args.silicon_id)
    silicon_type = int(data.atoms.loc[args.silicon_id, "type"])
    silicon_mass = float(data.masses.loc[silicon_type, "mass"])
    if abs(silicon_mass - ld.ELEMENT_MASSES["Si"]) > MASS_TOLERANCE:
        print(f"NOTE: atom {args.silicon_id} has type {silicon_type} with mass "
              f"{silicon_mass}, which is not silicon.", file=sys.stderr)

    start_position, start_notes = resolve_start(data, args)
    to_silicon = ld.displacement(data, start_position, silicon)
    start_distance = float(np.linalg.norm(to_silicon))
    if start_distance <= args.final_distance:
        raise ValueError(
            f"the starting position is {start_distance:.3f} A from silicon "
            f"{args.silicon_id}, which is already inside --final-distance "
            f"({args.final_distance:.3f} A); there is no path to walk"
        )
    approach = to_silicon / start_distance

    oxygen_type, hydrogen_type, type_notes = apply_water_types(data)
    h1_offset, h2_offset, geometry_note = resolve_offsets(args, args.atom_style, approach)
    molecule_id = ld.next_molecule_id(data)

    oxygen_charge = 0.0 if args.o_charge is None else args.o_charge
    hydrogen_charge = 0.0 if args.h_charge is None else args.h_charge
    charges = {"O": oxygen_charge, "H1": hydrogen_charge, "H2": hydrogen_charge}
    has_charge_column = "q" in data.atoms.columns
    if has_charge_column and args.o_charge is None and args.h_charge is None:
        print(f"NOTE: atom style {atom_style} carries charges and the water is being "
              f"written with q = 0; correct for a potential that computes charges itself, "
              f"otherwise pass --o-charge / --h-charge.", file=sys.stderr)
    elif not has_charge_column and (args.o_charge is not None or args.h_charge is not None):
        raise ValueError(
            f"--o-charge/--h-charge were given but atom style {atom_style} has no charge "
            f"column, so the charges cannot be written"
        )

    for note in start_notes + type_notes:
        print(f"  {note}")
    print(f"  water geometry {geometry_note}")
    print(f"  O-H {np.linalg.norm(h1_offset):.4f} / {np.linalg.norm(h2_offset):.4f} A, "
          f"H-O-H {np.degrees(np.arccos(np.dot(normalize(h1_offset), normalize(h2_offset)))):.2f} deg")
    print(f"  walking the oxygen from {start_distance:.3f} A to "
          f"{args.final_distance:.3f} A of silicon {args.silicon_id} "
          f"in {args.n_steps} steps, writing 1.data through {args.n_steps + 1}.data")

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
    print(f"{'frame':>5}  {'O-Si':>7}  {'frac':>5}  {'closest atom other than the target':<34}  file")

    # Frames are numbered from 1, so the file names line up with the frame count
    # rather than with the step count.
    for index, target_distance in enumerate(distances, start=1):
        # Measured back from the silicon so the last frame sits at exactly
        # --final-distance, whatever rounding the start position carries.
        oxygen = silicon - target_distance * approach
        positions = {"O": oxygen, "H1": oxygen + h1_offset, "H2": oxygen + h2_offset}

        data_file = outdir / f"{index}.data"
        where = f"frame {index}"

        frame, ids, wrapped, wrap_messages = build_frame(
            data, positions, oxygen_type, hydrogen_type, charges, molecule_id, where)
        messages.extend(wrap_messages)
        wrapped_atoms += wrapped

        # Measured against the original structure, so the water is not compared
        # with itself, and skipping the target silicon, which the oxygen is
        # deliberately closing in on and which would otherwise mask a real clash
        # on the last few frames.
        closest = min(
            ((label,) + ld.nearest_existing(data, positions[label],
                                            exclude={args.silicon_id})[0]
             for label in ("O", "H1", "H2")),
            key=lambda entry: entry[3],
        )
        label, closest_id, closest_type, closest_distance = closest
        clash = closest_distance < args.min_separation
        if clash:
            clashes += 1
            messages.append(
                f"{where}: {label} is {closest_distance:.3f} A from atom {closest_id} "
                f"(type {closest_type}), below --min-separation {args.min_separation}"
            )

        wrapped_oxygen = ld.position_of(frame, ids["O"])
        fraction = (start_distance - target_distance) / span if span else 0.0

        if not args.dry_run:
            ld.write_data_file(
                frame, str(data_file),
                f"LAMMPS data file via generate_water_path.py, from {args.input} - "
                f"frame {index}/{len(distances)}, O-Si {target_distance:.4f} A",
                atom_style)

        rows.append({
            "frame": index,
            "o_si_distance": round(target_distance, 4),
            "fraction": round(fraction, 4),
            "o_x": round(float(wrapped_oxygen[0]), 4),
            "o_y": round(float(wrapped_oxygen[1]), 4),
            "o_z": round(float(wrapped_oxygen[2]), 4),
            "o_id": ids["O"], "h1_id": ids["H1"], "h2_id": ids["H2"],
            "o_type": oxygen_type, "h_type": hydrogen_type,
            "closest_atom_id": closest_id,
            "closest_atom_type": closest_type,
            "closest_distance": round(closest_distance, 4),
            "closest_to": label,
            "clash": int(clash),
            "data_file": data_file.name,
        })

        marker = " CLASH" if clash else ""
        print(f"{index:>5}  {target_distance:>7.3f}  {fraction:>5.3f}  "
              f"{f'{closest_distance:.3f} A to id {closest_id} (type {closest_type}, {label})':<34}"
              f"{marker}  {data_file.name}")

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
        print(f"wrote 1.data through {len(rows)}.data and manifest.csv to {outdir}/")
    print("Every frame is a seed, not a relaxed geometry: relax each one with the O-Si "
          "distance restrained and the hydrogens free before reading any energy off it.")

    if clashes and args.fail_on_clash:
        print(f"{clashes} frame(s) below --min-separation", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, OSError) as error:
        sys.exit(f"generate_water_path.py: {error}")
