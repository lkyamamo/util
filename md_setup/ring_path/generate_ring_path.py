#!/usr/bin/env python3
"""
generate_ring_path.py

Seed a series of LAMMPS data files for a ring-permeation study: a rigid water
molecule is walked along a straight line through the middle of a silica ring,
one data file per step, entering on one side and leaving on the other.

The sibling of generate_water_path.py. Building, aiming and rolling the molecule
is shared with it through water_path/water_geometry.py; putting it into a file is
not, because that script writes POSCARs as well as data files and goes through
structure_io, while this one is LAMMPS-only and goes through lammps_data and
water_placement.

What differs in the physics is the coordinate. There the water walks *to* a
target silicon and stops at a bond length; here it starts outside the ring,
passes through the opening, and ends outside on the far side, so the coordinate
is a signed displacement along the ring's axis, running from -half_length through
0 to +half_length.

What this script is and is not
------------------------------
The frames are rigid: the cluster never moves and the water keeps one fixed
internal geometry, so the energy difference between two frames is a pure
interaction energy. That makes these frames usable for single points directly,
which is the opposite of the situation generate_water_path.py warns about - there
the hydrogens have to rearrange along a reaction coordinate, here nothing is
supposed to react.

What it does not capture is the ring dilating to let the water through, which is
most of the real permeation barrier. A frozen ring gives a steric barrier past a
fixed obstacle. That is a well-defined quantity and comparable across rings; it
is not the activated barrier.

The reference frame
-------------------
The energy of a path is referenced to its own first frame, which means that frame
has to be genuinely non-interacting: far enough out that no water atom is inside
the potential's cutoff of any cluster atom. The *hydrogens* set that distance,
not the oxygen, because under --orientation away they trail about an angstrom
behind it. Both ends are checked against --min-endpoint-separation and the run
fails if either is too close, because a silently interacting reference quietly
biases every energy on the path.

The roll
--------
The roll of the water about the axis is a free parameter that nothing else fixes,
and it matters more here than in an open approach: threading an opening, the
hydrogens hit the rim first. --h-rotation takes several angles at once and writes
each into its own subdirectory.

Run all of them and take the best. Screening at the ring centre alone - one frame
per roll, discard whatever clashes - looks like it should work and does not:
measured over 36 rings it keeps the best roll for only 22, discards it for 4, and
for 10 leaves no roll standing at all, because in a tight ring every roll clashes
at the centre. The constriction is not at the centre; it sits around |s| = 0.8 A
and reaches 2 A, since the centroid is not the middle of the opening. If a screen
is needed, screen on the minimum over a coarse span of s rather than on s = 0.

Output
------
  <outdir>/roll_<deg>/<i>.data   the frames, numbered from 1
  <outdir>/manifest.csv          one row per frame per roll: the axis position,
                                 the ids and types of the three new atoms, and
                                 the closest approach of each to a cluster atom
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "water_path"))
import lammps_data as ld  # noqa: E402
from ring_geometry import nearest_distance_fn  # noqa: E402
from water_geometry import (  # noqa: E402
    WATER_LABELS, describe_offsets, normalize, resolve_offsets,
)
from water_placement import (  # noqa: E402
    NEAREST_COLUMNS, apply_water_types, build_frame, nearest_columns,
    nearest_report, template_offsets,
)

MANIFEST_COLUMNS = [
    "roll", "frame", "axis_position", "fraction", "o_x", "o_y", "o_z",
    "o_id", "h1_id", "h2_id", "o_type", "h_type",
    *NEAREST_COLUMNS,
    "min_separation", "clash", "data_file",
]


def axis_positions(half_length: float, n_steps: int) -> list[float]:
    """
    The signed axis positions of the n_steps + 1 frames, evenly spaced from
    -half_length to +half_length.

    Evenly spaced and covering the full range rather than half of it: a ring in
    a glass has no symmetry, so the way in and the way out are different, and
    the barrier is not necessarily at zero.
    """
    return list(np.linspace(-half_length, half_length, n_steps + 1))


def required_half_length(nearest, center: np.ndarray, axis: np.ndarray,
                         offsets: list[tuple[np.ndarray, np.ndarray]],
                         threshold: float, margin: float = 0.5,
                         step: float = 0.25, limit: float = 40.0) -> float:
    """
    The smallest half-length at which both ends of the path are non-interacting:
    every water atom, at every roll, further than `threshold` from every cluster
    atom, plus `margin`, rounded up to a multiple of `step`.

    Half-length, not length: the water runs from -h through the ring centre to
    +h, so it travels 2h in total.

    Worth solving, because the intuitive estimate is wrong in both directions.

    It is not "cutoff plus an O-H bond". The atom that reaches furthest toward
    the departing water is neither at the ring centre nor on the axis: measured
    across this glass it sits a median 2.31 A off-axis, so it needs only
    sqrt(rc^2 - 2.31^2) = 4.99 A of axial separation, but it is also displaced a
    median 1.75 A along the axis, because the charge-balancing oxygens hang off
    the ring plane. With the trailing hydrogen that comes to ~6.9 A against a
    5.5 A cutoff. The hydrogens are the small term; the cluster's axial extent is
    the large one.

    Nor is a single fixed number the simpler option. Over 1,722 rings the exact
    minimum ranges 6.02 to 8.91 A, so any fixed value safe for all of them is
    ~9 A - longer than the median ring needs, which costs either resolution or
    frames on every short path. Solving per ring is what keeps the paths short.

    Clearance grows monotonically with |s| once the water is clear of the
    cluster, so a scan outward from zero finds the first value that works.
    """
    position = step
    while position <= limit:
        probes = []
        for h1, h2 in offsets:
            for sign in (-1.0, 1.0):
                oxygen = center + sign * position * axis
                probes += [oxygen, oxygen + h1, oxygen + h2]
        if float(nearest(np.array(probes)).min()) >= threshold:
            return float(np.ceil((position + margin) / step) * step)
        position += step
    raise ValueError(
        f"no half-length below {limit:g} A puts both ends further than {threshold:g} A "
        f"from the cluster; check the axis and the box")


def resolve_axis(args) -> tuple[np.ndarray, np.ndarray, str]:
    """The ring centre and axis, from the cluster manifest or given outright."""
    if args.center is not None:
        return (np.array(args.center, dtype=float),
                normalize(np.array(args.axis, dtype=float)),
                "given explicitly")

    with open(args.manifest) as handle:
        matches = [row for row in csv.DictReader(handle) if row["ring_id"] == args.ring]
    if not matches:
        raise ValueError(f"{args.manifest}: no ring {args.ring!r}")
    if len(matches) > 1:
        # Never silently take one of them: the whole path depends on this axis,
        # and picking the wrong ring's axis produces plausible-looking nonsense.
        raise ValueError(
            f"{args.manifest}: ring {args.ring!r} appears {len(matches)} times, so the "
            f"axis is ambiguous; ring ids must be unique")
    row = matches[0]
    return (np.array([float(row[f"center_{a}"]) for a in "xyz"]),
            normalize(np.array([float(row[f"axis_{a}"]) for a in "xyz"])),
            f"ring {args.ring} in {args.manifest}")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("--input", required=True,
                        help="the ring cluster, from build_ring_cluster.py")
    parser.add_argument("--outdir", required=True, help="directory to write the frames into")
    parser.add_argument("--atom-style", default=None,
                        help="override the style in the file's 'Atoms # <style>' comment")

    where = parser.add_argument_group("where the ring is")
    where.add_argument("--manifest", default=None,
                       help="cluster_manifest.csv from build_ring_cluster.py")
    where.add_argument("--ring", default=None, metavar="RING_ID",
                       help="which row of --manifest to use")
    where.add_argument("--center", nargs=3, type=float, default=None,
                       metavar=("X", "Y", "Z"),
                       help="the point the path passes through, instead of --manifest")
    where.add_argument("--axis", nargs=3, type=float, default=None,
                       metavar=("X", "Y", "Z"), help="the direction of travel")

    path = parser.add_argument_group("the path")
    path.add_argument("--half-length", type=float, default=None,
                      help="how far either side of the ring's centre the path reaches, in "
                           "A; the water travels twice this. Default is to solve for the "
                           "smallest value that leaves both ends outside "
                           "--min-endpoint-separation, which depends on the ring and is "
                           "not worth guessing")
    path.add_argument("--half-length-margin", type=float, default=0.5,
                      help="how far past the minimum the solved half-length goes "
                           "(default: 0.5)")
    path.add_argument("--n-steps", type=int, default=14,
                      help="steps along the path; writes n+1 frames (default: 14)")
    path.add_argument("--positions", nargs="+", type=float, default=None,
                      metavar="S",
                      help="explicit axis positions instead of a range - '--positions 0' "
                           "writes the single frame at the ring centre, which is the "
                           "cheap way to screen rolls")

    water = parser.add_argument_group("the water")
    water.add_argument("--orientation", choices=("away", "toward", "perpendicular"),
                       default="away",
                       help="which way the hydrogens face along the path. 'away' "
                            "(default) trails them behind the oxygen, so the oxygen "
                            "leads the whole way through")
    water.add_argument("--h-rotation", nargs="+", type=float, default=[0.0],
                       metavar="DEGREES",
                       help="roll the hydrogens about the H-O-H bisector, in degrees. "
                            "Several values are allowed and each gets its own "
                            "subdirectory (default: 0)")
    water.add_argument("--water-template", default=None, metavar="FILE",
                       help="take the O-H lengths and H-O-H angle from a data file "
                            "holding a single water molecule")
    water.add_argument("--oh-length", type=float, default=0.9572,
                       help="O-H bond length in A (default: 0.9572)")
    water.add_argument("--hoh-angle", type=float, default=104.52,
                       help="H-O-H angle in degrees (default: 104.52)")
    water.add_argument("--o-charge", type=float, default=None,
                       help="charge on the water oxygen, for atom styles with a q column")
    water.add_argument("--h-charge", type=float, default=None,
                       help="charge on each water hydrogen; see --o-charge")

    checks = parser.add_argument_group("checks")
    checks.add_argument("--min-separation", type=float, default=1.5,
                        help="warn when a water atom lands closer than this to a cluster "
                             "atom (default: 1.5). Expect this to fire in the middle of a "
                             "tight ring; that is the physics, not a setup error")
    checks.add_argument("--min-endpoint-separation", type=float, default=5.5,
                        help="the first and last frame must have every water atom at "
                             "least this far from every cluster atom, so they are outside "
                             "the potential's cutoff and usable as the energy reference "
                             "(default: 5.5, the USC two-body cutoff)")
    checks.add_argument("--no-endpoint-check", action="store_true",
                        help="skip the endpoint check; the energies are then not "
                             "referenced to a non-interacting frame")
    checks.add_argument("--dry-run", action="store_true",
                        help="report the path and the checks without writing any files")

    args = parser.parse_args(argv)

    if (args.center is None) != (args.axis is None):
        parser.error("--center and --axis must be given together")
    if args.center is None and not (args.manifest and args.ring):
        parser.error("give either --center/--axis or --manifest/--ring")
    if args.positions is None and args.n_steps < 1:
        parser.error("--n-steps must be at least 1")
    if args.half_length is not None and args.half_length <= 0:
        parser.error("--half-length must be positive")
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

    center, axis, axis_note = resolve_axis(args)

    # Offsets for every roll up front: the half-length has to clear the cluster
    # for all of them, since one manifest covers the whole set.
    # The template molecule, if any, is read once here rather than inside
    # resolve_offsets: that function is format-independent and does not open
    # files, so reading the LAMMPS one is this script's job.
    template = (template_offsets(args.water_template, args.atom_style)
                if args.water_template else None)
    geometry = {}
    for roll in args.h_rotation:
        geometry[roll] = resolve_offsets(
            axis, orientation=args.orientation, h_rotation=roll,
            oh_length=args.oh_length, hoh_angle=args.hoh_angle,
            offsets=template)

    half_length, half_length_note = args.half_length, "given"
    if args.positions is None and half_length is None:
        nearest = nearest_distance_fn(
            data.atoms[["x", "y", "z"]].to_numpy(dtype=float), ld.box_lengths(data))
        half_length = required_half_length(
            nearest, center, axis, [(h1, h2) for h1, h2, _ in geometry.values()],
            args.min_endpoint_separation, args.half_length_margin)
        half_length_note = (f"solved: the smallest that clears "
                            f"{args.min_endpoint_separation:g} A at both ends")

    positions_along = (args.positions if args.positions is not None
                       else axis_positions(half_length, args.n_steps))

    oxygen_type, hydrogen_type, type_notes = apply_water_types(data)
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

    for note in type_notes:
        print(f"  {note}")
    print(f"  axis {axis_note}: through "
          f"({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}) along "
          f"({axis[0]:.4f}, {axis[1]:.4f}, {axis[2]:.4f})")
    if half_length is not None:
        print(f"  half-length {half_length:.3f} A ({half_length_note})")
    print(f"  {len(positions_along)} frame(s) from {positions_along[0]:+.3f} to "
          f"{positions_along[-1]:+.3f} A along the axis, at "
          f"{len(args.h_rotation)} roll(s): "
          f"{', '.join(f'{r:g}' for r in args.h_rotation)} deg")

    outdir = Path(args.outdir)
    rows: list[dict] = []
    messages: list[str] = []
    clashes = 0
    wrapped_atoms = 0
    endpoint_failures: list[str] = []

    for roll in args.h_rotation:
        h1_offset, h2_offset, geometry_note = geometry[roll]
        roll_dir = outdir / f"roll_{int(round(roll)) % 360:03d}"
        if not args.dry_run:
            roll_dir.mkdir(parents=True, exist_ok=True)

        print(f"\nroll {roll:g} deg: water geometry {geometry_note}")
        print(f"  {describe_offsets(h1_offset, h2_offset)}")
        print(f"{'':>5}  {'':>8}  {'':>5}  "
              f"{'nearest cluster atom to each, in A':^70}")
        print(f"{'frame':>5}  {'axis':>8}  {'frac':>5}  "
              f"{'O':>22}  {'H1':>22}  {'H2':>22}  file")

        span = positions_along[-1] - positions_along[0]
        for index, along in enumerate(positions_along, start=1):
            oxygen = center + along * axis
            positions = {"O": oxygen, "H1": oxygen + h1_offset, "H2": oxygen + h2_offset}
            data_file = roll_dir / f"{index}.data"
            where = f"roll {roll:g} frame {index}"

            frame, ids, wrapped, wrap_messages = build_frame(
                data, positions, oxygen_type, hydrogen_type, charges, molecule_id, where)
            messages.extend(wrap_messages)
            wrapped_atoms += wrapped

            # No exclusions: unlike an approach to a named silicon, nothing here
            # is a deliberate target, so every cluster atom is an honest
            # neighbour and every close contact is a real one.
            nearest = nearest_report(data, positions)
            closest = min(nearest[label][2] for label in WATER_LABELS)

            clashing = [label for label in WATER_LABELS
                        if nearest[label][2] < args.min_separation]
            if clashing:
                clashes += 1
                for label in clashing:
                    near_id, near_type, near_distance = nearest[label]
                    messages.append(
                        f"{where}: {label} is {near_distance:.3f} A from atom {near_id} "
                        f"(type {near_type}), below --min-separation {args.min_separation}")

            # Only a full path has a reference frame to protect. With explicit
            # --positions the caller is sampling, not building a path, and the
            # ends are wherever they asked for.
            if (not args.no_endpoint_check and args.positions is None
                    and index in (1, len(positions_along))
                    and closest < args.min_endpoint_separation):
                endpoint_failures.append(
                    f"roll {roll:g} frame {index} (axis {along:+.3f} A): closest water "
                    f"atom is {closest:.3f} A from the cluster, inside the "
                    f"{args.min_endpoint_separation:.3f} A reference cutoff")

            if not args.dry_run:
                ld.write_data_file(
                    frame, str(data_file),
                    f"LAMMPS data file via generate_ring_path.py, from {args.input} - "
                    f"roll {roll:g} deg, frame {index}/{len(positions_along)}, "
                    f"axis {along:+.4f} A",
                    atom_style)

            wrapped_oxygen = ld.position_of(frame, ids["O"])
            rows.append({
                "roll": roll,
                "frame": index,
                "axis_position": round(float(along), 4),
                "fraction": round((along - positions_along[0]) / span, 4) if span else 0.0,
                "o_x": round(float(wrapped_oxygen[0]), 4),
                "o_y": round(float(wrapped_oxygen[1]), 4),
                "o_z": round(float(wrapped_oxygen[2]), 4),
                "o_id": ids["O"], "h1_id": ids["H1"], "h2_id": ids["H2"],
                "o_type": oxygen_type, "h_type": hydrogen_type,
                "min_separation": round(float(closest), 4),
                "clash": int(bool(clashing)),
                "data_file": str(data_file.relative_to(outdir)),
                **nearest_columns(nearest),
            })

            cells = "  ".join(
                f"{f'{nearest[label][2]:.3f} id {nearest[label][0]} (t{nearest[label][1]})':>22}"
                for label in WATER_LABELS
            )
            marker = "  CLASH " + ",".join(clashing) if clashing else ""
            print(f"{index:>5}  {along:>+8.3f}  {rows[-1]['fraction']:>5.3f}  {cells}  "
                  f"{data_file.name}{marker}")

    if not args.dry_run:
        outdir.mkdir(parents=True, exist_ok=True)
        with open(outdir / "manifest.csv", "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS)
            writer.writeheader()
            writer.writerows(rows)

    # Per-roll summary: the tightest squeeze anywhere on the path. This is what
    # the roll screen reads - a roll that puts a hydrogen in the rim shows up as
    # a much smaller number than its neighbours.
    print(f"\n{'roll':>6}  {'tightest approach (A)':>21}  {'at axis':>9}  {'clashing frames':>15}")
    for roll in args.h_rotation:
        subset = [r for r in rows if r["roll"] == roll]
        tightest = min(subset, key=lambda r: r["min_separation"])
        print(f"{roll:>6g}  {tightest['min_separation']:>21.3f}  "
              f"{tightest['axis_position']:>+9.3f}  "
              f"{sum(r['clash'] for r in subset):>15d}")

    if wrapped_atoms:
        print(f"\n{wrapped_atoms} of {3 * len(rows)} water atoms crossed a periodic "
              f"boundary and were wrapped back into the box; the structures are unchanged.")

    sys.stdout.flush()  # keep the warnings below the table when output is piped
    if messages:
        print(f"\n{len(messages)} close-contact warning(s):", file=sys.stderr)
        for message in messages:
            print(f"  {message}", file=sys.stderr)

    if endpoint_failures:
        print(f"\n{len(endpoint_failures)} endpoint(s) are not far enough out to serve as "
              f"the energy reference:", file=sys.stderr)
        for failure in endpoint_failures:
            print(f"  {failure}", file=sys.stderr)
        print(f"  Raise --half-length (currently {half_length:.3f} A) until both ends "
              f"clear {args.min_endpoint_separation:.3f} A, or drop it entirely and let "
              f"the half-length be solved for. --no-endpoint-check skips this, and then "
              f"the energies are not referenced to a non-interacting frame.",
              file=sys.stderr)
        return 1

    if args.dry_run:
        print(f"\ndry run: nothing written to {outdir}/")
    else:
        print(f"\nwrote {len(rows)} frames across {len(args.h_rotation)} roll(s) and "
              f"manifest.csv to {outdir}/")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, OSError) as error:
        sys.exit(f"generate_ring_path.py: {error}")
