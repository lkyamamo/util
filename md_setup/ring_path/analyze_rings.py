#!/usr/bin/env python3
"""
analyze_rings.py

Measure every ring written by silica_rings_to_lammps.py and write one row per
ring to a CSV, so that representatives can be chosen from a distribution rather
than by eye.

Per ring: size, pucker, eccentricity, mean radius, the aperture and the path
axis (see ring_geometry), plus two covariates that a single ring cannot show on
its own -

  smaller_ring_overlap  how many *smaller* rings share at least two atoms with
                        this one. A large ring fused to a 3- or 4-ring has its
                        geometry dictated by the strained neighbour. In a dense
                        glass this is close to universal for n >= 6, so it is
                        recorded as something to control for, not filtered on.
  undercoordinated_si   how many of the ring's silicons are already 3-fold in
                        the source structure, before any carving. Needs
                        --source; skipped without it.

Usage
-----
  python analyze_rings.py --ring-output ring_output --out rings.csv \\
      --source final.data

The printed summary is the table to check a run against: the aperture and pucker
medians per size are reproducible numbers that pin down the ring set. Apertures
are distances between the axis and atom centres and depend on no radius; the
contact radii only decide which close approaches get flagged.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from ring_geometry import CONTACT_RADII, O_TYPE, SI_TYPE, read_ring, ring_metrics  # noqa: E402

FIELDS = [
    "ring_id", "n", "file", "n_atoms", "axis_choice", "aperture",
    "planarity_rmsd", "eccentricity", "mean_radius", "mean_radius_inplane",
    "mean_radius_si", "mean_radius_o",
    "smaller_ring_overlap", "undercoordinated_si",
    "center_x", "center_y", "center_z",
    "axis_x", "axis_y", "axis_z",
    "aperture_centroid", "aperture_slid", "aperture_tilted",
    "hole_offset", "tilt_deg", "axis_at_cone_wall",
    "n_inside_contact_radius", "n_si_inside", "n_o_inside",
    "worst_contact_overlap", "nearest_atom_type",
    "centroid_x", "centroid_y", "centroid_z",
    "normal_x", "normal_y", "normal_z",
    "atom_ids",
]


def read_summary(path: Path) -> list[dict]:
    """
    Parse ring_summary.txt: one line per ring, "size file id id id ...", with the
    ids alternating Si, O in ring order.
    """
    rings = []
    for line in path.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        words = line.split()
        ids = [int(w) for w in words[2:]]
        rings.append({"n": int(words[0]), "file": words[1], "ids": ids,
                      "silicons": ids[0::2], "atom_set": set(ids)})
    return rings


def count_smaller_overlaps(rings: list[dict]) -> list[int]:
    """
    For each ring, how many rings of strictly smaller size share >= 2 atoms.

    Done through an atom -> rings index rather than by comparing every pair, so
    this stays linear in the number of shared atoms instead of quadratic in the
    number of rings.
    """
    holders = defaultdict(list)
    for index, ring in enumerate(rings):
        for atom in ring["atom_set"]:
            holders[atom].append(index)

    counts = []
    for index, ring in enumerate(rings):
        shared = Counter()
        for atom in ring["atom_set"]:
            for other in holders[atom]:
                if other != index and rings[other]["n"] < ring["n"]:
                    shared[other] += 1
        counts.append(sum(1 for c in shared.values() if c >= 2))
    return counts


def silicon_coordination(source: Path, cutoff: float = 1.9) -> dict[int, int]:
    """
    Si -> number of bonded oxygens in the source structure. Used only to flag
    rings that contain a silicon which was already under-coordinated before the
    ring was carved out.
    """
    from scipy.spatial import cKDTree

    box: dict[str, float] = {}
    rows, section = [], None
    for raw in source.open():
        line = raw.split("#")[0].strip()
        if not line:
            continue
        if line in ("Atoms", "Masses", "Velocities"):
            section = line
            continue
        words = line.split()
        if section == "Atoms":
            rows.append(words)
        elif section is None and words[-2:] in (["xlo", "xhi"], ["ylo", "yhi"], ["zlo", "zhi"]):
            box[words[-2]], box[words[-1]] = float(words[0]), float(words[1])

    ids = np.array([int(r[0]) for r in rows])
    types = np.array([int(r[1]) for r in rows])
    pos = np.array([[float(r[2]), float(r[3]), float(r[4])] for r in rows])
    lengths = np.array([box["xhi"] - box["xlo"], box["yhi"] - box["ylo"],
                        box["zhi"] - box["zlo"]])
    pos = np.mod(pos, lengths)

    si = np.where(types == SI_TYPE)[0]
    ox = np.where(types == O_TYPE)[0]
    tree = cKDTree(pos[ox], boxsize=lengths)
    return {int(ids[a]): len(hits)
            for a, hits in zip(si, tree.query_ball_point(pos[si], cutoff))}


def summarise(rows: list[dict]) -> None:
    """Print the per-size table that a run is checked against."""
    by_size = defaultdict(list)
    for row in rows:
        by_size[row["n"]].append(row)

    print(f"\n{'n':>2} {'N':>5} | {'aperture (A)':^32} | {'pucker':>7} | {'ecc':>6}")
    print(f"{'':>2} {'':>5} | {'p10':>6} {'med':>6} {'p90':>6} {'min':>6} {'max':>6} | "
          f"{'med':>7} | {'med':>6}")
    print("-" * 74)
    for n in sorted(by_size):
        group = by_size[n]
        ap = np.array([r["aperture"] for r in group])
        pl = np.array([r["planarity_rmsd"] for r in group])
        ecc = np.array([r["eccentricity"] for r in group])
        print(f"{n:2d} {len(group):5d} | {np.percentile(ap, 10):6.2f} {np.median(ap):6.2f} "
              f"{np.percentile(ap, 90):6.2f} {ap.min():6.2f} {ap.max():6.2f} | "
              f"{np.median(pl):7.2f} | {np.median(ecc):6.2f}")

    print("\nwhat the three candidate axes give (median aperture, and how far the "
          "opening sits\nfrom the centroid). The centroid drifts off the hole as rings "
          "get larger, which\nbiases any trend measured against ring size:")
    print(f"  {'n':>2} | {'centroid':>9} {'slid':>7} {'tilted':>7} | {'hole offset':>11} "
          f"| {'gained by sliding':>17}")
    for n in sorted(by_size):
        group = by_size[n]
        med = lambda key: np.median([r[key] for r in group])  # noqa: E731
        gain = np.median([r["aperture_slid"] - r["aperture_centroid"] for r in group])
        print(f"  {n:2d} | {med('aperture_centroid'):9.2f} {med('aperture_slid'):7.2f} "
              f"{med('aperture_tilted'):7.2f} | {med('hole_offset'):9.2f} A "
              f"| {gain:15.2f} A")
    wall = np.array([r["axis_at_cone_wall"] for r in rows])
    print(f"  {wall.sum()} of {len(rows)} rings ({100 * wall.mean():.1f}%) would tilt "
          f"further than the cone allows")

    print("\ncontact flags: atoms lying closer to the axis than their own contact\n"
          "radius (Si %.2f, O %.2f). These do not enter the aperture - they are the\n"
          "close approaches a hard-sphere picture would have called an overlap:"
          % (CONTACT_RADII[SI_TYPE], CONTACT_RADII[O_TYPE]))
    print(f"  {'n':>2} | {'rings with >=1':>14} | {'median atoms':>12} | {'worst overlap':>13} "
          f"| {'nearest atom':>12}")
    for n in sorted(by_size):
        group = by_size[n]
        inside = np.array([r["n_inside_contact_radius"] for r in group])
        worst = np.array([r["worst_contact_overlap"] for r in group])
        o_nearest = 100 * np.mean([r["nearest_atom_type"] == O_TYPE for r in group])
        print(f"  {n:2d} | {100 * (inside > 0).mean():13.0f}% | {np.median(inside):12.1f} "
              f"| {worst.max():11.2f} A | {o_nearest:11.0f}% O")

    overlap = np.array([r["smaller_ring_overlap"] for r in rows])
    print(f"\n{100 * (overlap > 0).mean():.1f}% of all rings share >= 2 atoms with a "
          f"smaller ring")
    if any(r["undercoordinated_si"] != "" for r in rows):
        bad = sum(1 for r in rows if r["undercoordinated_si"] not in ("", 0))
        print(f"{bad} ring(s) contain a silicon that is under-coordinated in the source")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--ring-output", required=True,
                        help="the OUTDIR written by silica_rings_to_lammps.py")
    parser.add_argument("--out", default="rings.csv", help="CSV to write")
    parser.add_argument("--source", default=None,
                        help="the structure the rings came from; enables the "
                             "under-coordinated silicon check")
    parser.add_argument("--si-radius", type=float, default=CONTACT_RADII[SI_TYPE],
                        help="contact radius used only to flag close approaches; "
                             "no aperture depends on it")
    parser.add_argument("--o-radius", type=float, default=CONTACT_RADII[O_TYPE],
                        help="see --si-radius")
    parser.add_argument("--axis", choices=("centroid", "slid", "tilted"), default="centroid",
                        help="which line to use as the path axis. 'centroid' (default) is "
                             "the ring centroid along the best-fit plane normal - simple "
                             "and reproducible. 'slid' keeps that direction but moves the "
                             "line to the widest point of the opening, which matters "
                             "because the centroid drifts off the hole as rings get "
                             "larger. 'tilted' also lets the direction tilt within a cone. "
                             "All three apertures are reported whichever is chosen")
    parser.add_argument("--limit", type=int, default=None,
                        help="only measure the first N rings, for a quick check")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    outdir = Path(args.ring_output)
    radii = {SI_TYPE: args.si_radius, O_TYPE: args.o_radius}
    print(f"axis: {args.axis}")
    print(f"contact radii (flags only): Si {radii[SI_TYPE]}, O {radii[O_TYPE]}")

    rings = read_summary(outdir / "ring_summary.txt")
    print(f"{len(rings)} rings in {outdir}/ring_summary.txt")
    overlaps = count_smaller_overlaps(rings)

    coordination = {}
    if args.source:
        coordination = silicon_coordination(Path(args.source))
        low = sum(1 for c in coordination.values() if c < 4)
        print(f"source: {len(coordination)} Si, {low} of them under-coordinated")

    if args.limit:
        rings, overlaps = rings[:args.limit], overlaps[:args.limit]

    rows = []
    for index, (entry, overlap) in enumerate(zip(rings, overlaps), start=1):
        ring = read_ring(str(outdir / entry["file"]))
        if ring.n != entry["n"]:
            raise ValueError(
                f"{entry['file']}: summary says {entry['n']}-membered but the file "
                f"holds {ring.n} silicons")

        # The size directory is part of the identity: ring_00018 exists in every
        # size_NN directory, so the bare stem collides and a lookup by it picks
        # an arbitrary one of them.
        row = {"ring_id": str(Path(entry["file"]).with_suffix("")),
               "file": entry["file"],
               "smaller_ring_overlap": overlap,
               "atom_ids": " ".join(str(i) for i in entry["ids"])}
        row.update(ring_metrics(ring, radii, axis_choice=args.axis))
        row["undercoordinated_si"] = (
            sum(1 for s in entry["silicons"] if coordination.get(s, 4) < 4)
            if coordination else "")
        rows.append(row)

        if index % 100 == 0 or index == len(rings):
            print(f"\r  measured {index}/{len(rings)}", end="", flush=True)
    print()

    identifiers = [row["ring_id"] for row in rows]
    if len(set(identifiers)) != len(identifiers):
        duplicated = [i for i, c in Counter(identifiers).items() if c > 1]
        raise ValueError(
            f"ring_id is not unique ({len(duplicated)} repeated, e.g. {duplicated[:3]}); "
            f"everything downstream looks rings up by it")

    with open(args.out, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    summarise(rows)
    print(f"\nwrote {len(rows)} rows to {args.out}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, OSError) as error:
        sys.exit(f"analyze_rings.py: {error}")
