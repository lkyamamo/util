#!/usr/bin/env python3
"""
build_ring_cluster.py

Turn a ring carved out of a silica network into a cluster that can be put in
vacuum and computed on.

The problem this solves
-----------------------
A Guttman ring is n silicons and n oxygens, and that is not a usable model of
anything. Under a fixed-charge potential like USC (Si +1.2, O -0.6) it carries a
net charge of +0.6n, so its energies are dominated by the net charge rather than
by its shape. The fix is stoichiometry: add n more oxygens, taken from the
silicons' own bonds in the source structure, to reach Si_n O_2n. That is SiO2,
and 1.2n - 0.6(2n) = 0 exactly.

Which oxygens get added is decided round-robin: walk the ring's silicons in ring
order, give each one an oxygen it is actually bonded to, skip any silicon with
nothing left, and repeat until n have been added. In a dense glass this finishes
in a single pass and every silicon ends up 3-fold. A silicon with two candidates
picks by --tie-break.

What this does not do
---------------------
The added oxygens are singly coordinated. That is what balances the charge, and
it is not the same as a chemically sensible termination: there are no hydrogens
here, and in a DFT calculation each of these is a bare oxygen with an unpaired
electron unless it pairs with its silicon's dangling bond. The cluster is a
charge-balanced fragment for a frozen single-point scan, not a molecule.

The aperture is re-measured with the added oxygens in place, because an added
oxygen can sit in the channel and narrow the opening the ring appeared to have.
Compare `aperture_cluster` with `aperture_ring` in the manifest; a large drop
means the tie-break put an oxygen in the way.

Output
------
  <outdir>/size_NN/ring_XXXXX.data   one cluster per ring, centred in its box
  <outdir>/cluster_manifest.csv      one row per cluster: composition, the
                                     charge sum, the path axis in the cluster's
                                     own coordinates, and both apertures
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

import numpy as np
from pymatgen.io.lammps.data import LammpsBox

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import lammps_data as ld  # noqa: E402
from ring_geometry import (  # noqa: E402
    CONTACT_RADII, O_TYPE, SI_TYPE, clearance_along, contact_violations,
    minimum_image, read_ring, read_structure, si_o_neighbours,
)

# Formal charges the USC potential assigns each species. They live in the
# potential file, not in the data file, so nothing downstream can check this for
# us - the stoichiometry assertion here is the only guard.
USC_CHARGE = {SI_TYPE: 1.2, O_TYPE: -0.6}

MANIFEST_COLUMNS = [
    "ring_id", "n", "file", "n_atoms", "n_si", "n_o", "charge_sum", "passes",
    "aperture_ring", "aperture_cluster", "aperture_drop",
    "n_inside_contact_radius", "worst_contact_overlap",
    "added_ids", "min_added_axis_distance",
    "cluster_radius", "box_side",
    "center_x", "center_y", "center_z", "axis_x", "axis_y", "axis_z",
]


def choose_oxygens(ring, neighbours: dict[int, list[int]], structure,
                   center: np.ndarray, axis: np.ndarray,
                   tie_break: str) -> tuple[list[int], int]:
    """
    Pick the n oxygens that balance the ring, round-robin over its silicons.

    Returns (oxygen ids in the order they were added, number of passes taken).
    Raises if the ring's silicons between them do not have n spare oxygens,
    which would leave the cluster charged; that is a real property of the
    structure, not a bug, so it is reported rather than worked around.
    """
    index = structure.index_of()
    ring_atoms = set(int(i) for i in ring.ids)
    silicons = [int(s) for s in ring.silicon_ids]
    added: list[int] = []
    passes = 0

    def rank(oxygen_id: int) -> float:
        if tie_break == "lowest-id":
            return float(oxygen_id)
        # Prefer the oxygen furthest from the path axis, so balancing the charge
        # does not narrow the channel the water has to pass through.
        position = structure.positions[index[oxygen_id]]
        offset = position - center
        return -float(np.linalg.norm(offset - np.dot(offset, axis) * axis))

    while len(added) < len(silicons):
        passes += 1
        progressed = False
        for silicon in silicons:
            if len(added) >= len(silicons):
                break
            candidates = [o for o in neighbours.get(silicon, [])
                          if o not in ring_atoms and o not in added]
            if not candidates:
                continue
            added.append(min(candidates, key=rank))
            progressed = True
        if not progressed:
            raise ValueError(
                f"ring has {len(silicons)} silicons but only {len(added)} spare "
                f"oxygens between them, so it cannot be balanced to SiO2")
    return added, passes


def build(ring, added: list[int], structure) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Assemble the cluster's atoms, unwrapped and contiguous.

    The ring already comes unwrapped. Each added oxygen is placed relative to the
    silicon it is bonded to, through the minimum image, so the bond stays intact
    even when the pair straddles a periodic boundary in the source.
    """
    index = structure.index_of()
    ring_atoms = {int(a): p for a, p in zip(ring.ids, ring.positions)}
    silicons = [int(s) for s in ring.silicon_ids]

    positions = [ring.positions]
    ids = [ring.ids]
    types = [ring.types]

    extra = []
    for oxygen in added:
        # The silicon this oxygen came from: the ring silicon it is closest to
        # under the minimum image, which is the one it is bonded to.
        source_o = structure.positions[index[oxygen]]
        best, best_distance = None, np.inf
        for silicon in silicons:
            delta = minimum_image(source_o - structure.positions[index[silicon]],
                                  structure.lengths)
            distance = float(np.linalg.norm(delta))
            if distance < best_distance:
                best, best_distance = ring_atoms[silicon] + delta, distance
        extra.append(best)

    positions.append(np.array(extra))
    ids.append(np.array(added))
    types.append(np.full(len(added), O_TYPE))
    return np.concatenate(ids), np.concatenate(types), np.concatenate(positions)


def write_cluster(ring, ids, types, positions, out_path: Path,
                  box_side: float | None, header: str) -> tuple[float, np.ndarray]:
    """
    Write the cluster, centred in its box.

    Built by loading the ring's own data file - which already has the right atom
    style, types and masses - adding the extra oxygens through lammps_data, then
    overwriting every coordinate with the centred ones and clearing the image
    flags, since the cluster is contiguous and no longer wrapped.

    Returns (box side actually used, the translation applied).
    """
    data, style = ld.load(ring.source, quiet=True)
    for atom_id, atom_type in zip(ids, types):
        if int(atom_id) not in data.atoms.index:
            ld.add_atom(data, int(atom_type), np.zeros(3), atom_id=int(atom_id))

    side = float(box_side) if box_side else float(data.box.bounds[0][1] - data.box.bounds[0][0])
    shift = np.full(3, side / 2.0) - positions.mean(axis=0)

    for atom_id, position in zip(ids, positions + shift):
        data.atoms.loc[int(atom_id), ["x", "y", "z"]] = position
        data.atoms.loc[int(atom_id), ["nx", "ny", "nz"]] = 0

    data.box = LammpsBox([[0.0, side]] * 3)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ld.write_data_file(data, str(out_path), header, style)
    return side, shift


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--source", required=True,
                        help="the structure the rings were carved from")
    parser.add_argument("--ring-output", required=True,
                        help="the OUTDIR written by silica_rings_to_lammps.py")
    parser.add_argument("--rings", required=True,
                        help="rings.csv from analyze_rings.py, for the path axis")
    parser.add_argument("--outdir", required=True, help="directory to write clusters into")
    parser.add_argument("--ring", action="append", default=None, metavar="RING_ID",
                        help="build only this ring (repeatable); default is all of them")
    parser.add_argument("--size", type=int, action="append", default=None,
                        help="build only rings of this size (repeatable)")
    parser.add_argument("--cutoff", type=float, default=1.9,
                        help="Si-O bond cutoff; must match the ring finder (default: 1.9)")
    parser.add_argument("--tie-break", choices=("axis", "lowest-id"), default="axis",
                        help="which oxygen a silicon contributes when it has more than "
                             "one spare. 'axis' (default) takes the one furthest from the "
                             "path axis, so balancing the charge does not block the "
                             "channel; 'lowest-id' is arbitrary but free of that bias")
    parser.add_argument("--box", type=float, default=None, metavar="SIDE",
                        help="cubic box side in A for the output clusters. Default keeps "
                             "the source box, which is generous for LAMMPS and wasteful "
                             "for a plane-wave code; pass something smaller for VASP")
    parser.add_argument("--path-half-length", type=float, default=7.0,
                        help="how far the probe will travel either side of the ring; the "
                             "box has to hold it (default: 7.0)")
    parser.add_argument("--cutoff-radius", type=float, default=5.5,
                        help="the potential's interaction cutoff, used to check the box is "
                             "big enough that the cluster cannot see its own image "
                             "(default: 5.5, the USC two-body cutoff)")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    ring_dir = Path(args.ring_output)
    outdir = Path(args.outdir)

    structure = read_structure(args.source)
    print(f"source: {len(structure.ids)} atoms, box "
          f"{' x '.join(f'{v:.3f}' for v in structure.lengths)}")
    neighbours = si_o_neighbours(structure, args.cutoff)
    print(f"{len(neighbours)} Si, Si-O cutoff {args.cutoff}")

    with open(args.rings) as handle:
        catalogue = list(csv.DictReader(handle))
    identifiers = [row["ring_id"] for row in catalogue]
    if len(set(identifiers)) != len(identifiers):
        duplicated = [i for i, c in Counter(identifiers).items() if c > 1]
        raise ValueError(
            f"{args.rings}: ring_id is not unique ({len(duplicated)} repeated, e.g. "
            f"{duplicated[:3]}); rebuild it with a current analyze_rings.py")
    if args.ring:
        wanted = set(args.ring)
        catalogue = [r for r in catalogue if r["ring_id"] in wanted]
    if args.size:
        sizes = set(args.size)
        catalogue = [r for r in catalogue if int(r["n"]) in sizes]
    if not catalogue:
        raise ValueError("no rings selected")
    print(f"building {len(catalogue)} cluster(s) into {outdir}/")

    rows, pass_counts, drops = [], Counter(), []
    for count, entry in enumerate(catalogue, start=1):
        ring = read_ring(str(ring_dir / entry["file"]))
        center = np.array([float(entry[f"center_{a}"]) for a in "xyz"])
        axis = np.array([float(entry[f"axis_{a}"]) for a in "xyz"])
        axis = axis / np.linalg.norm(axis)

        added, passes = choose_oxygens(ring, neighbours, structure, center, axis,
                                       args.tie_break)
        pass_counts[passes] += 1
        ids, types, positions = build(ring, added, structure)

        n_si = int((types == SI_TYPE).sum())
        n_o = int((types == O_TYPE).sum())
        charge = n_si * USC_CHARGE[SI_TYPE] + n_o * USC_CHARGE[O_TYPE]
        if n_o != 2 * n_si or abs(charge) > 1e-9:
            raise ValueError(
                f"{entry['ring_id']}: composition Si{n_si}O{n_o} has USC charge "
                f"{charge:+.3f}, expected Si_n O_2n at zero")

        aperture = clearance_along(positions, center, axis)
        flags = contact_violations(positions, types, center, axis)
        ring_aperture = float(entry["aperture"])
        drops.append(ring_aperture - aperture)

        offsets = positions[len(ring.ids):] - center
        perpendicular = offsets - np.outer(offsets @ axis, axis)
        cluster_radius = float(np.linalg.norm(positions - positions.mean(axis=0),
                                              axis=1).max())

        out_path = outdir / entry["file"]
        side, shift = write_cluster(ring, ids, types, positions, out_path,
                                    args.box, f"{entry['n']}-membered ring "
                                    f"{entry['ring_id']} balanced to SiO2, from "
                                    f"{Path(args.source).name}")

        # The box has to hold the cluster and the probe's full travel, with the
        # potential's cutoff to spare, or the cluster interacts with its image.
        needed = 2 * (max(cluster_radius, args.path_half_length) + args.cutoff_radius)
        if side < needed:
            raise ValueError(
                f"{entry['ring_id']}: box side {side:.2f} A is too small; the cluster "
                f"({cluster_radius:.2f} A radius) plus a path reaching "
                f"{args.path_half_length:.2f} A plus a {args.cutoff_radius:.2f} A cutoff "
                f"needs at least {needed:.2f} A")

        rows.append({
            "ring_id": entry["ring_id"], "n": entry["n"], "file": entry["file"],
            "n_atoms": len(ids), "n_si": n_si, "n_o": n_o,
            "charge_sum": round(charge, 12), "passes": passes,
            "aperture_ring": round(ring_aperture, 4),
            "aperture_cluster": round(aperture, 4),
            "aperture_drop": round(ring_aperture - aperture, 4),
            "n_inside_contact_radius": flags["n_inside_contact_radius"],
            "worst_contact_overlap": round(flags["worst_contact_overlap"], 4),
            "added_ids": " ".join(str(i) for i in added),
            "min_added_axis_distance": round(float(np.linalg.norm(perpendicular, axis=1).min()), 4),
            "cluster_radius": round(cluster_radius, 4),
            "box_side": round(side, 4),
            **{f"center_{a}": round(float(v), 6) for a, v in zip("xyz", center + shift)},
            **{f"axis_{a}": round(float(v), 6) for a, v in zip("xyz", axis)},
        })
        if count % 100 == 0 or count == len(catalogue):
            print(f"\r  built {count}/{len(catalogue)}", end="", flush=True)
    print()

    outdir.mkdir(parents=True, exist_ok=True)
    with open(outdir / "cluster_manifest.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=MANIFEST_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)

    drops = np.array(drops)
    print(f"\nall {len(rows)} clusters are Si_n O_2n with USC charge 0")
    print("passes needed to balance: " +
          ", ".join(f"{p}: {c}" for p, c in sorted(pass_counts.items())))
    print(f"aperture lost to the added oxygens: median {np.median(drops):.3f} A, "
          f"90th pct {np.percentile(drops, 90):.3f} A, worst {drops.max():.3f} A")
    flagged = sum(1 for r in rows if r["n_inside_contact_radius"] > 0)
    print(f"{flagged} of {len(rows)} clusters ({100 * flagged / len(rows):.1f}%) have at "
          f"least one atom closer to the axis than its contact radius "
          f"(Si {CONTACT_RADII[SI_TYPE]}, O {CONTACT_RADII[O_TYPE]}); this is a flag on "
          f"the geometry, not an aperture")
    print(f"wrote {len(rows)} clusters and cluster_manifest.csv to {outdir}/")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, OSError) as error:
        sys.exit(f"build_ring_cluster.py: {error}")
