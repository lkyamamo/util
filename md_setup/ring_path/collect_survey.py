#!/usr/bin/env python3
"""
collect_survey.py

Join the energies LAMMPS printed back onto the cases, turn each path into a
barrier, and pick a representative ring for every size.

The energy of a path is referenced to its own first frame:

    dE(s) = E(s) - E(s_first)

The water is rigid and identical everywhere and the cluster never moves, so the
water's self-energy is an additive constant and the cluster's is constant along a
path; both cancel in that difference and neither has to be computed. The first
frame is far enough out that no water atom is inside the potential's cutoff of
any cluster atom, which setup_survey.py solved for and asserted, so the
subtraction is exact rather than approximate under USC.

    E_barrier = max_s dE(s)

Each ring is run at several rolls, and the roll is a free parameter nothing else
fixes, so the ring's barrier is the *minimum* over its rolls - the easiest way
the molecule can be turned as it goes through.

Choosing the representative
---------------------------
Per size, the ring whose barrier is closest to the median of its size class:
representative in the quantity actually being measured, rather than in a
geometric proxy that may not map onto it. The median, not an extremum - the most
open ring of a size is an outlier by construction and flattens the size trend.

--per-size takes more than one, at evenly spaced percentiles, so the DFT points
come with a spread instead of as singletons.

Output
------
  <outdir>/frames.csv         one row per frame: the axis position, the raw pe
                              and dE, with the ring's size, aperture, radius and
                              pucker on every row, so it can be grouped or
                              plotted by any of them without a further join
  <outdir>/paths.csv          one row per case: the barrier, where it peaks, and
                              how well the far end closed
  <outdir>/barriers.csv       one row per ring: the best roll and its barrier,
                              joined to the ring's geometry
  <outdir>/representatives.csv the chosen rings
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

# Energies are in eV (LAMMPS metal units).
TAIL_TOLERANCE = 0.005   # eV; how close to zero the far end of a path must land


def resolve_energy_files(survey: Path, given: str | None,
                         run_dir: Path | None) -> list[Path]:
    """
    The energy files to read: what was asked for, or the concatenated
    energies.txt, or failing that the per-shard energies_*.txt directly.

    Looked for in the run directory as well as the survey directory, because
    under this project's layout LAMMPS runs in <case>/run while the cases and
    structures live in <case>/setup/survey, so the energies land beside the
    submit script rather than beside cases.csv.

    Reading the shards directly means the `cat` step in the submit script is a
    convenience rather than a requirement, and a survey where one shard is still
    running can be collected for what has finished.
    """
    if given:
        return sorted(Path().glob(given)) or [Path(given)]
    looked = [d for d in (run_dir, survey) if d is not None]
    for directory in looked:
        if (directory / "energies.txt").exists():
            return [directory / "energies.txt"]
    for directory in looked:
        shards = sorted(directory.glob("energies_*.txt"))
        if shards:
            return shards
    where = " or ".join(str(d) for d in looked)
    raise ValueError(
        f"no energies found: looked for energies.txt and energies_*.txt in {where}")


def read_energies(paths: list[Path]) -> dict[str, dict[int, float]]:
    """
    Parse what in.survey printed: "<label> <frame> <pe>" per line.

    Every line carries its own label, so the file is a tidy table and the order
    it was concatenated in does not matter. Later lines win, so re-running a
    shard over an existing file is harmless rather than silently mixing two runs
    together.
    """
    energies: dict[str, dict[int, float]] = defaultdict(dict)
    for path in paths:
        for line in path.read_text().splitlines():
            words = line.split()
            if len(words) != 3:
                continue
            label, frame, value = words
            try:
                energies[label][int(frame)] = float(value)
            except ValueError:
                continue
    return energies


def summarise_path(case: dict, frames: dict[int, float]) -> dict:
    """Turn one case's energies into a barrier, or say why it cannot be used."""
    # Everything arrives from the CSV as a string; coerce once here so that no
    # downstream comparison or format silently works on text.
    case = {**case,
            "n": int(case["n"]),
            "roll": float(case["roll"]),
            "n_frames": int(case["n_frames"]),
            "half_length": float(case["half_length"]),
            "step": float(case["step"]),
            "min_separation": float(case["min_separation"]),
            "reference_separation": float(case["reference_separation"])}

    n_frames = case["n_frames"]
    missing = [f for f in range(1, n_frames + 1) if f not in frames]
    if missing:
        return {**case, "status": f"missing {len(missing)} of {n_frames} frames",
                "e_barrier": "", "barrier_position": "", "tail": ""}

    half = case["half_length"]
    step = case["step"]
    positions = np.array([-half + step * (f - 1) for f in range(1, n_frames + 1)])
    energies = np.array([frames[f] for f in range(1, n_frames + 1)])
    delta = energies - energies[0]

    peak = int(np.argmax(delta))
    return {
        **case,
        "status": "ok",
        "e_barrier": round(float(delta[peak]), 6),
        "barrier_position": round(float(positions[peak]), 4),
        # The entry end is zero by construction, so the exit end is the real
        # test that the path started and finished outside the interaction.
        "tail": round(float(delta[-1]), 6),
        # Kept so the per-frame curve can be written out; stripped before the
        # summary CSVs, which are one row per path.
        "_positions": positions,
        "_energies": energies,
        "_delta": delta,
    }


def choose(rows: list[dict], per_size: int) -> list[dict]:
    """Per size, the ring(s) nearest the class median barrier."""
    by_size = defaultdict(list)
    for row in rows:
        by_size[int(row["n"])].append(row)

    chosen = []
    for n in sorted(by_size):
        group = sorted(by_size[n], key=lambda r: r["e_barrier"])
        if not group:
            continue
        # Evenly spaced percentiles centred on the median: one pick is the
        # median, three are the quartiles, and so on.
        if per_size == 1:
            fractions = [0.5]
        else:
            fractions = list(np.linspace(0.25, 0.75, per_size))
        for fraction in fractions:
            index = int(round(fraction * (len(group) - 1)))
            row = dict(group[index])
            row["percentile"] = round(100 * index / max(len(group) - 1, 1), 1)
            chosen.append(row)
    return chosen


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--survey", required=True, help="the setup_survey.py output directory")
    parser.add_argument("--energies", default=None,
                        help="the concatenated energies file (default: <survey>/energies.txt)")
    parser.add_argument("--rings", default=None,
                        help="rings.csv, to join the ring geometry onto the barriers")
    parser.add_argument("--run-dir", default=None,
                        help="the directory LAMMPS ran in (<case>/run). Searched for the "
                             "energies, and used as the default place to write, since that "
                             "is where this project keeps a run's outputs")
    parser.add_argument("--outdir", default=None,
                        help="where to write (default: --run-dir if given, else <survey>)")
    parser.add_argument("--per-size", type=int, default=1,
                        help="how many representatives per size (default: 1)")
    parser.add_argument("--no-frames", action="store_true",
                        help="skip frames.csv, the per-frame dE(s) table. It is one row "
                             "per frame rather than per path, so it is the big one")
    parser.add_argument("--tail-tolerance", type=float, default=TAIL_TOLERANCE,
                        help=f"how close to zero the far end of a path must land, in eV "
                             f"(default: {TAIL_TOLERANCE})")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    survey = Path(args.survey)
    run_dir = Path(args.run_dir) if args.run_dir else None
    outdir = Path(args.outdir) if args.outdir else (run_dir or survey)
    outdir.mkdir(parents=True, exist_ok=True)
    energy_files = resolve_energy_files(survey, args.energies, run_dir)
    print("reading " + (f"{len(energy_files)} shard files" if len(energy_files) > 1
                        else str(energy_files[0])))

    with open(survey / "cases.csv") as handle:
        cases = list(csv.DictReader(handle))
    energies = read_energies(energy_files)
    print(f"{len(cases)} cases, {len(energies)} of them with energies")

    paths = [summarise_path(case, energies.get(case["case_id"], {})) for case in cases]
    incomplete = [p for p in paths if p["status"] != "ok"]
    if incomplete:
        print(f"WARNING: {len(incomplete)} case(s) are incomplete, e.g. "
              f"{incomplete[0]['case_id']}: {incomplete[0]['status']}")

    usable = [p for p in paths if p["status"] == "ok"]
    if not usable:
        raise ValueError("no complete paths; check the energies file and the shard logs")

    # Everything below is one row per path or per ring, so drop the per-frame
    # arrays that summarise_path carried along for frames.csv.

    # The far end of every path should come back to the reference. A tail that
    # does not close means the path was too short, or something is on the axis
    # beyond the ring, and the barrier is measured against the wrong zero.
    tails = np.array([p["tail"] for p in usable])
    open_tails = int((np.abs(tails) > args.tail_tolerance).sum())
    print(f"far-end tail: median {np.median(np.abs(tails)):.2e} eV, "
          f"{open_tails} of {len(usable)} above {args.tail_tolerance} eV")

    geometry = {}
    if args.rings:
        with open(args.rings) as handle:
            geometry = {row["ring_id"]: row for row in csv.DictReader(handle)}

    carry = ["aperture", "aperture_centroid", "aperture_slid", "planarity_rmsd",
             "eccentricity", "mean_radius", "mean_radius_inplane", "smaller_ring_overlap",
             "undercoordinated_si", "n_inside_contact_radius", "nearest_atom_type"]

    fields = [k for k in usable[0] if not k.startswith("_")]
    with open(outdir / "paths.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(paths)

    # The per-frame curve, in long form with the ring's geometry on every row, so
    # it can be grouped or plotted by size, aperture, radius or pucker without a
    # further join. This is the only place dE(s) is written out; the other two
    # files are one row per path and one row per ring.
    if not args.no_frames:
        frame_fields = (["ring_id", "n", "roll", "frame", "axis_position", "pe", "delta_e"]
                        + (carry if geometry else []))
        with open(outdir / "frames.csv", "w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=frame_fields, extrasaction="ignore")
            writer.writeheader()
            for path in usable:
                shared = {"ring_id": path["ring_id"], "n": path["n"], "roll": path["roll"]}
                shared.update({k: geometry.get(path["ring_id"], {}).get(k, "")
                               for k in (carry if geometry else [])})
                for index, (position, energy, delta) in enumerate(
                        zip(path["_positions"], path["_energies"], path["_delta"]), start=1):
                    writer.writerow({**shared, "frame": index,
                                     "axis_position": round(float(position), 4),
                                     "pe": round(float(energy), 6),
                                     "delta_e": round(float(delta), 6)})
        written = sum(len(p["_delta"]) for p in usable)
        print(f"frames.csv: {written} rows, one per frame"
              + (", with the ring geometry joined on" if geometry else
                 " (pass --rings to join the geometry on)"))

    # One row per ring: the roll that makes the crossing easiest.
    by_ring: dict[str, list[dict]] = defaultdict(list)
    for path in usable:
        by_ring[path["ring_id"]].append(path)

    barriers = []
    for ring_id, group in by_ring.items():
        best = min(group, key=lambda p: p["e_barrier"])
        row = {
            "ring_id": ring_id, "n": best["n"],
            "best_roll": best["roll"], "e_barrier": best["e_barrier"],
            "barrier_position": best["barrier_position"],
            "n_rolls": len(group),
            "roll_spread": round(max(p["e_barrier"] for p in group)
                                 - min(p["e_barrier"] for p in group), 6),
            "half_length": best["half_length"],
            "min_separation": best["min_separation"],
            "worst_tail": round(max(abs(p["tail"]) for p in group), 6),
        }
        row.update({k: geometry.get(ring_id, {}).get(k, "") for k in carry})
        barriers.append(row)
    barriers.sort(key=lambda r: (r["n"], r["e_barrier"]))

    with open(outdir / "barriers.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(barriers[0].keys()))
        writer.writeheader()
        writer.writerows(barriers)

    print(f"\n{'n':>2} {'N':>5} | {'E_barrier (eV)':^30} | {'roll spread':>11} | "
          f"{'peak at':>8}")
    print(f"{'':>2} {'':>5} | {'p10':>6} {'med':>6} {'p90':>6} {'min':>6} | {'median':>11} | "
          f"{'median':>8}")
    print("-" * 74)
    by_size = defaultdict(list)
    for row in barriers:
        by_size[row["n"]].append(row)
    for n in sorted(by_size):
        e = np.array([r["e_barrier"] for r in by_size[n]])
        spread = np.array([r["roll_spread"] for r in by_size[n]])
        peak = np.array([abs(r["barrier_position"]) for r in by_size[n]])
        print(f"{n:2d} {len(e):5d} | {np.percentile(e, 10):6.2f} {np.median(e):6.2f} "
              f"{np.percentile(e, 90):6.2f} {e.min():6.2f} | {np.median(spread):9.2f} eV | "
              f"{np.median(peak):6.2f} A")

    if geometry:
        energy = np.array([r["e_barrier"] for r in barriers])
        sizes = np.array([r["n"] for r in barriers], dtype=float)
        apertures = np.array([float(r["aperture"]) if r["aperture"] != "" else np.nan
                              for r in barriers])

        def correlate(x, y):
            """r, or None when one of the two does not vary and r is undefined."""
            good = ~(np.isnan(x) | np.isnan(y))
            if good.sum() < 3 or x[good].std() == 0 or y[good].std() == 0:
                return None
            return float(np.corrcoef(x[good], y[good])[0, 1])

        print("\ncorrelation of the barrier with")
        for name, values, note in [("ring size n", sizes, ""),
                                   ("aperture   ", apertures,
                                    "   <- expect this to be the tighter one")]:
            r = correlate(values, energy)
            print(f"  {name} : " + ("undefined (no variation in this subset)"
                                    if r is None else f"r = {r:+.3f}{note}"))

    representatives = choose(barriers, args.per_size)
    with open(outdir / "representatives.csv", "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(representatives[0].keys()))
        writer.writeheader()
        writer.writerows(representatives)

    print(f"\nrepresentatives ({args.per_size} per size, by barrier percentile):")
    print(f"  {'n':>2}  {'ring_id':<22} {'pct':>5} {'E_barrier':>10} {'roll':>5} "
          f"{'aperture':>9}")
    for row in representatives:
        print(f"  {row['n']:2d}  {row['ring_id']:<22} {row['percentile']:>5.1f} "
              f"{row['e_barrier']:>10.3f} {row['best_roll']:>5g} "
              f"{str(row['aperture'])[:8]:>9}")

    print(f"\nwrote paths.csv, barriers.csv and representatives.csv to {outdir}/")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, OSError) as error:
        sys.exit(f"collect_survey.py: {error}")
