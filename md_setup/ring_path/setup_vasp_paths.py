#!/usr/bin/env python3
"""
setup_vasp_paths.py

Convert the representative rings of a LAMMPS ring-permeation survey into VASP
POSCARs, one per frame of the path, so the same single points can be run in DFT.

The survey stores one data file per case, with the water parked at the start of
its path; in.survey makes every later frame by moving the three water atoms by
(dx, dy, dz) between `run 0` calls. This does the same displacement in python, so
frame f here is exactly frame f there: the start structure from
<survey>/data/<case_id>.data, with the water moved (f - 1) steps, where the step
is read from the case's shard lists rather than recomputed.

Each representative is converted at its best roll, the one its barrier in
representatives.csv came from.

The cell
--------
Not the survey's by default. The clusters kept the 42.8 A box of the amorphous
cell they were carved from, which leaves 31-35 A of vacuum around every frame,
and plane waves fill all of it. --box puts every frame in a cubic cell of that
edge instead (default 25 A, about a fifth of the volume), with the ring centre -
where the oxygen is halfway along the path - at the middle of the cell.

The move is a rigid translation shared by every frame of a path, after unwrapping
each atom around the ring centre, so the structures are unchanged and the cell is
identical along a path, which is what lets frame energies be subtracted. Every
frame is checked to leave at least --min-vacuum between it and its periodic image
along each cell vector, or nothing more is written. --box 0 keeps the survey's
cell and coordinates exactly.

Output
------
  <outdir>/<case_id>/<frame>/POSCAR   frames numbered from 1 as in frames.csv,
                                      species in Si O H order (the POTCAR's)
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import lammps_data as ld  # noqa: E402

# The survey's atom types, fixed by pair_coeff * * SiOH.usc Si O H.
SPECIES = {1: "Si", 2: "O", 3: "H"}


def shard_step(shards: Path, shard: str, case_id: str) -> np.ndarray:
    """The (dx, dy, dz) in.survey moved this case's water by between frames."""
    labels = (shards / shard / "label.txt").read_text().split()
    if labels.count(case_id) != 1:
        raise ValueError(f"shard {shard}: {case_id} appears {labels.count(case_id)} times")
    index = labels.index(case_id)
    return np.array([float((shards / shard / f"{c}.txt").read_text().split()[index])
                     for c in ("dx", "dy", "dz")])


def write_poscar(path: Path, comment: str, box: np.ndarray, types: np.ndarray,
                 positions: np.ndarray) -> None:
    order = [t for t in SPECIES if t in set(types)]
    lines = [comment, "1.0",
             f"  {box[0]:.10f}  0.0000000000  0.0000000000",
             f"  0.0000000000  {box[1]:.10f}  0.0000000000",
             f"  0.0000000000  0.0000000000  {box[2]:.10f}",
             "  " + "  ".join(SPECIES[t] for t in order),
             "  " + "  ".join(str(int((types == t).sum())) for t in order),
             "Cartesian"]
    for t in order:
        lines += [f"  {x:16.10f}  {y:16.10f}  {z:16.10f}" for x, y, z in positions[types == t]]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--survey", required=True,
                        help="the LAMMPS survey's setup directory (<run>/setup/survey)")
    parser.add_argument("--representatives", required=True,
                        help="representatives.csv from collect_survey.py")
    parser.add_argument("--outdir", required=True, help="where the POSCARs go")
    parser.add_argument("--box", type=float, default=25.0,
                        help="cubic cell edge in A, the same for every frame; 0 keeps the "
                             "survey's cell and coordinates (default: 25)")
    parser.add_argument("--min-vacuum", type=float, default=10.0,
                        help="refuse a --box that leaves less than this between a frame "
                             "and its periodic image, in A (default: 10)")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    survey = Path(args.survey)
    outdir = Path(args.outdir)

    with open(survey / "cases.csv") as handle:
        cases = {(row["ring_id"], float(row["roll"])): row for row in csv.DictReader(handle)}
    with open(args.representatives) as handle:
        representatives = list(csv.DictReader(handle))

    written = 0
    for rep in representatives:
        key = (rep["ring_id"], float(rep["best_roll"]))
        if key not in cases:
            raise ValueError(f"no survey case for {rep['ring_id']} at roll {rep['best_roll']}")
        case = cases[key]
        case_id, n_frames = case["case_id"], int(case["n_frames"])

        data, _ = ld.load(str(survey / case["data_file"]), quiet=True)
        if ld.is_triclinic(data):
            raise ValueError(f"{case['data_file']}: triclinic box")
        atoms = data.atoms.sort_index()
        types = atoms["type"].to_numpy(dtype=int)
        if set(types) - set(SPECIES):
            raise ValueError(f"{case['data_file']}: types {sorted(set(types))}, expected "
                             f"{sorted(SPECIES)}")
        start = atoms[["x", "y", "z"]].to_numpy(dtype=float)
        # The same check in.survey makes: the water is the three ids above nfix.
        water = atoms.index.to_numpy() > int(case["nfix"])
        if water.sum() != 3:
            raise ValueError(f"{case_id}: {water.sum()} atoms above nfix, expected 3")

        step = shard_step(survey / "shards", case["shard"], case_id)
        if abs(np.linalg.norm(step) - float(case["step"])) > 1e-4:
            raise ValueError(f"{case_id}: the shard step and cases.csv disagree")

        bounds = np.asarray(data.box.bounds, dtype=float)
        survey_box = bounds[:, 1] - bounds[:, 0]
        half_length, step_length = float(case["half_length"]), float(case["step"])

        if args.box:
            # The oxygen is at the ring centre halfway along the path. Unwrapping
            # around that point undoes any periodic wrap of the survey box, so the
            # cluster and the water are contiguous before they are moved.
            oxygen = start[atoms.index.to_numpy() == int(case["o_id"])][0]
            ring_centre = oxygen + 0.5 * (n_frames - 1) * step
            start = ring_centre + ld.minimum_image(start - ring_centre, survey_box)
            shift = np.full(3, args.box / 2) - ring_centre
            box = np.full(3, args.box)
        else:
            shift, box = -bounds[:, 0], survey_box

        least_vacuum = np.inf
        for frame in range(1, n_frames + 1):
            positions = start.copy()
            positions[water] += (frame - 1) * step
            positions += shift
            vacuum = float((box - (positions.max(axis=0) - positions.min(axis=0))).min())
            if args.box and (vacuum < args.min_vacuum or positions.min() < 0
                             or (positions >= box).any()):
                raise ValueError(
                    f"{case_id} frame {frame}: {vacuum:.2f} A of vacuum in a {args.box:g} A "
                    f"cell (need {args.min_vacuum:g}), or atoms outside it; raise --box")
            least_vacuum = min(least_vacuum, vacuum)
            axis_position = -half_length + step_length * (frame - 1)
            write_poscar(outdir / case_id / str(frame) / "POSCAR",
                         f"{rep['ring_id']} roll {float(case['roll']):g} frame "
                         f"{frame}/{n_frames} axis {axis_position:+.4f} A",
                         box, types, positions)
        written += n_frames
        print(f"  {case_id:<32} {len(types):>3} atoms  {n_frames} frames  "
              f"vacuum >= {least_vacuum:.1f} A")

    print(f"\nwrote {written} POSCARs for {len(representatives)} representatives to {outdir}/")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, OSError) as error:
        sys.exit(f"setup_vasp_paths.py: {error}")
