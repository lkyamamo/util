#!/usr/bin/env python3
"""
calibrate_cutoffs.py — SKETCH, NOT VALIDATED, NOT WIRED INTO ANYTHING.

Parked here for possible future development. bond_events.py does not import or
call this; its cutoffs are explicit CLI arguments with stated defaults. Read the
STATUS section before trusting any number this prints.

WHAT IT IS FOR
--------------
The O-H "bond" cutoff in this repo is 1.1, 1.15, 1.2, 1.4, 2.3 or 2.47 A
depending on which script you open:

    md_setup/fix_broken_bonds/fix_broken_edge_bonds.py:29   1.1
    analysis/distributions/lammps_GrBa_ovito/1.ovito_GrBa.py:16   1.15
    analysis/structure/silanol.py:162                       1.2
    analysis/shock/voxel_analysis.py:36                     1.2
    visualization/ovito_rotation_config.py:55               1.2
    analysis/distributions/.../bad_freud.py:128             1.4
    analysis/distributions/water_kenichi_GrBa/stat/stat.f90:34   2.3
    analysis/distributions/kenichi_GrNrB/stat/stat.f90      2.47

None of them records where the number came from. The idea here is to stop
adding unsourced constants and instead derive the cutoff from the first minimum
of the O-H g(r) of the very trajectory being analysed — the same g(r)
rdf_freud.py already computes.

STATUS — why this is parked rather than shipped
-----------------------------------------------
The minimum-finding is the hard part and is not trustworthy yet:

  * Taking the first local minimum after the main peak picks up any single-bin
    noise dip on the peak's flank, giving a cutoff far below the real one.
  * Smoothing plus a prominence threshold helps, but a genuinely structured
    bond-length distribution still produces a dip *inside* the first shell that
    is mistaken for the shell edge. SHELL_FLOOR_FRACTION below is an attempt to
    exclude that by refusing to look until g(r) has fallen to 25% of the peak,
    but 25% is itself an unsourced constant — exactly the disease being cured.
  * It has only ever been exercised on synthetic trajectories. It has never
    been run against a production dynamics.lammpstrj, and the suggested
    +/-10% bracket around the minimum has no justification beyond being round.

TO PROMOTE THIS
---------------
  1. Run it on several real production trajectories at different temperatures
     and check the minimum lands where rdf_freud.py's published g(r) says it
     should.
  2. Justify or replace SHELL_FLOOR_FRACTION and the +/-10% bracket. The
     bracket in particular should probably come from the width of the g(r)
     minimum, not from a fixed percentage.
  3. Check the answer is stable against --calibrate-frames and against the
     bin count.
  4. Only then consider wiring it into bond_events.py.

USAGE (if you are developing it)
--------------------------------
    python calibrate_cutoffs.py --traj dynamics.lammpstrj --frames-used 20
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import freud

# bond_events.py owns the trajectory reader; reuse it rather than forking a
# fourth copy of the LAMMPS dump parser.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from bond_events import _freud_box_and_points, _parse_type_map, read_frames  # noqa: E402


#: A candidate first minimum must lie where g(r) has already fallen to this
#: fraction of the peak, so a dip within the covalent peak itself is not
#: mistaken for the shell edge. UNSOURCED — see STATUS above.
SHELL_FLOOR_FRACTION = 0.25

#: Fractional bracket placed either side of the minimum to get r_form/r_break.
#: UNSOURCED — see STATUS above.
BRACKET = 0.10


def _smooth(y: np.ndarray, window: int = 5) -> np.ndarray:
    """Moving average, used only to keep minimum-finding off single-bin noise."""
    if window < 2 or len(y) < window:
        return y
    return np.convolve(y, np.ones(window) / window, mode="same")


def first_minimum(r: np.ndarray, g: np.ndarray) -> Tuple[int, Optional[int], str]:
    """
    Locate the boundary of the first coordination shell in g(r).

    Returns (peak index, minimum index or None, how it was found). The search
    starts only once g has dropped to SHELL_FLOOR_FRACTION of the peak; from
    there the first prominent local minimum is taken, falling back to the global
    minimum of the gap region.
    """
    smooth = _smooth(g)
    peak = int(np.argmax(smooth))
    if smooth[peak] <= 0:
        return peak, None, "g(r) is empty"

    below = np.flatnonzero(smooth[peak:] < SHELL_FLOOR_FRACTION * smooth[peak])
    if len(below) == 0:
        return peak, None, (
            f"g(r) never falls below {SHELL_FLOOR_FRACTION:.0%} of its peak within "
            f"r_max, so the first shell has no visible edge"
        )

    start = peak + int(below[0])
    tail = smooth[start:]
    if len(tail) < 3:
        return peak, start, "shell edge at the end of the sampled range"

    try:
        from scipy.signal import find_peaks

        minima, _ = find_peaks(-tail, prominence=0.02 * float(smooth[peak]))
        if len(minima):
            return peak, start + int(minima[0]), "prominent local minimum"
    except ImportError:
        pass

    return peak, start + int(np.argmin(tail)), "global minimum of the gap region"


def second_peak(r: np.ndarray, g: np.ndarray, after: int) -> Optional[Tuple[float, float]]:
    """The hydrogen-bonded O-H shell, if there is one. Its presence is what makes
    the first minimum a real boundary between two populations rather than an
    arbitrary point on a decaying tail."""
    smooth = _smooth(g)
    tail = smooth[after:]
    if len(tail) < 3:
        return None
    idx = after + int(np.argmax(tail))
    if smooth[idx] <= smooth[after] * 1.5:
        return None
    return float(r[idx]), float(g[idx])


def compute_oh_rdf(args: argparse.Namespace) -> Tuple[np.ndarray, np.ndarray, int]:
    rdf = freud.density.RDF(bins=args.bins, r_max=args.r_max, r_min=0.4)
    n_used = 0
    for frame, topology in read_frames(args.traj, type_map=_parse_type_map(args.type_map)):
        box, points = _freud_box_and_points(frame)
        rdf.compute(
            system=(box, points[topology.h_rows]),
            query_points=points[topology.o_rows],
            reset=False,
        )
        n_used += 1
        if n_used >= args.frames_used:
            break
    return rdf.bin_centers, rdf.rdf, n_used


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="EXPERIMENTAL sketch: suggest O-H cutoffs from a trajectory's own g(r).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--traj", type=Path, required=True)
    p.add_argument("--frames-used", type=int, default=20)
    p.add_argument("--bins", type=int, default=400)
    p.add_argument("--r-max", type=float, default=2.5)
    p.add_argument("--type-map", default=None,
                   help="e.g. '1:O;2:H' for dumps with no element column")
    args = p.parse_args(argv)

    print("*** EXPERIMENTAL — see the STATUS section of this file before using ***\n")

    r, g, n_used = compute_oh_rdf(args)
    peak, minimum, note = first_minimum(r, g)

    print(f"O-H g(r) from {n_used} frame(s) of {args.traj}")
    print(f"  first peak     r = {r[peak]:.3f} A   g = {g[peak]:.2f}   (covalent O-H)")
    if minimum is None:
        print(f"  first minimum  not found — {note}")
        print("\nNo cutoff can be suggested from this g(r).")
        return 1

    r_min, g_min = float(r[minimum]), float(g[minimum])
    print(f"  first minimum  r = {r_min:.3f} A   g = {g_min:.3f}   [{note}]")

    second = second_peak(r, g, minimum)
    if second is not None:
        print(f"  second peak    r = {second[0]:.3f} A   g = {second[1]:.2f}   (hydrogen bonded)")
    else:
        print("  second peak    not found — see the warning below")

    print(f"\nSuggested cutoffs (bracketing the minimum by +/-{BRACKET:.0%}):")
    print(f"  --r-form  {r_min * (1 - BRACKET):.2f}")
    print(f"  --r-break {r_min * (1 + BRACKET):.2f}")

    if second is None:
        print("\nWARNING: no hydrogen-bonded second peak was found out to "
              f"r = {args.r_max} A, so the minimum above is just where g(r) happens to be")
        print("lowest rather than a boundary between two populations.")
    elif g_min > 0.1 * g[peak]:
        print("\nWARNING: the minimum is shallow relative to the peak, so the covalent and")
        print("hydrogen-bonded populations overlap. No distance criterion separates them")
        print("cleanly.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
