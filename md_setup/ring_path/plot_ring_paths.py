#!/usr/bin/env python3
"""
plot_ring_paths.py

Figures for the ring-permeation survey, from the CSVs collect_survey.py writes.

Three things get drawn:

  size_NN_paths.png  dE against distance from the ring centre, every ring of
                     that size overlaid, with the median over them picked out.
                     One figure per ring size.
  surface_3d.png     ring radius and distance from the centre on the x and y
                     axes, energy on z. Binned by radius and taken as a median,
                     because the rings do not share a radius.
  surface_heatmap.png the same array with magnitude as colour - the one to read
                     values off, since a surface hides its own far side.
  paths/<ring>.png   dE against position along the axis, for one ring. All six
                     rolls are on the plot, but only the best one is drawn as
                     data - the rest are recessive context, because the question
                     the figure answers is "how hard is the easiest crossing",
                     not "what did roll 90 do".
  summary.png        the barrier distribution per ring size, beside the barrier
                     against aperture. Two panels, because they are two
                     different measures - never one plot with two y scales.
  by_size.png        barrier against aperture again, faceted one panel per ring
                     size. This is the figure that shows whether aperture still
                     predicts the barrier *within* a size class, which is the
                     claim that ring size is the weaker variable.

Why symlog
----------
Nothing is relaxed, so a water pushed into a sealed ring is hugely repulsive:
the single-point runs in this project reach +45,000 eV while the interesting
structure is a few eV. A linear axis shows one spike and six flat lines. Every
energy axis here is symlog with a linear region around zero, and the path figure
carries a second linear panel zoomed on that region - the same treatment
plot_sp_pe.py uses for the same reason.

Colour
------
One accent (the project's #1f77e0) for the series being read, one recessive grey
for context. Ring size is never encoded in colour: it is ordinal with six levels
and the blue ramp cannot give six steps that stay apart, so size is carried by
position and by faceting instead.
"""

from __future__ import annotations

import argparse
import csv
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

# Ink and chrome, matching the other plotting scripts in this project.
ACCENT = "#1f77e0"      # the series being read
CONTEXT = "#b0b0b0"     # rolls that are not the best one: present, not read
TEXT = "#444444"
MUTED = "#888888"
GRID = "#dddddd"
SPINE = "#bbbbbb"

# The project's sequential blue ramp, light to dark, for continuous magnitude.
# One hue - never a rainbow - so the reader gets "more" from the darkness alone.
BLUE_RAMP = ["#cde2fb", "#b7d3f6", "#9ec5f4", "#86b6ef", "#6da7ec", "#5598e7",
             "#3987e5", "#2a78d6", "#256abf", "#1c5cab", "#184f95", "#104281",
             "#0d366b"]
BLUE_CMAP = LinearSegmentedColormap.from_list("ring_blue", BLUE_RAMP)

# "Ring radius" means the mean distance of the 2n ring atoms from the centroid,
# measured in the ring's own best-fit plane. The plain 3-D distance
# ("mean_radius") runs larger by an amount that grows with ring size - 0.005 A at
# n=3, 0.085 A at n=8 - which would put pucker on an axis meant to carry size.
# Silicon and oxygen are averaged together; they sit 0.26-0.38 A apart, so this
# is between the two shells rather than on either.
RADIUS_KEY = "mean_radius_inplane"

LINTHRESH = 1.0         # eV; half-width of the linear region of the symlog axis
ZOOM_EV = 5.0           # eV; half-width of the linear zoom panel


def style_axes(ax):
    ax.grid(axis="y", color=GRID, lw=1)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(SPINE)
    ax.tick_params(colors=TEXT, labelsize=11)


def read_csv(path: Path) -> list[dict]:
    with open(path) as handle:
        return list(csv.DictReader(handle))


def load_frames(path: Path) -> dict[str, dict[float, tuple[np.ndarray, np.ndarray]]]:
    """ring_id -> roll -> (axis positions, dE), both sorted by position."""
    buckets: dict[str, dict[float, list[tuple[float, float]]]] = defaultdict(
        lambda: defaultdict(list))
    for row in read_csv(path):
        buckets[row["ring_id"]][float(row["roll"])].append(
            (float(row["axis_position"]), float(row["delta_e"])))
    out = {}
    for ring, rolls in buckets.items():
        out[ring] = {}
        for roll, points in rolls.items():
            points.sort()
            out[ring][roll] = (np.array([p for p, _ in points]),
                               np.array([e for _, e in points]))
    return out


def plot_path(ring_id: str, rolls: dict, meta: dict, out_png: Path) -> None:
    """One ring: the full symlog scan beside a linear zoom on the low-energy region."""
    best_roll = float(meta["best_roll"])
    barrier = float(meta["e_barrier"])
    peak_at = float(meta["barrier_position"])

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 5.5), dpi=120)

    for panel, linear in ((ax, False), (ax2, True)):
        for roll, (positions, delta) in sorted(rolls.items()):
            if roll == best_roll:
                continue
            panel.plot(positions, delta, "-", color=CONTEXT, lw=1.4, zorder=2)
        positions, delta = rolls[best_roll]
        panel.plot(positions, delta, "-o", color=ACCENT, lw=2.5, ms=6, zorder=3)
        panel.axhline(0.0, color=SPINE, lw=1, zorder=1)
        panel.axvline(0.0, color=GRID, lw=1, ls=":", zorder=1)
        style_axes(panel)
        panel.set_xlabel("position along the ring axis (Å)", fontsize=13, color=TEXT)

    ax.set_yscale("symlog", linthresh=LINTHRESH)
    ax.set_ylabel(f"$\\Delta E$ vs the first frame (eV, symlog ±{LINTHRESH:g})",
                  fontsize=13, color=TEXT)
    ax.set_title("full scan, every roll", fontsize=12, color=TEXT)
    ax.plot([peak_at], [barrier], "o", ms=16, mfc="none", mec=ACCENT, mew=2.5, zorder=4)
    ax.annotate(f"barrier {barrier:,.1f} eV\n@ {peak_at:+.2f} Å, roll {best_roll:g}°",
                xy=(peak_at, barrier), xytext=(8, -4), textcoords="offset points",
                fontsize=11, color=TEXT)

    ax2.set_ylim(-ZOOM_EV, ZOOM_EV)
    ax2.set_ylabel("$\\Delta E$ vs the first frame (eV)", fontsize=13, color=TEXT)
    ax2.set_title(f"linear zoom, |$\\Delta E$| ≤ {ZOOM_EV:g} eV", fontsize=12, color=TEXT)
    inside = int(np.sum(np.abs(rolls[best_roll][1]) <= ZOOM_EV))
    ax2.annotate(f"{inside} of {len(rolls[best_roll][1])} frames of the best roll "
                 f"are in view",
                 xy=(0.97, 0.05), xycoords="axes fraction", ha="right", va="bottom",
                 fontsize=10, color=MUTED)

    # Two series on the plot, so both are named; the accent one is also direct-
    # labelled above, which is the relief the recessive grey needs.
    ax.plot([], [], "-", color=CONTEXT, lw=1.4, label="the other rolls")
    ax.plot([], [], "-o", color=ACCENT, lw=2.5, ms=6,
            label=f"best roll ({best_roll:g}°)")
    ax.legend(frameon=False, fontsize=11, labelcolor=TEXT, loc="upper left")

    aperture = meta.get("aperture", "")
    detail = f"aperture {float(aperture):.2f} Å" if aperture else ""
    fig.suptitle(f"{ring_id} — water walked through the ring, rigid single points\n"
                 f"{meta['n']}-membered ring{', ' + detail if detail else ''}",
                 fontsize=15, color="#222222")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def plot_summary(barriers: list[dict], out_png: Path) -> None:
    """Barrier by ring size, beside barrier against aperture."""
    by_size = defaultdict(list)
    for row in barriers:
        by_size[int(row["n"])].append(float(row["e_barrier"]))
    sizes = sorted(by_size)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12.5, 5.5), dpi=120)

    parts = ax.boxplot([by_size[n] for n in sizes], positions=range(len(sizes)),
                       widths=0.6, showfliers=False, patch_artist=True,
                       medianprops=dict(color=ACCENT, lw=2.5),
                       boxprops=dict(facecolor="#eef4fd", edgecolor=SPINE, lw=1.2),
                       whiskerprops=dict(color=SPINE, lw=1.2),
                       capprops=dict(color=SPINE, lw=1.2))
    del parts
    for index, n in enumerate(sizes):
        ax.annotate(f"n={len(by_size[n])}", xy=(index, 0), xycoords=("data", "axes fraction"),
                    xytext=(0, 4), textcoords="offset points", ha="center",
                    fontsize=9, color=MUTED)
    ax.set_xticks(range(len(sizes)), [str(n) for n in sizes])
    ax.set_yscale("symlog", linthresh=LINTHRESH)
    ax.set_xlabel("ring size n (silicons)", fontsize=13, color=TEXT)
    ax.set_ylabel("barrier (eV, symlog)", fontsize=13, color=TEXT)
    ax.set_title("barrier by ring size", fontsize=12, color=TEXT)
    style_axes(ax)

    aperture, energy = [], []
    for row in barriers:
        if row.get("aperture"):
            aperture.append(float(row["aperture"]))
            energy.append(float(row["e_barrier"]))
    if aperture:
        ax2.plot(aperture, energy, "o", color=ACCENT, ms=4, alpha=0.45,
                 mew=0, zorder=3)
        ax2.set_yscale("symlog", linthresh=LINTHRESH)
        ax2.set_xlabel("aperture (Å)", fontsize=13, color=TEXT)
        ax2.set_ylabel("barrier (eV, symlog)", fontsize=13, color=TEXT)
        ax2.set_title("barrier against the opening", fontsize=12, color=TEXT)
        r_ap = np.corrcoef(aperture, np.log10(np.clip(energy, 1e-3, None)))[0, 1]
        r_n = np.corrcoef([float(r["n"]) for r in barriers],
                          np.log10(np.clip([float(r["e_barrier"]) for r in barriers],
                                           1e-3, None)))[0, 1]
        ax2.annotate(f"r(log barrier, aperture) = {r_ap:+.2f}\n"
                     f"r(log barrier, n)        = {r_n:+.2f}",
                     xy=(0.97, 0.95), xycoords="axes fraction", ha="right", va="top",
                     fontsize=11, color=TEXT, family="monospace")
        style_axes(ax2)

    fig.suptitle("ring permeation barriers — rigid single points, best roll per ring",
                 fontsize=15, color="#222222")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def plot_by_size(barriers: list[dict], out_png: Path) -> None:
    """Barrier against aperture, one panel per ring size, shared axes."""
    by_size = defaultdict(list)
    for row in barriers:
        if row.get("aperture"):
            by_size[int(row["n"])].append((float(row["aperture"]),
                                           float(row["e_barrier"])))
    sizes = sorted(by_size)
    if not sizes:
        return

    fig, axes = plt.subplots(1, len(sizes), figsize=(3.1 * len(sizes), 4.2), dpi=120,
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    for panel, n in zip(axes, sizes):
        x = [a for a, _ in by_size[n]]
        y = [e for _, e in by_size[n]]
        panel.plot(x, y, "o", color=ACCENT, ms=4, alpha=0.5, mew=0, zorder=3)
        panel.set_yscale("symlog", linthresh=LINTHRESH)
        panel.set_title(f"n = {n}   ({len(x)} rings)", fontsize=12, color=TEXT)
        panel.set_xlabel("aperture (Å)", fontsize=12, color=TEXT)
        style_axes(panel)
    axes[0].set_ylabel("barrier (eV, symlog)", fontsize=12, color=TEXT)

    fig.suptitle("does the opening still predict the barrier within one ring size?",
                 fontsize=15, color="#222222")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def symlog(values, linthresh: float = LINTHRESH):
    """
    Signed log for the surface's z axis: linear-ish inside +/-linthresh, log
    beyond it. A plain log cannot show the shallow negative wells, and a linear
    axis cannot show anything else once one ring reaches a keV.
    """
    values = np.asarray(values, dtype=float)
    return np.sign(values) * np.log10(1.0 + np.abs(values) / linthresh)


def symlog_ticks(lo: float, hi: float):
    """Tick positions in transformed space, labelled with real eV values."""
    decades = [0.0]
    step = linthresh = LINTHRESH
    while step <= max(abs(lo), abs(hi)) * 10:
        decades += [step, -step]
        step *= 10
    del linthresh
    marks = sorted(v for v in decades if lo <= symlog(v) <= hi)
    return [symlog(v) for v in marks], [f"{v:,.0f}" if abs(v) >= 1 else "0" for v in marks]


def best_roll_curves(frames: dict, barriers: list[dict]):
    """
    ring_id -> (positions, dE) for the best roll only, with the ring's radius.

    The aggregate figures are about rings, not rolls, and the ring's barrier is
    defined as the easiest roll, so the other five would just blur the picture.
    """
    curves = {}
    for row in barriers:
        ring = row["ring_id"]
        if ring not in frames or not row.get(RADIUS_KEY):
            continue
        roll = float(row["best_roll"])
        if roll not in frames[ring]:
            continue
        positions, delta = frames[ring][roll]
        curves[ring] = (positions, delta, float(row[RADIUS_KEY]), int(row["n"]))
    return curves


def common_grid(curves, step: float = 0.25):
    """
    Put every ring's curve on one position grid so they can be compared and
    summarised. Outside a ring's own path the water is beyond the cutoff and the
    interaction is zero by construction, so padding with zero is the physics
    rather than a convenience.
    """
    reach = max(float(np.abs(p).max()) for p, _, _, _ in curves.values())
    grid = np.arange(-reach, reach + step / 2, step)
    stacked = {}
    for ring, (positions, delta, radius, n) in curves.items():
        stacked[ring] = (np.interp(grid, positions, delta, left=0.0, right=0.0),
                         radius, n)
    return grid, stacked


def plot_size_overlays(curves, outdir: Path) -> list[Path]:
    """One figure per ring size: every ring of that size, plus the median curve."""
    grid, stacked = common_grid(curves)
    by_size = defaultdict(list)
    for ring, (delta, _, n) in stacked.items():
        by_size[n].append(delta)

    written = []
    for n in sorted(by_size):
        block = np.array(by_size[n])
        fig, ax = plt.subplots(figsize=(8.5, 5.5), dpi=120)
        for row in block:
            ax.plot(grid, row, "-", color=CONTEXT, lw=0.7, alpha=0.35, zorder=2)
        median = np.median(block, axis=0)
        ax.plot(grid, median, "-", color=ACCENT, lw=2.5, zorder=4)
        ax.axhline(0.0, color=SPINE, lw=1, zorder=1)
        ax.axvline(0.0, color=GRID, lw=1, ls=":", zorder=1)
        ax.set_yscale("symlog", linthresh=LINTHRESH)
        ax.set_xlabel("distance from the ring centre along the axis (Å)",
                      fontsize=13, color=TEXT)
        ax.set_ylabel(f"$\\Delta E$ vs the first frame (eV, symlog ±{LINTHRESH:g})",
                      fontsize=13, color=TEXT)
        peak = float(np.max(median))
        ax.annotate(f"median peak {peak:,.1f} eV @ {grid[int(np.argmax(median))]:+.2f} Å",
                    xy=(0.97, 0.95), xycoords="axes fraction", ha="right", va="top",
                    fontsize=11, color=TEXT)
        ax.plot([], [], "-", color=CONTEXT, lw=1.4, label=f"each of the {len(block)} rings")
        ax.plot([], [], "-", color=ACCENT, lw=2.5, label="median over the rings")
        ax.legend(frameon=False, fontsize=11, labelcolor=TEXT, loc="upper left")
        style_axes(ax)
        fig.suptitle(f"{n}-membered rings — water through the ring, best roll each",
                     fontsize=15, color="#222222")
        fig.tight_layout()
        out = outdir / f"size_{n:02d}_paths.png"
        fig.savefig(out)
        plt.close(fig)
        written.append(out)
    return written


def bin_surface(curves, n_radius: int = 26, step: float = 0.25,
                min_count: int = 3):
    """
    Median dE on a (ring radius) x (distance from centre) grid.

    Binned rather than scattered: the rings do not share a radius, so the raw
    cloud is 1,722 unevenly spaced curves. The median in each cell is what makes
    it a surface, and it is a median so one sealed ring cannot drag a whole cell
    into the keV. Bins holding fewer than `min_count` rings are dropped rather
    than drawn: a "median" of one ring is that ring, and at the sparse ends of
    the radius range that would read as structure.
    """
    grid, stacked = common_grid(curves, step)
    radii = np.array([radius for _, radius, _ in stacked.values()])
    edges = np.linspace(radii.min(), radii.max(), n_radius + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])

    surface = np.full((n_radius, len(grid)), np.nan)
    counts = np.zeros(n_radius, dtype=int)
    for index in range(n_radius):
        lo, hi = edges[index], edges[index + 1]
        members = [delta for delta, radius, _ in stacked.values()
                   if (lo <= radius < hi) or (index == n_radius - 1 and radius == hi)]
        counts[index] = len(members)
        if len(members) >= min_count:
            surface[index] = np.median(np.array(members), axis=0)
        else:
            counts[index] = 0          # too thin to draw
    return centres, grid, surface, counts


def plot_surface_3d(curves, out_png: Path, min_count: int = 3) -> None:
    """Ring radius and distance from the centre against energy, as a surface."""
    centres, grid, surface, counts = bin_surface(curves, min_count=min_count)
    keep = counts > 0
    centres, surface = centres[keep], surface[keep]

    x, y = np.meshgrid(grid, centres)
    z = symlog(np.nan_to_num(surface, nan=0.0))

    fig = plt.figure(figsize=(11, 8), dpi=120)
    ax = fig.add_subplot(111, projection="3d")
    ax.plot_surface(y, x, z, cmap=BLUE_CMAP, linewidth=0, antialiased=True,
                    rcount=len(centres), ccount=len(grid))

    ax.set_xlabel("ring radius (Å)", fontsize=12, color=TEXT, labelpad=10)
    ax.set_ylabel("distance from ring centre (Å)", fontsize=12, color=TEXT, labelpad=10)
    ax.set_zlabel("$\\Delta E$ (eV)", fontsize=12, color=TEXT, labelpad=12)
    ticks, labels = symlog_ticks(float(np.nanmin(z)), float(np.nanmax(z)))
    ax.set_zticks(ticks)
    ax.set_zticklabels(labels)
    ax.tick_params(colors=TEXT, labelsize=10)
    ax.view_init(elev=26, azim=-128)
    ax.set_box_aspect((1.3, 1.5, 0.9))
    for pane in (ax.xaxis, ax.yaxis, ax.zaxis):
        pane.pane.set_facecolor("#fcfcfb")
        pane.pane.set_edgecolor(GRID)

    fig.suptitle("permeation barrier against ring radius and how far the water has gone\n"
                 "median over rings in each radius bin, best roll, z on a signed log scale",
                 fontsize=14, color="#222222", y=0.97)
    # tight_layout cannot measure a 3-D axes' decorations, so it leaves a band of
    # dead space under the title; set the margins directly instead.
    fig.subplots_adjust(left=0.0, right=1.0, top=1.04, bottom=0.02)
    fig.savefig(out_png)
    plt.close(fig)


def plot_surface_heatmap(curves, out_png: Path, min_count: int = 3) -> None:
    """
    The same surface as a heatmap.

    A 3-D surface shows the shape but values cannot be read off it - the near
    ridge hides what is behind it and the perspective distorts the scale. This
    is the same array with magnitude as colour, which is the one to read numbers
    from; the 3-D figure is the one to look at.
    """
    centres, grid, surface, counts = bin_surface(curves, min_count=min_count)
    keep = counts > 0
    centres, surface, counts = centres[keep], surface[keep], counts[keep]

    fig, ax = plt.subplots(figsize=(10, 6), dpi=120)
    z = symlog(np.ma.masked_invalid(surface))
    mesh = ax.pcolormesh(grid, centres, z, cmap=BLUE_CMAP, shading="nearest")
    bar = fig.colorbar(mesh, ax=ax, pad=0.02)
    ticks, labels = symlog_ticks(float(np.min(z)), float(np.max(z)))
    bar.set_ticks(ticks)
    bar.set_ticklabels(labels)
    bar.set_label("median $\\Delta E$ (eV, signed log scale)", fontsize=12, color=TEXT)
    bar.ax.tick_params(colors=TEXT, labelsize=10)

    ax.axvline(0.0, color="#ffffff", lw=1, ls=":", zorder=3)
    ax.set_xlabel("distance from the ring centre along the axis (Å)",
                  fontsize=13, color=TEXT)
    ax.set_ylabel("ring radius (Å)", fontsize=13, color=TEXT)
    ax.tick_params(colors=TEXT, labelsize=11)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(SPINE)

    fig.suptitle("permeation barrier against ring radius and distance travelled\n"
                 f"median over the rings in each radius bin "
                 f"({counts.min()}–{counts.max()} rings per bin), best roll",
                 fontsize=14, color="#222222")
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--results", required=True,
                        help="the directory collect_survey.py wrote into - its barriers.csv "
                             "and frames.csv. Under this project's layout that is "
                             "<case>/run, NOT <case>/setup/survey, which is what "
                             "collect_survey.py's own --survey points at")
    parser.add_argument("--outdir", default=None, help="default: <results>/plots")
    parser.add_argument("--ring", action="append", default=None, metavar="RING_ID",
                        help="draw a path figure for this ring (repeatable)")
    parser.add_argument("--all", action="store_true",
                        help="draw a path figure for every ring, not just the "
                             "representatives; this is one PNG per ring")
    parser.add_argument("--no-paths", action="store_true",
                        help="skip the per-ring path figures")
    parser.add_argument("--min-bin-count", type=int, default=3,
                        help="radius bins holding fewer rings than this are left "
                             "blank on the surface figures (default: 3)")
    parser.add_argument("--no-surface", action="store_true",
                        help="skip the radius x distance surface and its heatmap")
    return parser.parse_args(argv)


def main(argv=None) -> int:
    args = parse_args(argv)
    results = Path(args.results)
    outdir = Path(args.outdir) if args.outdir else results / "plots"
    outdir.mkdir(parents=True, exist_ok=True)

    barriers = read_csv(results / "barriers.csv")
    plot_summary(barriers, outdir / "summary.png")
    plot_by_size(barriers, outdir / "by_size.png")
    print(f"wrote {outdir}/summary.png and {outdir}/by_size.png")

    frames_csv = results / "frames.csv"
    if not frames_csv.exists():
        print("no frames.csv, so no path or surface figures; re-run "
              "collect_survey.py without --no-frames", file=sys.stderr)
        return 0
    frames = load_frames(frames_csv)

    curves = best_roll_curves(frames, barriers)
    if curves:
        written = plot_size_overlays(curves, outdir)
        print(f"wrote {len(written)} per-size overlay(s): "
              + ", ".join(w.name for w in written))
        if not args.no_surface:
            plot_surface_3d(curves, outdir / "surface_3d.png", args.min_bin_count)
            plot_surface_heatmap(curves, outdir / "surface_heatmap.png",
                                 args.min_bin_count)
            print(f"wrote {outdir}/surface_3d.png and {outdir}/surface_heatmap.png")
    else:
        print(f"no {RADIUS_KEY} in barriers.csv, so no per-size or surface figures; "
              "pass --rings to collect_survey.py", file=sys.stderr)

    if args.no_paths:
        return 0

    meta = {row["ring_id"]: row for row in barriers}
    if args.ring:
        wanted = list(args.ring)
    elif args.all:
        wanted = [row["ring_id"] for row in barriers]
    else:
        representatives = results / "representatives.csv"
        wanted = ([row["ring_id"] for row in read_csv(representatives)]
                  if representatives.exists() else [])

    if not wanted:
        return 0
    paths_dir = outdir / "paths"
    paths_dir.mkdir(exist_ok=True)
    drawn = 0
    for ring_id in wanted:
        if ring_id not in frames or ring_id not in meta:
            print(f"  skipping {ring_id}: not in frames.csv", file=sys.stderr)
            continue
        plot_path(ring_id, frames[ring_id], meta[ring_id],
                  paths_dir / f"{ring_id.replace('/', '__')}.png")
        drawn += 1
    print(f"wrote {drawn} path figure(s) to {paths_dir}/")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, KeyError, OSError) as error:
        sys.exit(f"plot_ring_paths.py: {error}")
