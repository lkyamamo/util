#!/usr/bin/env python3
"""Plot the LAMMPS potential-energy scans written to output_atom*/pe*.txt.

Each reaction-energy run writes one pe{i}.txt per walked water molecule, holding
one potential energy (eV, metal units) per scan frame, preceded by a non-numeric
header line ("start") that is skipped. The matching setup log (setup/output.txt)
holds one block per atom with the O-Si distance of every frame, so the default
x axis is that distance in Angstrom. This script writes ONE PNG PER FILE.

Usage:
    python plot_pe.py [--run-dir PATH] [options]

Options:
    --run-dir PATH    Directory holding the output_atom*/ subdirectories
                      (default: current directory).
    --glob PAT        Glob for the energy files, relative to --run-dir
                      (default: output_atom*/pe*.txt).
    --setup-log PATH  Setup log with the per-frame distances (default: the
                      first of <run-dir>/../setup/output.txt,
                      <run-dir>/setup/output.txt, <run-dir>/output.txt
                      that exists).
    --x {distance,frame}
                      X axis (default: distance; frame needs no setup log).
    --outdir PATH     Where the PNGs go (default: <run-dir>/plots).
    --absolute        Plot the raw potential energy instead of the default
                      dE = PE(frame) - PE(reference frame).
    --reference-frame N
                      Frame the energies are measured from (default: 1).
    --include-clash   Keep the frames the setup log flagged CLASH (dropped by
                      default: an added atom sitting on top of an existing one,
                      so the energy is an overlap artifact of tens of keV).
    --frames A:B      Plot only frames A through B, 1-based and inclusive.
                      Either side may be left empty (e.g. :18).
    --full-range      Keep the repulsive part of the curve. By default the y
                      window stops at the reference level, showing only the
                      well (negative dE).
    --ymin V / --ymax V
                      Hard y limits, in the plotted units; override the above.
    --summary PATH    CSV of the per-atom minima
                      (default: <outdir>/minimum_energies.csv).
    --figsize W H     Figure size in inches (default: 5 3.5).
    --dpi N           Output resolution (default: 200).
    --ext EXT         Image extension: png, pdf, svg, ... (default: png).

Each plot's minimum is ringed and labelled with dE_min = PE(min) - PE(reference),
the depth of the well; the same numbers are printed as a table and written to the
summary CSV. The reference stays fixed when --frames trims it out of the plotted
window, so curves from different trims stay comparable.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import MaxNLocator

# Categorical slot 1 (blue) — one series per figure, so one color throughout.
SERIES_COLOR = "#2a78d6"
TEXT_PRIMARY = "#0b0b0b"
TEXT_SECONDARY = "#52514e"
GRID_COLOR = "#d9d8d4"

ATOM_INDEX_RE = re.compile(r"pe(\d+)\.txt$")
# "walking the oxygen from 5.997 A to 0.500 A of silicon 5096 in 20 steps, ..."
BLOCK_HEAD_RE = re.compile(r"walking the oxygen from .* of silicon (\d+)")
# "    1    5.997  0.000      5.680 id 4005 (t2) ...  1.data  CLASH H1,H2"
TABLE_ROW_RE = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+.*?(\d+)\.data\b(.*)$")


@dataclass
class Result:
    """What one scan yielded: the minimum, and its offset from the reference."""

    label: str
    pe_path: Path
    png_path: Path
    silicon_id: str
    reference_frame: int
    reference_pe: float
    min_frame: int
    min_pe: float
    min_distance: float
    delta_min: float
    n_frames: int
    n_plotted: int


@dataclass
class ScanPath:
    """One walked water molecule, as described by a setup-log block."""

    silicon_id: str
    frames: List[int] = field(default_factory=list)
    distances: List[float] = field(default_factory=list)
    clashes: List[bool] = field(default_factory=list)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Plot output_atom*/pe*.txt potential-energy scans, one PNG per file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--run-dir", type=Path, default=Path.cwd())
    p.add_argument("--glob", default="output_atom*/pe*.txt")
    p.add_argument("--setup-log", type=Path, default=None)
    p.add_argument("--x", choices=("distance", "frame"), default="distance")
    p.add_argument("--outdir", type=Path, default=None)
    p.add_argument("--absolute", action="store_true")
    p.add_argument("--reference-frame", type=int, default=1, metavar="N")
    p.add_argument("--include-clash", action="store_true")
    p.add_argument("--summary", type=Path, default=None)
    p.add_argument("--frames", default=None, metavar="A:B")
    p.add_argument("--ymin", type=float, default=None)
    p.add_argument("--ymax", type=float, default=None)
    p.add_argument("--full-range", action="store_true")
    p.add_argument("--figsize", type=float, nargs=2, default=(5.0, 3.5), metavar=("W", "H"))
    p.add_argument("--dpi", type=int, default=200)
    p.add_argument("--ext", default="png")
    return p.parse_args()


def parse_frame_range(spec: Optional[str]) -> Tuple[Optional[int], Optional[int]]:
    """Turn "A:B" (either side optional) into 1-based inclusive bounds."""
    if spec is None:
        return None, None
    if ":" not in spec:
        raise SystemExit(f"--frames expects A:B, got {spec!r}")
    lo_s, hi_s = spec.split(":", 1)
    lo = int(lo_s) if lo_s.strip() else None
    hi = int(hi_s) if hi_s.strip() else None
    if lo is not None and lo < 1:
        raise SystemExit("--frames is 1-based; A must be >= 1")
    if lo is not None and hi is not None and hi < lo:
        raise SystemExit(f"--frames A must be <= B, got {spec!r}")
    return lo, hi


def find_setup_log(run_dir: Path, explicit: Optional[Path]) -> Optional[Path]:
    if explicit is not None:
        path = explicit.expanduser().resolve()
        if not path.is_file():
            raise SystemExit(f"--setup-log not found: {path}")
        return path
    for candidate in (run_dir.parent / "setup" / "output.txt", run_dir / "setup" / "output.txt", run_dir / "output.txt"):
        if candidate.is_file():
            return candidate
    return None


def parse_setup_log(path: Path) -> Dict[str, ScanPath]:
    """Read the setup log's per-atom frame tables.

    The log concatenates one block per walked molecule in atom1..atomN order,
    each opening with a "walking the oxygen ... of silicon <id>" line, so the
    k-th block belongs to atom{k} — the same order input.sh generated them in.
    """
    paths: Dict[str, ScanPath] = {}
    current: Optional[ScanPath] = None
    for line in path.read_text().splitlines():
        head = BLOCK_HEAD_RE.search(line)
        if head:
            current = ScanPath(silicon_id=head.group(1))
            paths[f"atom{len(paths) + 1}"] = current
            continue
        if current is None:
            continue
        row = TABLE_ROW_RE.match(line)
        if row and row.group(1) == row.group(4):  # frame number matches <n>.data
            current.frames.append(int(row.group(1)))
            current.distances.append(float(row.group(2)))
            current.clashes.append("CLASH" in row.group(5))
    if not paths:
        raise SystemExit(f"no scan blocks found in {path}")
    return paths


def read_pe(path: Path) -> np.ndarray:
    """Read one pe{i}.txt, skipping the "start" header and any blank lines."""
    values: List[float] = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        try:
            values.append(float(line))
        except ValueError:
            continue  # header line such as "start"
    if not values:
        raise SystemExit(f"no numeric energies found in {path}")
    return np.asarray(values, dtype=float)


def atom_label(path: Path) -> str:
    """Prefer the pe{i}.txt index, fall back to the parent directory name."""
    m = ATOM_INDEX_RE.search(path.name)
    return f"atom{m.group(1)}" if m else path.parent.name


def plot_one(
    path: Path,
    outdir: Path,
    args: argparse.Namespace,
    scan: Optional[ScanPath],
    frame_lo: Optional[int],
    frame_hi: Optional[int],
) -> Path:
    label = atom_label(path)
    pe = read_pe(path)
    frames = np.arange(1, pe.size + 1)

    if scan is not None:
        if len(scan.distances) != pe.size:
            raise SystemExit(
                f"{label}: {pe.size} energies in {path.name} but "
                f"{len(scan.distances)} frames in the setup log — cannot pair them"
            )
        x_all = np.asarray(scan.distances, dtype=float)
        clash_all = np.asarray(scan.clashes, dtype=bool)
        xlabel = "O–Si distance (Å)"
        subtitle = f"water O walked toward Si {scan.silicon_id}"
    else:
        x_all = frames.astype(float)
        clash_all = np.zeros(pe.size, dtype=bool)
        xlabel = "Frame"
        subtitle = None

    reference_frame = args.reference_frame
    if not 1 <= reference_frame <= pe.size:
        raise SystemExit(f"{label}: --reference-frame {reference_frame} outside 1..{pe.size}")
    if args.absolute:
        y_all = pe
        ylabel = "Potential energy (eV)"
    else:
        # The reference stays fixed even when --frames trims it out of the window.
        y_all = pe - pe[reference_frame - 1]
        ylabel = rf"$\Delta E$ relative to frame {reference_frame} (eV)"

    lo = 1 if frame_lo is None else frame_lo
    hi = pe.size if frame_hi is None else min(frame_hi, pe.size)
    keep = (frames >= lo) & (frames <= hi)
    if not keep.any():
        raise SystemExit(f"--frames {lo}:{hi} selects no frames of {path} ({pe.size} frames)")
    if not args.include_clash:
        # CLASH frames are overlap artifacts (tens of keV); they are not part of
        # the reaction path and would swamp the y scale.
        keep &= ~clash_all
        if not keep.any():
            raise SystemExit(f"{label}: every selected frame is flagged CLASH; pass --include-clash")
    frames_kept = frames[keep]
    order = np.argsort(x_all[keep])  # draw left to right along the x axis
    x, y, clash = x_all[keep][order], y_all[keep][order], clash_all[keep][order]
    frames_plot = frames_kept[order]

    fig, ax = plt.subplots(figsize=tuple(args.figsize), dpi=args.dpi)
    ax.plot(x, y, color=SERIES_COLOR, linewidth=2.0, zorder=2)
    ax.plot(
        x[~clash],
        y[~clash],
        linestyle="none",
        marker="o",
        markersize=4.5,
        color=SERIES_COLOR,
        markeredgecolor="white",
        markeredgewidth=0.8,
        zorder=3,
    )
    if clash.any():
        # Hollow markers: overlap artifacts, not points on the reaction path.
        ax.plot(
            x[clash],
            y[clash],
            linestyle="none",
            marker="o",
            markersize=5.5,
            markerfacecolor="white",
            markeredgecolor=SERIES_COLOR,
            markeredgewidth=1.6,
            zorder=3,
        )
        ax.legend(
            handles=[
                Line2D([], [], color=SERIES_COLOR, marker="o", markersize=4.5,
                       markeredgecolor="white", linewidth=2.0, label="scan"),
                Line2D([], [], color=SERIES_COLOR, marker="o", markersize=5.5,
                       markerfacecolor="white", markeredgewidth=1.6, linestyle="none",
                       label="clash (overlapping atoms)"),
            ],
            fontsize=8,
            frameon=False,
            labelcolor=TEXT_SECONDARY,
            loc="best",
        )

    # The result: the well depth relative to the reference frame. Direct-label
    # that one point rather than every point.
    i_min = int(np.argmin(y))
    delta_min = pe[frames_plot[i_min] - 1] - pe[reference_frame - 1]
    ax.plot(
        x[i_min],
        y[i_min],
        linestyle="none",
        marker="o",
        markersize=8,
        markerfacecolor="none",
        markeredgecolor=SERIES_COLOR,
        markeredgewidth=1.2,
        zorder=4,
    )
    ax.annotate(
        f"min {delta_min:+.4g} eV",
        xy=(x[i_min], y[i_min]),
        xytext=(11, 0),  # to the right of the well, where the curve leaves room
        textcoords="offset points",
        ha="left",
        va="center",
        fontsize=8,
        color=TEXT_PRIMARY,
        zorder=4,
    )

    if not args.absolute:
        # y = 0 is the reference frame; the drop to the labelled point is the answer.
        ax.axhline(0.0, color=GRID_COLOR, linewidth=1.0, zorder=1)
        ax.annotate(
            "",
            xy=(x[i_min], y[i_min]),
            xytext=(x[i_min], 0.0),
            arrowprops=dict(arrowstyle="-|>", color=GRID_COLOR, linewidth=0.9,
                            shrinkA=0, shrinkB=4),
            zorder=2,
        )

    title = f"{label} — potential energy scan"
    if subtitle:
        title += f"\n{subtitle}"
    ax.set_title(title, fontsize=11, color=TEXT_PRIMARY, pad=10)
    ax.set_xlabel(xlabel, fontsize=10, color=TEXT_SECONDARY)
    ax.set_ylabel(ylabel, fontsize=10, color=TEXT_SECONDARY)
    ax.grid(axis="y", color=GRID_COLOR, linewidth=0.6, alpha=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(GRID_COLOR)
    ax.tick_params(labelsize=9, colors=TEXT_SECONDARY, length=4, width=0.8)
    if scan is None:
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))  # frames are counts, not reals
        ax.set_xlim(x.min() - 0.5, x.max() + 0.5)
    else:
        pad = 0.03 * (x.max() - x.min())
        ax.set_xlim(x.min() - pad, x.max() + pad)
    depth = -min(float(y.min()), 0.0)
    if not args.absolute and not args.full_range and depth > 0:
        # Show only the well: the y window stops at the reference level (y = 0),
        # so the repulsive climb above it is cropped rather than compressing it.
        ax.set_ylim(bottom=y.min() - 0.12 * depth, top=0.06 * depth)
    elif depth == 0:
        # Purely repulsive scan — clipping to negative y would show nothing.
        print(f"  {label}: no frame falls below the reference; showing the full range")
    if args.ymin is not None:
        ax.set_ylim(bottom=args.ymin)
    if args.ymax is not None:
        ax.set_ylim(top=args.ymax)

    outdir.mkdir(parents=True, exist_ok=True)
    out_path = outdir / f"{label}_pe.{args.ext}"
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)

    return Result(
        label=label,
        pe_path=path,
        png_path=out_path,
        silicon_id=scan.silicon_id if scan is not None else "",
        reference_frame=reference_frame,
        reference_pe=float(pe[reference_frame - 1]),
        min_frame=int(frames_plot[i_min]),
        min_pe=float(pe[frames_plot[i_min] - 1]),
        min_distance=float(x[i_min]) if scan is not None else float("nan"),
        delta_min=float(delta_min),
        n_frames=int(pe.size),
        n_plotted=int(x.size),
    )


def main() -> None:
    args = parse_args()
    run_dir = args.run_dir.expanduser().resolve()
    if not run_dir.is_dir():
        raise SystemExit(f"--run-dir is not a directory: {run_dir}")

    paths = sorted(
        run_dir.glob(args.glob),
        key=lambda p: (int(m.group(1)) if (m := ATOM_INDEX_RE.search(p.name)) else 0, p.name),
    )
    if not paths:
        raise SystemExit(f"no files match {args.glob!r} under {run_dir}")

    scans: Dict[str, ScanPath] = {}
    if args.x == "distance":
        log_path = find_setup_log(run_dir, args.setup_log)
        if log_path is None:
            raise SystemExit(
                "no setup log found for the distance axis — pass --setup-log PATH, "
                "or use --x frame to plot against the frame index instead"
            )
        scans = parse_setup_log(log_path)
        print(f"distances from {log_path} ({len(scans)} scan block(s))")

    outdir = (args.outdir or run_dir / "plots").expanduser().resolve()
    frame_lo, frame_hi = parse_frame_range(args.frames)
    results: List[Result] = []
    for path in paths:
        label = atom_label(path)
        scan = scans.get(label) if scans else None
        if args.x == "distance" and scan is None:
            raise SystemExit(f"{label}: no matching block in the setup log (found {sorted(scans)})")
        results.append(plot_one(path, outdir, args, scan, frame_lo, frame_hi))

    print(f"\n{len(results)} plot(s) written to {outdir}\n")
    print(f"{'atom':6} {'Si':>6} {'ref frame':>9} {'min frame':>9} {'r_min (A)':>9} "
          f"{'E_ref (eV)':>15} {'E_min (eV)':>15} {'dE_min (eV)':>12}")
    for r in results:
        r_min = f"{r.min_distance:9.3f}" if np.isfinite(r.min_distance) else f"{'-':>9}"
        print(f"{r.label:6} {r.silicon_id:>6} {r.reference_frame:9d} {r.min_frame:9d} {r_min} "
              f"{r.reference_pe:15.6f} {r.min_pe:15.6f} {r.delta_min:12.4f}")

    summary_path = (args.summary or outdir / "minimum_energies.csv").expanduser().resolve()
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with summary_path.open("w") as fh:
        fh.write("atom,silicon_id,reference_frame,reference_pe_ev,min_frame,min_distance_ang,"
                 "min_pe_ev,delta_min_ev,n_frames,n_frames_plotted,pe_file,png_file\n")
        for r in results:
            fh.write(f"{r.label},{r.silicon_id},{r.reference_frame},{r.reference_pe:.10g},"
                     f"{r.min_frame},{r.min_distance:.6g},{r.min_pe:.10g},{r.delta_min:.10g},"
                     f"{r.n_frames},{r.n_plotted},{r.pe_path},{r.png_path}\n")
    print(f"\nsummary -> {summary_path}")


if __name__ == "__main__":
    main()
