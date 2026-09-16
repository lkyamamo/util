#!/usr/bin/env python3
"""Plot single-point potential-energy scans of a LAMMPS reaction-energy run.

plot_min_pe.py reads CG-minimized energies out of log.lammps.  A single-point
run never minimizes: each frame is read, a `run 0` is issued, and the potential
energy of the as-built structure is printed.  The block ordering is identical --
one block per frame, in loop order (atom1 frames 1-N, atom2 frames 1-N, ...),
each opened by the fully resolved

    read_data ../setup/atomN/frames/J.data

echo -- so only the line carrying the energy changes.  Here it is the single
row of the `run 0` thermo table,

    Step PotEng TotEng Fmax Fnorm
           0  -20002.2115429595  -19801.8872397809   0.0042003727   0.0110931577

whose PotEng column is the single-point energy.  That value is cross-checked
against the `print ${mype} append output_atomN/peN.txt` dump the run also
writes (output_atomN/peN.txt), which must agree frame for frame.

Because nothing is relaxed, the frames that push the water O inside the Si
first-neighbour shell are hugely repulsive (~+45000 eV at the last frame),
so the per-atom figure carries two panels: the full scan on a symlog axis and
a linear zoom on the bound well.

Usage:
    python plot_sp_pe.py [RUN_DIR] [--setup-dir PATH] [--zoom-ev V]

    RUN_DIR       Directory holding log.lammps and output_atom*/ (default: cwd).
    --setup-dir   Directory holding input.sh and atom*/frames/manifest.csv
                  (default: RUN_DIR/../setup).
    --zoom-ev     Half-width, in eV, of the zoom panel around dE = 0 (default: 5).

Writes into RUN_DIR/plots/.
"""

import argparse
import csv
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

READ_DATA = re.compile(r"^read_data\s+\.\./setup/(atom\d+)/frames/(\d+)\.data\s*$")
SILICON_ARG = re.compile(r"--silicon-id\s+(\d+)")
THERMO_HEADER = re.compile(r"^Step(\s+\S+)+\s*$")


def parse_args():
    p = argparse.ArgumentParser(
        description="Plot single-point potential-energy scans from log.lammps.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("run_dir", nargs="?", type=Path, default=Path.cwd())
    p.add_argument("--setup-dir", type=Path, default=None,
                   help="defaults to RUN_DIR/../setup")
    p.add_argument("--zoom-ev", type=float, default=5.0,
                   help="half-width, in eV, of the zoom panel around dE = 0")
    return p.parse_args()


def parse_silicon_ids(setup):
    """Map atomN -> silicon id, from the --silicon-id flags in setup/input.sh."""
    text = (setup / "input.sh").read_text()
    ids = {}
    current = None
    for line in text.splitlines():
        s = line.strip()
        m = re.match(r"^cd\s+(atom\d+)$", s)
        if m:
            current = m.group(1)
            continue
        m = SILICON_ARG.search(s)
        if m and current:
            ids[current] = int(m.group(1))
            current = None
    if not ids:
        raise RuntimeError(f"no --silicon-id flags found in {setup / 'input.sh'}")
    return ids


def parse_log(log_path):
    """Return ({(atom, frame): single_point_energy_eV}, ordered keys).

    Same block walk as plot_min_pe.py: the resolved `read_data` line opens a
    frame, and the first thermo row that follows closes it.  A `run 0` emits
    exactly one row (Step 0), so the frame's energy is that row's PotEng.
    """
    lines = log_path.read_text().splitlines()
    results = {}
    order = []
    current = None
    columns = None
    for idx, line in enumerate(lines):
        s = line.strip()
        m = READ_DATA.match(s)
        if m:
            if current is not None:
                raise RuntimeError(
                    f"{current[0]} frame {current[1]} has no thermo output before "
                    f"the next read_data at line {idx+1}"
                )
            current = (m.group(1), int(m.group(2)))
            columns = None
            continue
        if current is None:
            continue
        if THERMO_HEADER.match(s):
            columns = s.split()
            if "PotEng" not in columns:
                raise RuntimeError(f"thermo table at line {idx+1} has no PotEng column")
            continue
        if columns is None:
            continue
        parts = s.split()
        if len(parts) != len(columns):
            continue  # WARNING lines and the like between header and data
        try:
            vals = [float(p) for p in parts]
        except ValueError:
            continue
        if current in results:
            raise RuntimeError(f"duplicate energy for {current} at line {idx+1}")
        results[current] = vals[columns.index("PotEng")]
        order.append(current)
        current = None
        columns = None
    if current is not None:
        raise RuntimeError(f"log ends before {current[0]} frame {current[1]} reported an energy")
    return results, order


def read_pe_dump(run, atom):
    """Return the per-frame energies the run appended to output_atomN/peN.txt.

    The file is a "start" marker followed by one energy per frame, in loop
    order, so index i is frame i+1.  Returns None if the file is missing.
    """
    n = atom[len("atom"):]
    path = run / f"output_{atom}" / f"pe{n}.txt"
    if not path.exists():
        return None
    vals = []
    for tok in path.read_text().split():
        try:
            vals.append(float(tok))
        except ValueError:
            continue  # the leading "start" line
    return {i + 1: v for i, v in enumerate(vals)}


def check_against_pe_dumps(run, atoms, energies):
    """Cross-check every log energy against the peN.txt the run wrote itself."""
    checked = 0
    for atom in sorted(atoms):
        dump = read_pe_dump(run, atom)
        if dump is None:
            print(f"note: no pe dump for {atom}; skipping cross-check")
            continue
        for (a, frame), e in energies.items():
            if a != atom:
                continue
            if frame not in dump:
                raise RuntimeError(f"{atom} frame {frame} missing from its pe dump")
            if abs(dump[frame] - e) > 1e-6:
                raise RuntimeError(
                    f"{atom} frame {frame}: log PotEng {e:.10f} does not match "
                    f"pe dump {dump[frame]:.10f}; the by-block pairing is wrong"
                )
            checked += 1
    return checked


def read_manifest(setup, atom):
    """Return {frame: nominal O-Si distance} from the frame-generation manifest."""
    path = setup / atom / "frames" / "manifest.csv"
    with path.open() as fh:
        return {int(r["frame"]): float(r["o_si_distance"]) for r in csv.DictReader(fh)}


def style_axes(ax):
    ax.grid(axis="y", color="#dddddd", lw=1)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#bbbbbb")
    ax.tick_params(colors="#444444", labelsize=11)


def plot_atom(run_name, atom, si_id, frames, dists, energies, zoom_ev, out_png):
    """Full scan (symlog) beside a linear zoom on the bound well."""
    ref = energies[0]
    delta = [e - ref for e in energies]
    i_min = min(range(len(delta)), key=lambda k: delta[k])

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 5.5), dpi=120)

    ax.plot(dists, delta, "-o", color="#1f77e0", lw=2.5, ms=5, zorder=3)
    ax.axhline(0.0, color="#bbbbbb", lw=1, zorder=1)
    # the unrelaxed wall runs to ~+4.5e4 eV; the well is ~1 eV deep
    ax.set_yscale("symlog", linthresh=1.0)
    ax.set_xlabel("O–Si distance (Å)", fontsize=13, color="#444444")
    ax.set_ylabel("$\\Delta E$ relative to frame 1 (eV, symlog)", fontsize=13, color="#444444")
    ax.set_title("full scan", fontsize=12, color="#444444")
    style_axes(ax)

    zoom = [(d, v) for d, v in zip(dists, delta) if abs(v) <= zoom_ev]
    if zoom:
        ax2.plot([d for d, _ in zoom], [v for _, v in zoom],
                 "-o", color="#1f77e0", lw=2.5, ms=6, zorder=3)
        ax2.axhline(0.0, color="#bbbbbb", lw=1, zorder=1)
        if abs(delta[i_min]) <= zoom_ev:
            ax2.plot(dists[i_min], delta[i_min], "o", ms=16, mfc="none",
                     mec="#1f77e0", mew=2.5, zorder=4)
            ax2.annotate(
                f"min {delta[i_min]:+.3f} eV\n@ {dists[i_min]:.3f} Å (frame {frames[i_min]})",
                xy=(dists[i_min], delta[i_min]),
                xytext=(8, 2),
                textcoords="offset points",
                va="center",
                fontsize=11,
            )
        n_hidden = len(dists) - len(zoom)
        if n_hidden:
            ax2.annotate(
                f"{n_hidden} of {len(dists)} frames outside ±{zoom_ev:g} eV",
                xy=(0.97, 0.95), xycoords="axes fraction", ha="right", va="top",
                fontsize=10, color="#888888",
            )
    else:
        ax2.annotate(f"no frame within ±{zoom_ev:g} eV of frame 1", xy=(0.5, 0.5),
                     xycoords="axes fraction", ha="center", va="center",
                     fontsize=11, color="#888888")
    ax2.set_xlabel("O–Si distance (Å)", fontsize=13, color="#444444")
    ax2.set_ylabel("$\\Delta E$ relative to frame 1 (eV)", fontsize=13, color="#444444")
    ax2.set_title(f"zoom, |$\\Delta E$| ≤ {zoom_ev:g} eV", fontsize=12, color="#444444")
    style_axes(ax2)

    fig.suptitle(
        f"{run_name} {atom} — single-point potential energy scan\n"
        f"water O walked toward Si {si_id}",
        fontsize=15,
    )
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def main():
    args = parse_args()
    run = args.run_dir.expanduser().resolve()
    setup = (args.setup_dir or run.parent / "setup").expanduser().resolve()
    run_name = run.parent.name
    plots = run / "plots"
    plots.mkdir(exist_ok=True)

    silicon_ids = parse_silicon_ids(setup)
    log_path = run / "log.lammps"
    energies, order = parse_log(log_path)

    n_checked = check_against_pe_dumps(run, silicon_ids, energies)
    print(f"cross-checked {n_checked}/{len(energies)} log energies against the pe dumps\n")

    rows = []
    per_frame_rows = []
    for atom in sorted(silicon_ids):
        manifest = read_manifest(setup, atom)
        frames = sorted(f for (a, f) in energies if a == atom)
        if not frames:
            raise RuntimeError(f"no single-point energies found for {atom}")
        dists = [manifest[f] for f in frames]
        pes = [energies[(atom, f)] for f in frames]

        out_png = plots / f"{atom}_sp_pe.png"
        plot_atom(run_name, atom, silicon_ids[atom], frames, dists, pes, args.zoom_ev, out_png)

        ref = pes[0]
        i_min = min(range(len(pes)), key=lambda k: pes[k])
        rows.append(
            {
                "atom": atom,
                "silicon_id": silicon_ids[atom],
                "reference_frame": frames[0],
                "reference_pe_ev": f"{ref:.6f}",
                "min_frame": frames[i_min],
                "min_distance_ang": f"{dists[i_min]:.4f}",
                "min_pe_ev": f"{pes[i_min]:.6f}",
                "delta_min_ev": f"{pes[i_min] - ref:.6f}",
                "n_frames": len(frames),
                "log_file": str(log_path),
                "png_file": str(out_png),
            }
        )
        for f, d, e in zip(frames, dists, pes):
            per_frame_rows.append([atom, silicon_ids[atom], f, f"{d:.4f}",
                                   f"{e:.6f}", f"{e - ref:.6f}"])

    csv_path = plots / "single_point_minimum_energies.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    per_frame = plots / "single_point_energies_by_frame.csv"
    with per_frame.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["atom", "silicon_id", "frame", "o_si_distance_ang",
                    "pe_ev", "delta_pe_ev"])
        w.writerows(per_frame_rows)

    for r in rows:
        print(f"{r['atom']}  Si {r['silicon_id']:>5}  min frame {r['min_frame']} "
              f"@ {r['min_distance_ang']} A  dE = {r['delta_min_ev']} eV")
    print(f"\nwrote {csv_path}")
    print(f"wrote {per_frame}")
    print(f"wrote {len(rows)} figures under {plots}")


if __name__ == "__main__":
    main()
