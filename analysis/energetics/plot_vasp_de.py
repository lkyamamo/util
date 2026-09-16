#!/usr/bin/env python3
"""Plot VASP relaxed reaction energies along a water path.

The DFT counterpart of plot_min_pe.py.  Each frame of the water path is a full
VASP ionic relaxation living in its own numbered directory, and the reaction
energy is absolute rather than relative to frame 1:

    dE = E_frame - E_slab - E_water

with E_slab (the isolated silica slab) and E_water (the isolated water) passed
as --slab-energy and --water-energy.  Zero is therefore physically meaningful --
it is water infinitely far from the slab -- so these plots carry a zero line
where the LAMMPS ones do not.

Both references must come from the same functional, ENCUT, POTCAR set and cell
as the frames, or dE is meaningless.  If the slab is frozen in every frame, its
energy is a single constant; the water relaxes, so --water-energy should be
relaxed gas-phase H2O.  Use the same energy definition as the frames:
E0 = energy(sigma->0), the last "E0=" on the last ionic line of the reference
run's OSZICAR.

E_frame is the last "E0=" value on the last ionic line of OSZICAR, i.e.

     168 F= -.18981613E+03 E0= -.18971626E+03  d E =-.671262E-07

Those lines are written one per ionic step, so the same parse also yields the
energy-vs-ionic-step trace of every relaxation.  E0 is energy(sigma->0), the
smearing-extrapolated energy, which is not the free energy F: with ISMEAR=0 and
SIGMA=0.2 a run carries EENTRO ~ -0.2 eV and E0 = F - TS/2.  That entropy
does not cancel against a gas-phase water reference, so dE is written out three
ways -- from E0, from F, and from "energy without entropy" -- and only the E0
one is plotted.

The x axis is the frame index; no O-Si distance is read or computed here.

Usage:
    python plot_vasp_de.py [RUN_DIR] --slab-energy EV --water-energy EV

    RUN_DIR       Directory holding the numbered frame directories
                  (default: cwd).

Writes into RUN_DIR/plots/.
"""

import argparse
import csv
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors as mcolors
from matplotlib.ticker import MaxNLocator

# "  168 F= -.18981613E+03 E0= -.18971626E+03  d E =-.671262E-07"
# note there is no space after "d E =", so the pattern must not require one
IONIC_LINE = re.compile(
    r"^\s*(\d+)\s+F=\s*(\S+)\s+E0=\s*(\S+)\s+d\s*E\s*=\s*(\S+)"
)
# OUTCAR spells these with doubled spaces
OUTCAR_TOTEN = re.compile(r"^\s*free\s+energy\s+TOTEN\s*=\s*(\S+)")
OUTCAR_SIGMA0 = re.compile(
    r"^\s*energy\s+without\s+entropy\s*=\s*(\S+)\s+energy\(sigma->0\)\s*=\s*(\S+)"
)


def parse_args():
    p = argparse.ArgumentParser(
        description="Plot VASP relaxed reaction energies, dE = E_frame - E_slab - E_water.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("run_dir", nargs="?", type=Path, default=Path.cwd())
    p.add_argument("--slab-energy", type=float, required=True, metavar="EV",
                   help="E0 of the isolated silica slab (eV)")
    p.add_argument("--water-energy", type=float, required=True, metavar="EV",
                   help="E0 of the relaxed isolated water (eV)")
    return p.parse_args()


def parse_oszicar(path):
    """Return one dict per ionic step: {step, f_ev, e0_ev, d_e}.

    VASP writes an ionic line after each ionic step, so the list is the
    relaxation trace and its last entry is the final relaxed energy.
    """
    steps = []
    for line in path.read_text().splitlines():
        m = IONIC_LINE.match(line)
        if not m:
            continue
        steps.append(
            {
                "step": int(m.group(1)),
                "f_ev": float(m.group(2)),
                "e0_ev": float(m.group(3)),
                "d_e": float(m.group(4)),
            }
        )
    return steps


def parse_outcar_final(path):
    """Return (toten, e_without_entropy, e_sigma0, converged, finished).

    Single streaming pass keeping only the last match of each energy line --
    OUTCAR runs to ~21 MB per frame, so it is never slurped whole.  `converged`
    means the relaxation stopped on EDIFFG rather than exhausting NSW;
    `finished` means VASP wrote its closing timing block.
    """
    toten = sigma0 = noentropy = None
    converged = finished = False
    with path.open() as fh:
        for line in fh:
            m = OUTCAR_TOTEN.match(line)
            if m:
                toten = float(m.group(1))
                continue
            m = OUTCAR_SIGMA0.match(line)
            if m:
                noentropy = float(m.group(1))
                sigma0 = float(m.group(2))
                continue
            if "reached required accuracy" in line:
                converged = True
            elif "General timing and accounting informations" in line:
                finished = True
    return toten, noentropy, sigma0, converged, finished


def count_xdatcar_configs(path):
    """Number of ionic configurations in an XDATCAR, or None if absent."""
    if not path.exists():
        return None
    n = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith("Direct configuration="):
                n += 1
    return n


def collect_frames(run):
    """Parse every numbered frame directory under `run`, in numeric order.

    Frames without an OSZICAR are skipped and reported; frames that are
    unconverged or still running are kept and flagged, since a half-relaxed
    energy is still worth seeing on the scan.
    """
    frames = {}
    skipped = []
    for d in sorted(
        (p for p in run.iterdir() if p.is_dir() and p.name.isdigit()),
        key=lambda p: int(p.name),
    ):
        frame = int(d.name)
        oszicar = d / "OSZICAR"
        if not oszicar.exists():
            skipped.append((frame, "no OSZICAR"))
            continue
        steps = parse_oszicar(oszicar)
        if not steps:
            skipped.append((frame, "no ionic lines in OSZICAR"))
            continue

        f_ev = steps[-1]["f_ev"]
        e0_ev = steps[-1]["e0_ev"]

        # OSZICAR stays authoritative for f_ev/e0_ev so the summary row and the last
        # point of the trace are the same number; OUTCAR carries two more digits but
        # mixing the two makes them disagree at ~1e-6 eV for no benefit.
        outcar = d / "OUTCAR"
        if outcar.exists():
            toten, noentropy, sigma0, converged, finished = parse_outcar_final(outcar)
            if sigma0 is not None and abs(sigma0 - e0_ev) > 1e-3:
                raise RuntimeError(
                    f"frame {frame}: OUTCAR energy(sigma->0) {sigma0:.6f} does not match "
                    f"OSZICAR E0 {e0_ev:.6f}"
                )
            # entropy from OUTCAR's own TOTEN/no-entropy pair, so it is self-consistent
            eentro = None if (toten is None or noentropy is None) else toten - noentropy
        else:
            # exact for ISMEAR=0, where E0 = F - TS/2
            noentropy = 2 * e0_ev - f_ev
            eentro = f_ev - noentropy
            converged = finished = False
            print(f"note: frame {frame} has no OUTCAR; "
                  f"energy without entropy derived as 2*E0 - F, convergence unknown")

        n_config = count_xdatcar_configs(d / "XDATCAR")
        if n_config is not None and n_config != len(steps):
            print(f"warning: frame {frame} has {n_config} XDATCAR configurations but "
                  f"{len(steps)} ionic steps; the final energy may not be CONTCAR's")

        frames[frame] = {
            "steps": steps,
            "f_ev": f_ev,
            "e0_ev": e0_ev,
            "e_without_entropy_ev": noentropy,
            "eentro_ev": eentro,
            "converged": converged,
            "finished": finished,
            "has_outcar": outcar.exists(),
            "oszicar": oszicar,
        }
    return frames, skipped


def run_label(run):
    """Name the run for plot titles, e.g. 'vasp_0279 atom1'.

    Frame directories sit either directly under <run_id>/run/ or one level
    deeper at <run_id>/run/<atom>/, so the run id is not always run.parent.
    """
    if run.name == "run":
        return run.parent.name
    if run.parent.name == "run":
        return f"{run.parent.parent.name} {run.name}"
    return run.name


def style_axes(ax):
    ax.grid(axis="y", color="#dddddd", lw=1)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#bbbbbb")
    ax.tick_params(colors="#444444", labelsize=11)


def status(rec):
    """Convergence in words.  Absent OUTCAR is unknown, not failed."""
    if not rec["has_outcar"]:
        return "convergence unknown (no OUTCAR)"
    if not rec["finished"]:
        return "still running / killed"
    return "converged" if rec["converged"] else "stopped without reaching EDIFFG"


def plot_reaction_energies(run, deltas, out_png):
    """dE vs frame index.  Zero is water infinitely far from the slab.

    Convergence is not marked here -- it is carried by the `converged` and
    `finished` columns of reaction_energies.csv and by the printed summary.
    """
    keys = sorted(deltas)
    vals = [deltas[f] for f in keys]

    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=120)
    ax.axhline(0.0, color="#bbbbbb", lw=1, zorder=1)
    ax.plot(keys, vals, "-o", color="#1f77e0", lw=3, ms=6, zorder=3)

    i_min = min(range(len(vals)), key=lambda k: vals[k])
    ax.plot(keys[i_min], vals[i_min], "o", ms=16, mfc="none", mec="#1f77e0", mew=2.5, zorder=4)

    ax.set_xlabel("frame", fontsize=13, color="#444444")
    ax.set_ylabel("$\\Delta E$ (eV)", fontsize=13, color="#444444")
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))  # frames are integers
    ax.set_title(
        f"{run_label(run)} — relaxed reaction energy\n"
        "$\\Delta E = E_{frame} - E_{slab} - E_{water}$",
        fontsize=15,
    )
    style_axes(ax)

    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def plot_trace(run, frame, rec, out_png):
    """Energy vs ionic step for a single relaxation.

    Left panel is the drop from the starting energy; right panel is the same
    trace as distance above the final energy on a log axis, which is the only
    way to see the tail once the first few steps have dominated.
    """
    steps = [s["step"] for s in rec["steps"]]
    pe = [s["e0_ev"] for s in rec["steps"]]
    delta = [e - pe[0] for e in pe]
    floor = min(pe)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 5.5), dpi=120)

    ax.plot(steps, delta, "-o", color="#1f77e0", lw=2, ms=4, zorder=3)
    ax.axhline(0.0, color="#bbbbbb", lw=1, zorder=1)
    drop = f"{delta[-1]:+.4f}" if abs(delta[-1]) >= 1e-4 else f"{delta[-1]:+.2e}"
    nstep = len(steps) - 1
    ax.annotate(
        f"$E_0$ = {pe[0]:.4f} eV\n$E_f$ = {pe[-1]:.4f} eV\n"
        f"$\\Delta E$ = {drop} eV over {nstep} step{'' if nstep == 1 else 's'}\n"
        f"{status(rec)}",
        xy=(0.97, 0.95),
        xycoords="axes fraction",
        ha="right",
        va="top",
        fontsize=10,
        color="#444444",
    )
    ax.set_xlabel("ionic step", fontsize=13, color="#444444")
    ax.set_ylabel("$E_0 - E_0(\\mathrm{step\\ 1})$ (eV)", fontsize=13, color="#444444")
    style_axes(ax)

    tail = [(s, e - floor) for s, e in zip(steps, pe) if e - floor > 0]
    if tail:
        ax2.semilogy([s for s, _ in tail], [v for _, v in tail],
                     "-o", color="#e0721f", lw=2, ms=4, zorder=3)
    else:
        ax2.annotate("energy flat to machine precision", xy=(0.5, 0.5),
                     xycoords="axes fraction", ha="center", va="center",
                     fontsize=11, color="#888888")
    ax2.set_xlabel("ionic step", fontsize=13, color="#444444")
    ax2.set_ylabel("$E_0 - E_{min}$ (eV)", fontsize=13, color="#444444")
    style_axes(ax2)

    fig.suptitle(
        f"{run_label(run)} frame {frame} — VASP ionic relaxation",
        fontsize=15,
    )
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def plot_traces_overlay(run, frames, out_png):
    """Every frame's relaxation trace, coloured by frame index."""
    keys = sorted(frames)
    # a single frame would give Normalize(vmin == vmax); widen it so cmap() is defined
    lo, hi = min(keys), max(keys)
    norm = mcolors.Normalize(vmin=lo, vmax=hi if hi > lo else lo + 1)
    cmap = cm.viridis

    fig, ax = plt.subplots(figsize=(8.5, 5.5), dpi=120)
    for frame in keys:
        pe = [s["e0_ev"] for s in frames[frame]["steps"]]
        steps = [s["step"] for s in frames[frame]["steps"]]
        ax.plot(steps, [e - pe[0] for e in pe], "-", lw=1.8,
                color=cmap(norm(frame)), zorder=3)
    ax.axhline(0.0, color="#bbbbbb", lw=1, zorder=1)

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("frame", fontsize=12, color="#444444")
    cbar.ax.tick_params(colors="#444444", labelsize=10)

    # relaxations span many orders of magnitude; symlog keeps the tail readable
    ax.set_yscale("symlog", linthresh=1e-3)
    ax.set_ylim(top=1e-3)  # every trace goes downhill; drop the empty positive half
    ax.set_xlabel("ionic step", fontsize=13, color="#444444")
    ax.set_ylabel("$E_0 - E_0(\\mathrm{step\\ 1})$ (eV, symlog)",
                  fontsize=13, color="#444444")
    ax.set_title(f"{run_label(run)} — VASP relaxation traces", fontsize=15)
    style_axes(ax)

    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def write_traces_csv(frames, path):
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "step", "f_ev", "e0_ev", "de_e0_from_step0_ev"])
        for frame in sorted(frames):
            steps = frames[frame]["steps"]
            e0_first = steps[0]["e0_ev"]
            for s in steps:
                w.writerow([
                    frame, s["step"],
                    f"{s['f_ev']:.6f}",
                    f"{s['e0_ev']:.6f}",
                    f"{s['e0_ev'] - e0_first:.6f}",
                ])


def main():
    args = parse_args()
    run = args.run_dir.expanduser().resolve()
    plots = run / "plots"
    plots.mkdir(exist_ok=True)

    frames, skipped = collect_frames(run)
    if not frames:
        sys.exit(f"no frame directories with a parsable OSZICAR under {run}")

    if skipped:
        print(f"note: skipped {len(skipped)} frame director{'y' if len(skipped) == 1 else 'ies'}:")
        for frame, why in skipped:
            print(f"  frame {frame}: {why}")
        print()

    odd = [f for f in sorted(frames)
           if not (frames[f]["converged"] and frames[f]["finished"])]
    if odd:
        print(f"note: {len(odd)}/{len(frames)} relaxations are not confirmed converged:")
        for f in odd:
            print(f"  frame {f}: {status(frames[f])}")
        print()

    ref = args.slab_energy + args.water_energy
    deltas = {f: rec["e0_ev"] - ref for f, rec in frames.items()}

    png = plots / "reaction_energies.png"
    plot_reaction_energies(run, deltas, png)

    traces_dir = plots / "traces"
    traces_dir.mkdir(exist_ok=True)
    for frame, rec in sorted(frames.items()):
        plot_trace(run, frame, rec, traces_dir / f"{frame}_trace.png")
    plot_traces_overlay(run, frames, traces_dir / "all_traces.png")

    csv_path = plots / "reaction_energies.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["frame", "n_ionic_steps", "converged", "finished",
                    "e0_ev", "f_ev", "e_without_entropy_ev", "eentro_ev",
                    "slab_ev", "water_ev",
                    "de_e0_ev", "de_f_ev", "de_without_entropy_ev", "oszicar"])
        for frame in sorted(frames):
            rec = frames[frame]
            ne = rec["e_without_entropy_ev"]
            w.writerow([
                # blank rather than False when there is no OUTCAR to judge by
                frame, len(rec["steps"]),
                rec["converged"] if rec["has_outcar"] else "",
                rec["finished"] if rec["has_outcar"] else "",
                f"{rec['e0_ev']:.6f}",
                f"{rec['f_ev']:.6f}",
                f"{ne:.6f}" if ne is not None else "",
                f"{rec['eentro_ev']:.6f}" if rec["eentro_ev"] is not None else "",
                f"{args.slab_energy:.6f}",
                f"{args.water_energy:.6f}",
                f"{deltas[frame]:.6f}",
                f"{rec['f_ev'] - ref:.6f}",
                f"{ne - ref:.6f}" if ne is not None else "",
                str(rec["oszicar"]),
            ])

    traces_csv = plots / "relaxation_traces.csv"
    write_traces_csv(frames, traces_csv)

    for frame in sorted(frames):
        rec = frames[frame]
        flag = "" if (rec["converged"] and rec["finished"]) else f"  [{status(rec)}]"
        print(f"frame {frame:>3}  {len(rec['steps']):>4} steps  "
              f"E0 = {rec['e0_ev']:12.6f} eV  dE = {deltas[frame]:+9.4f} eV{flag}")

    best = min(deltas, key=lambda f: deltas[f])
    print(f"\nminimum: frame {best} at dE = {deltas[best]:+.6f} eV")
    print(f"wrote {png}")
    print(f"wrote {csv_path}")
    print(f"wrote {traces_csv}")
    print(f"wrote {len(frames)} relaxation traces under {traces_dir}")


if __name__ == "__main__":
    main()
