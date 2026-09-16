#!/usr/bin/env python3
"""Plot CG-minimized potential-energy scans of a LAMMPS reaction-energy run.

The minimized counterpart of plot_sp_pe.py: dE relative to frame 1 vs O-Si
distance, except the energies are the CG-minimized final energies parsed out of
log.lammps rather than single-point energies.

Each "Minimization stats:" block in log.lammps carries

    Stopping criterion = <...>
    Energy initial, next-to-last, final =
          <E_initial>   <E_next_to_last>   <E_final>

so the line two after "Stopping criterion" holds the three energies and the
third column is the minimized energy. Blocks appear in loop order
(atom1 frames 1-N, atom2 frames 1-N, ...), and are cross-checked against the
resolved `read_data ../setup/atomN/frames/J.data` line of each block.

The same log also carries the per-iteration thermo table of every CG
minimization ("Step PotEng TotEng ..." through "Loop time of"). Those tables
are read in file order and matched one-to-one, by order, with the minimization
blocks above, giving an energy-vs-iteration trace per frame.

Usage:
    python plot_min_pe.py [RUN_DIR] [--setup-dir PATH]

    RUN_DIR       Directory holding log.lammps and output_atom*/ (default: cwd).
    --setup-dir   Directory holding input.sh and atom*/frames/manifest.csv
                  (default: RUN_DIR/../setup).

Writes into RUN_DIR/plots/.
"""

import argparse
import csv
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm, colors as mcolors

READ_DATA = re.compile(r"^read_data\s+\.\./setup/(atom\d+)/frames/(\d+)\.data\s*$")
SILICON_ARG = re.compile(r"--silicon-id\s+(\d+)")
THERMO_HEADER = re.compile(r"^Step(\s+\S+)+\s*$")


def parse_args():
    p = argparse.ArgumentParser(
        description="Plot CG-minimized potential-energy scans from log.lammps.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("run_dir", nargs="?", type=Path, default=Path.cwd())
    p.add_argument("--setup-dir", type=Path, default=None,
                   help="defaults to RUN_DIR/../setup")
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
    """Return ({(atom, frame): (minimized_energy_eV, criterion)}, ordered keys)."""
    lines = log_path.read_text().splitlines()
    results = {}
    order = []
    current = None
    for idx, line in enumerate(lines):
        s = line.strip()
        m = READ_DATA.match(s)
        if m:
            current = (m.group(1), int(m.group(2)))
            continue
        if s.startswith("Stopping criterion"):
            if current is None:
                raise RuntimeError(f"minimization block at line {idx+1} has no read_data")
            criterion = s.split("=", 1)[1].strip()
            # line + 1 is "Energy initial, next-to-last, final ="
            # line + 2 holds the three values
            vals = lines[idx + 2].split()
            if len(vals) != 3:
                raise RuntimeError(f"unexpected energy line at {idx+3}: {lines[idx+2]!r}")
            if current in results:
                raise RuntimeError(f"duplicate minimization for {current} at line {idx+1}")
            results[current] = (float(vals[2]), criterion)
            order.append(current)
            current = None
    return results, order


def read_thermo_blocks(log_path):
    """Read every thermo table out of a LAMMPS log, in file order.

    A table runs from its "Step ..." header to the "Loop time of" line that
    closes the run/minimize.  Returns a list of
    {"columns": [name, ...], "rows": {name: [float, ...]}} — one entry per
    run or minimize, whatever the thermo_style happened to be.
    """
    blocks = []
    cur = None
    for line in log_path.read_text().splitlines():
        s = line.strip()
        if THERMO_HEADER.match(s):
            if cur is not None:  # header with no closing "Loop time of"
                blocks.append(cur)
            cols = s.split()
            cur = {"columns": cols, "rows": {c: [] for c in cols}}
            continue
        if cur is None:
            continue
        if s.startswith("Loop time of"):
            blocks.append(cur)
            cur = None
            continue
        parts = s.split()
        if len(parts) != len(cur["columns"]):
            continue  # WARNING lines and the like inside a run
        try:
            vals = [float(p) for p in parts]
        except ValueError:
            continue
        for c, v in zip(cur["columns"], vals):
            cur["rows"][c].append(v)
    if cur is not None:  # log truncated mid-run
        blocks.append(cur)
    return blocks


def match_thermo_to_minimizations(order, blocks):
    """Pair thermo tables with minimizations by file order, one-to-one."""
    if len(blocks) != len(order):
        raise RuntimeError(
            f"{len(blocks)} thermo tables but {len(order)} minimizations; "
            "cannot pair them by order"
        )
    return dict(zip(order, blocks))


def read_manifest(setup, atom):
    """Return ({frame: nominal O-Si distance}, walked water O id or None).

    The O id comes from the manifest's o_id column; it is the same on every row.
    """
    path = setup / atom / "frames" / "manifest.csv"
    with path.open() as fh:
        rows = list(csv.DictReader(fh))
    dists = {int(r["frame"]): float(r["o_si_distance"]) for r in rows}
    o_ids = {r["o_id"] for r in rows if r.get("o_id")}
    if len(o_ids) > 1:
        raise RuntimeError(f"{path}: o_id varies across frames ({sorted(o_ids)})")
    return dists, (int(o_ids.pop()) if o_ids else None)


def read_data_coords(path, wanted_ids):
    """Return {id: (x, y, z)} plus box lengths from a LAMMPS data file."""
    box = {}
    coords = {}
    in_atoms = False
    with path.open() as fh:
        for line in fh:
            s = line.strip()
            if not in_atoms:
                for lo, hi, key in (("xlo", "xhi", "x"), ("ylo", "yhi", "y"), ("zlo", "zhi", "z")):
                    if s.endswith(f"{lo} {hi}"):
                        parts = s.split()
                        box[key] = float(parts[1]) - float(parts[0])
                if s.startswith("Atoms"):
                    in_atoms = True
                continue
            if not s:
                continue
            parts = s.split()
            aid = int(parts[0])
            if aid in wanted_ids:
                coords[aid] = tuple(float(v) for v in parts[2:5])
                if len(coords) == len(wanted_ids):
                    break
    return coords, box


def mic_distance(a, b, box):
    """Minimum-image distance in an orthogonal periodic box."""
    total = 0.0
    for ai, bi, L in zip(a, b, (box["x"], box["y"], box["z"])):
        d = ai - bi
        d -= L * round(d / L)
        total += d * d
    return math.sqrt(total)


def relaxed_distance(run, atom, frame, o_id, si_id):
    """O-Si distance measured in the minimized structure, or None if absent."""
    path = run / f"output_{atom}" / f"{frame}_min.data"
    if o_id is None or not path.exists():
        return None
    coords, box = read_data_coords(path, {o_id, si_id})
    if o_id not in coords or si_id not in coords:
        return None
    return mic_distance(coords[o_id], coords[si_id], box)


def plot_atom(run_name, atom, si_id, dists, energies, out_png):
    ref = energies[0]
    delta = [e - ref for e in energies]

    fig, ax = plt.subplots(figsize=(8, 5.5), dpi=120)
    ax.plot(dists, delta, "-o", color="#1f77e0", lw=3, ms=6, zorder=3)

    i_min = min(range(len(delta)), key=lambda k: delta[k])
    ax.plot(dists[i_min], delta[i_min], "o", ms=16, mfc="none", mec="#1f77e0", mew=2.5, zorder=4)
    ax.annotate(
        f"min {delta[i_min]:+.3f} eV" if delta[i_min] else "min +0 eV",
        xy=(dists[i_min], delta[i_min]),
        xytext=(8, 2),
        textcoords="offset points",
        va="center",
        fontsize=11,
    )

    ax.set_xlabel("O–Si distance (Å)", fontsize=13, color="#444444")
    ax.set_ylabel("$\\Delta E$ relative to frame 1 (eV)", fontsize=13, color="#444444")
    ax.set_title(
        f"{run_name} {atom} — minimized potential energy scan\n"
        f"water O walked toward Si {si_id}",
        fontsize=15,
    )
    style_axes(ax)

    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def style_axes(ax):
    ax.grid(axis="y", color="#dddddd", lw=1)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#bbbbbb")
    ax.tick_params(colors="#444444", labelsize=11)


def plot_trace(run_name, atom, si_id, frame, dist, thermo, criterion, out_png):
    """Energy vs CG iteration for a single minimization.

    Left panel is the drop from the starting energy; right panel is the same
    trace as distance above the final energy on a log axis, which is the only
    way to see the tail once the first few iterations have dominated.
    """
    steps = thermo["rows"]["Step"]
    pe = thermo["rows"]["PotEng"]
    delta = [e - pe[0] for e in pe]
    floor = min(pe)

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 5.5), dpi=120)

    ax.plot(steps, delta, "-o", color="#1f77e0", lw=2, ms=4, zorder=3)
    ax.axhline(0.0, color="#bbbbbb", lw=1, zorder=1)
    drop = f"{delta[-1]:+.4f}" if abs(delta[-1]) >= 1e-4 else f"{delta[-1]:+.2e}"
    niter = len(steps) - 1
    ax.annotate(
        f"$E_0$ = {pe[0]:.4f} eV\n$E_f$ = {pe[-1]:.4f} eV\n"
        f"$\\Delta E$ = {drop} eV over {niter} iteration{'' if niter == 1 else 's'}\n"
        f"stopped on {criterion}",
        xy=(0.97, 0.95),
        xycoords="axes fraction",
        ha="right",
        va="top",
        fontsize=10,
        color="#444444",
    )
    ax.set_xlabel("CG iteration", fontsize=13, color="#444444")
    ax.set_ylabel("$E - E_0$ (eV)", fontsize=13, color="#444444")
    style_axes(ax)

    tail = [(s, e - floor) for s, e in zip(steps, pe) if e - floor > 0]
    if tail:
        ax2.semilogy([s for s, _ in tail], [v for _, v in tail],
                     "-o", color="#e0721f", lw=2, ms=4, zorder=3)
    else:
        ax2.annotate("energy flat to machine precision", xy=(0.5, 0.5),
                     xycoords="axes fraction", ha="center", va="center",
                     fontsize=11, color="#888888")
    ax2.set_xlabel("CG iteration", fontsize=13, color="#444444")
    ax2.set_ylabel("$E - E_{min}$ (eV)", fontsize=13, color="#444444")
    style_axes(ax2)

    fig.suptitle(
        f"{run_name} {atom} frame {frame} — CG minimization\n"
        f"O–Si {si_id} target distance {dist:.3f} Å",
        fontsize=15,
    )
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def plot_traces_overlay(run_name, atom, si_id, frames, dists, thermos, out_png):
    """Every frame's minimization trace for one atom, coloured by O-Si distance."""
    norm = mcolors.Normalize(vmin=min(dists), vmax=max(dists))
    cmap = cm.viridis

    fig, ax = plt.subplots(figsize=(8.5, 5.5), dpi=120)
    for frame, dist in zip(frames, dists):
        pe = thermos[(atom, frame)]["rows"]["PotEng"]
        steps = thermos[(atom, frame)]["rows"]["Step"]
        ax.plot(steps, [e - pe[0] for e in pe], "-", lw=1.8, color=cmap(norm(dist)), zorder=3)
    ax.axhline(0.0, color="#bbbbbb", lw=1, zorder=1)

    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_label("O–Si target distance (Å)", fontsize=12, color="#444444")
    cbar.ax.tick_params(colors="#444444", labelsize=10)

    # relaxations span ~1e-6 eV (far frames) to ~1e4 eV (overlapping frames)
    ax.set_yscale("symlog", linthresh=1e-3)
    ax.set_ylim(top=1e-3)  # every trace goes downhill; drop the empty positive half
    ax.set_xlabel("CG iteration", fontsize=13, color="#444444")
    ax.set_ylabel("$E - E_0$ (eV, symlog)", fontsize=13, color="#444444")
    ax.set_title(
        f"{run_name} {atom} — CG minimization traces\n"
        f"water O walked toward Si {si_id}",
        fontsize=15,
    )
    style_axes(ax)

    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def write_traces_csv(thermos, silicon_ids, path):
    with path.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["atom", "silicon_id", "frame", "step", "pe_ev",
                    "delta_pe_from_step0_ev", "fmax", "fnorm"])
        for (atom, frame), block in thermos.items():
            rows = block["rows"]
            pe0 = rows["PotEng"][0]
            for i, step in enumerate(rows["Step"]):
                w.writerow([
                    atom, silicon_ids[atom], frame, int(step),
                    f"{rows['PotEng'][i]:.6f}",
                    f"{rows['PotEng'][i] - pe0:.6f}",
                    f"{rows['Fmax'][i]:.8f}" if "Fmax" in rows else "",
                    f"{rows['Fnorm'][i]:.8f}" if "Fnorm" in rows else "",
                ])


def main():
    args = parse_args()
    run = args.run_dir.expanduser().resolve()
    setup = (args.setup_dir or run.parent / "setup").expanduser().resolve()
    run_name = run.parent.name
    plots = run / "plots"
    plots.mkdir(exist_ok=True)

    silicon_ids = parse_silicon_ids(setup)
    log_path = run / "log.lammps"
    parsed, order = parse_log(log_path)
    energies = {k: v[0] for k, v in parsed.items()}
    criteria = {k: v[1] for k, v in parsed.items()}

    # thermo tables come out of the log in the same order the minimizations do
    thermos = match_thermo_to_minimizations(order, read_thermo_blocks(log_path))
    for key, block in thermos.items():
        final = block["rows"]["PotEng"][-1]
        if abs(final - energies[key]) > 1e-4:
            raise RuntimeError(
                f"{key[0]} frame {key[1]}: thermo final PE {final:.6f} does not match "
                f"minimization stats energy {energies[key]:.6f}; the by-order pairing "
                "of thermo tables to minimizations is wrong"
            )

    traces_dir = plots / "traces"
    traces_dir.mkdir(exist_ok=True)

    odd = sorted(k for k, c in criteria.items() if c != "energy tolerance")
    if odd:
        print(f"note: {len(odd)}/{len(criteria)} minimizations did not stop on energy tolerance:")
        for atom, frame in odd:
            print(f"  {atom} frame {frame}: {criteria[(atom, frame)]}")
        print()

    rows = []
    per_frame_rows = []
    for atom in sorted(silicon_ids):
        si_id = silicon_ids[atom]
        manifest, o_id = read_manifest(setup, atom)
        frames = sorted(f for (a, f) in energies if a == atom)
        if not frames:
            raise RuntimeError(f"no minimized energies found for {atom}")
        dists = [manifest[f] for f in frames]
        pes = [energies[(atom, f)] for f in frames]
        relaxed = {f: relaxed_distance(run, atom, f, o_id, si_id) for f in frames}

        out_png = plots / f"{atom}_pe.png"
        plot_atom(run_name, atom, si_id, dists, pes, out_png)

        # one energy-vs-iteration plot per minimization, plus an overlay
        atom_dir = traces_dir / atom
        atom_dir.mkdir(exist_ok=True)
        for frame, dist in zip(frames, dists):
            plot_trace(
                run_name, atom, si_id, frame, dist, thermos[(atom, frame)],
                criteria[(atom, frame)], atom_dir / f"{frame}_trace.png",
            )
        plot_traces_overlay(run_name, atom, si_id, frames, dists, thermos,
                            traces_dir / f"{atom}_traces.png")

        ref = pes[0]
        i_min = min(range(len(pes)), key=lambda k: pes[k])
        d_min = relaxed[frames[i_min]]
        rows.append(
            {
                "atom": atom,
                "silicon_id": si_id,
                "reference_frame": frames[0],
                "reference_pe_ev": f"{ref:.6f}",
                "min_frame": frames[i_min],
                "min_distance_ang": f"{dists[i_min]:.4f}",
                "min_distance_relaxed_ang": f"{d_min:.4f}" if d_min is not None else "",
                "min_pe_ev": f"{pes[i_min]:.6f}",
                "delta_min_ev": f"{pes[i_min] - ref:.6f}",
                "min_stopping_criterion": criteria[(atom, frames[i_min])],
                "n_frames": len(frames),
                "log_file": str(log_path),
                "png_file": str(out_png),
            }
        )

        for f in frames:
            rel = relaxed[f]
            trace = thermos[(atom, f)]["rows"]["PotEng"]
            per_frame_rows.append(
                [atom, si_id, f, f"{manifest[f]:.4f}",
                 f"{rel:.4f}" if rel is not None else "",
                 f"{trace[0]:.6f}",
                 f"{energies[(atom, f)]:.6f}",
                 f"{energies[(atom, f)] - ref:.6f}",
                 f"{trace[-1] - trace[0]:.6f}",
                 len(trace) - 1,
                 criteria[(atom, f)]]
            )

    csv_path = plots / "minimum_energies.csv"
    with csv_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)

    # per-frame dump, useful for checking how far the O relaxed from its target
    per_frame = plots / "minimized_energies_by_frame.csv"
    with per_frame.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            ["atom", "silicon_id", "frame", "o_si_distance_ang",
             "o_si_distance_relaxed_ang", "initial_pe_ev", "min_pe_ev", "delta_pe_ev",
             "relaxation_ev", "n_iterations", "stopping_criterion"]
        )
        w.writerows(per_frame_rows)

    traces_csv = plots / "minimization_traces.csv"
    write_traces_csv(thermos, silicon_ids, traces_csv)

    for r in rows:
        print(f"{r['atom']}  Si {r['silicon_id']:>5}  min frame {r['min_frame']} "
              f"@ {r['min_distance_ang']} A  dE = {r['delta_min_ev']} eV")
    print(f"\nwrote {csv_path}")
    print(f"wrote {per_frame}")
    print(f"wrote {traces_csv}")
    print(f"wrote {len(thermos)} minimization traces under {traces_dir}")


if __name__ == "__main__":
    main()
