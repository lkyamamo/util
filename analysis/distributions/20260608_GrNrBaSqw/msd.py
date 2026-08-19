"""
msd.py — Mean Square Displacement from LAMMPS dump trajectories.

Same multiple-time-origin approach as analysis/dynamics/src/msd.cpp (and this
script's own vdos.py, which computes VACF the same way): no globally unwrapped
coordinates needed. See METHOD below.

QUICK START
-----------
1. Set DUMP_FILE to your LAMMPS custom dump trajectory (element, x, y, z
   columns required). Defaults to dynamics.lammpstrj (env var DYNAMICS_TRAJ)
   — the same higher-frequency trajectory vdos.py reads; it already has
   wrapped x y z, which is all this needs.
2. Set the required environment variables listed under CONFIGURATION below.
3. Run:  python msd.py

CONFIGURATION
-------------
Each setting comes from exactly one environment variable — no fallback names,
and no silent defaults for anything that changes the result.

REQUIRED (unset or empty is a hard error; the script will not guess):
  DYNAMICS_DT        fs between dumped frames — the ONLY variable msd.py and
                     vdos.py share, since it describes dynamics.lammpstrj
                     rather than either analysis
  MSD_CORR_LENGTH    fs; max time lag
  MSD_CORR_INTERVAL  fs; spacing between reference frames
  MSD_FIT_FRACTION   tail fraction of the window used for the D fit

Optional (the default is the identity choice — read the whole trajectory —
so leaving it unset cannot distort a result):
  MSD_N_FRAMES       max frames to read; 0 = all      (default 0)
  MSD_STRIDE         read every Nth frame             (default 1)

The MSD_ prefix is the point: msd.py and vdos.py read the same
dynamics.lammpstrj but want different settings — MSD needs a long lag to reach
the diffusive regime, VDOS a short one for frequency resolution — so a name
that reached both would make one of them silently wrong. DYNAMICS_DT is the
sole shared name: dt is a fact about the trajectory, so one value must reach
both scripts.

CORR_LENGTH previously defaulted to 75% of the trajectory when unset. It no
longer does: that default silently set both the lag range and the diffusion
fit window, so D depended on a number nobody chose.

METHOD
------
MSD(t) is built from multiple time origins, same accumulation as msd.cpp's
get_msd(): a reference frame is taken every CORR_INTERVAL fs, and compared
against every frame up to CORR_LENGTH fs later.

No unwrapped (xu yu zu) coordinates are needed. Instead, each reference-to-
current displacement is corrected by the minimum-image convention — same
trick msd.cpp uses (apply_pbc: shift by one box length if the raw difference
exceeds half the box) — which is valid because CORR_LENGTH is short enough
that no atom diffuses more than half a box length within one reference-to-
current comparison. This script's version, disp -= L*round(disp/L), is
vectorized over every atom and lag at once per reference (np.round handles
any number of wraps in one op, a minor generalization of msd.cpp's single
if/else), so the only Python-level loop is over the (few) reference frames —
same O(n_refs) looping vdos.py's compute_vacf_multi_origin uses, not
msd.cpp's O(n_refs x corr_length x n_atoms) nested loop.

Elements are combined into a total via a mole-fraction-weighted average
(Total(t) = sum_el (N_el/N_total)*MSD_el(t)) — MSD is a plain average
quantity, unlike vdos.py's phonon-DOS convention, so no extra normalization
is needed.

As a convenience, a self-diffusion coefficient is estimated per element (and
total) via the 3D Einstein relation MSD(t) = 6*D*t, linear-fit over the last
FIT_FRACTION of the correlation window (the early, non-diffusive/ballistic
part of MSD(t) is excluded from the fit) and printed — not written to the
CSV, since it's a derived summary rather than part of the MSD(t) curve.

OUTPUT
------
- msd.csv — time_fs, then MSD_<element> per element and MSD_total (Å²)
- <date>_msd/ — ONE PNG PER CURVE: <element>.png per species and total.png,
  plus all_species.png, the overlay kept because comparing species on one axes
  is how a diffusion run is read.  The directory is created by the runner, not
  by this script; set MSD_PLOT_DIR= (empty) to write only the CSV.

DEPENDENCIES
------------
  pip install numpy matplotlib

COLUMN LAYOUT
-------------
Expects LAMMPS custom dump format with at minimum columns:
id element x y z (velocities are not read/needed). Column positions are read
automatically from the ITEM: ATOMS header line. Atom order must be
consistent frame-to-frame (dump_modify ... sort id in the LAMMPS input,
already used throughout this pipeline). Positions must be wrapped x y z (the
dump style this pipeline already uses) — NOT unwrapped xu yu zu, which this
script does not look for and does not need.
"""

import os
import re
from datetime import date

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# =============================================================================
# CONFIGURATION — edit these variables between runs
# =============================================================================

# Input trajectory file — same one vdos.py reads (already has wrapped x y z).
DUMP_FILE = os.environ.get("DYNAMICS_TRAJ", "dynamics.lammpstrj")

# Every setting below is read from exactly one environment variable, all
# MSD_-prefixed except DYNAMICS_DT — no fallback names. msd.py and vdos.py read
# the same dynamics.lammpstrj but want different settings (MSD needs a long lag
# to reach the diffusive regime, VDOS a short one for frequency resolution), so
# a name that reached both would make one of them silently wrong.
# Required settings have no default: an unset or empty value is a hard error,
# never a silently substituted number. Missing ones are collected so a single
# run reports every one of them at once instead of failing one at a time.
_MISSING = []

# =============================================================================
# Sampling reporting — see "How much data is enough" in the README
# =============================================================================

# Relative statistical error of an averaged quantity is ~1/sqrt(M), where M is
# the number of INDEPENDENT contributions.  Ladder: 1e2 = 10%, 1e3 = 3% (the
# pass mark), 1e4 = 1%, 1e6 = 0.1%.
SAMPLING_TARGET = 1e3


def _verdict(m):
    if m <= 0:
        return float('inf'), 'EMPTY'
    return 1.0 / np.sqrt(m), ('ok' if m >= SAMPLING_TARGET else 'LOW')


def report_sampling_origins(counts, n_origins, n_independent, span_fs, corr_length_fs):
    """
    M = atoms of that species x time origins — the atoms x samples product.

    Two numbers per species, because they answer different questions.  M_raw
    counts every origin; M_indep counts only origins far enough apart to be
    statistically independent.  Overlapping windows are the severe case in this
    pipeline: with CORR_LENGTH=2000 fs and CORR_INTERVAL=100 fs consecutive
    origins share 95% of the same trajectory, so averaging them does NOT buy
    sqrt(n_origins).  The number of genuinely independent windows is capped at
    span/CORR_LENGTH, a ceiling set by trajectory length that shrinking
    CORR_INTERVAL cannot raise — past it you pay linearly in runtime for
    nothing.  To actually reduce noise, lengthen the trajectory or add atoms.
    """
    print(f"\nSampling achieved (relative error ~ 1/sqrt(M), target M >= {SAMPLING_TARGET:.0e}):")
    print(f"  {n_origins} time origins over {span_fs/1000:.1f} ps; only "
          f"{n_independent} are independent (span / CORR_LENGTH = "
          f"{span_fs/1000:.1f} / {corr_length_fs/1000:.1f} ps)")
    width = max((len(el) for el in counts), default=4)
    worst = None
    for el in sorted(counts):
        m_raw = counts[el] * n_origins
        m_ind = counts[el] * n_independent
        err_raw, _ = _verdict(m_raw)
        err_ind, verdict = _verdict(m_ind)
        print(f"  {el.ljust(width)}  {counts[el]:6d} atoms   "
              f"M_raw = {m_raw:9.3g} ({100*err_raw:5.2f}%)   "
              f"M_indep = {m_ind:9.3g} ({100*err_ind:5.2f}%)  {verdict}")
        if worst is None or m_ind < worst[1]:
            worst = (el, m_ind)
    if worst is not None:
        err, verdict = _verdict(worst[1])
        print(f"  limiting: {worst[0]} at {100*err:.2f}% on independent windows ({verdict})")
        if n_origins > n_independent:
            factor = np.sqrt(n_origins / max(n_independent, 1))
            print(f"  -> M_raw overstates precision by {factor:.1f}x here; "
                  f"CORR_INTERVAL below span/CORR_LENGTH costs runtime and buys nothing")


def _require(name, description):
    """Value of `name`, or None after recording it as missing."""
    value = os.environ.get(name, "")
    if value == "":
        _MISSING.append(f"  {name:<22} {description}")
        return None
    return value

def _check_required(script):
    if _MISSING:
        raise SystemExit(
            f"{script}: required environment variable(s) not set:\n"
            + "\n".join(_MISSING)
            + "\n\nThese determine the numbers this script produces, so it will not "
              "guess them.\nSet them in submit_pipeline.conf (driven by "
              "submit_pipeline.sh / submit_pipeline_local.sh)\nor in the Analysis "
              "parameters block of distribution_run.sh / distribution_submit.slurm."
        )

def _env(name, default):
    """Optional setting: value of `name` if non-empty, else `default`. Used only
    where the default cannot silently distort the result — reading every frame,
    or an algorithm choice documented in this file's header."""
    value = os.environ.get(name, "")
    return value if value != "" else default


# Trajectory sampling. Optional: the defaults are the identity choice (read
# every frame of the trajectory), so leaving them unset cannot distort a result.
N_FRAMES = int(_env("MSD_N_FRAMES", "0"))   # max frames to read; 0 = all
STRIDE   = int(_env("MSD_STRIDE", "1"))     # read every Nth frame

# --- Required. DYNAMICS_DT is the one variable msd.py and vdos.py share: it is
# a property of dynamics.lammpstrj itself, not of either analysis, so both
# scripts must read the same value. The rest are MSD's alone.
_DYNAMICS_DT_ENV   = _require("DYNAMICS_DT", "fs between dumped frames")
_CORR_LENGTH_ENV   = _require("MSD_CORR_LENGTH", "fs; max time lag")
_CORR_INTERVAL_ENV = _require("MSD_CORR_INTERVAL", "fs; spacing between reference frames")
_FIT_FRACTION_ENV  = _require("MSD_FIT_FRACTION", "tail fraction of the window used for the D fit")
_check_required("msd.py")

TIME_UNIT = float(_DYNAMICS_DT_ENV)
# Fraction of the tail of the correlation window used for the diffusion-
# coefficient linear fit (excludes the early ballistic/non-diffusive regime).
FIT_FRACTION = float(_FIT_FRACTION_ENV)

# Output files (set to None to skip writing)
OUTPUT_CSV  = "msd.csv"

# ---- Plot output ---------------------------------------------------------
# Directory every PNG is written into, one curve per file.  A BASENAME; the
# YYYYMMDD_ prefix is added below like every other output.  Empty disables
# plotting and the CSV is still written.  Created by the runner
# (distribution_run.sh / distribution_submit.slurm), not by this script.
MSD_PLOT_DIR   = os.environ.get("MSD_PLOT_DIR", "msd")
# 'analysis' (default) or 'publication'; see PLOT_STYLES after the config block.
PLOT_STYLE     = os.environ.get("MSD_PLOT_STYLE", "analysis")
PLOT_STYLE_KEY = "MSD_PLOT_STYLE"
PLOT_DIR_KEY   = "MSD_PLOT_DIR"
SCRIPT_NAME    = "msd.py"

# Plot appearance
PLOT_DPI = 150

# =============================================================================
# END CONFIGURATION
# =============================================================================

# Prepend today's date (YYYYMMDD_) to every output filename.
def _dated(filename):
    return None if filename is None else f"{date.today():%Y%m%d}_{filename}"

OUTPUT_CSV  = _dated(OUTPUT_CSV)
PLOT_DIR    = _dated(MSD_PLOT_DIR) if MSD_PLOT_DIR else None

# =============================================================================
# Plot output — one quantity per file
# =============================================================================
# Every curve is written as its own PNG into PLOT_DIR; nothing is packed into a
# subplot grid.  A composite survives only where the combination *is* the
# result — Wright's T(r) against the baseline it oscillates about, the mode
# character decomposition, the species overlays — and those are named so they
# read as composites rather than as one more quantity.
#
# PLOT_DIR is created by the runner (distribution_run.sh /
# distribution_submit.slurm), never by this script: a missing directory means
# the run was wired wrong, and quietly creating it would hide that.

PLOT_STYLES = {
    # 'analysis'    — titled and fully labelled, for reading a run.
    # 'publication' — heavy lines and spines, large bold labels, no y ticks.
    #                 This is the styling the separate *_plot.py pass used to
    #                 apply; it is a config choice now, not a second script.
    'analysis':    dict(figsize=(7.0, 4.5), linewidth=1.5, color='C0', spine_lw=0.8,
                        weight='normal', label_fs=12, tick_fs=10,
                        tick_len=4, tick_w=1.0, yticks=True, titles=True),
    'publication': dict(figsize=(4.0, 3.0), linewidth=3.0, color='steelblue', spine_lw=2.0,
                        weight='bold', label_fs=20, tick_fs=14,
                        tick_len=6, tick_w=2.0, yticks=False, titles=False),
}


def plot_style():
    if PLOT_STYLE not in PLOT_STYLES:
        raise ValueError(
            f"Unknown {PLOT_STYLE_KEY}={PLOT_STYLE!r}; use one of {list(PLOT_STYLES)}.")
    return PLOT_STYLES[PLOT_STYLE]


def check_plot_dir():
    """
    True if plots should be written.  PLOT_DIR must already exist — see the note
    above on why this refuses to create it.
    """
    if PLOT_DIR is None:
        return False
    if not os.path.isdir(PLOT_DIR):
        raise SystemExit(
            f"{SCRIPT_NAME}: plot directory {PLOT_DIR!r} does not exist.\n"
            f"  The runner creates it; running this script by hand, create it first:\n"
            f"      mkdir -p {PLOT_DIR}\n"
            f"  Or set {PLOT_DIR_KEY}= (empty) to skip plotting and write only the CSVs.")
    return True


def new_plot():
    """A single-axes figure in the configured style."""
    st = plot_style()
    return plt.subplots(figsize=st['figsize']) + (st,)


def save_plot(fig, ax, name, xlabel, ylabel, title=None, legend=False):
    """
    Finish one figure and write it as PLOT_DIR/<name>.png.

    `name` becomes the filename, so anything a path cannot carry is substituted
    rather than left to mangle the path silently.
    """
    st = plot_style()
    ax.set_xlabel(xlabel, fontsize=st['label_fs'], fontweight=st['weight'])
    ax.set_ylabel(ylabel, fontsize=st['label_fs'], fontweight=st['weight'])
    if title and st['titles']:
        ax.set_title(title)
    ax.tick_params(axis='x', labelsize=st['tick_fs'], length=st['tick_len'], width=st['tick_w'])
    if st['yticks']:
        ax.tick_params(axis='y', labelsize=st['tick_fs'], length=st['tick_len'], width=st['tick_w'])
    else:
        ax.yaxis.set_ticks([])
    if legend:
        ax.legend(fontsize=8 if st['titles'] else 10)
    for spine in ax.spines.values():
        spine.set_linewidth(st['spine_lw'])
    fig.tight_layout()
    path = os.path.join(PLOT_DIR, f"{re.sub(r'[^A-Za-z0-9._+-]', '_', name)}.png")
    fig.savefig(path, dpi=PLOT_DPI)
    plt.close(fig)
    return path

# 1 Angstrom^2/fs = 1e-16 cm^2 / 1e-15 s = 0.1 cm^2/s = 1e4 x(1e-5 cm^2/s)
ANG2_FS_TO_1E5_CM2_S = 1.0e4


def read_lammps_dump(filename, n_frames=0, stride=1):
    """
    Read only what MSD needs: per-frame element labels, wrapped positions,
    and box dimensions. Velocities are not parsed. Returns (elements,
    positions, box_L) where elements is (n_atoms,), positions is
    (n_kept_frames, n_atoms, 3), and box_L is (n_kept_frames, 3) — the box
    edge lengths (Lx, Ly, Lz) used for minimum-image correction.
    """
    elements = None
    pos_frames = []
    box_frames = []
    frame_idx = 0

    with open(filename) as f:
        while True:
            line = f.readline()
            if not line:
                break  # EOF

            # TIMESTEP
            f.readline()

            # NUMBER OF ATOMS
            f.readline()
            n_atoms = int(f.readline().strip())

            # BOX BOUNDS (header + 3 dims)
            f.readline()
            xlo, xhi = map(float, f.readline().split())
            ylo, yhi = map(float, f.readline().split())
            zlo, zhi = map(float, f.readline().split())

            # ATOMS header — parse column positions dynamically
            header = f.readline().split()  # ['ITEM:', 'ATOMS', 'id', 'type', 'element', ...]
            cols = header[2:]
            missing = [c for c in ('element', 'x', 'y', 'z') if c not in cols]
            if missing:
                raise ValueError(
                    f"Dump is missing required column(s) {missing} — MSD needs "
                    f"wrapped positions (dump ... x y z) in addition to element."
                )
            col_element = cols.index('element')
            col_x       = cols.index('x')
            col_y       = cols.index('y')
            col_z       = cols.index('z')

            keep = (frame_idx % stride == 0)
            if keep:
                frame_elements = []
                positions = np.empty((n_atoms, 3))
                for i in range(n_atoms):
                    parts = f.readline().split()
                    frame_elements.append(parts[col_element])
                    positions[i, 0] = float(parts[col_x])
                    positions[i, 1] = float(parts[col_y])
                    positions[i, 2] = float(parts[col_z])

                frame_elements = np.array(frame_elements)
                if elements is None:
                    elements = frame_elements
                elif not np.array_equal(elements, frame_elements):
                    raise ValueError(
                        "Atom element order changed between frames — MSD requires "
                        "consistent atom identity/order across the trajectory "
                        "(check 'dump_modify ... sort id' is set in the LAMMPS input)."
                    )
                pos_frames.append(positions)
                box_frames.append((xhi - xlo, yhi - ylo, zhi - zlo))
            else:
                for _ in range(n_atoms):
                    f.readline()

            frame_idx += 1
            if n_frames and len(pos_frames) >= n_frames:
                break

    if not pos_frames:
        raise ValueError(f"No frames read from {filename} (empty file or STRIDE too large).")

    return elements, np.stack(pos_frames, axis=0), np.array(box_frames)


def compute_msd_multi_origin(positions, box_L, elements, corr_length_frames, corr_interval_frames):
    """
    Per-element MSD via multiple time origins, minimum-image corrected
    without ever unwrapping coordinates globally — see METHOD in the module
    docstring. A reference frame every corr_interval_frames is compared
    against every frame up to corr_length_frames later; the whole
    (lag x atom) displacement array for one reference is computed in a
    single batched op, so only the (few) references are looped in Python.

    Returns (msd, n_refs, counts) where msd is {element: (corr_length_frames,)
    array} in Å², and counts is {element: n_atoms_of_that_element}.
    """
    n_frames = positions.shape[0]
    unique_els = sorted(set(elements.tolist()))
    masks = {el: (elements == el) for el in unique_els}
    counts = {el: int(masks[el].sum()) for el in unique_els}

    msd_sum = {el: np.zeros(corr_length_frames) for el in unique_els}

    refs = range(0, n_frames - corr_length_frames + 1, corr_interval_frames)
    n_refs = 0
    for r in refs:
        pos_ref = positions[r]                                       # (n_atoms, 3)
        segment = positions[r:r + corr_length_frames]                 # (corr_length_frames, n_atoms, 3)
        disp    = segment - pos_ref[None, :, :]
        L       = box_L[r:r + corr_length_frames][:, None, :]         # (corr_length_frames, 1, 3)
        disp   -= L * np.round(disp / L)                              # minimum-image correction
        dr2     = (disp ** 2).sum(axis=2)                             # (corr_length_frames, n_atoms)
        for el in unique_els:
            msd_sum[el] += dr2[:, masks[el]].sum(axis=1)
        n_refs += 1

    if n_refs == 0:
        raise ValueError(
            "No valid MSD reference frames: CORR_LENGTH is longer than the trajectory "
            "— shorten CORR_LENGTH or provide more frames."
        )

    msd = {el: msd_sum[el] / (counts[el] * n_refs) for el in unique_els}
    return msd, n_refs, counts


def combine_total(msd, counts):
    """Mole-fraction-weighted average across elements — MSD is a plain
    average quantity, so the total is just the atom-count-weighted mean."""
    n_total = sum(counts.values())
    total = np.zeros_like(next(iter(msd.values())))
    for el, curve in msd.items():
        total += (counts[el] / n_total) * curve
    results = dict(msd)
    results['total'] = total
    return results


def estimate_diffusion_coefficients(results, time_fs, fit_fraction):
    """
    Self-diffusion coefficient per curve via the 3D Einstein relation
    MSD(t) = 6*D*t, linear-fit over the last fit_fraction of the window
    (excludes the early ballistic/non-diffusive regime).

    Returns {label: D in units of 1e-5 cm^2/s} — the standard unit for
    reporting liquid/solid self-diffusion coefficients in the literature
    (e.g. water at room temperature: D ~ 2.3 x1e-5 cm^2/s).
    """
    n = len(time_fs)
    start = int(n * (1 - fit_fraction))
    D = {}
    for label, curve in results.items():
        slope, _intercept = np.polyfit(time_fs[start:], curve[start:], 1)
        D_ang2_per_fs = slope / 6.0
        D[label] = D_ang2_per_fs * ANG2_FS_TO_1E5_CM2_S
    return D


def save_csv(results, time_fs, filename):
    """Save results dict {element_or_'total': array} to CSV, one column per label."""
    order = [el for el in results if el != 'total'] + ['total']
    header_parts = ['time_fs']
    columns = [time_fs]
    for label in order:
        header_parts.append('MSD_total' if label == 'total' else f'MSD_{label}')
        columns.append(results[label])
    header = ','.join(header_parts)
    data = np.column_stack(columns)
    np.savetxt(filename, data, delimiter=',', header=header, comments='', fmt='%.6f')
    print(f"Data table saved to {filename}")


def plot_msd(results, time_fs):
    """
    One PNG per species, plus the species overlay.

    The overlay is kept because the comparison between species IS the result
    here — a light element outrunning a heavy one is the thing being read —
    whereas the individual curves are what a single number gets measured off.
    """
    for label, curve in results.items():
        fig, ax, st = new_plot()
        ax.plot(time_fs, curve, color=st['color'], linewidth=st['linewidth'])
        save_plot(fig, ax, label, 't (fs)', 'MSD (Å²)', title=f'MSD  {label}')

    fig, ax, st = new_plot()
    for label, curve in results.items():
        ax.plot(time_fs, curve, label=label, linewidth=1.5 if label == 'total' else 1.0)
    save_plot(fig, ax, 'all_species', 't (fs)', 'MSD (Å²)',
              title='MSD, all species', legend=True)
    print(f"{len(results) + 1} MSD plot(s) written to {PLOT_DIR}/")


if __name__ == '__main__':
    # A missing plot directory is a wiring error: catch it before the parse.
    plots = check_plot_dir()
    print(f"Reading trajectory: {DUMP_FILE}")
    print(f"  N_FRAMES={N_FRAMES or 'all'}, STRIDE={STRIDE}, TIME_UNIT={TIME_UNIT} fs")

    import time as _time
    t0 = _time.time()
    elements, positions, box_L = read_lammps_dump(DUMP_FILE, n_frames=N_FRAMES, stride=STRIDE)
    t1 = _time.time()
    n_frames_read = positions.shape[0]
    print(f"  Frames read: {n_frames_read}, atoms: {positions.shape[1]} ({t1 - t0:.2f}s)")

    unique_els = sorted(set(elements.tolist()))
    print(f"  Elements: {unique_els}")

    CORR_LENGTH = float(_CORR_LENGTH_ENV)
    CORR_INTERVAL = float(_CORR_INTERVAL_ENV)
    print(f"  CORR_LENGTH={CORR_LENGTH:.1f} fs, CORR_INTERVAL={CORR_INTERVAL:.1f} fs")

    corr_length_frames = min(round(CORR_LENGTH / TIME_UNIT), n_frames_read)
    corr_interval_frames = max(1, round(CORR_INTERVAL / TIME_UNIT))
    if corr_length_frames < round(CORR_LENGTH / TIME_UNIT):
        print(f"  Warning: CORR_LENGTH clipped to the trajectory length "
              f"({corr_length_frames} frames, ~{corr_length_frames * TIME_UNIT:.1f} fs).")

    print("Computing MSD (multiple time origins, minimum-image corrected)...")
    t2 = _time.time()
    msd, n_refs, counts = compute_msd_multi_origin(
        positions, box_L, elements, corr_length_frames, corr_interval_frames
    )
    results = combine_total(msd, counts)
    time_fs = np.arange(corr_length_frames) * TIME_UNIT
    t3 = _time.time()
    print(f"  Reference frames used: {n_refs}")
    _span_fs = n_frames_read * TIME_UNIT
    _n_indep = max(1, int(_span_fs // CORR_LENGTH))
    report_sampling_origins(counts, n_refs, _n_indep, _span_fs, CORR_LENGTH)
    print(f"  Compute: {t3 - t2:.2f}s")
    print(f"  Total: {t3 - t0:.2f}s")

    D = estimate_diffusion_coefficients(results, time_fs, FIT_FRACTION)
    print(f"  Self-diffusion coefficient (Einstein relation, "
          f"last {FIT_FRACTION:.0%} of window):")
    for label, D_val in D.items():
        print(f"    D({label}) = {D_val:.6f} x1e-5 cm²/s")

    if OUTPUT_CSV is not None:
        save_csv(results, time_fs, OUTPUT_CSV)
    if plots:
        plot_msd(results, time_fs)
