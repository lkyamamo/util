"""
dsf.py — Dynamic and Static Structure Factor from LAMMPS dump trajectories.

QUICK START
-----------
1. Set DUMP_FILE to your LAMMPS custom dump trajectory (element column required).
   Defaults to dynamics.lammpstrj (env var DYNAMICS_TRAJ) — the higher-frequency
   trajectory dumped alongside dump.lammpstrj specifically for dsf.py/vdos.py;
   see OH-therm.input/b-SiO-therm.input's second dump block.
2. Set DT to the time between consecutive dumped frames in femtoseconds.
3. Set N_FRAMES and WINDOW_SIZE.
4. Run:  python dsf.py

OUTPUT
------
COMPUTE_STATIC=True
  sq.csv    — q (Å⁻¹), partial S_AB(q), total S(q), neutron-weighted S(q)
  sq.png    — line plot of all S(q) curves

COMPUTE_DYNAMIC=True
  dsf.csv   — q (Å⁻¹), ω (THz), partial S_AB(q,ω), total, neutron-weighted
  dsf.png   — 2D heatmap S(q,ω) for total and neutron-weighted

DEPENDENCIES
------------
  pip install dynasor matplotlib
  pip install icc_rt          # optional: 5–10× numba speedup

PARALLELIZATION
---------------
  dynasor uses numba internally; set N_THREADS here or OMP_NUM_THREADS in the
  shell (the shell value takes precedence when already set, e.g. from SLURM).
"""

import os
from datetime import date

# =============================================================================
# CONFIGURATION — edit these variables between runs
# =============================================================================

DUMP_FILE       = os.environ.get("DYNAMICS_TRAJ", "dynamics.lammpstrj")

# Trajectory sampling.
# N_FRAMES is frame_stop, an index into the dump — NOT a count.  Frames actually
# used is N_FRAMES / STRIDE, so raising STRIDE uses fewer frames over the same
# span rather than the same number over a longer span.  Skipped frames are still
# parsed (measured: iterating a 100-frame span costs the same at STRIDE=1, 10 and
# 50), so set N_FRAMES to the whole trajectory and pick STRIDE for the frame count
# you can afford — spanning more time is free, computing more frames is not.
# Defaults: 30000 dumped frames at DT=2 fs = 60 ps, sampled every 1.2 ps.
N_FRAMES        = int(os.environ.get("N_FRAMES", "30000"))  # frame_stop in Trajectory
STRIDE          = int(os.environ.get("STRIDE", "600"))      # read every Nth frame (frame_step)

# Threading — 0 = use all available cores
# Only applied when OMP_NUM_THREADS is not already set in the environment.
N_THREADS       = 0

# Time axis
DT              = float(os.environ.get("DT", "2.0"))        # fs between consecutive dumped frames

# q-space, static S(q).
# Q_MAX=20 Å⁻¹ is the range needed to Fourier transform S(q) into G(r) without bad
# termination ripples (ripple period 2π/Q_MAX = 0.31 Å), and matches the range of
# neutron diffraction measurements on vitreous silica.  It costs 10.6M q-vectors on
# a ~43 Å cell versus 2.3M at Q_MAX=12, and cost is linear in N_q × frames.
# N_Q_BINS must stay ≤ Q_MAX / (2π/L) = 136 here, or low-q bins come back empty.
Q_MAX           = float(os.environ.get("Q_MAX", "20.0"))    # Å⁻¹, passed to get_spherical_qpoints
N_Q_BINS        = int(os.environ.get("N_Q_BINS", "130"))    # radial q-bins after spherical averaging

# --- Dynamic S(q,ω) parameters — NOT YET TUNED, see DYNAMIC NOTES below --------
WINDOW_SIZE     = int(os.environ.get("WINDOW_SIZE", "2000"))  # time lags; Δν = 1/(2 × WINDOW_SIZE × DT × STRIDE)
WINDOW_STEP     = int(os.environ.get("WINDOW_STEP", "1"))     # frames between window origins; 1 = maximum averaging
Q_MAX_DYN       = float(os.environ.get("Q_MAX_DYN", "4.0"))   # Å⁻¹, separate from Q_MAX — see notes
N_Q_BINS_DYN    = int(os.environ.get("N_Q_BINS_DYN", "25"))   # radial q-bins, dynamic
MAX_Q_POINTS_DYN = int(os.environ.get("MAX_Q_POINTS_DYN", "25000"))  # prune target; 0 = no pruning

# What to compute
COMPUTE_STATIC  = True      # S(q)
COMPUTE_DYNAMIC = False     # S(q,ω) and F(q,t) — off until the notes below are worked through
COMPUTE_SELF    = False     # incoherent/self part — can be slow

# Output files (set to None to skip writing)
OUTPUT_SQ_CSV   = "sq.csv"
OUTPUT_SQ_PLOT  = "sq.png"
OUTPUT_DSF_CSV  = "dsf.csv"
OUTPUT_DSF_PLOT = "dsf.png"

# Plot layout
PLOT_NCOLS = 2
PLOT_DPI   = 150

# =============================================================================
# END CONFIGURATION
# =============================================================================

# =============================================================================
# DYNAMIC NOTES — open items for when COMPUTE_DYNAMIC is turned back on
# =============================================================================
#
# The dynamic path runs and writes correct output, but its parameters have not
# been validated against physics.  Findings from benchmarking a 5184-atom,
# L = 42.7934 Å cell (q spacing 2π/L = 0.1468 Å⁻¹) on a laptop:
#
# MEMORY is the binding constraint, not time.  The raw DynamicSample holds
# 2 × (n_pairs + 1) arrays of shape (N_q, WINDOW_SIZE+1) — Fqt and Sqw, partials
# plus total — before spherical averaging:
#
#     memory ≈ 2 × (n_pairs + 1) × N_q × (WINDOW_SIZE + 1) × 8 bytes
#
#     Q_MAX_DYN   max_points   N_q          memory at WINDOW_SIZE=2000
#     20          none         10,587,961   2373 GB   <- the old default; would OOM
#     12          none          2,286,779    512 GB
#      4          none             84,823     19 GB
#      4          25,000           25,010    5.6 GB
#
# This is why Q_MAX_DYN and MAX_Q_POINTS_DYN exist separately from Q_MAX: the
# static run wants a large q range, the dynamic run cannot afford one.  Pruning
# keeps the low-q mesh intact and thins only above a cutoff dynasor chooses
# (max_points=25000 at Q_MAX_DYN=4 prunes only |q| > 1.44 Å⁻¹).
#
# TIME is dominated by computing ρ(q,t) once per frame, ≈ 7 µs × N_frames × N_q
# (laptop, all cores; a 32-thread node should be ~2-3× faster).  WINDOW_SIZE and
# WINDOW_STEP are nearly free — measured 20.5 s at W=30, 21.0 s at W=60, and
# 20.5 s at WINDOW_STEP=30 for the same 140 frames.  So spend on WINDOW_SIZE and
# leave WINDOW_STEP=1; economize on N_q instead.
#
# PHYSICS still to check before trusting the output:
#   - WINDOW_SIZE must cover several periods of the slowest mode of interest.
#     The lowest accessible q is 2π/L = 0.1468 Å⁻¹; a longitudinal acoustic mode
#     there sits near v_L·q/2π ≈ 1.4 THz for v_L ≈ 5900 m/s, i.e. a ~0.7 ps
#     period.  Five to ten periods needs a 3.5-7 ps window, so WINDOW_SIZE ≈
#     2000-3500 at DT = 2 fs.  WINDOW_SIZE=2000 (4 ps, Δν = 0.125 THz) is the
#     starting guess baked in above and has NOT been verified against a dispersion.
#   - STRIDE must stay 1.  Nyquist is 1/(2 × DT × STRIDE); at DT = 2 fs that is
#     250 THz, comfortably above the O-H stretch (~3700 cm⁻¹ ≈ 111 THz), but
#     STRIDE=4 would drop it to 62.5 THz and alias the O-H band down into the
#     Si-O region.
#   - N_FRAMES should be several times WINDOW_SIZE so many windows are averaged.
#
# UNVERIFIED: the COMPUTE_SELF=True branch.  Its column names match dynasor's
# correlation_functions.py (Sqw_incoh_<species>), but dynasor's incoherent kernel
# crashes numba locally ("workqueue threading layer is terminating") at any
# thread count, so that path has never actually been executed here.
# =============================================================================

# Prepend today's date (YYYYMMDD_) to every output filename.
def _dated(filename):
    return None if filename is None else f"{date.today():%Y%m%d}_{filename}"

OUTPUT_SQ_CSV   = _dated(OUTPUT_SQ_CSV)
OUTPUT_SQ_PLOT  = _dated(OUTPUT_SQ_PLOT)
OUTPUT_DSF_CSV  = _dated(OUTPUT_DSF_CSV)
OUTPUT_DSF_PLOT = _dated(OUTPUT_DSF_PLOT)

# Must be set before importing dynasor/numba — numba reads thread count at JIT time.
# os.environ.setdefault only writes when OMP_NUM_THREADS is not already present
# (e.g. already set by SLURM via --cpus-per-task).
if N_THREADS > 0:
    os.environ.setdefault('OMP_NUM_THREADS', str(N_THREADS))

import warnings

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# MDAnalysis' LAMMPSDUMP topology parser emits both of these for any dump that has
# an `element` column but no `type`/`mass` column — which is exactly the layout
# dsf.py wants.  Neither affects the results: species come from the element column
# via _element_indices() below, and structure factors are weighted by neutron
# scattering lengths, never by mass.  dynasor filters its own equivalent
# ('Guessed all Masses to 1.0') for the same reason.
warnings.filterwarnings('ignore', category=UserWarning, message='No mass column found')
warnings.filterwarnings('ignore', category=UserWarning, message='Set all atom types to 1')

from dynasor import (
    Trajectory,
    compute_static_structure_factors,
    compute_dynamic_structure_factors,
)
from dynasor.qpoints import get_spherical_qpoints
from dynasor.post_processing import (
    NeutronScatteringLengths,
    get_weighted_sample,
    get_spherically_averaged_sample_binned,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _element_indices():
    """Map element symbol -> atom indices, read from the first frame's element column.

    MDAnalysis' LAMMPSDUMP parser sets every atom's *type* to 1 when the dump has an
    `element` column but no `type` column ("Set all atom types to 1"), so neither
    dynasor's default (one lumped species 'X') nor atomic_indices='read_from_trajectory'
    gives species-resolved partials.  Parsing the element column here is what makes
    Sq_Si_O etc. and the neutron scattering length weighting work.
    """
    ids, elements = [], []
    with open(DUMP_FILE) as fh:
        for line in fh:
            if line.startswith('ITEM: ATOMS'):
                columns = line.split()[2:]
                break
        else:
            raise ValueError(f'No "ITEM: ATOMS" header found in {DUMP_FILE}')
        for required in ('id', 'element'):
            if required not in columns:
                raise ValueError(f'{DUMP_FILE} has no `{required}` column (found: {columns}); '
                                 'dsf.py needs `id` and element symbols, not numeric types')
        id_col, el_col = columns.index('id'), columns.index('element')
        for line in fh:
            if line.startswith('ITEM:'):     # start of the next frame
                break
            fields = line.split()
            ids.append(int(fields[id_col]))
            elements.append(fields[el_col])

    # MDAnalysis returns atoms sorted by id, which need not be the dump's line order
    # (LAMMPS only sorts with `dump_modify sort id`), so index by id rank, not file order.
    elements = np.array(elements)[np.argsort(ids)]
    return {str(el): np.flatnonzero(elements == el) for el in np.unique(elements)}


def _make_traj():
    """Return a fresh Trajectory (Trajectory is an exhaustible iterable)."""
    return Trajectory(
        DUMP_FILE,
        trajectory_format='LAMMPSDUMP',
        atomic_indices=_element_indices(),
        frame_stop=N_FRAMES,
        frame_step=STRIDE,
    )


def _n_atoms(sample):
    return sum(sample.particle_counts.values())


def _try_neutron_weights(sample):
    """Apply NeutronScatteringLengths weighting; return None and warn on failure."""
    try:
        nsl = NeutronScatteringLengths(sample.atom_types)
        return get_weighted_sample(sample, nsl)
    except Exception as exc:
        print(f"  Warning: neutron weighting skipped — {exc}")
        return None


# ---------------------------------------------------------------------------
# Static S(q) — I/O
# ---------------------------------------------------------------------------

def save_csv_sq(sample, sample_neutron, filename):
    """Write S(q) partials + total + neutron-weighted total to CSV.

    dynasor already divides by the atom count (correlation_functions.py:
    ``Sq = 1 / traj.n_atoms * ...``), so nothing is normalized again here.
    Sanity check: for uncorrelated positions S_total(q) -> 1 at all q > 0.
    """
    q = sample.q_norms              # (N_Q_BINS,) Å⁻¹

    header_parts = ['q_Ang-1']
    columns      = [q]

    for (a, b) in sample.pairs:
        header_parts.append(f'Sq_{a}_{b}')
        columns.append(sample[f'Sq_{a}_{b}'])

    header_parts.append('Sq_total')
    columns.append(sample.Sq)

    if sample_neutron is not None:
        header_parts.append('Sq_neutron')
        columns.append(sample_neutron.Sq)

    header = ','.join(header_parts)
    np.savetxt(filename, np.column_stack(columns),
               delimiter=',', header=header, comments='', fmt='%.6f')
    print(f"  S(q) data saved to {filename}")


def plot_sq(sample, sample_neutron, filename):
    """Line plot: one panel per partial + total + neutron-weighted."""
    q = sample.q_norms

    curves = {f'{a}-{b}': sample[f'Sq_{a}_{b}'] for (a, b) in sample.pairs}
    curves['total'] = sample.Sq
    if sample_neutron is not None:
        curves['neutron'] = sample_neutron.Sq

    n     = len(curves)
    ncols = PLOT_NCOLS
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols,
                              figsize=(6 * ncols, 4 * nrows), squeeze=False)
    axes = axes.flatten()

    for ax, (label, sq) in zip(axes, curves.items()):
        ax.plot(q, sq)
        ax.set_xlabel('q (Å⁻¹)')
        ax.set_ylabel('S(q)')
        ax.set_title(label)

    for ax in axes[n:]:
        ax.set_visible(False)

    fig.tight_layout()
    fig.savefig(filename, dpi=PLOT_DPI)
    plt.close(fig)
    print(f"  S(q) plot saved to {filename}")


# ---------------------------------------------------------------------------
# Dynamic S(q,ω) — I/O
# ---------------------------------------------------------------------------

def _omega_THz(sample):
    """Convert sample.omega (rad/fs) to THz.  1 rad/fs = 1e3/(2π) THz."""
    return sample.omega * 1e3 / (2.0 * np.pi)


def save_csv_dsf(sample, sample_neutron, filename):
    """Write S(q,ω) as a long-format CSV: one row per (q, ω) pair."""
    q         = sample.q_norms          # (N_Q_BINS,)
    omega     = _omega_THz(sample)      # (N_omega,)

    # meshgrid with indexing='ij': QQ[i,j] = q[i], WW[i,j] = omega[j]
    QQ, WW = np.meshgrid(q, omega, indexing='ij')   # (N_Q_BINS, N_omega)

    header_parts = ['q_Ang-1', 'freq_THz']
    flat_cols    = [QQ.ravel(), WW.ravel()]

    # dynasor names the dynamic correlation functions Sqw_coh[_A_B] (and Sqw_incoh_A
    # when calculate_incoherent=True) — unlike the static sample's plain Sq_A_B.
    # Already per-atom normalized by dynasor; see save_csv_sq.
    for (a, b) in sample.pairs:
        header_parts.append(f'Sqw_{a}_{b}')
        flat_cols.append(sample[f'Sqw_coh_{a}_{b}'].ravel())

    header_parts.append('Sqw_total')
    flat_cols.append(sample.Sqw_coh.ravel())

    if sample_neutron is not None:
        header_parts.append('Sqw_neutron')
        flat_cols.append(sample_neutron.Sqw_coh.ravel())

    if COMPUTE_SELF:
        for a in sample.atom_types:
            header_parts.append(f'Sqw_self_{a}')
            flat_cols.append(sample[f'Sqw_incoh_{a}'].ravel())
        header_parts.append('Sqw_self_total')
        flat_cols.append(sample.Sqw_incoh.ravel())

    header = ','.join(header_parts)
    np.savetxt(filename, np.column_stack(flat_cols),
               delimiter=',', header=header, comments='', fmt='%.6f')
    print(f"  S(q,ω) data saved to {filename}")


def plot_dsf(sample, sample_neutron, filename):
    """2D imshow heatmaps: total S(q,ω) and neutron-weighted S(q,ω)."""
    q     = sample.q_norms          # (N_Q_BINS,)
    omega = _omega_THz(sample)      # (N_omega,)

    # sample.Sqw_coh shape: (N_Q_BINS, N_omega) — q on axis-0, omega on axis-1
    panels = {'S(q,ω) total': sample.Sqw_coh}
    if sample_neutron is not None:
        panels['S(q,ω) neutron'] = sample_neutron.Sqw_coh

    n_panels = len(panels)
    fig, axes = plt.subplots(1, n_panels,
                              figsize=(7 * n_panels, 5), squeeze=False)
    axes = axes.flatten()

    extent = [omega[0], omega[-1], q[0], q[-1]]

    for ax, (title, Sqw) in zip(axes, panels.items()):
        vmax = np.nanpercentile(Sqw, 99)
        im = ax.imshow(
            Sqw,
            origin='lower',
            aspect='auto',
            extent=extent,
            vmin=0,
            vmax=vmax,
            cmap='inferno',
        )
        ax.set_xlabel('ω (THz)')
        ax.set_ylabel('q (Å⁻¹)')
        ax.set_title(title)
        fig.colorbar(im, ax=ax, label='S(q,ω)')

    fig.tight_layout()
    fig.savefig(filename, dpi=PLOT_DPI)
    plt.close(fig)
    print(f"  S(q,ω) plot saved to {filename}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    print(f"Reading trajectory: {DUMP_FILE}")
    print(f"  N_FRAMES={N_FRAMES}, STRIDE={STRIDE}, DT={DT} fs")
    print(f"  OMP_NUM_THREADS={os.environ.get('OMP_NUM_THREADS', '(all cores)')}")

    # Report the species assignment before any long-running work, so an interrupted
    # run still shows what the element column was parsed into.
    for element, indices in sorted(_element_indices().items()):
        print(f"  Species {element}: {len(indices)} atoms")

    # Build q-points from the cell of the first frame.
    # get_spherical_qpoints only reads .cell — does not consume trajectory frames.
    traj_cell = _make_traj()
    q_min = 2.0 * np.pi / np.linalg.norm(traj_cell.cell, axis=1).max()
    print(f"  q spacing 2π/L = {q_min:.4f} Å⁻¹ (lowest accessible q)")

    def _bin_warning(label, num_bins, q_max):
        """Bins narrower than the q spacing come back empty (dynasor: 'No q-points for bin')."""
        if num_bins > q_max / q_min:
            print(f"  Warning: {label}={num_bins} gives bins narrower than {q_min:.4f} Å⁻¹; "
                  f"low-q bins will be empty. Use ≲ {int(q_max / q_min)}.")

    # ------------------------------------------------------------------
    # Static S(q)
    # ------------------------------------------------------------------
    if COMPUTE_STATIC:
        q_points = get_spherical_qpoints(traj_cell.cell, q_max=Q_MAX)
        print(f"  static q-points:  {len(q_points):,} vectors, |q| ≤ {Q_MAX} Å⁻¹")
        _bin_warning('N_Q_BINS', N_Q_BINS, Q_MAX)

        print("Computing S(q)...")
        static_raw = compute_static_structure_factors(
            _make_traj(), q_points
        )
        print(f"  Atom types: {static_raw.atom_types}")
        print(f"  N atoms:    {_n_atoms(static_raw)}")

        print("  Spherically averaging S(q)...")
        static_avg = get_spherically_averaged_sample_binned(
            static_raw, num_q_bins=N_Q_BINS
        )

        print("  Applying neutron scattering length weights...")
        static_neutron = _try_neutron_weights(static_avg)

        if OUTPUT_SQ_CSV:
            save_csv_sq(static_avg, static_neutron, OUTPUT_SQ_CSV)
        if OUTPUT_SQ_PLOT:
            plot_sq(static_avg, static_neutron, OUTPUT_SQ_PLOT)

    # ------------------------------------------------------------------
    # Dynamic S(q,ω)
    # ------------------------------------------------------------------
    if COMPUTE_DYNAMIC:
        q_points_dyn = get_spherical_qpoints(
            traj_cell.cell, q_max=Q_MAX_DYN,
            max_points=MAX_Q_POINTS_DYN if MAX_Q_POINTS_DYN > 0 else None,
        )
        # Stored arrays: 6 partials + total, for Fqt and Sqw alike -> 2 * (n_pairs + 1).
        n_types = len(traj_cell.atom_types)
        n_pairs = n_types * (n_types + 1) // 2
        gb = 2 * (n_pairs + 1) * len(q_points_dyn) * (WINDOW_SIZE + 1) * 8 / 1e9
        print(f"  dynamic q-points: {len(q_points_dyn):,} vectors, |q| ≤ {Q_MAX_DYN} Å⁻¹")
        print(f"  dynamic sample will need ~{gb:.1f} GB of memory")
        _bin_warning('N_Q_BINS_DYN', N_Q_BINS_DYN, Q_MAX_DYN)

        delta_nu = 1e3 / (2 * WINDOW_SIZE * DT * STRIDE)   # THz
        nyquist  = 1e3 / (2 * DT * STRIDE)                  # THz
        print(f"Computing S(q,ω)  [WINDOW_SIZE={WINDOW_SIZE}, WINDOW_STEP={WINDOW_STEP}, "
              f"DT={DT} fs]...")
        print(f"  time window {WINDOW_SIZE * DT * STRIDE / 1e3:.2f} ps, "
              f"Δν = {delta_nu:.3f} THz, Nyquist ν = {nyquist:.1f} THz")
        dynamic_raw = compute_dynamic_structure_factors(
            _make_traj(), q_points_dyn,
            dt=DT,
            window_size=WINDOW_SIZE,
            window_step=WINDOW_STEP,
            calculate_incoherent=COMPUTE_SELF,
        )

        print("  Spherically averaging S(q,ω)...")
        dynamic_avg = get_spherically_averaged_sample_binned(
            dynamic_raw, num_q_bins=N_Q_BINS_DYN
        )

        print("  Applying neutron scattering length weights...")
        dynamic_neutron = _try_neutron_weights(dynamic_avg)

        if OUTPUT_DSF_CSV:
            save_csv_dsf(dynamic_avg, dynamic_neutron, OUTPUT_DSF_CSV)
        if OUTPUT_DSF_PLOT:
            plot_dsf(dynamic_avg, dynamic_neutron, OUTPUT_DSF_PLOT)

    print("Done.")
