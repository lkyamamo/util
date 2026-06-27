"""
dsf.py — Dynamic and Static Structure Factor from LAMMPS dump trajectories.

QUICK START
-----------
1. Set DUMP_FILE to your LAMMPS custom dump trajectory (element column required).
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

# =============================================================================
# CONFIGURATION — edit these variables between runs
# =============================================================================

DUMP_FILE       = os.environ.get("TRAJ", "dump.lammpstrj")

# Trajectory sampling
N_FRAMES        = 500       # max frames to read (frame_stop in Trajectory)
STRIDE          = 1         # read every Nth frame (frame_step in Trajectory)

# Threading — 0 = use all available cores
# Only applied when OMP_NUM_THREADS is not already set in the environment.
N_THREADS       = 0

# Time axis
DT              = 1.0       # fs between consecutive dumped frames

# Dynamic S(q,ω) parameters
WINDOW_SIZE     = 500       # number of time lags; frequency resolution Δω ~ 1/(WINDOW_SIZE × DT)

# q-space
Q_MAX           = 20.0      # Å⁻¹, passed to get_spherical_qpoints
N_Q_BINS        = 200       # radial q-bins after spherical averaging

# What to compute
COMPUTE_STATIC  = True      # S(q)
COMPUTE_DYNAMIC = True      # S(q,ω) and F(q,t)
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

# Must be set before importing dynasor/numba — numba reads thread count at JIT time.
# os.environ.setdefault only writes when OMP_NUM_THREADS is not already present
# (e.g. already set by SLURM via --cpus-per-task).
if N_THREADS > 0:
    os.environ.setdefault('OMP_NUM_THREADS', str(N_THREADS))

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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

def _make_traj():
    """Return a fresh Trajectory (Trajectory is an exhaustible iterable)."""
    return Trajectory(
        DUMP_FILE,
        trajectory_format='LAMMPSDUMP',
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
    """Write S(q) partials + total + neutron-weighted total to CSV."""
    N   = _n_atoms(sample)
    q   = sample.q_norms            # (N_Q_BINS,) Å⁻¹

    header_parts = ['q_Ang-1']
    columns      = [q]

    for (a, b) in sample.pairs:
        header_parts.append(f'Sq_{a}_{b}')
        columns.append(sample[f'Sq_{a}_{b}'] / N)

    header_parts.append('Sq_total')
    columns.append(sample.Sq / N)

    if sample_neutron is not None:
        header_parts.append('Sq_neutron')
        columns.append(sample_neutron.Sq / N)

    header = ','.join(header_parts)
    np.savetxt(filename, np.column_stack(columns),
               delimiter=',', header=header, comments='', fmt='%.6f')
    print(f"  S(q) data saved to {filename}")


def plot_sq(sample, sample_neutron, filename):
    """Line plot: one panel per partial + total + neutron-weighted."""
    N = _n_atoms(sample)
    q = sample.q_norms

    curves = {f'{a}-{b}': sample[f'Sq_{a}_{b}'] / N for (a, b) in sample.pairs}
    curves['total'] = sample.Sq / N
    if sample_neutron is not None:
        curves['neutron'] = sample_neutron.Sq / N

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
    N         = _n_atoms(sample)
    q         = sample.q_norms          # (N_Q_BINS,)
    omega     = _omega_THz(sample)      # (N_omega,)

    # meshgrid with indexing='ij': QQ[i,j] = q[i], WW[i,j] = omega[j]
    QQ, WW = np.meshgrid(q, omega, indexing='ij')   # (N_Q_BINS, N_omega)

    header_parts = ['q_Ang-1', 'freq_THz']
    flat_cols    = [QQ.ravel(), WW.ravel()]

    for (a, b) in sample.pairs:
        header_parts.append(f'Sqw_{a}_{b}')
        flat_cols.append((sample[f'Sqw_{a}_{b}'] / N).ravel())

    header_parts.append('Sqw_total')
    flat_cols.append((sample.Sqw / N).ravel())

    if sample_neutron is not None:
        header_parts.append('Sqw_neutron')
        flat_cols.append((sample_neutron.Sqw / N).ravel())

    header = ','.join(header_parts)
    np.savetxt(filename, np.column_stack(flat_cols),
               delimiter=',', header=header, comments='', fmt='%.6f')
    print(f"  S(q,ω) data saved to {filename}")


def plot_dsf(sample, sample_neutron, filename):
    """2D imshow heatmaps: total S(q,ω) and neutron-weighted S(q,ω)."""
    N     = _n_atoms(sample)
    q     = sample.q_norms          # (N_Q_BINS,)
    omega = _omega_THz(sample)      # (N_omega,)

    # sample.Sqw shape: (N_Q_BINS, N_omega) — q on axis-0, omega on axis-1
    panels = {'S(q,ω) total': sample.Sqw / N}
    if sample_neutron is not None:
        panels['S(q,ω) neutron'] = sample_neutron.Sqw / N

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

    # Build q-points from the cell of the first frame.
    # get_spherical_qpoints only reads .cell — does not consume trajectory frames.
    traj_cell = _make_traj()
    q_points  = get_spherical_qpoints(traj_cell.cell, q_max=Q_MAX)
    print(f"  q-points: {len(q_points)} vectors, |q| ≤ {Q_MAX} Å⁻¹")

    # ------------------------------------------------------------------
    # Static S(q)
    # ------------------------------------------------------------------
    if COMPUTE_STATIC:
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
        print(f"Computing S(q,ω)  [WINDOW_SIZE={WINDOW_SIZE}, DT={DT} fs]...")
        dynamic_raw = compute_dynamic_structure_factors(
            _make_traj(), q_points,
            dt=DT,
            window_size=WINDOW_SIZE,
            calculate_incoherent=COMPUTE_SELF,
        )

        print("  Spherically averaging S(q,ω)...")
        dynamic_avg = get_spherically_averaged_sample_binned(
            dynamic_raw, num_q_bins=N_Q_BINS
        )

        print("  Applying neutron scattering length weights...")
        dynamic_neutron = _try_neutron_weights(dynamic_avg)

        if OUTPUT_DSF_CSV:
            save_csv_dsf(dynamic_avg, dynamic_neutron, OUTPUT_DSF_CSV)
        if OUTPUT_DSF_PLOT:
            plot_dsf(dynamic_avg, dynamic_neutron, OUTPUT_DSF_PLOT)

    print("Done.")
