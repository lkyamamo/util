"""
Plot scalar quantities vs time from trajectory.h5:
  - jet tip position    (jet_tip_x)
  - hydronium ion count (sum of hydronium_count over x bins)
  - penetration depth   (crater/depth, from sphere fit)
  - pit volume          (spherical cap volume from crater/sphere_radius and crater/depth)
  - hydronium ion count vs x position, overlaid at several timesteps

Time axis is real time in picoseconds, computed from the stored per-frame
timestep array and the simulation's MD step size (--timestep-size, ps/step).
The dump frequency (--dump-frequency, steps/frame) is used only as a sanity
check that frame spacing matches what's expected.

Usage:
    python plot_trajectory.py <trajectory_h5> --timestep-size 0.001 --dump-frequency 50 [--out-dir DIR]
"""

import argparse
import os
import sys

import numpy as np
import h5py
import matplotlib.pyplot as plt

# --- Publication style (APS-like: serif font, inward ticks on all sides) ---
FONT_SIZE  = 20
LINE_COLOR = 'tab:blue'   # plotted data color; axes/text stay black

plt.rcParams.update({
    'font.family':       'serif',
    'font.size':         FONT_SIZE,
    'axes.labelsize':    FONT_SIZE,
    'axes.titlesize':    FONT_SIZE,
    'xtick.labelsize':   FONT_SIZE,
    'ytick.labelsize':   FONT_SIZE,
    'legend.fontsize':   FONT_SIZE * 0.6,
    'axes.edgecolor':    'black',
    'axes.labelcolor':   'black',
    'text.color':        'black',
    'xtick.color':       'black',
    'ytick.color':       'black',
    'axes.linewidth':    1.5,
    'lines.linewidth':   2.5,
    'lines.markersize':  6,
    'xtick.direction':   'in',
    'ytick.direction':   'in',
    'xtick.top':         True,
    'ytick.right':       True,
    'xtick.major.size':  6,
    'xtick.major.width': 1.5,
    'xtick.minor.size':  3,
    'xtick.minor.width': 1.0,
    'ytick.major.size':  6,
    'ytick.major.width': 1.5,
    'ytick.minor.size':  3,
    'ytick.minor.width': 1.0,
    'axes.grid':         False,
    'savefig.dpi':       300,
    'figure.figsize':    (8, 6),
})


def spherical_cap_volume(depth, radius):
    """V = (pi*h^2/3)*(3R - h), elementwise. NaN-safe."""
    h = depth
    R = radius
    V = (np.pi * h ** 2 / 3.0) * (3.0 * R - h)
    return np.where(V < 0, np.nan, V)


def compute_time_ps(timesteps, timestep_size_ps, dump_frequency):
    """Real time in ps from the stored per-frame timestep array.

    timestep_size_ps: MD integration step size (ps/step), from the LAMMPS input.
    dump_frequency: expected steps between dumps, used only to sanity-check
    that frame spacing matches expectations — not used in the calculation
    itself, since the actual per-frame timestep is already known.
    """
    diffs = np.diff(timesteps)
    if len(diffs) > 0 and not np.allclose(diffs, dump_frequency):
        unique_diffs = np.unique(diffs)
        print(f"WARNING: frame spacing does not match --dump-frequency={dump_frequency}. "
              f"Observed step differences: {unique_diffs}", file=sys.stderr)
    return timesteps.astype(np.float64) * timestep_size_ps


def style_axes(ax):
    ax.minorticks_on()
    ax.tick_params(direction='in', which='both', top=True, right=True)


def plot_quantity(time_ps, y, ylabel, out_path):
    fig, ax = plt.subplots()

    n_markers = 30
    markevery = max(1, len(time_ps) // n_markers)
    ax.plot(time_ps, y, marker='o', markevery=markevery, color=LINE_COLOR)

    ax.set_xlabel('Time (ps)')
    ax.set_ylabel(ylabel)
    style_axes(ax)

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  wrote {out_path}")


def plot_hydronium_vs_position(f, time_ps, out_path, n_curves=6):
    hydronium_count = f['hydronium_count'][()]   # (T, n_hydronium_x)
    T, n_hx = hydronium_count.shape

    hydronium_voxel_x = float(f.attrs['hydronium_voxel_x'])
    xlo = float(f.attrs['xlo'])
    bin_centers = (np.arange(n_hx) + 0.5) * hydronium_voxel_x + xlo

    frame_indices = np.unique(np.linspace(0, T - 1, min(n_curves, T)).astype(int))

    fig, ax = plt.subplots()
    cmap = plt.colormaps['viridis']
    denom = max(1, len(frame_indices) - 1)
    for i, frame_idx in enumerate(frame_indices):
        color = cmap(i / denom)
        ax.plot(bin_centers, hydronium_count[frame_idx], color=color,
                label=f"t = {time_ps[frame_idx]:.2f} ps")

    ax.set_xlabel('x position (Å)')
    ax.set_ylabel('Hydronium ion count')
    style_axes(ax)

    legend = ax.legend(frameon=False)
    for text in legend.get_texts():
        text.set_color('black')

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  wrote {out_path}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument('trajectory_h5')
    p.add_argument('--timestep-size', type=float, required=True,
                    help="MD integration step size in ps/step (from the LAMMPS input)")
    p.add_argument('--dump-frequency', type=int, required=True,
                    help="Expected number of MD steps between dumps (used as a sanity check)")
    p.add_argument('--n-hydronium-curves', type=int, default=6,
                    help="Number of timesteps to overlay in the hydronium-vs-position plot")
    p.add_argument('--out-dir', default=None,
                    help="Directory for output plots (default: alongside trajectory_h5)")
    args = p.parse_args()

    out_dir = args.out_dir or os.path.dirname(os.path.abspath(args.trajectory_h5))
    os.makedirs(out_dir, exist_ok=True)

    with h5py.File(args.trajectory_h5, 'r') as f:
        t = f['timestep'][()]   # (T,)
        time_ps = compute_time_ps(t, args.timestep_size, args.dump_frequency)

        # --- jet tip position ---
        jet_tip_x = f['jet_tip_x'][()]   # (T,)
        plot_quantity(time_ps, jet_tip_x, 'Jet tip x position (Å)',
                      os.path.join(out_dir, 'jet_tip_position.png'))

        # --- hydronium ion count ---
        hydronium_total = f['hydronium_count'][()].sum(axis=1)   # (T,)
        plot_quantity(time_ps, hydronium_total, 'Hydronium ion count',
                      os.path.join(out_dir, 'hydronium_count.png'))

        # --- hydronium ion count vs position, overlaid at several timesteps ---
        plot_hydronium_vs_position(f, time_ps, os.path.join(out_dir, 'hydronium_vs_position.png'),
                                    n_curves=args.n_hydronium_curves)

        # --- penetration depth ---
        if 'crater' in f and 'depth' in f['crater']:
            depth = f['crater']['depth'][()]   # (T,)
            plot_quantity(time_ps, depth, 'Penetration depth (Å)',
                          os.path.join(out_dir, 'penetration_depth.png'))

            # --- pit volume (spherical cap) ---
            radius = f['crater']['sphere_radius'][()]   # (T,)
            volume = spherical_cap_volume(depth, radius)
            plot_quantity(time_ps, volume, 'Pit volume (Å³)',
                          os.path.join(out_dir, 'pit_volume.png'))
        else:
            print("  no crater group found — skipping penetration depth / pit volume plots")

    print("Done.")


if __name__ == '__main__':
    main()
