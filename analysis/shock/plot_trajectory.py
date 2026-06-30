"""
Plot scalar quantities vs time from trajectory.h5:
  - jet tip position    (jet_tip_x)
  - hydronium ion count (sum of hydronium_count over x bins)
  - penetration depth   (crater/depth, from sphere fit)
  - pit volume          (spherical cap volume from crater/sphere_radius and crater/depth)

Usage:
    python plot_trajectory.py <trajectory_h5> [--out-dir DIR]
"""

import argparse
import os

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
    'legend.fontsize':   FONT_SIZE,
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


def plot_quantity(t, y, ylabel, out_path):
    fig, ax = plt.subplots()

    n_markers = 30
    markevery = max(1, len(t) // n_markers)
    ax.plot(t, y, marker='o', markevery=markevery, color=LINE_COLOR)

    ax.set_xlabel('Timestep')
    ax.set_ylabel(ylabel)

    ax.minorticks_on()
    ax.tick_params(direction='in', which='both', top=True, right=True)

    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(f"  wrote {out_path}")


def main():
    p = argparse.ArgumentParser()
    p.add_argument('trajectory_h5')
    p.add_argument('--out-dir', default=None,
                    help="Directory for output plots (default: alongside trajectory_h5)")
    args = p.parse_args()

    out_dir = args.out_dir or os.path.dirname(os.path.abspath(args.trajectory_h5))
    os.makedirs(out_dir, exist_ok=True)

    with h5py.File(args.trajectory_h5, 'r') as f:
        t = f['timestep'][()]   # (T,)

        # --- jet tip position ---
        jet_tip_x = f['jet_tip_x'][()]   # (T,)
        plot_quantity(t, jet_tip_x, 'Jet tip x position (Å)',
                      os.path.join(out_dir, 'jet_tip_position.png'))

        # --- hydronium ion count ---
        hydronium_total = f['hydronium_count'][()].sum(axis=1)   # (T,)
        plot_quantity(t, hydronium_total, 'Hydronium ion count',
                      os.path.join(out_dir, 'hydronium_count.png'))

        # --- penetration depth ---
        if 'crater' in f and 'depth' in f['crater']:
            depth = f['crater']['depth'][()]   # (T,)
            plot_quantity(t, depth, 'Penetration depth (Å)',
                          os.path.join(out_dir, 'penetration_depth.png'))

            # --- pit volume (spherical cap) ---
            radius = f['crater']['sphere_radius'][()]   # (T,)
            volume = spherical_cap_volume(depth, radius)
            plot_quantity(t, volume, 'Pit volume (Å³)',
                          os.path.join(out_dir, 'pit_volume.png'))
        else:
            print("  no crater group found — skipping penetration depth / pit volume plots")

    print("Done.")


if __name__ == '__main__':
    main()
