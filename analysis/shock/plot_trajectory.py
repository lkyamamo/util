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


def spherical_cap_volume(depth, radius):
    """V = (pi*h^2/3)*(3R - h), elementwise. NaN-safe."""
    h = depth
    R = radius
    V = (np.pi * h ** 2 / 3.0) * (3.0 * R - h)
    return np.where(V < 0, np.nan, V)


def plot_quantity(t, y, ylabel, title, out_path):
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(t, y, marker='o', markersize=3, linewidth=1)
    ax.set_xlabel('Timestep')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_path, dpi=150)
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
                      'Jet Tip Position vs Time',
                      os.path.join(out_dir, 'jet_tip_position.png'))

        # --- hydronium ion count ---
        hydronium_total = f['hydronium_count'][()].sum(axis=1)   # (T,)
        plot_quantity(t, hydronium_total, 'Hydronium ion count',
                      'Hydronium Ion Count vs Time',
                      os.path.join(out_dir, 'hydronium_count.png'))

        # --- penetration depth ---
        if 'crater' in f and 'depth' in f['crater']:
            depth = f['crater']['depth'][()]   # (T,)
            plot_quantity(t, depth, 'Penetration depth (Å)',
                          'Penetration Depth vs Time',
                          os.path.join(out_dir, 'penetration_depth.png'))

            # --- pit volume (spherical cap) ---
            radius = f['crater']['sphere_radius'][()]   # (T,)
            volume = spherical_cap_volume(depth, radius)
            plot_quantity(t, volume, 'Pit volume (Å³)',
                          'Pit Volume vs Time',
                          os.path.join(out_dir, 'pit_volume.png'))
        else:
            print("  no crater group found — skipping penetration depth / pit volume plots")

    print("Done.")


if __name__ == '__main__':
    main()
