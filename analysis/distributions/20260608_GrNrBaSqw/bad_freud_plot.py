"""
bad_freud_plot.py — Individual BAD plot exporter

Reads the CSV produced by bad_freud.py and saves one PNG per triplet
into OUTPUT_DIR.

QUICK START
-----------
1. Run bad_freud.py first to generate the CSV.
2. Adjust the config below as desired.
3. Run:  python bad_freud_plot.py

OUTPUT
------
One PNG per triplet column in the CSV, written to OUTPUT_DIR/<label>.png.
"""

import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# =============================================================================
# CONFIGURATION — edit these variables between runs
# =============================================================================

# Input CSV produced by bad_freud.py
INPUT_CSV  = "bads.csv"

# Subdirectory where individual PNGs are written
OUTPUT_DIR = "bads"

# Plot appearance
PLOT_FIGSIZE  = (4, 3)       # (width, height) in inches
PLOT_DPI      = 150
PLOT_COLOR    = 'steelblue'  # any matplotlib color string
PLOT_LWIDTH   = 3
PLOT_LSTYLE   = '-'          # '-', '--', ':', '-.'
SPINE_LWIDTH  = 2.0          # border (axes frame) line width
FONT_FAMILY   = 'DejaVu Sans' # any font available to matplotlib
FONT_WEIGHT   = 'bold'        # 'normal' or 'bold'
FONT_LABEL    = 20
FONT_TICK_X   = 12           # x-axis tick label font size (pt)
FONT_TICK_Y   = 9            # y-axis tick label font size (pt); only used when YTICKS = True
TICK_LENGTH   = 6            # tick mark length in points
TICK_WIDTH    = 2.5          # tick mark width in points
XTICK_POSITIONS = [i for i in range(0, 181, 30)]       # e.g. [0, 60, 120, 180]; None = matplotlib auto ticks
YTICKS        = False        # set True to show y-axis tick marks and labels

# =============================================================================
# END CONFIGURATION
# =============================================================================


def load_csv(filename):
    """
    Read the BAD CSV.

    Returns
    -------
    angles : ndarray, shape (N,)
    labels : list of str
    hists  : ndarray, shape (N, n_triplets)
    """
    with open(filename) as f:
        header = f.readline().strip()

    cols   = header.split(',')
    labels = cols[1:]           # first column is angle_deg

    data   = np.loadtxt(filename, delimiter=',', skiprows=1)
    angles = data[:, 0]
    hists  = data[:, 1:]
    return angles, labels, hists


def plot_individual(angles, labels, hists):
    mpl.rcParams['font.family'] = FONT_FAMILY
    mpl.rcParams['font.weight'] = FONT_WEIGHT
    mpl.rcParams['axes.labelweight'] = FONT_WEIGHT
    mpl.rcParams['axes.titleweight'] = FONT_WEIGHT

    os.makedirs(OUTPUT_DIR, exist_ok=True)

    for i, label in enumerate(labels):
        fig, ax = plt.subplots(figsize=PLOT_FIGSIZE)

        ax.plot(angles, hists[:, i],
                color=PLOT_COLOR,
                linewidth=PLOT_LWIDTH,
                linestyle=PLOT_LSTYLE)

        ax.set_xlabel('Angle (deg)', fontsize=FONT_LABEL, fontweight=FONT_WEIGHT)
        ax.set_ylabel('P(θ)', fontsize=FONT_LABEL, fontweight=FONT_WEIGHT)
        ax.tick_params(axis='x', labelsize=FONT_TICK_X, length=TICK_LENGTH, width=TICK_WIDTH)
        ax.set_xlim(0, 180)
        if XTICK_POSITIONS is not None:
            ax.set_xticks(XTICK_POSITIONS)

        if YTICKS:
            ax.tick_params(axis='y', labelsize=FONT_TICK_Y, length=TICK_LENGTH, width=TICK_WIDTH)
        else:
            ax.yaxis.set_ticks([])

        for spine in ax.spines.values():
            spine.set_linewidth(SPINE_LWIDTH)

        plt.tight_layout()
        out_path = os.path.join(OUTPUT_DIR, f"{label}.png")
        plt.savefig(out_path, dpi=PLOT_DPI)
        plt.close(fig)
        print(f"Saved {out_path}")


if __name__ == '__main__':
    print(f"Reading CSV: {INPUT_CSV}")
    angles, labels, hists = load_csv(INPUT_CSV)
    print(f"Triplets found: {labels}")
    plot_individual(angles, labels, hists)
    print(f"Done. {len(labels)} plot(s) written to '{OUTPUT_DIR}/'")
