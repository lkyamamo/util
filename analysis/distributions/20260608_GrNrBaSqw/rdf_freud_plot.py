"""
rdf_freud_plot.py — Individual RDF / n(r) plot exporter

Reads the CSVs produced by rdf_freud.py and saves one PNG per pair column
into separate subdirectories for g(r) and n(r).

QUICK START
-----------
1. Run rdf_freud.py first to generate the CSVs.
2. Adjust the config below as desired.
3. Run:  python rdf_freud_plot.py

OUTPUT
------
- OUTPUT_RDF_DIR/<label>.png  — one g(r) plot per pair/total/neutron column
- OUTPUT_NR_DIR/<label>.png   — one n(r) plot per pair column
Set either INPUT_*_CSV to None to skip that set entirely.
"""

import os

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# =============================================================================
# CONFIGURATION — edit these variables between runs
# =============================================================================

# Input CSVs produced by rdf_freud.py; set to None to skip
INPUT_RDF_CSV = "rdfs.csv"
INPUT_NR_CSV  = "nrs.csv"

# Subdirectories for individual PNGs
OUTPUT_RDF_DIR = "rdfs"
OUTPUT_NR_DIR  = "nrs"

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
FONT_TICK_X   = 14           # x-axis tick label font size (pt)
FONT_TICK_Y   = 9            # y-axis tick label font size (pt); only used when YTICKS = True
TICK_LENGTH   = 6            # tick mark length in points
TICK_WIDTH    = 2.0          # tick mark width in points
XTICK_POSITIONS = None       # e.g. [0, 2, 4, 6, 8, 10]; None = matplotlib auto ticks
YTICKS        = False        # set True to show y-axis tick marks and labels

# g(r) reference line at y = 1 (ideal gas limit)
SHOW_REF_LINE    = False      # draw a horizontal dashed line at g(r) = 1
REF_LINE_COLOR   = 'gray'
REF_LINE_LSTYLE  = '--'
REF_LINE_LWIDTH  = 1.0

# =============================================================================
# END CONFIGURATION
# =============================================================================


def load_csv(filename):
    """
    Read a two-column-or-more CSV (first column = x axis, rest = data columns).

    Returns
    -------
    x      : ndarray, shape (N,)
    labels : list of str
    data   : ndarray, shape (N, n_cols)
    """
    with open(filename) as f:
        header = f.readline().strip()

    cols   = header.split(',')
    labels = cols[1:]

    arr  = np.loadtxt(filename, delimiter=',', skiprows=1)
    x    = arr[:, 0]
    data = arr[:, 1:]
    return x, labels, data


def _apply_style():
    mpl.rcParams['font.family']      = FONT_FAMILY
    mpl.rcParams['font.weight']      = FONT_WEIGHT
    mpl.rcParams['axes.labelweight'] = FONT_WEIGHT
    mpl.rcParams['axes.titleweight'] = FONT_WEIGHT


def _save_plots(x, labels, data, out_dir, xlabel, ylabel, ref_line=False):
    os.makedirs(out_dir, exist_ok=True)

    for i, label in enumerate(labels):
        fig, ax = plt.subplots(figsize=PLOT_FIGSIZE)

        ax.plot(x, data[:, i],
                color=PLOT_COLOR,
                linewidth=PLOT_LWIDTH,
                linestyle=PLOT_LSTYLE)

        if ref_line and SHOW_REF_LINE:
            ax.axhline(1.0,
                       color=REF_LINE_COLOR,
                       linestyle=REF_LINE_LSTYLE,
                       linewidth=REF_LINE_LWIDTH)

        ax.set_xlabel(xlabel, fontsize=FONT_LABEL, fontweight=FONT_WEIGHT)
        ax.set_ylabel(ylabel, fontsize=FONT_LABEL, fontweight=FONT_WEIGHT)
        ax.tick_params(axis='x', labelsize=FONT_TICK_X, length=TICK_LENGTH, width=TICK_WIDTH)

        if XTICK_POSITIONS is not None:
            ax.set_xticks(XTICK_POSITIONS)

        if YTICKS:
            ax.tick_params(axis='y', labelsize=FONT_TICK_Y, length=TICK_LENGTH, width=TICK_WIDTH)
        else:
            ax.yaxis.set_ticks([])

        for spine in ax.spines.values():
            spine.set_linewidth(SPINE_LWIDTH)

        plt.tight_layout()
        out_path = os.path.join(out_dir, f"{label}.png")
        plt.savefig(out_path, dpi=PLOT_DPI)
        plt.close(fig)
        print(f"Saved {out_path}")


if __name__ == '__main__':
    _apply_style()

    if INPUT_RDF_CSV is not None:
        print(f"Reading g(r) CSV: {INPUT_RDF_CSV}")
        r, rdf_labels, rdf_data = load_csv(INPUT_RDF_CSV)
        print(f"Pairs found: {rdf_labels}")
        _save_plots(r, rdf_labels, rdf_data,
                    out_dir=OUTPUT_RDF_DIR,
                    xlabel='r (Å)',
                    ylabel='g(r)',
                    ref_line=True)
        print(f"Done. {len(rdf_labels)} g(r) plot(s) written to '{OUTPUT_RDF_DIR}/'")

    if INPUT_NR_CSV is not None:
        print(f"Reading n(r) CSV: {INPUT_NR_CSV}")
        r, nr_labels, nr_data = load_csv(INPUT_NR_CSV)
        print(f"Pairs found: {nr_labels}")
        _save_plots(r, nr_labels, nr_data,
                    out_dir=OUTPUT_NR_DIR,
                    xlabel='r (Å)',
                    ylabel='n(r)',
                    ref_line=False)
        print(f"Done. {len(nr_labels)} n(r) plot(s) written to '{OUTPUT_NR_DIR}/'")
