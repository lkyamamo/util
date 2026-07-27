import argparse
import os

import log_lammps_reader
import polars
import matplotlib
matplotlib.use("Agg")  # no display on compute nodes
import matplotlib.pyplot as plt

# Validated against the categorical palette's slot 1: on white it clears 3:1
# contrast, and a single series needs no legend since the axis label names it.
SERIES_COLOR = "#2a78d6"

# Figures are drawn at final print size so the type lands at the requested
# point size in the journal column -- never drawn large and scaled down.
SINGLE_COLUMN_IN = 3.375

# thermo_style is "custom step temp pe ke etotal press vol", so the plotted
# quantities sit at positions 0, 1 and 4. Units are metal: eV and ps.
QUANTITIES = [
    ("Temp", 1, "temperature (K)", "temperature"),
    ("TotEng", 4, "total energy (eV)", "total_energy"),
]

parser = argparse.ArgumentParser(description="Plot LAMMPS thermodynamic output.")
parser.add_argument("--log", default="log.lammps", help="LAMMPS log file to read")
parser.add_argument("--outdir", default=".", help="directory for the csv and figures")
parser.add_argument("--width", type=float, default=SINGLE_COLUMN_IN,
                    help="figure width in inches (journal column width)")
parser.add_argument("--fontsize", type=float, default=8.0,
                    help="base font size in points, at final figure size")
parser.add_argument("--timestep", type=float, default=None,
                    help="ps per step (metal units); plots time instead of step number")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)


def pick(df, name, index):
    """Column by thermo header name, falling back to position."""
    for column in df.columns:
        if column.lower() == name.lower():
            return df[column]
    return df[:, index]


# 1 to avoid cg
# 0 if no cg

i = 0
blocks = [log_lammps_reader.parse(args.log, i)] # pyright: ignore[reportAttributeAccessIssue]

i += 1

while True:
    try:
        blocks.append(log_lammps_reader.parse(args.log, i)) # pyright: ignore[reportAttributeAccessIssue]
    except:
        print("there are {} instances".format(i))
        break
    i += 1

# A reset_timestep between runs restarts the step counter. Concatenating across
# one would fold the line back over itself, so each reset starts a new segment
# that gets its own csv and its own figures.
segments = [[blocks[0]]]
for previous, current in zip(blocks, blocks[1:]):
    if pick(current, "Step", 0)[0] <= pick(previous, "Step", 0)[-1]:
        segments.append([current])
    else:
        segments[-1].append(current)

if len(segments) > 1:
    print("step counter resets {} time(s); writing {} separate segments".format(
        len(segments) - 1, len(segments)))

plt.rcParams.update({
    "font.size": args.fontsize,
    "axes.labelsize": args.fontsize,
    "xtick.labelsize": args.fontsize - 1,
    "ytick.labelsize": args.fontsize - 1,
    "axes.linewidth": 0.6,
    "lines.linewidth": 1.0,
    # ticks inward on all four sides, minor ticks shown, no grid
    "xtick.direction": "in", "ytick.direction": "in",
    "xtick.top": True, "ytick.right": True,
    "xtick.minor.visible": True, "ytick.minor.visible": True,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
    "axes.grid": False,
    # opaque white: a transparent background composites unpredictably in proofs
    "figure.facecolor": "white", "savefig.facecolor": "white",
    "savefig.transparent": False,
    # Type 42 rather than Type 3 -- journals reject Type 3 fonts
    "pdf.fonttype": 42, "ps.fonttype": 42,
})

written = []

for number, segment in enumerate(segments, start=1):
    # single segment keeps the plain names; only split runs get numbered
    suffix = "" if len(segments) == 1 else "_{:02d}".format(number)
    output = polars.concat(segment)

    csvname = os.path.join(args.outdir, "thermodynamics{}.csv".format(suffix))
    output.write_csv(csvname)
    written.append(csvname)

    step = pick(output, "Step", 0)
    if (step.diff().drop_nulls() < 0).any():
        print("WARNING: step numbers still decrease within segment {} -- the x axis "
              "is not monotonic and the plotted line will fold back".format(number))

    if args.timestep is not None:
        x, xlabel = step * args.timestep, "time (ps)"
    else:
        x, xlabel = step, "step"

    for name, index, label, stem in QUANTITIES:
        y = pick(output, name, index)
        # constrained_layout rather than bbox_inches="tight": tight cropping changes
        # the saved width, which would defeat drawing at exact column width
        fig, ax = plt.subplots(figsize=(args.width, args.width * 0.75),
                               layout="constrained")
        # a dense vector line bloats the pdf; rasterize the line, keep text vector
        ax.plot(x, y, color=SERIES_COLOR, rasterized=len(y) > 20000)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(label)
        ax.margins(x=0.02)
        for extension in ("pdf", "png"):
            path = os.path.join(args.outdir, "{}{}.{}".format(stem, suffix, extension))
            fig.savefig(path, dpi=600)
            written.append(path)
        plt.close(fig)

    # Statistics over the final run block of the segment -- the production
    # portion, after any earlier equilibration blocks have been left behind.
    last_block = segment[-1]
    stats = polars.DataFrame({
        "quantity": last_block.columns,
        "count":    [last_block.height] * last_block.width,
        "mean":     [last_block[c].mean() for c in last_block.columns],
        "std":      [last_block[c].std() for c in last_block.columns],
        "min":      [last_block[c].min() for c in last_block.columns],
        "max":      [last_block[c].max() for c in last_block.columns],
    })
    statsname = os.path.join(args.outdir, "statistics{}.csv".format(suffix))
    stats.write_csv(statsname)
    written.append(statsname)

    print("segment {} of {}: {} block(s), {} rows; statistics over its last block "
          "({} rows):".format(number, len(segments), len(segment), output.height,
                              last_block.height))
    with polars.Config(tbl_rows=-1, tbl_hide_dataframe_shape=True):
        print(stats)

for path in written:
    print("wrote {}".format(path))
