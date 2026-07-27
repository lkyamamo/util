import argparse
import os

import log_lammps_reader
import polars
import matplotlib
matplotlib.use("Agg")  # no display on compute nodes
import matplotlib.pyplot as plt

parser = argparse.ArgumentParser(description="Plot LAMMPS thermodynamic output.")
parser.add_argument("--log", default="log.lammps", help="LAMMPS log file to read")
parser.add_argument("--outdir", default=".", help="directory for the csv and figures")
args = parser.parse_args()

os.makedirs(args.outdir, exist_ok=True)
filename = os.path.join(args.outdir, "thermodynamics.csv")

# 1 to avoid cg
# 0 if no cg

i = 0
output = log_lammps_reader.parse(args.log, i)

i += 1

while True:
    try:
        df = log_lammps_reader.parse(args.log, i)
    except:
        print("there are {} instances".format(i))
        break
    output = polars.concat([output, df])
    i += 1

output.write_csv(filename)


# thermo_style is "custom step temp pe ke etotal press vol", so the wanted
# quantities sit at positions 0, 1 and 4. Prefer the header names in case the
# style changes, but fall back to position if the reader names them otherwise.
def pick(df, name, index):
    for column in df.columns:
        if column.lower() == name.lower():
            return df[column]
    return df[:, index]


step = pick(output, "Step", 0)

written = [filename]

for name, index, label, plotname in [
    ("Temp", 1, "temperature (K)", "temperature.png"),
    ("TotEng", 4, "total energy (eV)", "total_energy.png"),
]:
    path = os.path.join(args.outdir, plotname)
    fig, ax = plt.subplots(figsize=(5, 4), dpi=300)
    ax.plot(step, pick(output, name, index))
    ax.set_xlabel("step")
    ax.set_ylabel(label)
    fig.savefig(path, dpi=300, bbox_inches="tight", transparent=True)
    plt.close(fig)
    written.append(path)

for path in written:
    print("wrote {}".format(path))
