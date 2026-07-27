import log_lammps_reader
import polars
import matplotlib
matplotlib.use("Agg")  # no display on compute nodes
import matplotlib.pyplot as plt

filename = "thermodynamics.csv"
plotname = "thermodynamics_plots.png"

# 1 to avoid cg
# 0 if no cg

i = 0
output = log_lammps_reader.parse('log.lammps', i)

i += 1

while True:
    try:
        df = log_lammps_reader.parse('log.lammps', i)
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

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 4), dpi=300)

axs[0].plot(step, pick(output, "Temp", 1))
axs[0].set_xlabel("step")
axs[0].set_ylabel("temperature (K)")

axs[1].plot(step, pick(output, "TotEng", 4))
axs[1].set_xlabel("step")
axs[1].set_ylabel("total energy (eV)")

fig.savefig(plotname, dpi=300, bbox_inches="tight", transparent=True)
plt.close(fig)

print("wrote {} and {}".format(filename, plotname))
