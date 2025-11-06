import log_lammps_reader
import polars

filename = "thermodynamics.csv"

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
