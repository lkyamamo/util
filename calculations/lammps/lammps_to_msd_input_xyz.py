from collections import Counter

a = 28.5400009156
b = 28.5400009156
c = 28.5400009156
alpha = 90
beta = 90
gamma = 90

start_timestep=10000

type_dict = {"1":"Si", "2":"O"}

filename = "all_lammps.xyz"
file_read = open(filename, "r")
lines = file_read.readlines()

box_dim_string = "{:.6f} {:.6f} {:.6f} {:.6f} {:.6f} {:.6f}\n".format(a, b, c, alpha, beta, gamma)
num_atoms=int(lines[3])
start_line = start_timestep*(num_atoms+9) + 9

with open("all.xyz", "w") as file:
    for j in range(start_line, len(lines), num_atoms+9):
        file.write("{}\n".format(num_atoms))
        file.write(box_dim_string)
        for line in lines[j:j+num_atoms]:
            file.write(type_dict[line[0]] + line[1:])
