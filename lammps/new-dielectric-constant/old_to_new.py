from fileinput import filename
import os

template = "303."
output_filename = "all_lammps_new.xyz"


with open(output_filename, "w") as f_out:
    for num in range(0, 10000001, 10):
        filename = f"{template}{num}"
        if not os.path.exists(filename):
            print(f"Warning: File {filename} not found, skipping...")
            continue

        with open(filename, "r") as f:
            #line 1: ITEM: TIMESTEP
            line = f.readline()

            #line 2: timestep number
            line = f.readline()
            timestep = int(line.strip())     

            #line 3: "ITEM: NUMBER OF ATOMS"
            line = f.readline()        

            #line 4: ITEM: NUMBER OF ATOMS
            line = f.readline().strip()
            num_atoms = int(line[4])

            f_out.write(f"{num_atoms}\n")
            f_out.write(f"Atoms. Timestep: {timestep}\n")

            #line 5: ITEM: BOX BOUNDS pp pp pp
            line = f.readline()

            #line 6: a bound
            line = f.readline()
            #line 7: b bound
            line = f.readline()
            #line 8: c bound
            line = f.readline()
            #line 9: "ITEM: ATOMS id mol type x y z q"
            line = f.readline()

            #line 10: atom data starts
            for i in range(num_atoms):
                line = f.readline()
                atom_id, mol_id, atom_type, x, y, z, q = line.split()
                f_out.write(f"{atom_type} {x} {y} {z}\n")