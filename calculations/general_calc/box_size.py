import numpy as np
import sys

def get_unit_cell_mass(data):
    """Compute unit cell mass in grams. data is (quantity, mass_amu) per species;
    quantity is the total count in the unit cell."""
    print("input masses and quantities (total in unit cell)")
    total_mass_amu = 0
    for pair in data:
        print("{0} of mass {1} amu".format(pair[0], pair[1]))
        total_mass_amu += pair[0] * pair[1]

    mass_g = total_mass_amu * 1.66054e-24
    print(f"total mass: {total_mass_amu} amu")
    print(f"total mass: {mass_g} g")

    return mass_g

def get_supercell_mass(data, repeat):
    u_mass_g = get_unit_cell_mass(data)
    total_mass_g = u_mass_g*repeat[0]*repeat[1]*repeat[2]
    total_mass_amu = total_mass_g/1.66054e-24
    print(f"supercell mass")
    print(f"total mass: {total_mass_amu} amu")
    print(f"total mass: {total_mass_g} g")
    total_atoms = 0
    for pair in data:
        current_atoms = pair[0]*repeat[0]*repeat[1]*repeat[2]
        total_atoms += current_atoms
        print(f"Number of atoms of mass {pair[1]} amu: {current_atoms}")
    print(f"total atoms: {total_atoms}")
    return total_mass_g

def get_volume(dimensions):
    return dimensions[0]*dimensions[1]*dimensions[2]

    

# side_length in angstroms
def get_density(data, dimensions):

    volume_A3 = get_volume(dimensions)
    mass = get_unit_cell_mass(data)
    density_g_cm3 = mass/volume_A3 * 10**24
    return density_g_cm3

# given 2 supercell dimensions return the last that would yield the target density
def get_third_dimension(data, known_dimensions, repeat, target_density):
    total_mass_g = get_supercell_mass(data, repeat) # total mass in the supercell (g)
    third_dim = total_mass_g/(target_density * known_dimensions[0] * known_dimensions[1]) * 10**24 # convert to A
    
    print(f"final supercell dimension: {third_dim} A")

    return third_dim


    


filename = 'masses_quantities.dat'
data = np.loadtxt(filename)

target_density = 1.0
current_unit_cell_dim = (30.1181,30.1194,38.2032)
repeat = (1,1,1)
known_dimensions = (10.0163, 10.0162)

current_density = get_density(data, current_unit_cell_dim)
print(f"current density: {current_density} g/cm^3")
print(f"target density: {target_density} g/cm^3")
print(f"Number of atoms of mass {data[0][1]} amu: {data[0][0]}")
print(f"Number of atoms of mass {data[1][1]} amu: {data[1][0]}")

get_third_dimension(data, known_dimensions, repeat, target_density)
