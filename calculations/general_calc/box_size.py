import numpy as np
import sys

def get_unit_cell_mass_grams(data, repeat):
    print("input masses and quantities")
    unit_cell_mass_amu = 0
    for pair in data:
        print("{0} of mass {1} amu".format(pair[0],pair[1]))
        unit_cell_mass_amu += pair[0] * pair[1]

    print("repeat dimensions: {}".format(repeat))

    unit_cell_mass_g = unit_cell_mass_amu * 1.66054e-24

    print("unit cell mass: {0} g".format(unit_cell_mass_g))
    print("unit cell mass: {0} amu".format(unit_cell_mass_amu))

    return unit_cell_mass_g

def get_volume(target_density, data, repeat):
    density = target_density

    print("target density: {0} g/cm^3".format(density))

    volume_cm3 = (get_unit_cell_mass_grams(data,repeat)/density)
    volume_A3 = volume_cm3*1.0e24

    print("target unit cell volume: {0} A^3".format(volume_A3))
    print("target unit cell cube side length: {0} A".format((volume_A3**(1/3))))
    print("target supercell cube side length: {0} A".format((volume_A3**(1/3))*repeat[0]))

    multiplier = repeat[0]*repeat[1]*repeat[2]

    print("target total volume: {} A^3".format(volume_A3*repeat[0]*repeat[1]*repeat[2]))
    print("target total volume: {} cm^3".format(volume_cm3*repeat[0]*repeat[1]*repeat[2]))

# side_length in angstroms
def get_density(unit_cell_dims, data, repeat):

    volume = unit_cell_dims[0]*unit_cell_dims[1]*unit_cell_dims[2] 
    print("current unit cell volume: {0} A^3".format(volume))
    supercell_volume = volume*repeat[0]*repeat[1]*repeat[2]
    print("current supercell dimensions (A): {0}x{1}x{2}".format(unit_cell_dims[0]*repeat[0],unit_cell_dims[1]*repeat[1],unit_cell_dims[2]*repeat[2])) 
    print("current supercell volume: {0} A^3".format(supercell_volume))

    volume_cm = volume * 1.0e-24

    mass = get_unit_cell_mass_grams(data, repeat)
    print("current density: {0} g/cm^3".format(mass/volume_cm))    


filename = 'masses_quantities.dat'
data = np.loadtxt(filename)

target_density = 2.2
current_unit_cell_dim = (30.1181, 30.1194201, 40.564142)
repeat = (1,1,1)

get_density(current_unit_cell_dim, data, repeat)
get_volume(target_density, data, repeat)
