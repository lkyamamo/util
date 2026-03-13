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


# given two supercell dimensions, return the third supercell dimension
# that achieves the target density
# known_dims: (dim1, dim2) in angstroms for the two known supercell axes
# axis: which supercell axis is unknown (0=x, 1=y, 2=z)
def get_third_dimension(target_density, known_dims, data, repeat, axis=2):
    unit_cell_mass_g = get_unit_cell_mass_grams(data, repeat)
    total_mass_g = unit_cell_mass_g * repeat[0] * repeat[1] * repeat[2]

    target_volume_A3 = (total_mass_g / target_density) * 1.0e24
    third_dim = target_volume_A3 / (known_dims[0] * known_dims[1])

    axis_labels = ['x', 'y', 'z']
    known_axes = [l for i, l in enumerate(axis_labels) if i != axis]

    supercell = [None, None, None]
    for i, val in zip([j for j in range(3) if j != axis], known_dims):
        supercell[i] = val
    supercell[axis] = third_dim

    print("target density: {0} g/cm^3".format(target_density))
    print("known supercell dims: {0}={1} A, {2}={3} A".format(known_axes[0], known_dims[0], known_axes[1], known_dims[1]))
    print("required supercell {0} dimension: {1} A".format(axis_labels[axis], third_dim))
    print("final supercell dimensions (A): {0}x{1}x{2}".format(*supercell))
    print("final supercell volume: {0} A^3".format(target_volume_A3))

    return third_dim


filename = 'masses_quantities.dat'
data = np.loadtxt(filename)

target_density = 2.2
current_unit_cell_dim = (30.1181, 30.1194201, 40.564142)
repeat = (1,1,1)

get_density(current_unit_cell_dim, data, repeat)
get_volume(target_density, data, repeat)
