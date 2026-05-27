import numpy as np
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_MASSES = HERE / "masses_quantities.dat"
DEFAULT_BOX_PARAMS = HERE / "box_parameters.dat"


def load_masses_quantities(path):
    """Load (quantity, mass_amu) per species; quantity is count per unit cell."""
    raw = np.loadtxt(path)
    raw = np.atleast_2d(raw)
    if raw.shape[1] < 2:
        raise ValueError(f"{path}: need at least two columns (quantity, mass_amu)")
    return raw[:, :2].astype(float)


def load_box_parameters(path):
    """
    Plain-text key/value lines (optional # comments).
    Keys: unit_cell_lengths_A (3), repeat (3), target_density_g_cm3 (1),
          known_supercell_lengths_A (2).
    """
    expected = {
        "unit_cell_lengths_A": 3,
        "repeat": 3,
        "target_density_g_cm3": 1,
        "known_supercell_lengths_A": 2,
    }
    found = {}
    with open(path) as f:
        for line in f:
            s = line.split("#", 1)[0].strip()
            if not s:
                continue
            parts = s.split()
            key = parts[0]
            if key not in expected:
                raise ValueError(f"{path}: unknown key {key!r}")
            n = expected[key]
            if len(parts) != n + 1:
                raise ValueError(
                    f"{path}: {key} expects {n} value(s), got {len(parts) - 1}"
                )
            vals = tuple(float(x) for x in parts[1 : n + 1])
            found[key] = vals[0] if n == 1 else vals
    missing = set(expected) - set(found)
    if missing:
        raise ValueError(f"{path}: missing keys: {sorted(missing)}")
    return {
        "unit_cell_lengths_A": found["unit_cell_lengths_A"],
        "repeat": found["repeat"],
        "target_density_g_cm3": found["target_density_g_cm3"],
        "known_supercell_lengths_A": found["known_supercell_lengths_A"],
    }


def get_unit_cell_mass(data):
    """Compute unit cell mass in grams. data is (quantity, mass_amu) per species;
    quantity is the total count in the unit cell."""
    print("input masses and quantities (total in unit cell)")
    total_mass_amu = 0.0
    for pair in data:
        print("{0} of mass {1} amu".format(pair[0], pair[1]))
        total_mass_amu += pair[0] * pair[1]

    mass_g = total_mass_amu * 1.66054e-24
    print(f"total mass: {total_mass_amu} amu")
    print(f"total mass: {mass_g} g")

    return mass_g


def get_supercell_mass(data, repeat):
    u_mass_g = get_unit_cell_mass(data)
    total_mass_g = u_mass_g * repeat[0] * repeat[1] * repeat[2]
    total_mass_amu = total_mass_g / 1.66054e-24
    print(f"supercell mass")
    print(f"total mass: {total_mass_amu} amu")
    print(f"total mass: {total_mass_g} g")
    total_atoms = 0
    for pair in data:
        current_atoms = int(pair[0] * repeat[0] * repeat[1] * repeat[2])
        total_atoms += current_atoms
        print(f"Number of atoms of mass {pair[1]} amu: {current_atoms}")
    print(f"total atoms: {total_atoms}")
    return total_mass_g


def get_volume(dimensions):
    return dimensions[0] * dimensions[1] * dimensions[2]


# side_length in angstroms
def get_density(data, dimensions):
    volume_A3 = get_volume(dimensions)
    mass = get_unit_cell_mass(data)
    density_g_cm3 = mass / volume_A3 * 10**24
    return density_g_cm3


# given 2 supercell dimensions return the last that would yield the target density
def get_third_dimension(data, known_dimensions, repeat, target_density):
    total_mass_g = get_supercell_mass(data, repeat)  # total mass in the supercell (g)
    third_dim = (
        total_mass_g
        / (target_density * known_dimensions[0] * known_dimensions[1])
        * 10**24
    )  # convert to A

    print(f"final supercell dimension: {third_dim} A")

    return third_dim


def main(
    masses_path=None,
    box_params_path=None,
):
    masses_path = Path(masses_path or DEFAULT_MASSES)
    box_params_path = Path(box_params_path or DEFAULT_BOX_PARAMS)

    data = load_masses_quantities(masses_path)
    cfg = load_box_parameters(box_params_path)

    current_unit_cell_dim = cfg["unit_cell_lengths_A"]
    repeat = cfg["repeat"]
    target_density = cfg["target_density_g_cm3"]
    known_dimensions = cfg["known_supercell_lengths_A"]

    current_density = get_density(data, current_unit_cell_dim)
    print(f"current density: {current_density} g/cm^3")
    print(f"target density: {target_density} g/cm^3")
    for row in data:
        print(f"Number of atoms of mass {row[1]} amu in unit cell: {int(row[0])}")

    get_third_dimension(data, known_dimensions, repeat, target_density)


if __name__ == "__main__":
    m = sys.argv[1] if len(sys.argv) > 1 else None
    b = sys.argv[2] if len(sys.argv) > 2 else None
    main(masses_path=m, box_params_path=b)
