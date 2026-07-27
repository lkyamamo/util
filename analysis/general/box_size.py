import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
DEFAULT_MASSES = HERE / "masses_quantities.dat"
DEFAULT_BOX_PARAMS = HERE / "box_parameters.dat"

AMU_TO_G = 1.66054e-24
A3_TO_CM3 = 1e-24


def load_masses_quantities(path):
    """Load (quantity, mass_amu) per species; quantity is count per unit cell."""
    raw = np.atleast_2d(np.loadtxt(path))
    if raw.shape[1] < 2:
        raise ValueError(f"{path}: need at least two columns (quantity, mass_amu)")
    return raw[:, :2].astype(float)


def load_box_parameters(path):
    """
    Plain-text key/value lines (optional # comments).
    Required keys: unit_cell_lengths_A (3), repeat (3), target_density_g_cm3 (1).
    Optional key: known_supercell_lengths_A (2) - if given, the third supercell
    edge is solved for; otherwise a cubic supercell edge is solved for.
    """
    field_sizes = {
        "unit_cell_lengths_A": 3,
        "repeat": 3,
        "target_density_g_cm3": 1,
        "known_supercell_lengths_A": 2,
    }
    required = {"unit_cell_lengths_A", "repeat", "target_density_g_cm3"}

    found = {}
    with open(path) as f:
        for line in f:
            s = line.split("#", 1)[0].strip()
            if not s:
                continue
            parts = s.split()
            key = parts[0]
            if key not in field_sizes:
                raise ValueError(f"{path}: unknown key {key!r}")
            n = field_sizes[key]
            if len(parts) != n + 1:
                raise ValueError(
                    f"{path}: {key} expects {n} value(s), got {len(parts) - 1}"
                )
            vals = tuple(float(x) for x in parts[1 : n + 1])
            found[key] = vals[0] if n == 1 else vals

    missing = required - set(found)
    if missing:
        raise ValueError(f"{path}: missing keys: {sorted(missing)}")
    found.setdefault("known_supercell_lengths_A", None)
    return found


def unit_cell_mass_g(data):
    """Total unit cell mass in grams. data is (quantity, mass_amu) per species."""
    total_mass_amu = sum(qty * mass for qty, mass in data)
    return total_mass_amu * AMU_TO_G


def supercell_mass_g(data, repeat):
    return unit_cell_mass_g(data) * repeat[0] * repeat[1] * repeat[2]


def volume_A3(dimensions):
    return dimensions[0] * dimensions[1] * dimensions[2]


def density_g_cm3(mass_g, volume_A3_value):
    return mass_g / volume_A3_value / A3_TO_CM3


def third_dimension_A(mass_g, target_density, dim_a, dim_b):
    """Given 2 supercell dimensions, the edge length that yields target_density."""
    return mass_g / (target_density * dim_a * dim_b) / A3_TO_CM3


def cubic_dimension_A(mass_g, target_density):
    """Cubic supercell edge length that yields target_density."""
    return (mass_g / target_density / A3_TO_CM3) ** (1 / 3)


def main(masses_path=None, box_params_path=None):
    masses_path = Path(masses_path or DEFAULT_MASSES)
    box_params_path = Path(box_params_path or DEFAULT_BOX_PARAMS)

    data = load_masses_quantities(masses_path)
    cfg = load_box_parameters(box_params_path)

    unit_cell_dim = cfg["unit_cell_lengths_A"]
    repeat = cfg["repeat"]
    target_density = cfg["target_density_g_cm3"]
    known_dims = cfg["known_supercell_lengths_A"]

    unit_mass_g = unit_cell_mass_g(data)
    total_mass_g = unit_mass_g * repeat[0] * repeat[1] * repeat[2]
    total_atoms = sum(int(qty * repeat[0] * repeat[1] * repeat[2]) for qty, _ in data)
    current_density = density_g_cm3(unit_mass_g, volume_A3(unit_cell_dim))

    print("unit cell composition:")
    for qty, mass in data:
        print(f"  {int(qty)} atoms of mass {mass} amu")
    print(f"unit cell mass: {unit_mass_g:.6e} g")
    print(f"supercell repeat: {repeat[0]:g} x {repeat[1]:g} x {repeat[2]:g}")
    print(f"supercell mass: {total_mass_g:.6e} g")
    print(f"supercell atom count: {total_atoms}")
    print(f"current density: {current_density:.6f} g/cm^3")
    print(f"target density: {target_density:.6f} g/cm^3")

    if known_dims is not None:
        dim_a, dim_b = known_dims
        third = third_dimension_A(total_mass_g, target_density, dim_a, dim_b)
        print(
            f"supercell dimensions for target density: "
            f"{dim_a:.4f} x {dim_b:.4f} x {third:.4f} A"
        )
        return third
    else:
        cubic = cubic_dimension_A(total_mass_g, target_density)
        print(f"cubic supercell dimension for target density: {cubic:.4f} A")
        return cubic


if __name__ == "__main__":
    m = sys.argv[1] if len(sys.argv) > 1 else None
    b = sys.argv[2] if len(sys.argv) > 2 else None
    main(masses_path=m, box_params_path=b)
