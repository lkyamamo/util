import argparse

from ase.io import read, write
from ase.constraints import FixAtoms

parser = argparse.ArgumentParser()
parser.add_argument("--direction", choices=["x", "y", "z"], default="x",
                     help="axis along which to freeze atoms below the cutoff")
parser.add_argument("--cutoff", type=float, default=2.0,
                     help="freeze atoms with coordinate below this value (Å)")
args = parser.parse_args()

axis = {"x": 0, "y": 1, "z": 2}[args.direction]

atoms = read(
    "2.int_nvt.data",
    format="lammps-data",
    Z_of_type={1: 14, 2: 8, 3: 1},  # Si, O, H
    read_image_flags=False,
)

mask = atoms.positions[:, axis] < args.cutoff
print(f"{mask.sum()} atoms frozen ({args.direction} < {args.cutoff} Å)")
atoms.set_constraint(FixAtoms(mask=mask))

write("POSCAR_sd", atoms, format="vasp", vasp5=True, direct=True, sort=True)
