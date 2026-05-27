from ase.io import read, write
from ase.constraints import FixAtoms

atoms = read("/Users/loganyamamoto/Desktop/Research/grants/geo_sciences/bubble_collapse/data/systems/sioh_cov/a-SiO/1_surface/N192/x-axis/1.nvt/0085/1.int_nvt.POSCAR")

# Freeze atoms below z = 2.5 Å
xmax = 2.0

mask = atoms.positions[:, 0] < xmax
print(f"{mask.sum()} atoms frozen")
atoms.set_constraint(FixAtoms(mask=mask))

write("POSCAR_sd", atoms, format="vasp", vasp5=True, direct=True)
