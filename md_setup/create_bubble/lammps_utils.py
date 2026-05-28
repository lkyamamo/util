"""
lammps_utils.py

Shared helpers for reading LAMMPS data files.
Imported by generate_lammps_script.py and cleanup_orphan_atoms.py.
"""

from __future__ import annotations


def read_header(data_file: str) -> dict[str, int]:
    """
    Stream just the header block of a LAMMPS data file.

    Returns a dict with integer values for any of:
        n_atoms, n_bonds, n_angles, n_dihedrals, n_impropers,
        n_atom_types, n_bond_types, n_angle_types,
        n_dihedral_types, n_improper_types
    """
    info: dict[str, int] = {}
    with open(data_file) as fh:
        next(fh)  # skip comment line
        for line in fh:
            stripped = line.split("#")[0].strip()
            if not stripped:
                continue
            # first alphabetic section keyword signals end of header
            if stripped[0].isalpha() and "types" not in stripped:
                break
            parts = stripped.split()
            if len(parts) >= 2:
                if parts[1] == "atoms":
                    info["n_atoms"] = int(parts[0])
                elif parts[1] == "bonds":
                    info["n_bonds"] = int(parts[0])
                elif parts[1] == "angles":
                    info["n_angles"] = int(parts[0])
                elif parts[1] == "dihedrals":
                    info["n_dihedrals"] = int(parts[0])
                elif parts[1] == "impropers":
                    info["n_impropers"] = int(parts[0])
            if "atom types" in stripped:
                info["n_atom_types"] = int(parts[0])
            elif "bond types" in stripped:
                info["n_bond_types"] = int(parts[0])
            elif "angle types" in stripped:
                info["n_angle_types"] = int(parts[0])
            elif "dihedral types" in stripped:
                info["n_dihedral_types"] = int(parts[0])
            elif "improper types" in stripped:
                info["n_improper_types"] = int(parts[0])
    return info


def count_atoms_by_type(data_file: str) -> dict[int, int]:
    """
    Stream the Atoms section of a LAMMPS data file written with
    atom_style atomic (columns: atom-ID  type  x  y  z).

    Returns {atom_type: count, ...} for every type present.
    """
    counts: dict[int, int] = {}
    in_atoms = False
    with open(data_file) as fh:
        for line in fh:
            stripped = line.split("#")[0].strip()
            if not stripped:
                if in_atoms:
                    # blank line after Atoms section ends it
                    in_atoms = False
                continue
            if stripped == "Atoms":
                in_atoms = True
                continue
            # any new section keyword ends Atoms
            if in_atoms and stripped[0].isalpha():
                break
            if in_atoms:
                parts = stripped.split()
                if len(parts) >= 2:
                    atom_type = int(parts[1])
                    counts[atom_type] = counts.get(atom_type, 0) + 1
    return counts
