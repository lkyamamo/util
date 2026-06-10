#!/usr/bin/env python3
"""
Single-frame dipole moment calculator for LAMMPS custom dump files.

Reads one LAMMPS custom dump file (one frame), detects 3-atom molecules
(O + 2 nearest H) via vectorized cutoff-based neighbor search with minimum-image
PBC, computes per-molecule dipole moments, and writes the total moment as a
single line: timestep Mx My Mz

Box dimensions are read directly from the dump ITEM: BOX BOUNDS header.
Orthogonal boxes only; triclinic (xy xz yz tilt) is not supported.

Usage:
    python 1.calc_dipole_from_dump.py <dump_file>
        --cutoff <cutoff>
        --type-o <oxygen_type>
        --type-h <hydrogen_type>
        [--output <output_file>]
        [--log <log_file>]
        [--charge-o <charge>]
        [--charge-h <charge>]

Example:
    python 1.calc_dipole_from_dump.py dump.0.lammpstrj \\
        --cutoff 1.2 --type-o 1 --type-h 2 --output dipole_0.txt --log warn_0.log
"""

import sys
import argparse
import numpy as np


DEFAULT_TYPE_TO_CHARGE = {
    1: -0.813976,  # oxygen
    2:  0.406988,  # hydrogen
}


def parse_args():
    parser = argparse.ArgumentParser(
        description="Compute total dipole moment from a single LAMMPS custom dump frame."
    )
    parser.add_argument("dump_file", help="Path to the LAMMPS custom dump file (one frame)")
    parser.add_argument("--cutoff", type=float, required=True,
                        help="O-H bond cutoff distance in Angstroms (e.g. 1.2)")
    parser.add_argument("--type-o", type=int, required=True,
                        help="LAMMPS atom type integer for oxygen")
    parser.add_argument("--type-h", type=int, required=True,
                        help="LAMMPS atom type integer for hydrogen")
    parser.add_argument("--output", default=None,
                        help="Output file path (default: print to stdout)")
    parser.add_argument("--log", default=None,
                        help="Log file path for unassigned H warnings (default: print to stderr)")
    parser.add_argument("--charge-o", type=float, default=None,
                        help=f"Override charge for oxygen type (default: {DEFAULT_TYPE_TO_CHARGE[1]})")
    parser.add_argument("--charge-h", type=float, default=None,
                        help=f"Override charge for hydrogen type (default: {DEFAULT_TYPE_TO_CHARGE[2]})")
    return parser.parse_args()


def read_lammps_dump(filepath):
    """
    Parse a single-frame LAMMPS custom dump file.

    Expected columns (at minimum): id type x y z
    Optional column: q (per-atom charge)

    Returns:
        timestep  : int
        box_dims  : (la, lb, lc) tuple of floats  (orthogonal only)
        atoms     : dict  { atom_id: (atom_type, x, y, z, q_or_None) }
        has_charge: bool  — True if 'q' column was present in dump
    """
    atoms = {}
    timestep = None
    box_dims = None
    has_charge = False
    col_id = col_type = col_x = col_y = col_z = col_q = None

    xlo = xhi = ylo = yhi = zlo = zhi = None

    with open(filepath, "r") as f:
        lines = f.readlines()

    i = 0
    n = len(lines)

    while i < n:
        line = lines[i].strip()

        if line == "ITEM: TIMESTEP":
            i += 1
            timestep = int(lines[i].strip())

        elif line == "ITEM: NUMBER OF ATOMS":
            i += 1
            # num_atoms = int(lines[i].strip())  # not strictly needed

        elif line.startswith("ITEM: BOX BOUNDS"):
            # Next 3 lines: xlo xhi, ylo yhi, zlo zhi
            i += 1
            parts = lines[i].strip().split()
            xlo, xhi = float(parts[0]), float(parts[1])
            i += 1
            parts = lines[i].strip().split()
            ylo, yhi = float(parts[0]), float(parts[1])
            i += 1
            parts = lines[i].strip().split()
            zlo, zhi = float(parts[0]), float(parts[1])
            box_dims = (xhi - xlo, yhi - ylo, zhi - zlo)

        elif line.startswith("ITEM: ATOMS"):
            # Parse column headers
            headers = line.split()[2:]  # drop "ITEM:" and "ATOMS"
            col_id   = headers.index("id")   if "id"   in headers else None
            col_type = headers.index("type") if "type" in headers else None
            col_x    = headers.index("x")    if "x"    in headers else None
            col_y    = headers.index("y")    if "y"    in headers else None
            col_z    = headers.index("z")    if "z"    in headers else None
            col_q    = headers.index("q")    if "q"    in headers else None
            has_charge = col_q is not None

            if any(c is None for c in [col_id, col_type, col_x, col_y, col_z]):
                raise ValueError(
                    f"Dump file '{filepath}' is missing required columns. "
                    f"Found: {headers}. Need: id type x y z"
                )

            # Read atom lines until next ITEM or EOF
            i += 1
            while i < n and not lines[i].startswith("ITEM:"):
                parts = lines[i].strip().split()
                if parts:
                    atom_id   = int(parts[col_id])
                    atom_type = int(parts[col_type])
                    x = float(parts[col_x])
                    y = float(parts[col_y])
                    z = float(parts[col_z])
                    q = float(parts[col_q]) if has_charge else None
                    atoms[atom_id] = (atom_type, x, y, z, q)
                i += 1
            continue  # skip the i += 1 at the bottom

        i += 1

    if timestep is None:
        raise ValueError(f"Could not find ITEM: TIMESTEP in '{filepath}'")
    if box_dims is None:
        raise ValueError(f"Could not find ITEM: BOX BOUNDS in '{filepath}'")
    if not atoms:
        raise ValueError(f"No atom data found in '{filepath}'")

    return timestep, box_dims, atoms, has_charge


def apply_pbc_vectorized(dx, dy, dz, la, lb, lc):
    """Minimum-image PBC correction for numpy arrays."""
    dx = np.select([dx >= la / 2.0, dx <= -la / 2.0], [dx - la, dx + la], default=dx)
    dy = np.select([dy >= lb / 2.0, dy <= -lb / 2.0], [dy - lb, dy + lb], default=dy)
    dz = np.select([dz >= lc / 2.0, dz <= -lc / 2.0], [dz - lc, dz + lc], default=dz)
    return dx, dy, dz


def apply_pbc_scalar(dx, dy, dz, la, lb, lc):
    """Minimum-image PBC correction for scalar values."""
    if   dx >=  la / 2.0: dx -= la
    elif dx <= -la / 2.0: dx += la
    if   dy >=  lb / 2.0: dy -= lb
    elif dy <= -lb / 2.0: dy += lb
    if   dz >=  lc / 2.0: dz -= lc
    elif dz <= -lc / 2.0: dz += lc
    return dx, dy, dz


def detect_bonds(atoms, type_o, type_h, cutoff, box_dims):
    """
    Find O + 2 nearest H molecules via vectorized cutoff neighbor search.

    Returns:
        bonds          : dict { bond_id: [o_atom_id, h_atom_id1, h_atom_id2] }
        unassigned_h   : list of H atom IDs not claimed by any oxygen
        bond_stats     : dict with diagnostic counts
    """
    la, lb, lc = box_dims
    cutoff_sq = cutoff * cutoff

    o_atoms = []   # list of (atom_id, x, y, z)
    h_ids = []     # list of atom_id
    h_pos = []     # list of [x, y, z]

    for atom_id, (atom_type, x, y, z, _q) in atoms.items():
        if atom_type == type_o:
            o_atoms.append((atom_id, x, y, z))
        elif atom_type == type_h:
            h_ids.append(atom_id)
            h_pos.append([x, y, z])

    if not h_pos:
        return {}, list(h_ids), {
            "bonds_created": 0, "bonds_missing": len(o_atoms),
            "type_o_count": len(o_atoms), "type_h_count": 0,
        }

    h_pos_arr = np.array(h_pos, dtype=np.float32)
    h_ids_arr = np.array(h_ids, dtype=np.int32)

    bonds = {}
    bond_id = 1
    assigned_h_indices = set()   # indices into h_ids_arr
    bonds_created = 0
    bonds_missing = 0

    for o_id, ox, oy, oz in o_atoms:
        dx = h_pos_arr[:, 0] - ox
        dy = h_pos_arr[:, 1] - oy
        dz = h_pos_arr[:, 2] - oz

        dx, dy, dz = apply_pbc_vectorized(dx, dy, dz, la, lb, lc)
        dist_sq = dx * dx + dy * dy + dz * dz

        within_mask = dist_sq <= cutoff_sq
        n_within = int(np.sum(within_mask))

        if n_within == 2:
            idxs = np.where(within_mask)[0]
            i1, i2 = int(idxs[0]), int(idxs[1])

        elif n_within > 2:
            within_idxs = np.where(within_mask)[0]
            top2 = np.argpartition(dist_sq[within_mask], 2)[:2]
            i1, i2 = int(within_idxs[top2[0]]), int(within_idxs[top2[1]])

        else:
            # Fewer than 2 within cutoff — fall back to 2 globally nearest
            if len(h_ids_arr) >= 2:
                sorted_idxs = np.argsort(dist_sq)
                i1, i2 = int(sorted_idxs[0]), int(sorted_idxs[1])
                # Only accept if both are actually within cutoff
                if dist_sq[i1] > cutoff_sq or dist_sq[i2] > cutoff_sq:
                    bonds_missing += 1
                    continue
            else:
                bonds_missing += 1
                continue

        h1 = int(h_ids_arr[i1])
        h2 = int(h_ids_arr[i2])
        bonds[bond_id] = [o_id, h1, h2]
        bond_id += 1
        bonds_created += 1
        assigned_h_indices.add(i1)
        assigned_h_indices.add(i2)

    all_h_indices = set(range(len(h_ids_arr)))
    unassigned_h_indices = all_h_indices - assigned_h_indices
    unassigned_h = [int(h_ids_arr[i]) for i in sorted(unassigned_h_indices)]

    bond_stats = {
        "bonds_created": bonds_created,
        "bonds_missing": bonds_missing,
        "type_o_count": len(o_atoms),
        "type_h_count": len(h_ids),
    }

    return bonds, unassigned_h, bond_stats


def calc_total_dipole(atoms, bonds, box_dims, type_to_charge):
    """
    Sum per-molecule dipole moments.

    Uses the same formula as the original pipeline:
        dipole_x = (x_O - x_H1) * q_H - (x_H2 - x_O) * q_H
    with minimum-image PBC on each bond vector.

    Returns:
        total  : [Mx, My, Mz]
    """
    la, lb, lc = box_dims
    total = [0.0, 0.0, 0.0]

    for _bond_id, (o_id, h1_id, h2_id) in bonds.items():
        type_o_local, x1, y1, z1, q1_raw = atoms[o_id]
        type_h1,      x2, y2, z2, q2_raw = atoms[h1_id]
        type_h2,      x3, y3, z3, q3_raw = atoms[h2_id]

        q2 = q2_raw if q2_raw is not None else type_to_charge.get(type_h1, 0.0)
        q3 = q3_raw if q3_raw is not None else type_to_charge.get(type_h2, 0.0)

        # dx12: O -> H1 displacement (atom1 - atom2 convention from original)
        dx12, dy12, dz12 = apply_pbc_scalar(x1 - x2, y1 - y2, z1 - z2, la, lb, lc)
        # dx31: H2 -> O displacement
        dx31, dy31, dz31 = apply_pbc_scalar(x3 - x1, y3 - y1, z3 - z1, la, lb, lc)

        total[0] += dx12 * q2 - dx31 * q3
        total[1] += dy12 * q2 - dy31 * q3
        total[2] += dz12 * q2 - dz31 * q3

    return total


def write_log(log_dest, timestep, unassigned_h, bond_stats):
    """Write unassigned H warning entry. log_dest may be a file path or None (stderr)."""
    if not unassigned_h:
        return

    n = len(unassigned_h)
    ids_str = " ".join(str(a) for a in unassigned_h)
    msg = (
        f"WARN [timestep={timestep}] "
        f"{n} unassigned H atom(s) (not within cutoff of any O): [{ids_str}]  "
        f"| bonds_created={bond_stats['bonds_created']} "
        f"bonds_missing={bond_stats['bonds_missing']} "
        f"n_O={bond_stats['type_o_count']} n_H={bond_stats['type_h_count']}\n"
    )

    if log_dest is None:
        sys.stderr.write(msg)
    else:
        with open(log_dest, "a") as lf:
            lf.write(msg)


def main():
    args = parse_args()

    type_to_charge = dict(DEFAULT_TYPE_TO_CHARGE)
    if args.charge_o is not None:
        type_to_charge[args.type_o] = args.charge_o
    if args.charge_h is not None:
        type_to_charge[args.type_h] = args.charge_h

    timestep, box_dims, atoms, has_charge = read_lammps_dump(args.dump_file)

    bonds, unassigned_h, bond_stats = detect_bonds(
        atoms, args.type_o, args.type_h, args.cutoff, box_dims
    )

    total = calc_total_dipole(atoms, bonds, box_dims, type_to_charge)

    result_line = f"{timestep}  {total[0]:14.6f}  {total[1]:14.6f}  {total[2]:14.6f}\n"

    if args.output is None:
        sys.stdout.write(result_line)
    else:
        with open(args.output, "w") as out:
            out.write(result_line)

    if unassigned_h:
        write_log(args.log, timestep, unassigned_h, bond_stats)


if __name__ == "__main__":
    main()
