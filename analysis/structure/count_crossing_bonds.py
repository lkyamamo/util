"""
count_crossing_bonds.py
-----------------------
Counts the number of bonds crossing each face of a LAMMPS simulation box.
Bonds are inferred from a distance cutoff (with full PBC support).

A bond "crosses" a face if the two atoms are connected via a PBC image
(i.e. the minimum-image vector wraps through that periodic boundary).

Usage:
    python count_crossing_bonds.py structure.lammps --cutoff 2.0
    python count_crossing_bonds.py structure.lammps --cutoff 2.0 --types 1 2

Supports orthogonal boxes. Triclinic (tilt) support is included but
assumes small tilt factors (standard LAMMPS convention).
"""

import argparse
import numpy as np
from itertools import combinations
from collections import defaultdict


# ---------------------------------------------------------------------------
# LAMMPS data file parser
# ---------------------------------------------------------------------------

def parse_lammps_data(filepath):
    """Parse a LAMMPS data file and return atoms and box info."""
    atoms = {}       # atom_id -> {'type': int, 'x': float, 'y': float, 'z': float}
    box = {}         # xlo, xhi, ylo, yhi, zlo, zhi, xy, xz, yz
    box.update({'xy': 0.0, 'xz': 0.0, 'yz': 0.0})  # default ortho

    with open(filepath, 'r') as f:
        lines = f.readlines()

    section = None
    for line in lines:
        line = line.strip()
        if not line or line.startswith('#'):
            continue

        # Detect box bounds
        if 'xlo xhi' in line:
            parts = line.split()
            box['xlo'], box['xhi'] = float(parts[0]), float(parts[1])
        elif 'ylo yhi' in line:
            parts = line.split()
            box['ylo'], box['yhi'] = float(parts[0]), float(parts[1])
        elif 'zlo zhi' in line:
            parts = line.split()
            box['zlo'], box['zhi'] = float(parts[0]), float(parts[1])
        elif 'xy xz yz' in line:
            parts = line.split()
            box['xy'], box['xz'], box['yz'] = float(parts[0]), float(parts[1]), float(parts[2])

        # Detect sections
        elif line == 'Atoms' or line.startswith('Atoms '):
            section = 'atoms'
            continue
        elif line in ('Velocities', 'Bonds', 'Angles', 'Dihedrals',
                      'Impropers', 'Masses', 'Pair Coeffs',
                      'Bond Coeffs', 'Angle Coeffs'):
            section = None
            continue

        elif section == 'atoms':
            parts = line.split()
            if len(parts) >= 5:
                try:
                    # Support both atom styles:
                    # full:   id mol type q x y z [ix iy iz]
                    # atomic: id type x y z [ix iy iz]
                    atom_id = int(parts[0])
                    # Detect style: if parts[3] looks like a charge (float) use full style
                    try:
                        _ = float(parts[3])
                        # Could be full style (id mol type charge x y z)
                        if len(parts) >= 7:
                            atype = int(parts[2])
                            x, y, z = float(parts[4]), float(parts[5]), float(parts[6])
                        else:
                            raise ValueError
                    except (ValueError, IndexError):
                        atype = int(parts[1])
                        x, y, z = float(parts[2]), float(parts[3]), float(parts[4])

                    atoms[atom_id] = {'type': atype, 'x': x, 'y': y, 'z': z}
                except (ValueError, IndexError):
                    continue  # skip malformed lines

    if not atoms:
        raise ValueError("No atoms found. Check that your file has an 'Atoms' section.")
    if 'xlo' not in box:
        raise ValueError("Box bounds not found in LAMMPS data file.")

    return atoms, box


# ---------------------------------------------------------------------------
# PBC minimum image & bond-crossing detection
# ---------------------------------------------------------------------------

def build_cell_matrix(box):
    """Build 3x3 cell matrix from box dict (handles triclinic tilt)."""
    lx = box['xhi'] - box['xlo']
    ly = box['yhi'] - box['ylo']
    lz = box['zhi'] - box['zlo']
    xy, xz, yz = box['xy'], box['xz'], box['yz']

    # LAMMPS triclinic convention
    cell = np.array([
        [lx,  0,   0],
        [xy,  ly,  0],
        [xz,  yz,  lz]
    ], dtype=float)
    return cell


def wrap_to_box(pos, box):
    """Wrap a position back into the primary box (fractional coords)."""
    lx = box['xhi'] - box['xlo']
    ly = box['yhi'] - box['ylo']
    lz = box['zhi'] - box['zlo']
    x = pos[0] - box['xlo']
    y = pos[1] - box['ylo']
    z = pos[2] - box['zlo']
    # For triclinic, simple wrapping per LAMMPS convention
    z -= lz * np.round(z / lz)
    y -= ly * np.round(y / ly)
    x -= lx * np.round(x / lx)
    return np.array([x + box['xlo'], y + box['ylo'], z + box['zlo']])


def minimum_image_vector(r1, r2, box):
    """
    Compute the minimum image vector from atom1 to atom2,
    also returning the PBC image shift (nx, ny, nz).
    Returns (dr_cartesian, (nx, ny, nz))
    """
    lx = box['xhi'] - box['xlo']
    ly = box['yhi'] - box['ylo']
    lz = box['zhi'] - box['zlo']
    xy, xz, yz = box['xy'], box['xz'], box['yz']

    dx = r2[0] - r1[0]
    dy = r2[1] - r1[1]
    dz = r2[2] - r1[2]

    # Apply minimum image in z, then y, then x (LAMMPS triclinic order)
    nz = round(dz / lz)
    dz -= nz * lz
    dy -= nz * yz
    dx -= nz * xz

    ny = round(dy / ly)
    dy -= ny * ly
    dx -= ny * xy

    nx = round(dx / lx)
    dx -= nx * lx

    return np.array([dx, dy, dz]), (int(nx), int(ny), int(nz))


# ---------------------------------------------------------------------------
# Main bond-counting logic
# ---------------------------------------------------------------------------

def find_crossing_bonds(atoms, box, cutoff, atom_types=None):
    """
    For every pair of atoms within cutoff distance (using PBC),
    check whether the bond crosses a periodic boundary.

    Returns a dict counting bonds crossing each face pair:
        {'+x/-x': n, '+y/-y': n, '+z/-z': n}

    If atom_types is specified, only consider atoms of those types.
    """
    ids = list(atoms.keys())
    if atom_types:
        ids = [i for i in ids if atoms[i]['type'] in atom_types]

    positions = np.array([[atoms[i]['x'], atoms[i]['y'], atoms[i]['z']] for i in ids])
    n = len(ids)

    lx = box['xhi'] - box['xlo']
    ly = box['yhi'] - box['ylo']
    lz = box['zhi'] - box['zlo']

    # Count bonds crossing each axis boundary
    cross_counts = defaultdict(int)
    cross_bonds = defaultdict(list)   # stores (id1, id2, type1, type2)

    cutoff2 = cutoff ** 2

    print(f"  Checking {n} atoms ({n*(n-1)//2:,} pairs)...")

    for i in range(n):
        for j in range(i + 1, n):
            dr, (nx, ny, nz) = minimum_image_vector(positions[i], positions[j], box)
            dist2 = dr[0]**2 + dr[1]**2 + dr[2]**2

            if dist2 <= cutoff2:
                # This is a bond. Does it cross a PBC boundary?
                if nx != 0:
                    label = 'x-face (±x)'
                    cross_counts[label] += abs(nx)
                    cross_bonds[label].append((ids[i], ids[j],
                                               atoms[ids[i]]['type'],
                                               atoms[ids[j]]['type']))
                if ny != 0:
                    label = 'y-face (±y)'
                    cross_counts[label] += abs(ny)
                    cross_bonds[label].append((ids[i], ids[j],
                                               atoms[ids[i]]['type'],
                                               atoms[ids[j]]['type']))
                if nz != 0:
                    label = 'z-face (±z)'
                    cross_counts[label] += abs(nz)
                    cross_bonds[label].append((ids[i], ids[j],
                                               atoms[ids[i]]['type'],
                                               atoms[ids[j]]['type']))

    return cross_counts, cross_bonds


def compute_face_areas(box):
    """Compute the area of each face of the simulation box."""
    lx = box['xhi'] - box['xlo']
    ly = box['yhi'] - box['ylo']
    lz = box['zhi'] - box['zlo']
    xy, xz, yz = box['xy'], box['xz'], box['yz']

    # Face normals (for orthogonal: just lx*ly etc.)
    # For triclinic, use cross products of cell vectors
    a = np.array([lx, 0, 0])
    b = np.array([xy, ly, 0])
    c = np.array([xz, yz, lz])

    area_xy = np.linalg.norm(np.cross(a, b))   # z-face area
    area_xz = np.linalg.norm(np.cross(a, c))   # y-face area
    area_yz = np.linalg.norm(np.cross(b, c))   # x-face area

    return {
        'x-face (±x)': area_yz,
        'y-face (±y)': area_xz,
        'z-face (±z)': area_xy,
    }


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Count bonds crossing PBC faces of a LAMMPS simulation box.')
    parser.add_argument('datafile', help='LAMMPS data file')
    parser.add_argument('--cutoff', type=float, required=True,
                        help='Bond cutoff distance (same units as data file)')
    parser.add_argument('--types', type=int, nargs='+', default=None,
                        help='Only consider atoms of these types (optional)')
    parser.add_argument('--per-area', action='store_true',
                        help='Also report bond density per unit area')
    args = parser.parse_args()

    print(f"\n{'='*55}")
    print(f"  Bond-crossing counter (PBC faces)")
    print(f"{'='*55}")
    print(f"  File   : {args.datafile}")
    print(f"  Cutoff : {args.cutoff}")
    print(f"  Types  : {args.types if args.types else 'all'}")
    print(f"{'='*55}\n")

    # Parse
    print("Parsing LAMMPS data file...")
    atoms, box = parse_lammps_data(args.datafile)
    print(f"  Found {len(atoms)} atoms")
    print(f"  Box: x=[{box['xlo']:.3f}, {box['xhi']:.3f}]  "
          f"y=[{box['ylo']:.3f}, {box['yhi']:.3f}]  "
          f"z=[{box['zlo']:.3f}, {box['zhi']:.3f}]")
    if any(box[k] != 0 for k in ('xy', 'xz', 'yz')):
        print(f"  Tilt: xy={box['xy']:.4f}  xz={box['xz']:.4f}  yz={box['yz']:.4f}")

    # Count
    print("\nFinding crossing bonds...")
    cross_counts, cross_bonds = find_crossing_bonds(atoms, box, args.cutoff, args.types)

    # Face areas
    face_areas = compute_face_areas(box)

    # Report
    print(f"\n{'='*55}")
    print(f"  Results")
    print(f"{'='*55}")

    all_labels = ['x-face (±x)', 'y-face (±y)', 'z-face (±z)']
    for label in all_labels:
        count = cross_counts.get(label, 0)
        area = face_areas[label]
        print(f"\n  {label}:")
        print(f"    Bonds crossing : {count}")
        print(f"    Face area      : {area:.4f} (length units)^2")
        if args.per_area and area > 0:
            print(f"    Bond density   : {count/area:.6f} bonds / (length unit)^2")

        # Type breakdown
        if label in cross_bonds:
            type_pairs = defaultdict(int)
            for _, _, t1, t2 in cross_bonds[label]:
                key = tuple(sorted([t1, t2]))
                type_pairs[key] += 1
            if type_pairs:
                print(f"    By type pair   :", end='')
                for (t1, t2), n in sorted(type_pairs.items()):
                    print(f"  ({t1}-{t2}): {n}", end='')
                print()

    total = sum(cross_counts.values())
    print(f"\n  Total crossing bonds (all faces): {total}")
    print(f"{'='*55}\n")


if __name__ == '__main__':
    main()