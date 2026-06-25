#!/usr/bin/env python3
"""
rotate_lammps.py

Rotate a LAMMPS data file 90° about the positive Y axis and recenter
the system so all box dimensions start at 0.

Rotation: +90° about Y
    x' =  z
    y' =  y
    z' = -x

Applied to atom positions, image flags, and velocity vectors.

Assumptions
-----------
- OVITO atomic format:  id type x y z ix iy iz
- Velocities section:   id vx vy vz
- Rectangular (orthogonal) box only
- Masses section may have inline  # comments

Optimized for large files (10M+ atoms) via numpy vectorized operations.
np.loadtxt/savetxt require numpy >= 1.23 for fast C-based I/O.

Usage
-----
    python rotate_lammps.py input.data output.data
"""

import re
import numpy as np
import argparse


def main():
    parser = argparse.ArgumentParser(
        description="Rotate a LAMMPS data file 90° about the +Y axis and recenter.")
    parser.add_argument("input",  help="Input LAMMPS data file")
    parser.add_argument("output", help="Output LAMMPS data file")
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Read
    # ------------------------------------------------------------------
    print(f"Reading: {args.input}")

    n_atoms      = None
    n_atom_types = None
    box          = {}
    masses       = {}   # type_id -> (mass_float, comment_str)
    header_comment = ""

    count_re = re.compile(r'^(\d+)\s+(atoms|atom types)$')
    box_re   = re.compile(
        r'^([+-]?[\d.eE+-]+)\s+([+-]?[\d.eE+-]+)\s+(xlo xhi|ylo yhi|zlo zhi)$')

    with open(args.input, 'r') as fin:
        header_comment = fin.readline().rstrip()

        # --- header + Masses ---
        in_masses = False
        for line in fin:
            s = line.strip()

            if not s:
                in_masses = False
                continue

            if in_masses:
                core    = s.split('#')[0].strip()
                comment = s[s.find('#'):] if '#' in s else ''
                if core:
                    parts = core.split()
                    masses[int(parts[0])] = (float(parts[1]), comment.strip())
                continue

            m = count_re.match(s)
            if m:
                if m.group(2) == 'atoms':
                    n_atoms = int(m.group(1))
                else:
                    n_atom_types = int(m.group(1))
                continue

            m = box_re.match(s)
            if m:
                lo, hi, label = float(m.group(1)), float(m.group(2)), m.group(3)
                axis = label[0]
                box[axis + 'lo'] = lo
                box[axis + 'hi'] = hi
                continue

            if s == 'Masses':
                in_masses = True
                continue

            if s.startswith('Atoms'):
                fin.readline()  # skip blank line after section header
                break

        assert n_atoms is not None, "Could not find atom count in header"

        # --- atom data: id type x y z ix iy iz ---
        print(f"Loading {n_atoms:,} atoms ...")
        atom_data = np.loadtxt(fin, max_rows=n_atoms)
        # shape: (n_atoms, 8)  cols: id type x y z ix iy iz

        # --- velocities (optional): id vx vy vz ---
        vel_data = None
        for line in fin:
            s = line.strip()
            if not s:
                continue
            if s.startswith('Velocities'):
                fin.readline()  # skip blank line after section header
                print(f"Loading {n_atoms:,} velocities ...")
                vel_data = np.loadtxt(fin, max_rows=n_atoms)
                # shape: (n_atoms, 4)  cols: id vx vy vz
                break

    # ------------------------------------------------------------------
    # Transform
    # ------------------------------------------------------------------
    xlo, xhi = box['xlo'], box['xhi']
    ylo, yhi = box['ylo'], box['yhi']
    zlo, zhi = box['zlo'], box['zhi']

    print(f"Box (in) : x [{xlo:.4f}, {xhi:.4f}]  "
          f"y [{ylo:.4f}, {yhi:.4f}]  z [{zlo:.4f}, {zhi:.4f}]")

    # After +90° Y rotation x'=z, y'=y, z'=-x, the axis ranges are:
    #   x' in [zlo, zhi]   → shift_x = -zlo
    #   y' in [ylo, yhi]   → shift_y = -ylo
    #   z' in [-xhi, -xlo] → shift_z = +xhi
    shift_x = -zlo
    shift_y = -ylo
    shift_z =  xhi

    new_box = {
        'xlo': 0.0, 'xhi': zhi - zlo,
        'ylo': 0.0, 'yhi': yhi - ylo,
        'zlo': 0.0, 'zhi': xhi - xlo,
    }

    print(f"Box (out): x [{new_box['xlo']:.4f}, {new_box['xhi']:.4f}]  "
          f"y [{new_box['ylo']:.4f}, {new_box['yhi']:.4f}]  "
          f"z [{new_box['zlo']:.4f}, {new_box['zhi']:.4f}]")

    # Rotate + shift positions (vectorized, in-place)
    x = atom_data[:, 2].copy()
    y = atom_data[:, 3].copy()
    z = atom_data[:, 4].copy()
    atom_data[:, 2] =  z + shift_x
    atom_data[:, 3] =  y + shift_y
    atom_data[:, 4] = -x + shift_z

    # Rotate image flags (same mapping, no shift)
    ix = atom_data[:, 5].copy()
    iy = atom_data[:, 6].copy()
    iz = atom_data[:, 7].copy()
    atom_data[:, 5] =  iz
    atom_data[:, 6] =  iy
    atom_data[:, 7] = -ix

    # Rotate velocity vectors (no shift — vectors, not positions)
    if vel_data is not None:
        vx = vel_data[:, 1].copy()
        vy = vel_data[:, 2].copy()
        vz = vel_data[:, 3].copy()
        vel_data[:, 1] =  vz
        vel_data[:, 2] =  vy
        vel_data[:, 3] = -vx

    # ------------------------------------------------------------------
    # Write
    # ------------------------------------------------------------------
    print(f"Writing: {args.output}")
    b = new_box
    with open(args.output, 'w') as fout:
        fout.write(f"{header_comment} (rotated 90deg about +Y, recentered)\n\n")
        fout.write(f"{n_atoms} atoms\n")
        fout.write(f"{n_atom_types} atom types\n\n")
        fout.write(f"{b['xlo']:.10g} {b['xhi']:.10g} xlo xhi\n")
        fout.write(f"{b['ylo']:.10g} {b['yhi']:.10g} ylo yhi\n")
        fout.write(f"{b['zlo']:.10g} {b['zhi']:.10g} zlo zhi\n\n")

        if masses:
            fout.write("Masses\n\n")
            for atype in sorted(masses):
                mass_val, comment = masses[atype]
                if comment:
                    fout.write(f"{atype} {mass_val}  {comment}\n")
                else:
                    fout.write(f"{atype} {mass_val}\n")
            fout.write("\n")

        fout.write("Atoms  # atomic\n\n")
        np.savetxt(fout, atom_data,
                   fmt=['%d', '%d', '%.10g', '%.10g', '%.10g', '%d', '%d', '%d'])
        fout.write("\n")

        if vel_data is not None:
            fout.write("Velocities\n\n")
            np.savetxt(fout, vel_data, fmt=['%d', '%.10g', '%.10g', '%.10g'])
            fout.write("\n")

    print("Done.")


if __name__ == "__main__":
    main()
