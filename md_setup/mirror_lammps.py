#!/usr/bin/env python3
"""
mirror_lammps.py

Mirror a LAMMPS data file across one or more planes and recenter
the system so all box dimensions start at 0.

Mirror planes:
    xy  ->  z = -z  (negate z component)
    xz  ->  y = -y  (negate y component)
    yz  ->  x = -x  (negate x component)

Applied to atom positions, image flags, and velocity vectors.

Assumptions
-----------
- Rectangular (orthogonal) box only
- Masses section may have inline # comments
- Image flags (ix iy iz) are present in the Atoms section

Optimized for large files (10M+ atoms) via numpy vectorized operations.
np.loadtxt/savetxt require numpy >= 1.23 for fast C-based I/O.

Usage
-----
    python mirror_lammps.py input.data output.data --planes yz
    python mirror_lammps.py input.data output.data --planes xy xz
    python mirror_lammps.py input.data output.data --planes xy xz yz

Arguments
---------
    input           Input LAMMPS data file
    output          Output LAMMPS data file
    --planes        One or more of: xy, xz, yz  (can specify multiple)
    --atom_style    LAMMPS atom style (default: atomic)
"""

import re
import numpy as np
import argparse


# Per atom style: (pos_cols, img_cols, fmt_list)
#   pos_cols  — (cx, cy, cz) 0-based column indices for x, y, z
#   img_cols  — (icx, icy, icz) 0-based column indices for ix, iy, iz
#   fmt_list  — savetxt format per column
STYLE_INFO = {
    # id type x y z ix iy iz
    'atomic': (
        (2, 3, 4), (5, 6, 7),
        ['%d', '%d', '%.10g', '%.10g', '%.10g', '%d', '%d', '%d'],
    ),
    # id type charge x y z ix iy iz
    'charge': (
        (3, 4, 5), (6, 7, 8),
        ['%d', '%d', '%.10g', '%.10g', '%.10g', '%.10g', '%d', '%d', '%d'],
    ),
    # id mol type x y z ix iy iz
    'molecular': (
        (3, 4, 5), (6, 7, 8),
        ['%d', '%d', '%d', '%.10g', '%.10g', '%.10g', '%d', '%d', '%d'],
    ),
    # id mol type charge x y z ix iy iz
    'full': (
        (4, 5, 6), (7, 8, 9),
        ['%d', '%d', '%d', '%.10g', '%.10g', '%.10g', '%.10g', '%d', '%d', '%d'],
    ),
}


def main():
    parser = argparse.ArgumentParser(
        description="Mirror a LAMMPS data file across one or more planes and recenter.")
    parser.add_argument("input",  help="Input LAMMPS data file")
    parser.add_argument("output", help="Output LAMMPS data file")
    parser.add_argument(
        "--planes", nargs="+", choices=["xy", "xz", "yz"], required=True,
        help="Mirror plane(s): xy, xz, yz (can specify multiple)")
    parser.add_argument(
        "--atom_style", default="atomic", choices=sorted(STYLE_INFO),
        help="LAMMPS atom style (default: atomic)")
    args = parser.parse_args()

    planes = list(set(args.planes))  # deduplicate
    print(f"Input        : {args.input}")
    print(f"Output       : {args.output}")
    print(f"Mirror planes: {', '.join(sorted(planes))}")
    print(f"Atom style   : {args.atom_style}")

    pos_cols, img_cols, fmt = STYLE_INFO[args.atom_style]
    cx,  cy,  cz  = pos_cols
    icx, icy, icz = img_cols

    # Mirror signs: +1 = unchanged, -1 = mirrored
    sx = -1 if 'yz' in planes else 1
    sy = -1 if 'xz' in planes else 1
    sz = -1 if 'xy' in planes else 1

    # ------------------------------------------------------------------
    # Read
    # ------------------------------------------------------------------
    print(f"Reading: {args.input}")

    n_atoms        = None
    n_atom_types   = None
    box            = {}
    masses         = {}    # type_id -> (mass_float, comment_str)
    header_comment = ""
    atoms_header   = "Atoms"  # preserved verbatim (may contain "# atomic" etc.)

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
                atoms_header = s
                fin.readline()  # skip blank line after section header
                break

        assert n_atoms is not None, "Could not find atom count in header"

        # --- atom data ---
        print(f"Loading {n_atoms:,} atoms ...")
        atom_data = np.loadtxt(fin, max_rows=n_atoms)

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

    # Recentering shift per axis so the output box starts at 0.
    #   sign=+1: new_i =  i - lo   → shift = -lo
    #   sign=-1: new_i = -i + hi   → shift = +hi  (maps [lo,hi] → [0, hi-lo])
    shift_x = xhi if sx < 0 else -xlo
    shift_y = yhi if sy < 0 else -ylo
    shift_z = zhi if sz < 0 else -zlo

    new_box = {
        'xlo': 0.0, 'xhi': xhi - xlo,
        'ylo': 0.0, 'yhi': yhi - ylo,
        'zlo': 0.0, 'zhi': zhi - zlo,
    }

    print(f"Box (out): x [{new_box['xlo']:.4f}, {new_box['xhi']:.4f}]  "
          f"y [{new_box['ylo']:.4f}, {new_box['yhi']:.4f}]  "
          f"z [{new_box['zlo']:.4f}, {new_box['zhi']:.4f}]")

    # Mirror + shift positions (vectorized, in-place)
    x = atom_data[:, cx].copy()
    y = atom_data[:, cy].copy()
    z = atom_data[:, cz].copy()
    atom_data[:, cx] = sx * x + shift_x
    atom_data[:, cy] = sy * y + shift_y
    atom_data[:, cz] = sz * z + shift_z

    # Mirror image flags (no shift — lattice-vector counts, not positions)
    atom_data[:, icx] = sx * atom_data[:, icx]
    atom_data[:, icy] = sy * atom_data[:, icy]
    atom_data[:, icz] = sz * atom_data[:, icz]

    # Mirror velocity vectors (no shift — vectors, not positions)
    if vel_data is not None:
        vel_data[:, 1] = sx * vel_data[:, 1]
        vel_data[:, 2] = sy * vel_data[:, 2]
        vel_data[:, 3] = sz * vel_data[:, 3]

    # ------------------------------------------------------------------
    # Write
    # ------------------------------------------------------------------
    print(f"Writing: {args.output}")
    plane_str = '+'.join(sorted(planes))
    b = new_box
    with open(args.output, 'w') as fout:
        fout.write(f"{header_comment} (mirrored {plane_str}, recentered)\n\n")
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

        fout.write(f"{atoms_header}\n\n")
        np.savetxt(fout, atom_data, fmt=fmt)
        fout.write("\n")

        if vel_data is not None:
            fout.write("Velocities\n\n")
            np.savetxt(fout, vel_data, fmt=['%d', '%.10g', '%.10g', '%.10g'])
            fout.write("\n")

    print("Done.")


if __name__ == "__main__":
    main()
