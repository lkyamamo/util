#!/usr/bin/env python3
"""
mirror_lammps.py

Mirror a LAMMPS data file across one or more planes using the OVITO Python API.
Mirrors are applied to atom positions, velocities, and the simulation cell.

Mirror planes:
    xy  ->  z = -z  (negate z component)
    xz  ->  y = -y  (negate y component)
    yz  ->  x = -x  (negate x component)

After mirroring the system is recentered so all box dimensions are positive.

Usage
-----
    python mirror_lammps.py input.data output.data --planes yz
    python mirror_lammps.py input.data output.data --planes xy xz
    python mirror_lammps.py input.data output.data --planes xy xz yz

Arguments
---------
    input           Input LAMMPS data file
    output          Output LAMMPS data file
    --planes        One or more of: xy, xz, yz  (order does not matter)
    --atom_style    LAMMPS atom style (default: full)
"""

import argparse
import numpy as np

try:
    from ovito.io import import_file, export_file
    from ovito.modifiers import AffineTransformationModifier
except ImportError:
    raise SystemExit(
        "OVITO is required: pip install ovito"
    )


def _fix_cell_header(filepath):
    """
    OVITO 3.9+ writes cell geometry using the newer LAMMPS avec/bvec/cvec/
    abc origin format even for orthogonal boxes, which LAMMPS treats as
    triclinic.  This rewrites that 4-line block as the standard 3-line
    xlo xhi / ylo yhi / zlo zhi format.
    """
    with open(filepath) as f:
        lines = f.readlines()

    avec = bvec = cvec = origin = None
    avec_idx = bvec_idx = cvec_idx = origin_idx = None

    for i, line in enumerate(lines):
        parts = line.split()
        if len(parts) >= 4 and parts[3] == 'avec':
            avec, avec_idx = list(map(float, parts[:3])), i
        elif len(parts) >= 4 and parts[3] == 'bvec':
            bvec, bvec_idx = list(map(float, parts[:3])), i
        elif len(parts) >= 4 and parts[3] == 'cvec':
            cvec, cvec_idx = list(map(float, parts[:3])), i
        elif len(parts) >= 5 and parts[3] == 'abc' and parts[4] == 'origin':
            origin, origin_idx = list(map(float, parts[:3])), i

    if avec is None:
        return  # already in xlo xhi format

    xlo, xhi = origin[0], origin[0] + avec[0]
    ylo, yhi = origin[1], origin[1] + bvec[1]
    zlo, zhi = origin[2], origin[2] + cvec[2]

    skip = {avec_idx, bvec_idx, cvec_idx}
    new_lines = []
    for i, line in enumerate(lines):
        if i in skip:
            continue
        elif i == origin_idx:
            new_lines.append(f"{xlo:.10g} {xhi:.10g} xlo xhi\n")
            new_lines.append(f"{ylo:.10g} {yhi:.10g} ylo yhi\n")
            new_lines.append(f"{zlo:.10g} {zhi:.10g} zlo zhi\n")
        else:
            new_lines.append(line)

    with open(filepath, 'w') as f:
        f.writelines(new_lines)



# ---------------------------------------------------------------------------
# Mirror matrices (4x3 affine — OVITO convention)
# Each row is [m00 m01 m02 translation] for x, y, z respectively
# ---------------------------------------------------------------------------

MIRROR_MATRICES = {
    # XY plane: negate z
    "xy": np.array([
        [ 1,  0,  0,  0],
        [ 0,  1,  0,  0],
        [ 0,  0, -1,  0],
    ], dtype=float),

    # XZ plane: negate y
    "xz": np.array([
        [ 1,  0,  0,  0],
        [ 0, -1,  0,  0],
        [ 0,  0,  1,  0],
    ], dtype=float),

    # YZ plane: negate x
    "yz": np.array([
        [-1,  0,  0,  0],
        [ 0,  1,  0,  0],
        [ 0,  0,  1,  0],
    ], dtype=float),
}


def combined_mirror_matrix(planes):
    """
    Compose mirror matrices for the requested planes into a single
    4x3 affine matrix. Order doesn't matter for axis-aligned mirrors
    since they commute (diagonal matrices).
    """
    # Build full 3x3 diagonal from selected planes
    scale = np.ones(3)
    for plane in planes:
        m = MIRROR_MATRICES[plane]
        # Extract the diagonal of the 3x3 part
        scale *= np.diag(m[:, :3])

    # Pack into OVITO's 3x4 affine format
    affine = np.zeros((3, 4))
    affine[0, 0] = scale[0]
    affine[1, 1] = scale[1]
    affine[2, 2] = scale[2]
    return affine


def recenter_matrix(pipeline, mirror_matrix):
    """
    Compute the translation needed to shift the box so that all
    dimensions are positive (xlo=ylo=zlo=0) after mirroring,
    and return an updated affine matrix with that translation baked in.
    """
    # Get the current simulation cell (3x3 column vectors + origin)
    data = pipeline.compute()
    cell = data.cell

    # The 3 cell vectors after mirroring
    R = mirror_matrix[:, :3]
    new_vectors = R @ cell.matrix[:, :3]        # rotate the 3 cell edge vectors
    new_origin  = R @ cell.matrix[:, 3]         # rotate the origin

    # For each axis find the minimum coordinate extent
    # (cell could be flipped so hi < lo after negation)
    mins = np.minimum(new_origin, new_origin + new_vectors.sum(axis=1))

    # Translation to push lo -> 0
    translation = -mins

    mat = mirror_matrix.copy()
    mat[:, 3] += translation
    return mat


def main():
    parser = argparse.ArgumentParser(
        description="Mirror a LAMMPS data file across one or more planes.")
    parser.add_argument("input",  help="Input LAMMPS data file")
    parser.add_argument("output", help="Output LAMMPS data file")
    parser.add_argument(
        "--planes", nargs="+", choices=["xy", "xz", "yz"], required=True,
        help="Mirror plane(s): xy, xz, yz (can specify multiple)")
    parser.add_argument(
        "--atom_style", default="atomic",
        help="LAMMPS atom style (default: atomic)")
    args = parser.parse_args()

    planes = list(set(args.planes))  # deduplicate
    print(f"Input      : {args.input}")
    print(f"Output     : {args.output}")
    print(f"Mirror planes: {', '.join(planes)}")
    print(f"Atom style : {args.atom_style}")

    # Load
    pipeline = import_file(args.input, atom_style=args.atom_style)

    # Build combined mirror + recenter matrix
    mirror_mat = combined_mirror_matrix(planes)
    affine_mat = recenter_matrix(pipeline, mirror_mat)

    print(f"\nAffine transformation matrix (3x4):")
    for row in affine_mat:
        print(f"  {row}")

    # Apply via OVITO modifier
    # operate_on replaces the old transform_particles / transform_velocities /
    # transform_box booleans removed in OVITO 3.x.
    # 'vector_properties' ensures velocities (and forces, etc.) are mirrored.
    pipeline.modifiers.append(
        AffineTransformationModifier(
            transformation=affine_mat,
            operate_on={"particles", "vector_properties", "cell"},
        )
    )

    # Normalize the cell to orthogonal (xlo xhi) format.
    #
    # After AffineTransformationModifier, a mirrored axis produces a negative
    # cell-vector component (e.g. a = (-lx, 0, 0)) with the origin shifted to
    # compensate.  OVITO's LAMMPS exporter cannot express a negative cell-
    # vector length in the xlo/xhi form, so it falls back to triclinic output.
    #
    # Fix: for each axis with a negative diagonal, absorb the vector into the
    # origin (shifting it to 0) and flip the sign, producing a valid positive
    # cell dimension.  Tilt factors are zeroed afterward as a noise guard.
    def _enforce_orthogonal(frame, data):
        cell = data.cell_
        m = cell.matrix.copy()
        for i in range(3):
            if m[i, i] < 0:
                m[i, 3] += m[i, i]   # origin_i += negative_vector_i → 0
                m[i, i]  = -m[i, i]  # flip to positive length
        m[0, 1] = m[0, 2] = m[1, 2] = 0.0  # zero xy, xz, yz tilt factors
        cell.matrix = m

    pipeline.modifiers.append(_enforce_orthogonal)

    # Export
    export_file(
        pipeline, args.output,
        format="lammps/data",
        atom_style=args.atom_style,
    )

    # OVITO 3.9+ defaults to avec/bvec/cvec cell format; rewrite as xlo xhi.
    _fix_cell_header(args.output)

    print(f"\nWritten: {args.output}")


if __name__ == "__main__":
    main()