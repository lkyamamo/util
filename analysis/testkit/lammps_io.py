"""Writers (and a minimal reader) for LAMMPS dump and data files."""

import numpy as np

INT_COLUMNS = {"id", "type", "mol"}


def write_dump(path, frames, box, columns, timesteps=None, bounds="pp pp pp", wrap=False):
    """Write a LAMMPS custom dump.

    frames     list of dicts, column name -> sequence of per-atom values
    box        (lo, hi), each length 3
    columns    column order for the ATOMS header, e.g. ("id", "element", "x", "y", "z")
    timesteps  one per frame (default 0, 1, 2, ...)
    wrap       wrap x y z into the box before writing
    """
    lo = np.asarray(box[0], dtype=float)
    hi = np.asarray(box[1], dtype=float)
    timesteps = list(range(len(frames))) if timesteps is None else list(timesteps)
    with open(path, "w") as f:
        for step, frame in zip(timesteps, frames):
            n_atoms = len(frame[columns[0]])
            data = {c: np.asarray(frame[c]) for c in columns}
            if wrap:
                for d, c in enumerate("xyz"):
                    if c in data:
                        data[c] = (data[c].astype(float) - lo[d]) % (hi[d] - lo[d]) + lo[d]
            f.write(f"ITEM: TIMESTEP\n{step}\nITEM: NUMBER OF ATOMS\n{n_atoms}\n")
            f.write(f"ITEM: BOX BOUNDS {bounds}\n")
            for d in range(3):
                f.write(f"{lo[d]:.10e} {hi[d]:.10e}\n")
            f.write("ITEM: ATOMS " + " ".join(columns) + "\n")
            for i in range(n_atoms):
                f.write(" ".join(_fmt(c, data[c][i]) for c in columns) + "\n")


def _fmt(column, value):
    if column in INT_COLUMNS:
        return str(int(value))
    if isinstance(value, str):
        return value
    return f"{float(value):.6f}"


def read_dump(path):
    """Read a LAMMPS custom dump into a list of (timestep, box, {column: array})."""
    frames = []
    with open(path) as f:
        while f.readline():
            step = int(f.readline())
            f.readline()
            n_atoms = int(f.readline())
            f.readline()
            box = np.array([[float(v) for v in f.readline().split()[:2]] for _ in range(3)])
            columns = f.readline().split()[2:]
            rows = [f.readline().split() for _ in range(n_atoms)]
            table = np.array(rows, dtype=object).reshape(n_atoms, len(columns))
            data = {}
            for j, c in enumerate(columns):
                col = table[:, j]
                if c == "element":
                    data[c] = col.astype(str)
                elif c in INT_COLUMNS:
                    data[c] = col.astype(int)
                else:
                    data[c] = col.astype(float)
            frames.append((step, box, data))
    return frames


def write_data(path, positions, types, masses, box, charges=None):
    """Write a LAMMPS data file, atom_style atomic (or charge when charges are given).

    masses  dict type -> mass (amu)
    box     (lo, hi), each length 3; orthogonal boxes only
    """
    positions = np.asarray(positions, dtype=float)
    types = np.asarray(types, dtype=int)
    lo, hi = np.asarray(box[0], dtype=float), np.asarray(box[1], dtype=float)
    style = "atomic" if charges is None else "charge"
    with open(path, "w") as f:
        f.write("LAMMPS data file written by analysis/testkit\n\n")
        f.write(f"{len(positions)} atoms\n{len(masses)} atom types\n\n")
        for d, axis in enumerate("xyz"):
            f.write(f"{lo[d]:.10f} {hi[d]:.10f} {axis}lo {axis}hi\n")
        f.write("\nMasses\n\n")
        for t in sorted(masses):
            f.write(f"{t} {masses[t]}\n")
        f.write(f"\nAtoms # {style}\n\n")
        for i, (p, t) in enumerate(zip(positions, types), start=1):
            q = "" if charges is None else f" {charges[i - 1]:.6f}"
            f.write(f"{i} {t}{q} {p[0]:.10f} {p[1]:.10f} {p[2]:.10f}\n")
