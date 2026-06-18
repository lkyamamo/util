#!/usr/bin/env python3
"""
MPI dielectric constant calculator for multi-trajectory LAMMPS custom dump files.

One MPI rank per trajectory file. Each rank reads its assigned file sequentially
from line 0 (or a checkpointed offset) with no byte-seeking. Dipole results are
appended directly to ranks/dipole_rank_N.txt; checkpointing is by line count so
restarts resume from where the rank left off.

Usage: launched via srun from 0.submit.slurm — not called directly.
"""

import argparse
import glob
import os
import re
import sys

import numpy as np
from mpi4py import MPI
from scipy.spatial import cKDTree


DEFAULT_TYPE_TO_CHARGE = {
    1: -0.813976,  # oxygen
    2:  0.406988,  # hydrogen
}


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--dump-dir",   required=True)
    p.add_argument("--dump-glob",  required=True)
    p.add_argument("--dump-every", type=int, required=True)
    p.add_argument("--cutoff",     type=float, required=True)
    p.add_argument("--type-o",     type=int, required=True)
    p.add_argument("--type-h",     type=int, required=True)
    p.add_argument("--output-dir", required=True)
    return p.parse_args()


# ---------------------------------------------------------------------------
# File discovery — sorted by integer index in filename
# ---------------------------------------------------------------------------

def discover_files(dump_dir, dump_glob):
    """
    Return trajectory files sorted by the integer index embedded in the filename.
    Expects filenames of the form dielectric.{i}.custom where i is an integer.
    """
    pattern = os.path.join(dump_dir, dump_glob)
    paths = glob.glob(pattern)
    if not paths:
        raise FileNotFoundError(f"No files matching {pattern}")

    def extract_index(path):
        m = re.search(r'\.(\d+)\.', os.path.basename(path))
        if m is None:
            raise ValueError(f"Cannot extract integer index from filename: {path}")
        return int(m.group(1))

    return sorted(paths, key=extract_index)


# ---------------------------------------------------------------------------
# Structure discovery — rank 0 peeks at the first file
# ---------------------------------------------------------------------------

def discover_structure(filepath):
    """
    Read the first frame of filepath to determine:
      - num_atoms
      - lines_per_frame  (num_atoms + 9)
      - column indices for type, x, y, z
      - frames_per_file  (total lines // lines_per_frame)

    Aborts if required columns (type, x, y, z) are missing or if the line
    count is not an exact multiple of lines_per_frame.
    """
    with open(filepath) as f:
        f.readline()                        # ITEM: TIMESTEP
        f.readline()                        # timestep value
        f.readline()                        # ITEM: NUMBER OF ATOMS
        num_atoms = int(f.readline())
        f.readline()                        # ITEM: BOX BOUNDS
        f.readline()                        # xlo xhi
        f.readline()                        # ylo yhi
        f.readline()                        # zlo zhi
        atoms_header = f.readline()         # ITEM: ATOMS col1 col2 ...

    headers = atoms_header.split()[2:]

    if "type" in headers and "x" in headers and "y" in headers and "z" in headers:
        col_type = headers.index("type")
        col_x    = headers.index("x")
        col_y    = headers.index("y")
        col_z    = headers.index("z")
    else:
        missing = [c for c in ("type", "x", "y", "z") if c not in headers]
        raise ValueError(
            f"Missing required columns {missing} in {filepath}. Found: {headers}"
        )

    lines_per_frame = num_atoms + 9

    total_lines = sum(1 for _ in open(filepath))
    if total_lines % lines_per_frame != 0:
        raise ValueError(
            f"{filepath}: {total_lines} lines is not a multiple of "
            f"lines_per_frame={lines_per_frame}"
        )
    frames_per_file = total_lines // lines_per_frame

    return {
        "num_atoms":       num_atoms,
        "lines_per_frame": lines_per_frame,
        "frames_per_file": frames_per_file,
        "col_type":        col_type,
        "col_x":           col_x,
        "col_y":           col_y,
        "col_z":           col_z,
    }


# ---------------------------------------------------------------------------
# Frame reader — streaming, one line at a time
# ---------------------------------------------------------------------------

def read_frame(f, num_atoms, col_type, col_x, col_y, col_z):
    """
    Read one frame from an already-open file handle positioned at the start
    of a frame. Returns (timestep, box_dims, o_pos, h_pos) where o_pos and
    h_pos are numpy arrays of shape (N, 3), or None at EOF.
    """
    line = f.readline()         # ITEM: TIMESTEP  (or empty at EOF)
    if not line:
        return None
    timestep = int(f.readline())

    f.readline()                            # ITEM: NUMBER OF ATOMS
    f.readline()                            # num_atoms (already known)

    f.readline()                            # ITEM: BOX BOUNDS ...
    xlo, xhi = map(float, f.readline().split())
    ylo, yhi = map(float, f.readline().split())
    zlo, zhi = map(float, f.readline().split())
    box_dims = np.array([xhi - xlo, yhi - ylo, zhi - zlo])

    f.readline()                            # ITEM: ATOMS header (columns known)

    o_pos = []
    h_pos = []
    for _ in range(num_atoms):
        parts = f.readline().split()
        atype = int(parts[col_type])
        x = float(parts[col_x])
        y = float(parts[col_y])
        z = float(parts[col_z])
        if atype == 1:
            o_pos.append((x, y, z))
        elif atype == 2:
            h_pos.append((x, y, z))

    return timestep, box_dims, np.array(o_pos), np.array(h_pos)


# ---------------------------------------------------------------------------
# Bond detection + dipole calculation
# ---------------------------------------------------------------------------

def calc_frame_dipole(o_pos, h_pos, cutoff, box_dims):
    box = box_dims

    o_pos_arr = o_pos % box
    h_pos_arr = h_pos % box

    tree = cKDTree(h_pos_arr, boxsize=box)
    dists, idxs = tree.query(o_pos_arr, k=2, distance_upper_bound=cutoff, workers=1)

    valid = ~np.isinf(dists).any(axis=1)
    bond_stats = {
        "bonds_created": int(np.sum(valid)),
        "bonds_missing": int(np.sum(~valid)),
        "type_o_count":  len(o_pos),
        "type_h_count":  len(h_pos),
    }

    assigned_h = set(idxs[valid].ravel().tolist())
    unassigned_h_indices = [i for i in range(len(h_pos)) if i not in assigned_h]

    if bond_stats["bonds_created"] == 0:
        return [0.0, 0.0, 0.0], unassigned_h_indices, bond_stats

    o_v  = o_pos_arr[valid]
    h1_v = h_pos_arr[idxs[valid, 0]]
    h2_v = h_pos_arr[idxs[valid, 1]]

    dr_oh1 = o_v - h1_v;  dr_oh1 -= box * np.round(dr_oh1 / box)
    dr_h2o = h2_v - o_v;  dr_h2o -= box * np.round(dr_h2o / box)

    q_H = DEFAULT_TYPE_TO_CHARGE[2]
    dipole = q_H * np.sum(dr_oh1 - dr_h2o, axis=0)

    return dipole.tolist(), unassigned_h_indices, bond_stats


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    args = parse_args()

    # --- File discovery and validation ---
    files = discover_files(args.dump_dir, args.dump_glob)
    N = len(files)

    if size != N:
        if rank == 0:
            print(
                f"ERROR: ntasks ({size}) must equal number of trajectory files ({N}).\n"
                f"Set --ntasks={N} in your SLURM script.",
                file=sys.stderr,
            )
        comm.Abort(1)

    # --- Structure discovery (rank 0, then broadcast) ---
    if rank == 0:
        try:
            info = discover_structure(files[0])
        except Exception as e:
            print(f"ERROR during structure discovery: {e}", file=sys.stderr)
            comm.Abort(1)
    else:
        info = None
    info = comm.bcast(info, root=0)

    num_atoms       = info["num_atoms"]
    lines_per_frame = info["lines_per_frame"]
    frames_per_file = info["frames_per_file"]
    col_type        = info["col_type"]
    col_x           = info["col_x"]
    col_y           = info["col_y"]
    col_z           = info["col_z"]

    # --- Output paths ---
    rank_output = os.path.join(args.output_dir, "ranks", f"dipole_rank_{rank}.txt")
    warn_log    = os.path.join(args.output_dir, "ranks", f"dipole_rank_{rank}_warn.log")
    os.makedirs(os.path.dirname(rank_output), exist_ok=True)

    # --- Checkpointing ---
    frames_done = 0
    if os.path.isfile(rank_output):
        with open(rank_output) as f:
            frames_done = sum(1 for _ in f)

    if frames_done >= frames_per_file:
        if rank == 0:
            print(f"Rank {rank}: already complete, skipping.")
        comm.Barrier()
        return

    frames_to_skip = frames_done % frames_per_file
    assigned_file  = files[rank]

    # --- Process frames ---
    errors = 0
    with open(assigned_file) as f_in, open(rank_output, "a") as out:
        # Skip already-completed frames by streaming past them.
        lines_to_skip = frames_to_skip * lines_per_frame
        for _ in range(lines_to_skip):
            f_in.readline()

        while True:
            result = read_frame(f_in, num_atoms, col_type, col_x, col_y, col_z)
            if result is None:
                break

            timestep, box_dims, o_pos, h_pos = result

            if len(o_pos) == 0 or len(h_pos) == 0:
                print(
                    f"ERROR: rank {rank} timestep {timestep} has no O or H atoms.",
                    file=sys.stderr,
                )
                errors += 1
                continue

            dipole, unassigned_h, bond_stats = calc_frame_dipole(
                o_pos, h_pos, args.cutoff, box_dims
            )

            out.write(
                f"{timestep}  {dipole[0]:14.6f}  {dipole[1]:14.6f}  {dipole[2]:14.6f}\n"
            )
            out.flush()

            if unassigned_h:
                with open(warn_log, "a") as wf:
                    wf.write(
                        f"WARN [timestep={timestep}] "
                        f"{len(unassigned_h)} unassigned H (indices): "
                        f"[{' '.join(str(i) for i in unassigned_h)}]  "
                        f"| bonds_created={bond_stats['bonds_created']} "
                        f"bonds_missing={bond_stats['bonds_missing']} "
                        f"n_O={bond_stats['type_o_count']} "
                        f"n_H={bond_stats['type_h_count']}\n"
                    )

    if errors:
        print(f"Rank {rank}: {errors} frame(s) failed", file=sys.stderr)
        comm.Abort(1)

    comm.Barrier()


if __name__ == "__main__":
    main()
