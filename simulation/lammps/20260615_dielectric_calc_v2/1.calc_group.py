#!/usr/bin/env python3
"""
Group dipole moment calculator for LAMMPS custom dump files.

Processes GROUP_SIZE consecutive frames for a single SLURM array task.
Frames are handled one at a time — only one frame's atom data is in memory
at any point. Each frame's result line is appended directly to the group
output file immediately after computation.

Checkpointing: on restart, the number of already-written lines is counted
and that many frames are skipped from the start of the group's range.
Consecutive timestep order is guaranteed by the loop structure, so the
output file is always sorted by timestep without any sorting step.

Usage:
    python 1.calc_group.py
        --group-id       <int>
        --group-size     <int>
        --num-frames     <int>
        --start-timestep <int>
        --dump-every     <int>
        --dump-dir       <path>
        --dump-glob      <glob-pattern>   (e.g. "dielectric.*.custom")
        --cutoff         <float>
        --type-o         <int>
        --type-h         <int>
        --output         <path>
"""

import argparse
import os
import sys
import numpy as np
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
    p.add_argument("--group-id",       type=int,   required=True)
    p.add_argument("--group-size",     type=int,   required=True)
    p.add_argument("--num-frames",     type=int,   required=True)
    p.add_argument("--start-timestep", type=int,   required=True)
    p.add_argument("--dump-every",     type=int,   required=True)
    p.add_argument("--dump-dir",                   required=True)
    p.add_argument("--dump-glob",                  required=True)
    p.add_argument("--cutoff",         type=float, required=True)
    p.add_argument("--type-o",         type=int,   required=True)
    p.add_argument("--type-h",         type=int,   required=True)
    p.add_argument("--output",                     required=True)
    return p.parse_args()


# ---------------------------------------------------------------------------
# Dump file parser — streaming, one line at a time
# ---------------------------------------------------------------------------

# LAMMPS custom dump frame structure (fixed 9-line header):
#   line 0 : ITEM: TIMESTEP
#   line 1 : <timestep>
#   line 2 : ITEM: NUMBER OF ATOMS
#   line 3 : <num_atoms>
#   line 4 : ITEM: BOX BOUNDS ...
#   line 5 : xlo xhi
#   line 6 : ylo yhi
#   line 7 : zlo zhi
#   line 8 : ITEM: ATOMS <col1> <col2> ...
#   lines 9 .. 9+num_atoms-1 : atom data
HEADER_LINES = 9


def _peek_num_atoms(f):
    """Read the first 4 lines to extract num_atoms, then seek back to 0."""
    f.readline()                        # ITEM: TIMESTEP
    f.readline()                        # timestep value
    f.readline()                        # ITEM: NUMBER OF ATOMS
    num_atoms = int(f.readline())       # num_atoms
    f.seek(0)
    return num_atoms


def read_lammps_dump(filepath):
    """
    Stream a single-frame LAMMPS custom dump file one line at a time.

    Peeks at lines 0-3 to get num_atoms, then reads exactly
    HEADER_LINES + num_atoms lines via iteration — no readlines().

    Returns:
        timestep  : int
        box_dims  : (la, lb, lc)
        atoms     : dict { atom_id: (atom_type, x, y, z) }
    """
    with open(filepath) as f:
        num_atoms = _peek_num_atoms(f)

        # --- Header (lines 0-8) ---
        f.readline()                                    # ITEM: TIMESTEP
        timestep  = int(f.readline())                   # timestep

        f.readline()                                    # ITEM: NUMBER OF ATOMS
        f.readline()                                    # num_atoms (already known)

        f.readline()                                    # ITEM: BOX BOUNDS ...
        xlo, xhi  = map(float, f.readline().split())
        ylo, yhi  = map(float, f.readline().split())
        zlo, zhi  = map(float, f.readline().split())
        box_dims  = (xhi - xlo, yhi - ylo, zhi - zlo)

        atoms_header = f.readline()                     # ITEM: ATOMS <cols>
        headers  = atoms_header.split()[2:]
        try:
            col_id   = headers.index("id")
            col_type = headers.index("type")
            col_x    = headers.index("x")
            col_y    = headers.index("y")
            col_z    = headers.index("z")
        except ValueError:
            raise ValueError(
                f"'{filepath}' missing required columns. Found: {headers}"
            )
        # --- Atom data (exactly num_atoms lines) ---
        atoms = {}
        for _ in range(num_atoms):
            parts     = f.readline().split()
            atom_id   = int(parts[col_id])
            atom_type = int(parts[col_type])
            x = float(parts[col_x])
            y = float(parts[col_y])
            z = float(parts[col_z])
            atoms[atom_id] = (atom_type, x, y, z)

    return timestep, box_dims, atoms



# ---------------------------------------------------------------------------
# Bond detection + dipole calculation (single frame, fully vectorized)
# ---------------------------------------------------------------------------

def calc_frame_dipole(atoms, type_o, type_h, cutoff, box_dims):
    """
    Find O-H2 molecules via KD-tree and compute the total dipole moment,
    all in vectorized numpy operations.

    Returns:
        total        : [Mx, My, Mz]
        unassigned_h : list of H atom IDs not bonded to any O
        bond_stats   : dict with diagnostic counts
    """
    la, lb, lc = box_dims
    box = np.array([la, lb, lc])

    o_ids, o_pos = [], []
    h_ids, h_pos = [], []
    for atom_id, (atom_type, x, y, z) in atoms.items():
        if atom_type == type_o:
            o_ids.append(atom_id)
            o_pos.append([x, y, z])
        elif atom_type == type_h:
            h_ids.append(atom_id)
            h_pos.append([x, y, z])

    bond_stats = {
        "bonds_created": 0, "bonds_missing": len(o_ids),
        "type_o_count": len(o_ids), "type_h_count": len(h_ids),
    }

    if not h_pos:
        return [0.0, 0.0, 0.0], list(h_ids), bond_stats

    # Wrap Cartesian coords into [0, L) as required by cKDTree PBC mode.
    o_pos_arr = np.array(o_pos) % box
    h_pos_arr = np.array(h_pos) % box
    h_ids_arr = np.array(h_ids, dtype=np.int32)

    # Find 2 nearest H per O within cutoff. dists == inf → no neighbor found.
    tree = cKDTree(h_pos_arr, boxsize=box)
    dists, idxs = tree.query(o_pos_arr, k=2, distance_upper_bound=cutoff, workers=1)

    valid = ~np.isinf(dists).any(axis=1)          # (N_O,) bool mask
    bonds_missing  = int(np.sum(~valid))
    bonds_created  = int(np.sum(valid))

    assigned_h = set(idxs[valid].ravel().tolist())
    unassigned_h = [int(h_ids_arr[i])
                    for i in range(len(h_ids_arr)) if i not in assigned_h]

    bond_stats = {
        "bonds_created": bonds_created,
        "bonds_missing": bonds_missing,
        "type_o_count":  len(o_ids),
        "type_h_count":  len(h_ids),
    }

    if bonds_created == 0:
        return [0.0, 0.0, 0.0], unassigned_h, bond_stats

    # Slice positions for bonded molecules only.
    o_v  = o_pos_arr[valid]                        # (N_bond, 3)
    h1_v = h_pos_arr[idxs[valid, 0]]               # (N_bond, 3)
    h2_v = h_pos_arr[idxs[valid, 1]]               # (N_bond, 3)

    # Minimum-image bond vectors (wrapped coords → displacements in (-L/2, L/2)).
    dr_oh1 = o_v - h1_v
    dr_oh1 -= box * np.round(dr_oh1 / box)         # O → H1 (original convention)
    dr_h2o = h2_v - o_v
    dr_h2o -= box * np.round(dr_h2o / box)         # H2 → O

    # All H are the same type so q_H is a scalar; it factors out of the sum.
    q_H = DEFAULT_TYPE_TO_CHARGE[type_h]
    dipole = q_H * np.sum(dr_oh1 - dr_h2o, axis=0)

    return dipole.tolist(), unassigned_h, bond_stats


# ---------------------------------------------------------------------------
# Per-frame processing — returns result, releases atom data
# ---------------------------------------------------------------------------

def process_frame(dump_path, cutoff, type_o, type_h):
    """
    Read one dump file, compute dipole, release all atom data.
    Returns (result_line, unassigned_h, bond_stats), or None if dump missing.
    """
    if not os.path.isfile(dump_path):
        print(f"ERROR: dump not found: {dump_path}", file=sys.stderr)
        return None

    timestep, box_dims, atoms = read_lammps_dump(dump_path)
    total, unassigned_h, bond_stats = calc_frame_dipole(
        atoms, type_o, type_h, cutoff, box_dims
    )
    del atoms

    result_line = f"{timestep}  {total[0]:14.6f}  {total[1]:14.6f}  {total[2]:14.6f}\n"
    return result_line, unassigned_h, bond_stats


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    args = parse_args()

    group_start  = args.group_id * args.group_size
    group_end    = min(group_start + args.group_size, args.num_frames)
    group_frames = group_end - group_start

    # Checkpoint: count complete lines already written to resume after failure.
    # Consecutive order is guaranteed by the loop below, so line N always
    # corresponds to the Nth frame of this group — no sorting needed.
    frames_done = 0
    if os.path.isfile(args.output):
        with open(args.output) as f:
            frames_done = sum(1 for _ in f)

    warn_log = args.output.replace(".txt", "_warn.log")

    errors = 0
    with open(args.output, "a") as out:
        for local_idx in range(frames_done, group_frames):
            frame_idx = group_start + local_idx
            timestep  = args.start_timestep + frame_idx * args.dump_every
            dump_name = args.dump_glob.replace("*", str(timestep))
            dump_path = os.path.join(args.dump_dir, dump_name)

            result = process_frame(dump_path, args.cutoff, args.type_o, args.type_h)
            if result is None:
                errors += 1
                continue

            result_line, unassigned_h, bond_stats = result
            out.write(result_line)
            out.flush()

            if unassigned_h:
                ids_str = " ".join(str(a) for a in unassigned_h)
                with open(warn_log, "a") as wf:
                    wf.write(
                        f"WARN [timestep={timestep}] "
                        f"{len(unassigned_h)} unassigned H: [{ids_str}]  "
                        f"| bonds_created={bond_stats['bonds_created']} "
                        f"bonds_missing={bond_stats['bonds_missing']} "
                        f"n_O={bond_stats['type_o_count']} "
                        f"n_H={bond_stats['type_h_count']}\n"
                    )

    if errors:
        print(f"ERROR: {errors} frame(s) failed in group {args.group_id}", file=sys.stderr)
        sys.exit(1)

    print(f"Group {args.group_id}: wrote {group_frames} frames → {args.output}")


if __name__ == "__main__":
    main()
