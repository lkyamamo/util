#!/usr/bin/env python3
"""
MPI dielectric constant calculator for multi-trajectory LAMMPS custom dump files.

Ranks may be fewer than trajectory files: each rank gets a contiguous block
of files (remainder files go one-each to the lowest-numbered ranks) and
processes them in order. Each file is read sequentially from line 0 (or a
checkpointed offset) with no byte-seeking. Dipole results from all of a
rank's files are appended, in order, to a single ranks/dipole_rank_N.txt;
checkpointing is by line count in that file so restarts resume from wherever
the rank left off, in whichever of its assigned files that falls in.

With --delete-processed-dumps, a trajectory file is removed once every frame
in it has become a line in the rank's output and those lines are fsynced —
each line carries everything the dielectric constant needs from its frame, so
the dump has no further use. This is irreversible: recovering a deleted dump
means re-running LAMMPS. It also means the file list has to come from
ranks/dumps_manifest.txt rather than a live glob, since assign_files slices by
position and a shrinking directory would silently repartition the work.

Usage: launched via srun from 0.submit.slurm — not called directly.
"""

import argparse
import glob
import json
import os
import re
import sys

import numpy as np
from mpi4py import MPI
from scipy.spatial import cKDTree


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
    p.add_argument("--charge-h",   type=float, required=True)
    p.add_argument("--output-dir", required=True)
    p.add_argument(
        "--delete-processed-dumps", action="store_true",
        help="Delete each trajectory file once every frame in it has been "
             "reduced to dipole lines and those lines are fsynced. IRREVERSIBLE: "
             "the dumps cannot be recomputed without re-running LAMMPS. Off by "
             "default; the sweep pipeline turns it on.",
    )
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


def resolve_files(dump_dir, dump_glob, ranks_dir):
    """The trajectory file list for this run, pinned across restarts.

    A live glob cannot be used once --delete-processed-dumps is in play.
    assign_files slices the list BY POSITION, so a resumed run that globbed a
    directory with the processed files removed would hand every rank a
    different block than the one its checkpointed frame count describes, and
    the dipole output would be silently wrong rather than merely incomplete.

    So the list is globbed once and written to a manifest, and every later run
    reads the manifest instead. It records what the run was launched against,
    not what happens to still be on disk.
    """
    manifest = os.path.join(ranks_dir, "dumps_manifest.txt")

    if os.path.isfile(manifest):
        with open(manifest) as f:
            files = [line.strip() for line in f if line.strip()]
        if not files:
            raise ValueError(f"{manifest} is empty; remove {ranks_dir} to start over")
        return files

    files = discover_files(dump_dir, dump_glob)
    tmp = f"{manifest}.partial"
    with open(tmp, "w") as f:
        f.write("".join(f"{p}\n" for p in files))
    os.replace(tmp, manifest)
    return files


def assign_files(files, rank, size):
    """
    Return this rank's contiguous block of files. N files over `size` ranks
    distributes evenly with any remainder (N % size) given one extra file
    each to the lowest-numbered ranks. When size == N (the historical
    default), every rank gets exactly one file.
    """
    n = len(files)
    base, rem = divmod(n, size)
    if rank < rem:
        start = rank * (base + 1)
        count = base + 1
    else:
        start = rem * (base + 1) + (rank - rem) * base
        count = base
    return files[start:start + count]


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

def read_frame(f, num_atoms, col_type, col_x, col_y, col_z, type_o, type_h):
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
        if atype == type_o:
            o_pos.append((x, y, z))
        elif atype == type_h:
            h_pos.append((x, y, z))

    return timestep, box_dims, np.array(o_pos), np.array(h_pos)


# ---------------------------------------------------------------------------
# Bond detection + dipole calculation
# ---------------------------------------------------------------------------

def calc_frame_dipole(o_pos, h_pos, cutoff, box_dims, q_h):
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

    dipole = q_h * np.sum(dr_oh1 - dr_h2o, axis=0)

    return dipole.tolist(), unassigned_h_indices, bond_stats


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    args = parse_args()

    # --- Discovery and restart safety (rank 0 only, then broadcast) ---
    #
    # Rank 0 owns the three files in ranks/ that have to agree with each other
    # — the rank-count marker, the trajectory manifest and the structure cache
    # — and hands the results to everyone else. Having every rank glob the same
    # directory and re-read the same first file bought nothing and, once the
    # manifest exists, would be a write race.
    ranks_dir = os.path.join(args.output_dir, "ranks")
    structure_cache = os.path.join(ranks_dir, "structure.json")

    payload = None
    if rank == 0:
        os.makedirs(ranks_dir, exist_ok=True)

        # Rank->file assignment is a function of rank count (size). If a partial
        # run resumes with a different size, ranks get reassigned different
        # files than whatever their checkpointed frame counts actually cover,
        # silently corrupting output. Record the rank count a run started with
        # and refuse to resume under a different one.
        nranks_marker = os.path.join(ranks_dir, "nranks_used.txt")
        if os.path.isfile(nranks_marker):
            with open(nranks_marker) as f:
                prev_size = int(f.read().strip())
            if prev_size != size:
                print(
                    f"ERROR: this run is using {size} ranks, but a previous run in "
                    f"{args.output_dir} used {prev_size} ranks. Rank->file assignment "
                    f"depends on rank count, so resuming with a different count would "
                    f"corrupt checkpoints. Re-run with {prev_size} ranks (matching "
                    f"--dielectric-ntasks), or remove {ranks_dir} to start over.",
                    file=sys.stderr,
                )
                comm.Abort(1)
        else:
            with open(nranks_marker, "w") as f:
                f.write(f"{size}\n")

        try:
            files = resolve_files(args.dump_dir, args.dump_glob, ranks_dir)
        except (OSError, ValueError) as e:
            print(f"ERROR resolving trajectory files: {e}", file=sys.stderr)
            comm.Abort(1)

        if size > len(files):
            print(
                f"ERROR: ntasks ({size}) must not exceed number of trajectory files "
                f"({len(files)}).\nSet --ntasks to at most {len(files)} in your "
                f"SLURM script.",
                file=sys.stderr,
            )
            comm.Abort(1)

        # Cached rather than re-derived, because with --delete-processed-dumps
        # a resumed run may find files[0] already gone. The layout is a property
        # of the trajectory as a whole, so reading it once is also just cheaper.
        if os.path.isfile(structure_cache):
            with open(structure_cache) as f:
                info = json.load(f)
        else:
            first_present = next((p for p in files if os.path.isfile(p)), None)
            if first_present is None:
                print(
                    f"ERROR: none of the {len(files)} files in the manifest exist, and "
                    f"there is no {structure_cache} to fall back on. If the dumps were "
                    f"deleted after processing, re-run LAMMPS to regenerate them.",
                    file=sys.stderr,
                )
                comm.Abort(1)
            try:
                info = discover_structure(first_present)
            except Exception as e:
                print(f"ERROR during structure discovery: {e}", file=sys.stderr)
                comm.Abort(1)
            tmp = f"{structure_cache}.partial"
            with open(tmp, "w") as f:
                json.dump(info, f)
            os.replace(tmp, structure_cache)

        payload = {"files": files, "info": info}

    payload = comm.bcast(payload, root=0)
    files = payload["files"]
    info = payload["info"]
    N = len(files)

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

    # --- This rank's assigned files (contiguous block; see assign_files) ---
    assigned_files = assign_files(files, rank, size)

    # --- Checkpointing ---
    frames_done = 0
    if os.path.isfile(rank_output):
        with open(rank_output) as f:
            frames_done = sum(1 for _ in f)

    # Frame 0 of each file is always skipped because consecutive trajectory
    # files share a boundary frame (the last frame of file N equals the
    # first frame of file N+1). Skipping the first frame of each file
    # removes the duplicate.
    processable_per_file = frames_per_file - 1
    total_processable = processable_per_file * len(assigned_files)

    done_flag = os.path.join(args.output_dir, "ranks", f"dipole_rank_{rank}.done")

    if frames_done >= total_processable:
        if rank == 0:
            print(f"Rank {rank}: already complete, skipping.")
        open(done_flag, "w").close()
        comm.Barrier()
        return

    # frames_done is a flat count across this rank's whole file list — map it
    # onto (which assigned file to resume in, how many of that file's frames
    # are already written). Files before resume_file_idx are fully done, so
    # only the first file processed below needs the extra partial-frame skip.
    resume_file_idx = frames_done // processable_per_file
    resume_frame_offset = frames_done % processable_per_file
    files_to_process = assigned_files[resume_file_idx:]

    # Files before resume_file_idx are fully reduced to dipole lines already.
    # If they survived an earlier run that had deletion off, drop them now —
    # the flag means "processed dumps do not stay on disk", and the usual
    # reason to turn it on mid-sweep is that the disk is already filling up.
    if args.delete_processed_dumps:
        for spent in assigned_files[:resume_file_idx]:
            try:
                os.remove(spent)
            except FileNotFoundError:
                pass
            except OSError as e:
                print(f"Rank {rank}: could not delete {spent}: {e}", file=sys.stderr)

    # --- Process frames ---
    errors = 0
    with open(rank_output, "a") as out:
        for i, assigned_file in enumerate(files_to_process):
            frame_offset = resume_frame_offset if i == 0 else 0

            with open(assigned_file) as f_in:
                # Skip frame 0 unconditionally, then skip any checkpointed frames.
                lines_to_skip = (1 + frame_offset) * lines_per_frame
                for _ in range(lines_to_skip):
                    f_in.readline()

                while True:
                    result = read_frame(
                        f_in, num_atoms, col_type, col_x, col_y, col_z,
                        args.type_o, args.type_h,
                    )
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
                        o_pos, h_pos, args.cutoff, box_dims, args.charge_h
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

            # This file is now fully reduced: every frame in it has become a
            # line in rank_output, and each line carries everything the
            # dielectric constant needs from the frame it came from.
            #
            # fsync before unlinking, not just flush. flush() only pushes the
            # lines out of Python's buffer into the page cache; a node that
            # dies before writeback would take the dipole lines with it while
            # the dump was already gone, and nothing could reconstruct either.
            # Once per file, so the cost is irrelevant next to reading it.
            if args.delete_processed_dumps:
                out.flush()
                os.fsync(out.fileno())
                try:
                    os.remove(assigned_file)
                except OSError as e:
                    print(
                        f"Rank {rank}: could not delete {assigned_file}: {e}",
                        file=sys.stderr,
                    )

    if errors:
        print(f"Rank {rank}: {errors} frame(s) failed", file=sys.stderr)
        comm.Abort(1)

    open(done_flag, "w").close()
    comm.Barrier()


if __name__ == "__main__":
    main()
