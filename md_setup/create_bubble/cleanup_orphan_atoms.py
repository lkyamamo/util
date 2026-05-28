#!/usr/bin/env python3
"""
cleanup_orphan_atoms.py

After LAMMPS has deleted atoms inside a spherical bubble with
atom_style atomic (no molecule IDs), this script removes any
orphaned water atoms left at the bubble boundary to restore
charge neutrality (1 O : 2 H ratio).

Algorithm
---------
1.  Count O and H atoms in the original file and in the LAMMPS
    output file to determine the deletion imbalance.
2a. deficit > 0 (too few H deleted): remove the most-isolated shell H
    atoms (sorted by nearest-O distance, descending).
2b. deficit < 0 (too few O deleted): remove the most-isolated shell O
    atoms (sorted by nearest-H distance, descending); if the surplus H
    count is odd, also remove the single most-isolated shell H after
    the O KDTree is rebuilt without the marked O atoms.
3.  Renumber atom IDs sequentially and overwrite the output file.

Usage:
    python3 cleanup_orphan_atoms.py \\
        --original system.data \\
        --data     system_bubble.data \\
        --cx 400.0 --cy 400.0 --cz 881.0 \\
        --radius 150.0 \\
        [--shell-thickness 2.0] \\
        [--o-type 2] \\
        [--h-type 3] \\
        [--si-type 1]
"""

from __future__ import annotations

import argparse
import tempfile
import warnings
from pathlib import Path

import numpy as np
from scipy.spatial import KDTree

from lammps_utils import count_atoms_by_type

# ── constants ─────────────────────────────────────────────────────────────────
_MAX_SHELL_EXPANSIONS = 5
_SHELL_EXPANSION_STEP = 0.1  # Å


# ── data file I/O ─────────────────────────────────────────────────────────────

class _Atom:
    __slots__ = ("atom_id", "atom_type", "x", "y", "z")

    def __init__(self, atom_id: int, atom_type: int, x: float, y: float, z: float):
        self.atom_id   = atom_id
        self.atom_type = atom_type
        self.x         = x
        self.y         = y
        self.z         = z


def _parse_atoms_section(data_file: str) -> list[_Atom]:
    """
    Parse the Atoms section of an atom_style atomic LAMMPS data file.
    Format per line: atom-ID  type  x  y  z
    Returns a list of _Atom objects in file order.
    """
    atoms: list[_Atom] = []
    in_atoms = False
    with open(data_file) as fh:
        for line in fh:
            stripped = line.split("#")[0].strip()
            if not stripped:
                if in_atoms:
                    in_atoms = False
                continue
            if stripped == "Atoms":
                in_atoms = True
                continue
            if in_atoms and stripped[0].isalpha():
                break
            if in_atoms:
                parts = stripped.split()
                if len(parts) >= 5:
                    atoms.append(_Atom(
                        atom_id   = int(parts[0]),
                        atom_type = int(parts[1]),
                        x         = float(parts[2]),
                        y         = float(parts[3]),
                        z         = float(parts[4]),
                    ))
    return atoms


def _rewrite_data_file(
    data_file: str,
    delete_ids: set[int],
) -> None:
    """
    Rewrite data_file in-place:
      - Skip atoms whose atom_id is in delete_ids (Atoms and Velocities sections).
      - Update the atom count in the header.
      - Renumber remaining atom IDs sequentially (1, 2, 3, …) in both sections.
      - Pass all other sections through unchanged.

    LAMMPS always writes Atoms before Velocities, so the old→new ID mapping
    built during the Atoms pass is fully available when Velocities is reached.
    """
    path = Path(data_file)
    original_text = path.read_text()
    lines = original_text.splitlines(keepends=True)

    n_delete = len(delete_ids)
    new_lines: list[str] = []
    in_atoms = False
    in_velocities = False
    new_id = 0
    id_map: dict[int, int] = {}  # old atom-ID → new sequential atom-ID

    for line in lines:
        stripped = line.split("#")[0].strip()

        # Update atom count in header (e.g. "1234567 atoms")
        if not in_atoms and not in_velocities and stripped.endswith(" atoms") and stripped.split()[0].isdigit():
            old_count = int(stripped.split()[0])
            new_count = old_count - n_delete
            new_lines.append(line.replace(str(old_count), str(new_count), 1))
            continue

        if stripped == "Atoms":
            in_atoms = True
            new_lines.append(line)
            continue

        if stripped == "Velocities":
            in_velocities = True
            new_lines.append(line)
            continue

        if in_atoms:
            if not stripped:
                new_lines.append(line)
                continue
            if stripped[0].isalpha():
                in_atoms = False
                new_lines.append(line)
                continue
            parts = stripped.split()
            if len(parts) >= 5:
                atom_id = int(parts[0])
                if atom_id in delete_ids:
                    continue
                new_id += 1
                id_map[atom_id] = new_id
                new_line = f"{new_id:>8}  " + "  ".join(parts[1:]) + "\n"
                new_lines.append(new_line)
                continue

        if in_velocities:
            if not stripped:
                new_lines.append(line)
                continue
            if stripped[0].isalpha():
                in_velocities = False
                new_lines.append(line)
                continue
            parts = stripped.split()
            if len(parts) >= 4:
                atom_id = int(parts[0])
                if atom_id in delete_ids:
                    continue
                mapped_id = id_map.get(atom_id, atom_id)
                new_line = f"{mapped_id:>8}  " + "  ".join(parts[1:]) + "\n"
                new_lines.append(new_line)
                continue

        new_lines.append(line)

    # Write atomically via a temp file in the same directory
    tmp = tempfile.NamedTemporaryFile(
        mode="w", dir=path.parent, delete=False, suffix=".tmp"
    )
    try:
        tmp.writelines(new_lines)
        tmp.flush()
        tmp.close()
        Path(tmp.name).replace(path)
    except Exception:
        Path(tmp.name).unlink(missing_ok=True)
        raise


# ── shell helpers ─────────────────────────────────────────────────────────────

def _distances_from_center(
    atoms: list[_Atom],
    cx: float,
    cy: float,
    cz: float,
) -> np.ndarray:
    if not atoms:
        return np.array([], dtype=float)
    coords = np.array([[a.x, a.y, a.z] for a in atoms])
    center = np.array([cx, cy, cz])
    return np.linalg.norm(coords - center, axis=1)


def _shell_atoms(
    atoms: list[_Atom],
    atom_type: int,
    cx: float,
    cy: float,
    cz: float,
    radius: float,
    shell_thickness: float,
) -> list[_Atom]:
    """Return atoms of atom_type whose distance from center is in (radius, radius+shell_thickness]."""
    result = []
    for a in atoms:
        if a.atom_type != atom_type:
            continue
        d = ((a.x - cx)**2 + (a.y - cy)**2 + (a.z - cz)**2) ** 0.5
        if radius < d <= radius + shell_thickness:
            result.append(a)
    return result


def _nearest_partner_distances(
    query_atoms: list[_Atom],
    partner_atoms: list[_Atom],
) -> np.ndarray:
    """
    For each atom in query_atoms, return the distance to its nearest
    atom in partner_atoms.  Returns inf if partner_atoms is empty.
    """
    if not partner_atoms or not query_atoms:
        return np.full(len(query_atoms), np.inf)
    partner_coords = np.array([[a.x, a.y, a.z] for a in partner_atoms])
    query_coords   = np.array([[a.x, a.y, a.z] for a in query_atoms])
    tree = KDTree(partner_coords)
    dists, _ = tree.query(query_coords)
    return dists


def _get_shell_with_expansion(
    atoms: list[_Atom],
    atom_type: int,
    cx: float,
    cy: float,
    cz: float,
    radius: float,
    shell_thickness: float,
    n_needed: int,
    label: str,
) -> tuple[list[_Atom], float]:
    """
    Return shell atoms of atom_type with enough members (>= n_needed).
    Expands shell_thickness by _SHELL_EXPANSION_STEP up to _MAX_SHELL_EXPANSIONS
    times before raising RuntimeError.
    """
    thickness = shell_thickness
    for attempt in range(_MAX_SHELL_EXPANSIONS + 1):
        shell = _shell_atoms(atoms, atom_type, cx, cy, cz, radius, thickness)
        if len(shell) >= n_needed:
            if attempt > 0:
                warnings.warn(
                    f"Shell expanded to {thickness:.1f} Å to find {n_needed} {label} atoms "
                    f"(found {len(shell)})."
                )
            return shell, thickness
        if attempt < _MAX_SHELL_EXPANSIONS:
            warnings.warn(
                f"Shell thickness {thickness:.1f} Å contains only {len(shell)} {label} atoms "
                f"(need {n_needed}); expanding by {_SHELL_EXPANSION_STEP} Å "
                f"(attempt {attempt + 1}/{_MAX_SHELL_EXPANSIONS})."
            )
            thickness += _SHELL_EXPANSION_STEP

    raise RuntimeError(
        f"Cannot find {n_needed} {label} atoms within shell thickness "
        f"{thickness:.1f} Å after {_MAX_SHELL_EXPANSIONS} expansions. "
        f"Only {len(shell)} found."
    )


# ── logging ───────────────────────────────────────────────────────────────────

def _log_shell_candidates(
    shell_atoms: list[_Atom],
    dists: np.ndarray,
    delete_ids: set[int],
    label: str,
) -> None:
    """Print a table of all shell candidates with their nearest-partner distance."""
    print(f"\n  {'atom-ID':>10}  {'nearest-' + label + ' dist (Å)':>24}  status")
    print(f"  {'-'*10}  {'-'*24}  {'-'*8}")
    order = np.argsort(dists)[::-1]
    for i in order:
        a = shell_atoms[i]
        tag = "[DELETE]" if a.atom_id in delete_ids else "[KEEP]  "
        print(f"  {a.atom_id:>10}  {dists[i]:>24.6f}  {tag}")


# ── main cleanup logic ────────────────────────────────────────────────────────

def cleanup_orphan_atoms(
    original_file: str,
    data_file: str,
    cx: float,
    cy: float,
    cz: float,
    radius: float,
    shell_thickness: float,
    o_type: int,
    h_type: int,
    si_type: int,
) -> None:
    # ── Step 1: compute imbalance ─────────────────────────────────────────────
    counts_before = count_atoms_by_type(original_file)
    counts_after  = count_atoms_by_type(data_file)

    n_si_before = counts_before.get(si_type, 0)
    n_si_after  = counts_after.get(si_type, 0)
    if n_si_after != n_si_before:
        warnings.warn(
            f"{n_si_before - n_si_after} Si atom(s) (type {si_type}) were deleted by LAMMPS. "
            "The bubble may intersect the silica slab. "
            "The O/H deficit calculation may be unreliable."
        )

    n_O_deleted = counts_before.get(o_type, 0) - counts_after.get(o_type, 0)
    n_H_deleted = counts_before.get(h_type, 0) - counts_after.get(h_type, 0)

    deficit = 2 * n_O_deleted - n_H_deleted

    print(f"\nOrphan atom cleanup")
    print(f"  O deleted : {n_O_deleted}")
    print(f"  H deleted : {n_H_deleted}")
    print(f"  deficit   : {deficit:+d}  "
          f"({'orphan H' if deficit > 0 else 'orphan O' if deficit < 0 else 'balanced'})")

    if deficit == 0:
        print("  No orphan atoms — renumbering IDs only.")
        _rewrite_data_file(data_file, delete_ids=set())
        return

    # Parse the full Atoms section of the output file for coordinate work
    atoms = _parse_atoms_section(data_file)
    delete_ids: set[int] = set()

    if deficit > 0:
        # ── Phase 2a: orphan H ────────────────────────────────────────────────
        n_h_extra = deficit
        print(f"\n  Orphan H cleanup: need to delete {n_h_extra} H atom(s) from shell.")

        shell_H, final_thickness = _get_shell_with_expansion(
            atoms, h_type, cx, cy, cz, radius, shell_thickness, n_h_extra, "H"
        )
        print(f"  Shell thickness used: {final_thickness:.1f} Å  "
              f"({len(shell_H)} H candidates)")

        all_O = [a for a in atoms if a.atom_type == o_type]
        dists_H = _nearest_partner_distances(shell_H, all_O)

        order = np.argsort(dists_H)[::-1]
        for i in order[:n_h_extra]:
            delete_ids.add(shell_H[i].atom_id)

        print(f"\n  Shell H candidates (sorted by nearest-O distance):")
        _log_shell_candidates(shell_H, dists_H, delete_ids, "O")

    else:
        # ── Phase 2b: orphan O ────────────────────────────────────────────────
        surplus       = -deficit
        n_o_extra     = surplus // 2
        odd_remainder = surplus % 2

        print(f"\n  Orphan O cleanup: need to delete {n_o_extra} O atom(s) from shell.")
        if odd_remainder:
            print(f"  Odd surplus — will also delete 1 additional shell H after O removal.")

        shell_O, final_thickness_O = _get_shell_with_expansion(
            atoms, o_type, cx, cy, cz, radius, shell_thickness, n_o_extra, "O"
        )
        print(f"  Shell thickness used: {final_thickness_O:.1f} Å  "
              f"({len(shell_O)} O candidates)")

        all_H = [a for a in atoms if a.atom_type == h_type]
        dists_O = _nearest_partner_distances(shell_O, all_H)

        order_O = np.argsort(dists_O)[::-1]
        for i in order_O[:n_o_extra]:
            delete_ids.add(shell_O[i].atom_id)

        print(f"\n  Shell O candidates (sorted by nearest-H distance):")
        _log_shell_candidates(shell_O, dists_O, delete_ids, "H")

        if odd_remainder:
            # Rebuild O set excluding atoms already marked for deletion
            remaining_O = [a for a in atoms
                           if a.atom_type == o_type and a.atom_id not in delete_ids]

            shell_H, final_thickness_H = _get_shell_with_expansion(
                atoms, h_type, cx, cy, cz, radius, shell_thickness, 1, "H"
            )
            print(f"\n  Odd-remainder H cleanup: shell thickness {final_thickness_H:.1f} Å  "
                  f"({len(shell_H)} H candidates, O KDTree rebuilt without marked O atoms)")

            dists_H = _nearest_partner_distances(shell_H, remaining_O)
            most_isolated_idx = int(np.argmax(dists_H))
            delete_ids.add(shell_H[most_isolated_idx].atom_id)

            print(f"\n  Shell H candidates (sorted by nearest-O distance, post-O deletion):")
            _log_shell_candidates(shell_H, dists_H, delete_ids, "O")

    print(f"\n  Total atoms to delete: {len(delete_ids)}")
    print(f"  Rewriting {data_file} ...")
    _rewrite_data_file(data_file, delete_ids)
    print("  Done.")


# ── entry point ───────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Remove orphan H or O atoms left at a bubble boundary."
    )
    parser.add_argument("--original",        required=True,
                        help="Original input data file (before LAMMPS deletion)")
    parser.add_argument("--data",            required=True,
                        help="LAMMPS output data file to clean up (modified in-place)")
    parser.add_argument("--cx",              required=True, type=float)
    parser.add_argument("--cy",              required=True, type=float)
    parser.add_argument("--cz",              required=True, type=float)
    parser.add_argument("--radius",          required=True, type=float,
                        help="Bubble radius (Å)")
    parser.add_argument("--shell-thickness", default=2.0,   type=float,
                        help="Initial shell thickness outside bubble to search (Å, default 2.0)")
    parser.add_argument("--o-type",          default=2,     type=int,
                        help="Atom type for oxygen (default 2)")
    parser.add_argument("--h-type",          default=3,     type=int,
                        help="Atom type for hydrogen (default 3)")
    parser.add_argument("--si-type",         default=1,     type=int,
                        help="Atom type for silicon (default 1)")
    args = parser.parse_args()

    cleanup_orphan_atoms(
        original_file   = args.original,
        data_file       = args.data,
        cx              = args.cx,
        cy              = args.cy,
        cz              = args.cz,
        radius          = args.radius,
        shell_thickness = args.shell_thickness,
        o_type          = args.o_type,
        h_type          = args.h_type,
        si_type         = args.si_type,
    )


if __name__ == "__main__":
    main()
