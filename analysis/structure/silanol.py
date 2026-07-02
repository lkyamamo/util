#!/usr/bin/env python3
"""Count surface silanols and export analysis artifacts with OVITO.

Usage:
    python silanol.py --local-path PATH [options]

Primary options:
    --local-path PATH            Local trajectory file, or a local directory
                                  containing one file per frame (files are
                                  sorted and loaded as a sequence).
    --frames START END           Optional frame window [START, END).
    --surface-axis {x,y,z}       Axis along which the Si slab extent is measured.
    --surface-thickness FLOAT    Depth (Angstrom) inward from the padded Si
                                 extent (lowest/highest Si position, each
                                 extended by 2 Angstrom); unrestricted beyond
                                 that padding.
    --no-surface-filter          Disable surface filtering and count all
                                  silanol candidates within the padded Si extent.
    --type-si/--type-o/--type-h  Optional explicit particle type IDs.
    --atoms-frame INT            Frame to export full silanol atom membership
                                  (O + bonded Si/H) for. Defaults to the first
                                  analysed frame.

A silanol is an O atom bonded to exactly 1 Si and at least 1 H (Si-O-H).

Output behavior:
    - Creates output directory next to trajectory source: <local parent>/output
    - Writes:
      boundary_si_ids.csv
      silanol_ids_per_frame.jsonl
      silanol_expression_selection.txt
      silanol_atoms_expression_selection.txt
      silanol_counts_per_frame.csv
      silanol_statistics.txt
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from ovito.io import import_file
from ovito.modifiers import CreateBondsModifier


# ---------------------------------------------------------------------------
# Statistics container
# ---------------------------------------------------------------------------

@dataclass
class SilanolStatistics:
    """Summary statistics over all analysed frames."""

    counts_per_frame: List[int] = field(default_factory=list)

    # Populated by finalize().
    n_frames: int = 0
    mean: float = 0.0
    median: float = 0.0
    std: float = 0.0
    sem: float = 0.0          # standard error of the mean
    variance: float = 0.0
    minimum: int = 0
    maximum: int = 0
    count_range: int = 0      # max - min
    q1: float = 0.0
    q3: float = 0.0
    iqr: float = 0.0

    def finalize(self) -> None:
        """Compute all derived statistics from counts_per_frame."""
        n = len(self.counts_per_frame)
        if n == 0:
            return

        self.n_frames = n
        self.minimum = min(self.counts_per_frame)
        self.maximum = max(self.counts_per_frame)
        self.count_range = self.maximum - self.minimum

        self.mean = sum(self.counts_per_frame) / n
        self.variance = sum((x - self.mean) ** 2 for x in self.counts_per_frame) / n
        self.std = math.sqrt(self.variance)
        self.sem = self.std / math.sqrt(n) if n > 1 else 0.0

        sorted_counts = sorted(self.counts_per_frame)
        self.median = _percentile(sorted_counts, 50)
        self.q1 = _percentile(sorted_counts, 25)
        self.q3 = _percentile(sorted_counts, 75)
        self.iqr = self.q3 - self.q1

    def as_text(self) -> str:
        lines = [
            "Silanol count statistics",
            "========================",
            f"  Frames analysed : {self.n_frames}",
            f"  Mean            : {self.mean:.4f}",
            f"  Median          : {self.median:.4f}",
            f"  Std dev         : {self.std:.4f}",
            f"  Variance        : {self.variance:.4f}",
            f"  Std error       : {self.sem:.4f}",
            f"  Minimum         : {self.minimum}",
            f"  Maximum         : {self.maximum}",
            f"  Range           : {self.count_range}",
            f"  Q1 (25th pct)   : {self.q1:.4f}",
            f"  Q3 (75th pct)   : {self.q3:.4f}",
            f"  IQR             : {self.iqr:.4f}",
        ]
        return "\n".join(lines) + "\n"


def _percentile(sorted_data: List[int], pct: float) -> float:
    """Linear-interpolation percentile on a pre-sorted list."""
    n = len(sorted_data)
    if n == 1:
        return float(sorted_data[0])
    idx = pct / 100.0 * (n - 1)
    lo = int(idx)
    hi = lo + 1
    if hi >= n:
        return float(sorted_data[-1])
    frac = idx - lo
    return sorted_data[lo] * (1 - frac) + sorted_data[hi] * frac


# ---------------------------------------------------------------------------
# Bond / type helpers
# ---------------------------------------------------------------------------

def _ensure_bonds(
    pipeline,
    *,
    si_o_cutoff: float = 1.9,
    o_h_cutoff: float = 1.2,
    type_si: Optional[int] = None,
    type_o: Optional[int] = None,
    type_h: Optional[int] = None,
) -> Tuple[int, int, int]:
    """Ensure the pipeline has a CreateBondsModifier with desired cutoffs."""
    cb = next((m for m in pipeline.modifiers if isinstance(m, CreateBondsModifier)), None)
    if cb is None:
        cb = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)
        pipeline.modifiers.append(cb)

    if type_si is None or type_o is None or type_h is None:
        data0 = pipeline.compute(0)
        available = [(int(t.id), str(t.name)) for t in data0.particles.particle_types.types]

        def _norm(s: str) -> str:
            return "".join(ch for ch in s.lower() if ch.isalpha())

        norm_to_id: dict[str, int] = {}
        for tid, tname in available:
            n = _norm(tname)
            if n:
                norm_to_id.setdefault(n, tid)

        if type_si is None:
            type_si = norm_to_id.get("si")
        if type_o is None:
            type_o = norm_to_id.get("o")
        if type_h is None:
            type_h = norm_to_id.get("h")

        if type_si is None or type_o is None or type_h is None:
            raise ValueError(
                "Could not infer particle type IDs for Si/O/H from the file. "
                "Pass them explicitly with --type-si/--type-o/--type-h. "
                f"Available types: {available}"
            )

    t_si = int(type_si)
    t_o = int(type_o)
    t_h = int(type_h)
    cb.set_pairwise_cutoff(t_si, t_o, si_o_cutoff)
    cb.set_pairwise_cutoff(t_o, t_h, o_h_cutoff)
    return t_si, t_o, t_h


# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------

SI_EDGE_PADDING_ANGSTROM = 2.0

def compute_silanols(
    pipeline,
    *,
    frames: Optional[Tuple[int, int]] = None,
    si_o_cutoff: float = 1.9,
    o_h_cutoff: float = 1.2,
    type_si: Optional[int] = None,
    type_o: Optional[int] = None,
    type_h: Optional[int] = None,
    surface_axis: str = "x",
    apply_surface_filter: bool = True,
    surface_thickness_angstrom: float = 5.0,
    atoms_frame: Optional[int] = None,
) -> Tuple[SilanolStatistics, List[List[int]], List[Tuple[int, int]], List[int]]:
    """Return (statistics, ids_per_frame, boundary_si_ids_per_frame, silanol_atom_ids).

    A silanol is an O atom bonded to exactly 1 Si and at least 1 H.
    ids are OVITO particle identifiers if present, otherwise particle indices (0-based).

    The surface region is defined by the lowest and highest Si position along
    surface_axis, each padded outward by SI_EDGE_PADDING_ANGSTROM. Surface
    filtering then keeps only silanols within surface_thickness_angstrom of
    those padded edges. Assumes the surface axis is non-periodic and all Si
    atoms sit strictly inside the box. The other two axes may still be
    periodic; OVITO's CreateBondsModifier handles that automatically via the
    cell's PBC flags.

    silanol_atom_ids collects every atom (the O plus its bonded Si and H
    neighbors) belonging to any silanol in a single frame: atoms_frame if
    given, otherwise the first analysed frame.
    """
    t_Si, t_O, t_H = _ensure_bonds(
        pipeline,
        si_o_cutoff=si_o_cutoff,
        o_h_cutoff=o_h_cutoff,
        type_si=type_si,
        type_o=type_o,
        type_h=type_h,
    )

    stats = SilanolStatistics()
    ids: List[List[int]] = []
    boundary_si_ids_per_frame: List[Tuple[int, int]] = []

    start_frame = frames[0] if frames else 0
    end_frame = frames[1] if frames else len(pipeline.frames)

    effective_atoms_frame = atoms_frame if atoms_frame is not None else start_frame
    if not (start_frame <= effective_atoms_frame < end_frame):
        raise ValueError(
            f"atoms_frame={effective_atoms_frame} is outside the analysed "
            f"frame range [{start_frame}, {end_frame})."
        )
    silanol_atom_ids: List[int] = []

    axis = surface_axis.lower()
    axis_to_idx = {"x": 0, "y": 1, "z": 2}
    if axis not in axis_to_idx:
        raise ValueError(f"surface_axis must be one of x/y/z, got: {surface_axis!r}")
    ax = axis_to_idx[axis]

    for frame in range(start_frame, end_frame):
        data = pipeline.compute(frame)
        ptype = data.particles["Particle Type"].array
        bonds = data.particles.bonds
        pos = data.particles.positions

        si_indices = [i for i in range(data.particles.count) if ptype[i] == t_Si]
        if len(si_indices) == 0:
            raise ValueError("No Si particles found; cannot define surface extent.")

        si_axis_values = [float(pos[i, ax]) for i in si_indices]
        has_pid = "Particle Identifier" in data.particles
        if has_pid:
            pid_arr = data.particles["Particle Identifier"].array
            si_ids = [int(pid_arr[i]) for i in si_indices]
        else:
            si_ids = [int(i) for i in si_indices]

        si_lowest, boundary_lo_id = min(zip(si_axis_values, si_ids))
        si_highest, boundary_hi_id = max(zip(si_axis_values, si_ids))
        boundary_si_ids_per_frame.append((boundary_lo_id, boundary_hi_id))

        lower_bound = si_lowest - SI_EDGE_PADDING_ANGSTROM
        upper_bound = si_highest + SI_EDGE_PADDING_ANGSTROM

        neighbors: List[List[int]] = [[] for _ in range(data.particles.count)]
        for i, j in bonds.topology:
            neighbors[i].append(j)
            neighbors[j].append(i)

        matched_indices: List[int] = []
        for i in range(data.particles.count):
            if ptype[i] != t_O:
                continue
            if not (lower_bound <= float(pos[i, ax]) <= upper_bound):
                continue

            nb = neighbors[i]
            n_si = sum(ptype[j] == t_Si for j in nb)
            n_h = sum(ptype[j] == t_H for j in nb)
            if n_si == 1 and n_h >= 1:
                matched_indices.append(i)

        if apply_surface_filter:
            if surface_thickness_angstrom <= 0:
                matched_indices = []
            else:
                thickness = float(surface_thickness_angstrom)
                matched_indices = [
                    i
                    for i in matched_indices
                    if (
                        float(pos[i, ax]) - lower_bound <= thickness
                        or upper_bound - float(pos[i, ax]) <= thickness
                    )
                ]

        if has_pid:
            matched_ids = [int(pid_arr[i]) for i in matched_indices]
        else:
            matched_ids = matched_indices

        if frame == effective_atoms_frame:
            atom_id_set: set[int] = set()
            for i in matched_indices:
                atom_id_set.add(int(pid_arr[i]) if has_pid else i)
                for j in neighbors[i]:
                    if ptype[j] == t_Si or ptype[j] == t_H:
                        atom_id_set.add(int(pid_arr[j]) if has_pid else j)
            silanol_atom_ids = sorted(atom_id_set)

        stats.counts_per_frame.append(len(matched_ids))
        ids.append(matched_ids)

    stats.finalize()
    return stats, ids, boundary_si_ids_per_frame, silanol_atom_ids


# ---------------------------------------------------------------------------
# Serialisation helpers
# ---------------------------------------------------------------------------

def _boundary_si_csv(boundary_si_ids_per_frame: List[Tuple[int, int]]) -> str:
    lines = ["frame,lower_boundary_si_id,upper_boundary_si_id"]
    lines.extend(
        f"{frame},{lo},{hi}"
        for frame, (lo, hi) in enumerate(boundary_si_ids_per_frame)
    )
    return "\n".join(lines) + "\n"


def _counts_per_frame_csv(stats: SilanolStatistics) -> str:
    lines = ["frame,silanol_count"]
    lines.extend(
        f"{frame},{count}"
        for frame, count in enumerate(stats.counts_per_frame)
    )
    return "\n".join(lines) + "\n"


def _silanol_ids_jsonl(ids_per_frame: List[List[int]]) -> str:
    lines = [
        f'{{"frame": {int(frame)}, "silanol_ids": [{", ".join(str(int(pid)) for pid in frame_ids)}]}}'
        for frame, frame_ids in enumerate(ids_per_frame)
    ]
    return "\n".join(lines) + "\n"


def _ovito_expression_for_particle_identifiers(frame_ids: List[int]) -> str:
    if not frame_ids:
        return "(ParticleIdentifier == -1) && (ParticleIdentifier == -2)"
    return " || ".join(f"(ParticleIdentifier == {int(pid)})" for pid in frame_ids)


def _silanol_expression_selection_text(ids_per_frame: List[List[int]]) -> str:
    lines: List[str] = []
    for frame, frame_ids in enumerate(ids_per_frame):
        lines.append(f"# frame {frame}")
        lines.append(_ovito_expression_for_particle_identifiers(frame_ids))
        lines.append("")
    return "\n".join(lines)


def _infer_frame_wildcard(names: List[str]) -> Optional[str]:
    """Infer a wildcard pattern like 'dump.*.lammpstrj' from per-frame filenames.

    Splits each name into alternating non-digit/digit tokens and requires
    that exactly one digit-run token varies across all files (the frame
    number), with everything else identical. This avoids mis-detecting the
    extension when the frame number is the last dot-separated component
    (e.g. "dump.100000", where Path.suffix would wrongly be ".100000").
    """
    tokenized = [re.split(r"(\d+)", n) for n in names]
    if len({len(t) for t in tokenized}) != 1:
        return None

    n_tokens = len(tokenized[0])
    varying = [i for i in range(n_tokens) if len({t[i] for t in tokenized}) > 1]
    if len(varying) != 1 or varying[0] % 2 == 0:
        return None

    template = list(tokenized[0])
    template[varying[0]] = "*"
    return "".join(template)


def _resolve_local_import_path(local_path: Path) -> str:
    """Return the path/pattern to hand to OVITO's import_file.

    If local_path is a directory containing one file per frame, build a
    wildcard pattern (e.g. "<dir>/dump.*.lammpstrj") so OVITO loads the
    sorted files as a multi-frame sequence. Otherwise return the file path
    unchanged.
    """
    if not local_path.is_dir():
        return str(local_path)

    files = sorted(p for p in local_path.iterdir() if p.is_file() and not p.name.startswith("."))
    if not files:
        raise ValueError(f"No frame files found in directory: {local_path}")
    if len(files) == 1:
        return str(files[0])

    pattern_name = _infer_frame_wildcard([p.name for p in files])
    if pattern_name is None:
        raise ValueError(
            "Could not infer a per-frame filename pattern from the files in "
            f"{local_path}. Expected filenames to share a common prefix/suffix "
            f"around a single varying frame number. Example filenames: "
            f"{[p.name for p in files[:5]]}"
        )
    return str(local_path / pattern_name)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def _parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute silanol count (OVITO).")
    p.add_argument("--local-path", type=Path, required=True,
                   help="Local trajectory file, or a directory with one file per "
                        "frame; creates '<local parent>/output'.")
    p.add_argument("--type-si", type=int, default=None, help="Particle type ID for Si.")
    p.add_argument("--type-o",  type=int, default=None, help="Particle type ID for O.")
    p.add_argument("--type-h",  type=int, default=None, help="Particle type ID for H.")
    p.add_argument("--surface-axis", default="x", choices=["x", "y", "z"],
                   help="Surface normal axis used to identify top/bottom surfaces.")
    p.add_argument("--surface-thickness", type=float, default=5.0,
                   help="Depth in Å inward from the padded Si extent "
                        "(lowest/highest Si position, each extended by "
                        f"{SI_EDGE_PADDING_ANGSTROM} Å).")
    p.add_argument("--no-surface-filter", action="store_true",
                   help="Count all silanols (do not restrict to surface slabs).")
    p.add_argument("--frames", nargs=2, type=int, default=None,
                   metavar=("START", "END"),
                   help="Analyse only frames [START, END) of a LAMMPS dump.")
    p.add_argument("--atoms-frame", type=int, default=None,
                   help="Frame index to export the full silanol atom "
                        "membership (O + bonded Si/H) for. Defaults to the "
                        "first analysed frame.")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Output writing
# ---------------------------------------------------------------------------

def _write_outputs_local(
    output_dir: Path,
    *,
    stats: SilanolStatistics,
    ids: List[List[int]],
    boundary_si_ids_per_frame: List[Tuple[int, int]],
    silanol_atom_ids: List[int],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "boundary_si_ids.csv").write_text(_boundary_si_csv(boundary_si_ids_per_frame))
    (output_dir / "silanol_counts_per_frame.csv").write_text(_counts_per_frame_csv(stats))
    (output_dir / "silanol_ids_per_frame.jsonl").write_text(_silanol_ids_jsonl(ids))
    (output_dir / "silanol_expression_selection.txt").write_text(
        _silanol_expression_selection_text(ids)
    )
    (output_dir / "silanol_statistics.txt").write_text(stats.as_text())
    (output_dir / "silanol_atoms_expression_selection.txt").write_text(
        _ovito_expression_for_particle_identifiers(silanol_atom_ids) + "\n"
    )
    print(f"output_dir={output_dir}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> None:
    args = _parse_args(argv)

    local_input_path = args.local_path.expanduser().resolve()
    output_dir = local_input_path.parent / "output"
    analysis_frames = (
        (int(args.frames[0]), int(args.frames[1])) if args.frames is not None else None
    )

    import_path = _resolve_local_import_path(local_input_path)
    pipeline = import_file(import_path)

    stats, ids, boundary_si_ids_per_frame, silanol_atom_ids = compute_silanols(
        pipeline,
        frames=analysis_frames,
        type_si=args.type_si,
        type_o=args.type_o,
        type_h=args.type_h,
        surface_axis=args.surface_axis,
        apply_surface_filter=not args.no_surface_filter,
        surface_thickness_angstrom=args.surface_thickness,
        atoms_frame=args.atoms_frame,
    )

    print(stats.as_text())
    print(f"boundary_si_ids[:20]={boundary_si_ids_per_frame[:20]}")

    _write_outputs_local(
        output_dir,
        stats=stats,
        ids=ids,
        boundary_si_ids_per_frame=boundary_si_ids_per_frame,
        silanol_atom_ids=silanol_atom_ids,
    )


if __name__ == "__main__":
    main()
