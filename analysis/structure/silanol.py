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
    --surface-method {padded-extent,alpha-shape}
                                  How to decide which silanols are "surface"
                                  silanols. padded-extent (default) uses two
                                  flat planes at the lowest/highest Si
                                  position; only correct for flat slabs.
                                  alpha-shape builds an alpha-shape surface
                                  mesh from the Si sublattice and uses region
                                  topology to find Si atoms bordering the
                                  surrounding medium (vacuum, water, etc.);
                                  shape-agnostic (also works for curved
                                  surfaces).
    --surface-thickness FLOAT    Depth (Angstrom) inward from the padded Si
                                 extent (lowest/highest Si position, each
                                 extended by 2 Angstrom); unrestricted beyond
                                 that padding. Only used by --surface-method
                                 padded-extent.
    --alpha-shape-radius FLOAT   Probe sphere radius (Angstrom) for
                                  --surface-method alpha-shape.
    --alpha-shape-smoothing INT  Mesh smoothing iterations for
                                  --surface-method alpha-shape.
    --no-surface-filter          Disable surface filtering and count all
                                  silanol candidates within the padded Si extent.
    --atoms-frame INT            Frame to export full silanol atom membership
                                  (O + bonded Si/H) for. Defaults to the first
                                  analysed frame.

A silanol is an O atom bonded to exactly 1 Si and at least 1 H (Si-O-H).
Si/O/H particle types are found by matching each type's name to its element
symbol (case-insensitive) — the file must label types by element.

Output behavior:
    - Creates an "output" directory next to this script (not next to the
      trajectory source) and writes:
      boundary_si_ids.csv
      silanol_ids_per_frame.jsonl
      silanol_expression_selection.txt
      silanol_atoms_expression_selection.txt
      silanol_counts_per_frame.csv
      silanol_statistics.txt
      surface_si_ids_per_frame.jsonl
      surface_si_expression_selection.txt
"""

from __future__ import annotations

import argparse
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from ovito.io import import_file
from ovito.modifiers import ConstructSurfaceModifier, CreateBondsModifier, SelectTypeModifier


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
) -> Tuple[int, int, int]:
    """Ensure the pipeline has a CreateBondsModifier with desired cutoffs.

    Si/O/H type IDs are found by matching each particle type's name against
    the element symbols Si/O/H (case-insensitive) — the file is expected to
    label types by element, not by opaque numeric type alone.
    """
    data0 = pipeline.compute(0)
    available = [(int(t.id), str(t.name)) for t in data0.particles.particle_types.types]
    type_by_symbol = {tname.strip().lower(): tid for tid, tname in available}

    t_si = type_by_symbol.get("si")
    t_o = type_by_symbol.get("o")
    t_h = type_by_symbol.get("h")
    if t_si is None or t_o is None or t_h is None:
        raise ValueError(
            "Could not find particle types named Si/O/H (case-insensitive) "
            f"in the file. Available types: {available}"
        )

    cb = next((m for m in pipeline.modifiers if isinstance(m, CreateBondsModifier)), None)
    if cb is None:
        cb = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)
        pipeline.modifiers.append(cb)
    cb.set_pairwise_cutoff(t_si, t_o, si_o_cutoff)
    cb.set_pairwise_cutoff(t_o, t_h, o_h_cutoff)
    return t_si, t_o, t_h


def _ensure_alpha_shape_surface(
    pipeline, *, type_si: int, radius: float, smoothing_level: int
) -> None:
    """Add a Si-only alpha-shape surface mesh (with region identification) to the pipeline."""
    select_si = SelectTypeModifier(operate_on="particles", property="Particle Type")
    select_si.types = {int(type_si)}
    pipeline.modifiers.append(select_si)
    pipeline.modifiers.append(ConstructSurfaceModifier(
        method=ConstructSurfaceModifier.Method.AlphaShape,
        radius=radius,
        smoothing_level=smoothing_level,
        identify_regions=True,
        only_selected=True,
    ))


def _alpha_shape_surface_si_indices(data) -> set[int]:
    """Particle indices of Si atoms on the true surrounding-medium surface, any shape.

    Distinguishes the real surface from internal voids/pores by region
    topology (OVITO's identify_regions), not by position, so it treats flat
    and curved surfaces the same way, and doesn't care whether the medium on
    the other side is vacuum or something else (e.g. water in a fully
    periodic cell, where nothing is ever "Exterior" since every region is
    topologically enclosed). The surrounding medium is simply the largest
    non-filled region by volume — internal voids/pores in an amorphous
    network are orders of magnitude smaller.

    A mesh face's Region tag is whichever region its normal points into; a
    boundary triangle is represented as two oppositely-oriented faces (one
    into the filled region, one into whatever region is on the other side),
    so filtering faces to those pointing into the largest non-filled region
    picks out exactly the true-surface triangles.
    """
    surf = data.surfaces["surface"]
    filled = surf.regions["Filled"].array
    volume = surf.regions["Volume"].array
    medium_region_id = max(
        (i for i, is_filled in enumerate(filled) if not is_filled),
        key=lambda i: volume[i],
    )

    face_region = surf.faces["Region"].array
    particle_index = surf.vertices["Particle Index"].array
    topo = surf.topology

    surface_si_indices: set[int] = set()
    for face in range(topo.face_count):
        if int(face_region[face]) != medium_region_id:
            continue
        e0 = topo.first_face_edge(face)
        e = e0
        while True:
            surface_si_indices.add(int(particle_index[topo.first_edge_vertex(e)]))
            e = topo.next_face_edge(e)
            if e == e0:
                break
    return surface_si_indices


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
    surface_axis: str = "x",
    apply_surface_filter: bool = True,
    surface_thickness_angstrom: float = 5.0,
    surface_method: str = "padded-extent",
    alpha_shape_radius: float = 4.2,
    alpha_shape_smoothing: int = 0,
    atoms_frame: Optional[int] = None,
) -> Tuple[
    SilanolStatistics, List[List[int]], List[Tuple[int, int]], List[int], List[List[int]]
]:
    """Return (statistics, ids_per_frame, boundary_si_ids_per_frame,
    silanol_atom_ids, surface_si_ids_per_frame).

    A silanol is an O atom bonded to exactly 1 Si and at least 1 H.
    ids are OVITO particle identifiers if present, otherwise particle indices (0-based).

    Two surface_method choices for deciding which silanols count as "surface":
    - "padded-extent" (default): the lowest and highest Si position along
      surface_axis, each padded outward by SI_EDGE_PADDING_ANGSTROM, define
      two flat planes; a silanol counts if its O is within
      surface_thickness_angstrom of either plane. Assumes a flat slab.
    - "alpha-shape": builds an alpha-shape surface mesh from the Si
      sublattice (radius=alpha_shape_radius, smoothing_level=
      alpha_shape_smoothing) and, via region topology, finds Si atoms
      bordering the surrounding medium (see _alpha_shape_surface_si_indices)
      — a silanol counts if its bonded Si is one of those atoms. Shape-agnostic:
      treats flat and curved surfaces identically, with no distance
      threshold. surface_thickness_angstrom is not used by this method.

    Assumes the surface axis is non-periodic and all Si atoms sit strictly
    inside the box. The other two axes may still be periodic; OVITO's
    CreateBondsModifier (and ConstructSurfaceModifier) handle that
    automatically via the cell's PBC flags.

    silanol_atom_ids collects every atom (the O plus its bonded Si and H
    neighbors) belonging to any silanol in a single frame: atoms_frame if
    given, otherwise the first analysed frame.

    surface_si_ids_per_frame lists, for every analysed frame, the Si atoms
    classified as "surface" by surface_method (empty per frame if
    apply_surface_filter is False, since no surface classification applies).
    """
    if surface_method not in ("padded-extent", "alpha-shape"):
        raise ValueError(
            f"surface_method must be 'padded-extent' or 'alpha-shape', got: {surface_method!r}"
        )

    t_Si, t_O, t_H = _ensure_bonds(
        pipeline,
        si_o_cutoff=si_o_cutoff,
        o_h_cutoff=o_h_cutoff,
    )

    if apply_surface_filter and surface_method == "alpha-shape":
        _ensure_alpha_shape_surface(
            pipeline,
            type_si=t_Si,
            radius=alpha_shape_radius,
            smoothing_level=alpha_shape_smoothing,
        )

    stats = SilanolStatistics()
    ids: List[List[int]] = []
    boundary_si_ids_per_frame: List[Tuple[int, int]] = []
    surface_si_ids_per_frame: List[List[int]] = []

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

        surface_si_indices: set[int] = set()
        if apply_surface_filter:
            if surface_method == "alpha-shape":
                surface_si_indices = _alpha_shape_surface_si_indices(data)
            elif surface_thickness_angstrom > 0:
                thickness = float(surface_thickness_angstrom)
                surface_si_indices = {
                    i
                    for i, v in zip(si_indices, si_axis_values)
                    if v - si_lowest <= thickness or si_highest - v <= thickness
                }

        if has_pid:
            surface_si_ids_per_frame.append(sorted(int(pid_arr[i]) for i in surface_si_indices))
        else:
            surface_si_ids_per_frame.append(sorted(surface_si_indices))

        matched_indices: List[int] = []
        for i in range(data.particles.count):
            if ptype[i] != t_O:
                continue
            if not (lower_bound <= float(pos[i, ax]) <= upper_bound):
                continue

            nb = neighbors[i]
            si_neighbors = [j for j in nb if ptype[j] == t_Si]
            n_h = sum(ptype[j] == t_H for j in nb)
            if len(si_neighbors) != 1 or n_h < 1:
                continue

            if apply_surface_filter:
                if surface_method == "padded-extent":
                    if surface_thickness_angstrom <= 0:
                        continue
                    p = float(pos[i, ax])
                    thickness = float(surface_thickness_angstrom)
                    if not (p - lower_bound <= thickness or upper_bound - p <= thickness):
                        continue
                else:  # alpha-shape
                    if si_neighbors[0] not in surface_si_indices:
                        continue

            matched_indices.append(i)

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
    return stats, ids, boundary_si_ids_per_frame, silanol_atom_ids, surface_si_ids_per_frame


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


def _ids_per_frame_jsonl(ids_per_frame: List[List[int]], key: str) -> str:
    lines = [
        f'{{"frame": {int(frame)}, "{key}": [{", ".join(str(int(pid)) for pid in frame_ids)}]}}'
        for frame, frame_ids in enumerate(ids_per_frame)
    ]
    return "\n".join(lines) + "\n"


def _ovito_expression_for_particle_identifiers(frame_ids: List[int]) -> str:
    if not frame_ids:
        return "(ParticleIdentifier == -1) && (ParticleIdentifier == -2)"
    return " || ".join(f"(ParticleIdentifier == {int(pid)})" for pid in frame_ids)


def _expression_selection_text_per_frame(ids_per_frame: List[List[int]]) -> str:
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
                        "frame.")
    p.add_argument("--surface-axis", default="x", choices=["x", "y", "z"],
                   help="Surface normal axis used to identify top/bottom surfaces.")
    p.add_argument("--surface-method", default="padded-extent",
                   choices=["padded-extent", "alpha-shape"],
                   help="How to decide which silanols are 'surface' silanols. "
                        "'padded-extent' (default) uses two flat planes at "
                        "the lowest/highest Si position; only correct for "
                        "flat slabs. 'alpha-shape' builds an alpha-shape "
                        "surface mesh from the Si sublattice and uses region "
                        "topology to find Si atoms bordering the surrounding "
                        "medium (vacuum, water, etc.); shape-agnostic (works "
                        "for curved surfaces too).")
    p.add_argument("--surface-thickness", type=float, default=5.0,
                   help="Depth in Å inward from the padded Si extent "
                        "(lowest/highest Si position, each extended by "
                        f"{SI_EDGE_PADDING_ANGSTROM} Å). Only used by "
                        "--surface-method padded-extent.")
    p.add_argument("--alpha-shape-radius", type=float, default=4.2,
                   help="Probe sphere radius in Å for --surface-method "
                        "alpha-shape.")
    p.add_argument("--alpha-shape-smoothing", type=int, default=0,
                   help="Mesh smoothing iterations for --surface-method "
                        "alpha-shape.")
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
    surface_si_ids_per_frame: List[List[int]],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "boundary_si_ids.csv").write_text(_boundary_si_csv(boundary_si_ids_per_frame))
    (output_dir / "silanol_counts_per_frame.csv").write_text(_counts_per_frame_csv(stats))
    (output_dir / "silanol_ids_per_frame.jsonl").write_text(
        _ids_per_frame_jsonl(ids, "silanol_ids")
    )
    (output_dir / "silanol_expression_selection.txt").write_text(
        _expression_selection_text_per_frame(ids)
    )
    (output_dir / "silanol_statistics.txt").write_text(stats.as_text())
    (output_dir / "silanol_atoms_expression_selection.txt").write_text(
        _ovito_expression_for_particle_identifiers(silanol_atom_ids) + "\n"
    )
    (output_dir / "surface_si_ids_per_frame.jsonl").write_text(
        _ids_per_frame_jsonl(surface_si_ids_per_frame, "surface_si_ids")
    )
    (output_dir / "surface_si_expression_selection.txt").write_text(
        _expression_selection_text_per_frame(surface_si_ids_per_frame)
    )
    print(f"output_dir={output_dir}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> None:
    args = _parse_args(argv)

    local_input_path = args.local_path.expanduser().resolve()
    output_dir = Path(__file__).resolve().parent / "output"
    analysis_frames = (
        (int(args.frames[0]), int(args.frames[1])) if args.frames is not None else None
    )

    import_path = _resolve_local_import_path(local_input_path)
    pipeline = import_file(import_path)

    stats, ids, boundary_si_ids_per_frame, silanol_atom_ids, surface_si_ids_per_frame = (
        compute_silanols(
            pipeline,
            frames=analysis_frames,
            surface_axis=args.surface_axis,
            apply_surface_filter=not args.no_surface_filter,
            surface_thickness_angstrom=args.surface_thickness,
            surface_method=args.surface_method,
            alpha_shape_radius=args.alpha_shape_radius,
            alpha_shape_smoothing=args.alpha_shape_smoothing,
            atoms_frame=args.atoms_frame,
        )
    )

    print(stats.as_text())
    print(f"boundary_si_ids[:20]={boundary_si_ids_per_frame[:20]}")

    _write_outputs_local(
        output_dir,
        stats=stats,
        ids=ids,
        boundary_si_ids_per_frame=boundary_si_ids_per_frame,
        silanol_atom_ids=silanol_atom_ids,
        surface_si_ids_per_frame=surface_si_ids_per_frame,
    )


if __name__ == "__main__":
    main()
