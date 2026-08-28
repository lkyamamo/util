#!/usr/bin/env python3
"""Count surface silanols and export analysis artifacts with OVITO.

Usage:
    python silanol.py --local-path PATH [options]

Primary options:
    --local-path PATH            Local trajectory file, or a local directory
                                  containing one file per frame (files are
                                  sorted and loaded as a sequence).
    --frames START END           Optional frame window [START, END).
    --type-si/--type-o/--type-h  Numeric particle type IDs for Si/O/H. Only
                                  needed if the file has no element names
                                  (opaque numeric types); ignored if it does.
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
If the file labels particle types by element, Si/O/H are found by matching
each type's name to its element symbol (case-insensitive) and --type-si/
--type-o/--type-h are ignored. If the file only has opaque numeric types (no
element names), --type-si/--type-o/--type-h are required.

With --surface-method alpha-shape, the area of the surface mesh separating
the slab from the surrounding medium (i.e. both exposed surfaces) is measured
per frame in nm^2, and a surface silanol density (silanols in that frame
divided by that frame's surface area, nm^-2) is reported alongside the raw
counts. The padded-extent method has no mesh, so it reports counts only.

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
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
from ovito.io import import_file
from ovito.modifiers import ConstructSurfaceModifier, CreateBondsModifier, SelectTypeModifier


# ---------------------------------------------------------------------------
# Statistics container
# ---------------------------------------------------------------------------

@dataclass
class SeriesStatistics:
    """Summary statistics for one per-frame series of values."""

    n_frames: int = 0
    mean: float = 0.0
    median: float = 0.0
    std: float = 0.0
    sem: float = 0.0          # standard error of the mean
    variance: float = 0.0
    minimum: float = 0.0
    maximum: float = 0.0
    value_range: float = 0.0  # max - min
    q1: float = 0.0
    q3: float = 0.0
    iqr: float = 0.0

    def as_text(self, title: str, *, fmt: str = ".4f",
                extreme_fmt: Optional[str] = None) -> str:
        ext = fmt if extreme_fmt is None else extreme_fmt
        lines = [
            title,
            "=" * len(title),
            f"  Frames analysed : {self.n_frames}",
            f"  Mean            : {self.mean:{fmt}}",
            f"  Median          : {self.median:{fmt}}",
            f"  Std dev         : {self.std:{fmt}}",
            f"  Variance        : {self.variance:{fmt}}",
            f"  Std error       : {self.sem:{fmt}}",
            f"  Minimum         : {self.minimum:{ext}}",
            f"  Maximum         : {self.maximum:{ext}}",
            f"  Range           : {self.value_range:{ext}}",
            f"  Q1 (25th pct)   : {self.q1:{fmt}}",
            f"  Q3 (75th pct)   : {self.q3:{fmt}}",
            f"  IQR             : {self.iqr:{fmt}}",
        ]
        return "\n".join(lines) + "\n"


def _summarize(values: Sequence[float]) -> SeriesStatistics:
    """Compute the full summary of one per-frame series."""
    summary = SeriesStatistics()
    n = len(values)
    if n == 0:
        return summary

    summary.n_frames = n
    summary.minimum = float(min(values))
    summary.maximum = float(max(values))
    summary.value_range = summary.maximum - summary.minimum

    summary.mean = sum(values) / n
    summary.variance = sum((x - summary.mean) ** 2 for x in values) / n
    summary.std = math.sqrt(summary.variance)
    summary.sem = summary.std / math.sqrt(n) if n > 1 else 0.0

    ordered = sorted(float(v) for v in values)
    summary.median = _percentile(ordered, 50)
    summary.q1 = _percentile(ordered, 25)
    summary.q3 = _percentile(ordered, 75)
    summary.iqr = summary.q3 - summary.q1
    return summary


@dataclass
class SilanolStatistics:
    """Per-frame series over all analysed frames, plus their summaries.

    surface_areas_per_frame (nm^2) and densities_per_frame (nm^-2) are only
    populated by the alpha-shape surface method, which is the only one that
    builds a mesh to measure an area with; they stay empty otherwise.
    """

    counts_per_frame: List[int] = field(default_factory=list)
    surface_areas_per_frame: List[float] = field(default_factory=list)
    densities_per_frame: List[float] = field(default_factory=list)

    # Populated by finalize(); area/density stay None without a surface mesh.
    counts: SeriesStatistics = field(default_factory=SeriesStatistics)
    area: Optional[SeriesStatistics] = None
    density: Optional[SeriesStatistics] = None

    def finalize(self) -> None:
        """Compute all derived statistics from the per-frame series."""
        self.counts = _summarize(self.counts_per_frame)
        if self.densities_per_frame:
            self.area = _summarize(self.surface_areas_per_frame)
            self.density = _summarize(self.densities_per_frame)

    def as_text(self) -> str:
        text = self.counts.as_text("Silanol count statistics", extreme_fmt=".0f")
        if self.area is not None and self.density is not None:
            text += "\n" + self.area.as_text("Surface area statistics (nm^2)")
            text += "\n" + self.density.as_text(
                "Surface silanol density statistics (nm^-2)"
            )
        return text


def _percentile(sorted_data: Sequence[float], pct: float) -> float:
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
    """Ensure the pipeline has a CreateBondsModifier with desired cutoffs.

    Reads the file's particle types to decide how Si/O/H are identified:
    - If types are labelled by element (e.g. a LAMMPS dump "element"
      column), Si/O/H are found by matching each type's name to its element
      symbol (case-insensitive) — type_si/type_o/type_h are ignored.
    - If types are only opaque numeric IDs (no element names), type_si/
      type_o/type_h must be passed explicitly (--type-si/--type-o/
      --type-h).
    - If the file has no particle types at all, raises ValueError.
    """
    data0 = pipeline.compute(0)
    if "Particle Type" not in data0.particles:
        raise ValueError(
            "The input file provides no particle types (neither element "
            "names nor numeric types); cannot identify Si/O/H atoms."
        )

    available = [(int(t.id), str(t.name)) for t in data0.particles.particle_types.types]
    has_element_names = any(tname.strip() for _, tname in available)

    if has_element_names:
        type_by_symbol = {tname.strip().lower(): tid for tid, tname in available}
        t_si = type_by_symbol.get("si")
        t_o = type_by_symbol.get("o")
        t_h = type_by_symbol.get("h")
        if t_si is None or t_o is None or t_h is None:
            raise ValueError(
                "File labels types by element, but Si/O/H were not all "
                f"found among them. Available types: {available}"
            )
    else:
        if type_si is None or type_o is None or type_h is None:
            raise ValueError(
                "File provides only numeric particle types (no element "
                "names); pass --type-si/--type-o/--type-h explicitly. "
                f"Available types: {available}"
            )
        t_si, t_o, t_h = int(type_si), int(type_o), int(type_h)

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


def _minimum_image_delta(cell):
    """Return a delta(a, b) -> a - b that applies the minimum image convention.

    Mesh vertices are stored wrapped into the cell, so a face straddling a
    periodic boundary has edge vectors that look box-sized unless they are
    folded back. Faces are far smaller than the cell, so the minimum image is
    always the intended edge.
    """
    matrix = np.asarray(cell[...], dtype=float)[:, :3]
    inverse = np.linalg.inv(matrix)
    pbc = [bool(flag) for flag in cell.pbc]

    def delta(a, b):
        reduced = inverse @ (np.asarray(a, dtype=float) - np.asarray(b, dtype=float))
        for k, periodic in enumerate(pbc):
            if periodic:
                reduced[k] -= round(reduced[k])
        return matrix @ reduced

    return delta


def _alpha_shape_surface_si_and_area(data) -> Tuple[set[int], float]:
    """Si atoms on the true surrounding-medium surface, and that surface's area.

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

    The returned area is the summed area of exactly those triangles, i.e. the
    whole slab/medium interface. For a slab that is periodic in the two
    in-plane directions this is the combined area of both exposed surfaces
    (internal void surfaces are excluded, since their faces point into other
    regions). Each boundary triangle contributes once, because only one of
    its two oppositely-oriented faces is tagged with the medium region.

    The area is in A^2, the units of the mesh vertex coordinates; callers
    convert (see SQ_ANGSTROM_PER_SQ_NM).
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
    vertex_pos = surf.vertices["Position"].array
    topo = surf.topology
    delta = _minimum_image_delta(surf.domain)

    surface_si_indices: set[int] = set()
    surface_area = 0.0
    for face in range(topo.face_count):
        if int(face_region[face]) != medium_region_id:
            continue

        face_vertices: List[int] = []
        e0 = topo.first_face_edge(face)
        e = e0
        while True:
            v = topo.first_edge_vertex(e)
            face_vertices.append(v)
            surface_si_indices.add(int(particle_index[v]))
            e = topo.next_face_edge(e)
            if e == e0:
                break

        # Fan-triangulate from the first vertex; alpha-shape faces are already
        # triangles, so this is a single cross product in practice.
        origin = vertex_pos[face_vertices[0]]
        for k in range(1, len(face_vertices) - 1):
            edge_a = delta(vertex_pos[face_vertices[k]], origin)
            edge_b = delta(vertex_pos[face_vertices[k + 1]], origin)
            surface_area += 0.5 * float(np.linalg.norm(np.cross(edge_a, edge_b)))

    return surface_si_indices, surface_area


# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------

SI_EDGE_PADDING_ANGSTROM = 2.0
SQ_ANGSTROM_PER_SQ_NM = 100.0

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

    type_si/type_o/type_h are only used if the file's particle types are
    opaque numeric IDs with no element names; if the file labels types by
    element, Si/O/H are found by matching type names to element symbols and
    these arguments are ignored. See _ensure_bonds.

    Two surface_method choices for deciding which silanols count as "surface":
    - "padded-extent" (default): the lowest and highest Si position along
      surface_axis, each padded outward by SI_EDGE_PADDING_ANGSTROM, define
      two flat planes; a silanol counts if its O is within
      surface_thickness_angstrom of either plane. Assumes a flat slab.
    - "alpha-shape": builds an alpha-shape surface mesh from the Si
      sublattice (radius=alpha_shape_radius, smoothing_level=
      alpha_shape_smoothing) and, via region topology, finds Si atoms
      bordering the surrounding medium (see _alpha_shape_surface_si_and_area)
      — a silanol counts if its bonded Si is one of those atoms. Shape-agnostic:
      treats flat and curved surfaces identically, with no distance
      threshold. surface_thickness_angstrom is not used by this method.
      This method also measures the area of that mesh surface per frame
      (nm^2) and records a surface silanol density in nm^-2 (that frame's
      silanol count divided by that frame's surface area) on the returned
      statistics.

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
        type_si=type_si,
        type_o=type_o,
        type_h=type_h,
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
        surface_area: Optional[float] = None
        if apply_surface_filter:
            if surface_method == "alpha-shape":
                surface_si_indices, surface_area = _alpha_shape_surface_si_and_area(data)
                if surface_area <= 0.0:
                    raise ValueError(
                        f"Frame {frame}: the alpha-shape mesh has no surface "
                        "bordering the surrounding medium (zero area); try a "
                        "different --alpha-shape-radius."
                    )
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
        if surface_area is not None:
            area_nm2 = surface_area / SQ_ANGSTROM_PER_SQ_NM
            stats.surface_areas_per_frame.append(area_nm2)
            stats.densities_per_frame.append(len(matched_ids) / area_nm2)
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
    """Per-frame counts, plus area and density columns when a mesh was built."""
    if not stats.densities_per_frame:
        lines = ["frame,silanol_count"]
        lines.extend(
            f"{frame},{count}"
            for frame, count in enumerate(stats.counts_per_frame)
        )
        return "\n".join(lines) + "\n"

    lines = ["frame,silanol_count,surface_area_nm2,surface_density_per_nm2"]
    lines.extend(
        f"{frame},{count},{area:.6f},{density:.6f}"
        for frame, (count, area, density) in enumerate(
            zip(stats.counts_per_frame, stats.surface_areas_per_frame,
                stats.densities_per_frame)
        )
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
    p.add_argument("--type-si", type=int, default=None,
                   help="Numeric particle type ID for Si. Only needed if the "
                        "file has no element names (opaque numeric types).")
    p.add_argument("--type-o", type=int, default=None,
                   help="Numeric particle type ID for O. Only needed if the "
                        "file has no element names (opaque numeric types).")
    p.add_argument("--type-h", type=int, default=None,
                   help="Numeric particle type ID for H. Only needed if the "
                        "file has no element names (opaque numeric types).")
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
                        "for curved surfaces too), and additionally reports "
                        "the per-frame mesh surface area and surface silanol "
                        "density.")
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
            type_si=args.type_si,
            type_o=args.type_o,
            type_h=args.type_h,
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
