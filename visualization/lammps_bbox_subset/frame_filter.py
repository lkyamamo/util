"""Format-agnostic bbox filtering for trajectory frames.

Assumptions (enforced by callers, not checked here):
  - natoms is constant across all source frames (NVT/NVE; grand-canonical not supported)
  - Every bbox-filtered frame contains at least one atom

Performance note:
  filter_by_bbox uses a fully vectorized NumPy boolean test when NumPy is
  available — a single C-level operation on the whole position array:

      mask = np.all((xyz >= lo) & (xyz <= hi), axis=1)

  This is equivalent to the reference pattern from filter_lammps.py and is
  ~10-50× faster than a per-atom Python comparison for large atom counts.
  Falls back to a pure-Python loop if NumPy is not installed.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

try:
    import numpy as np
    _HAS_NUMPY = True
except ImportError:
    _HAS_NUMPY = False


@dataclass
class BBox:
    xmin: float
    xmax: float
    ymin: float
    ymax: float
    zmin: float
    zmax: float

    def validate(self) -> None:
        for lo, hi, axis in (
            (self.xmin, self.xmax, "x"),
            (self.ymin, self.ymax, "y"),
            (self.zmin, self.zmax, "z"),
        ):
            if lo >= hi:
                raise ValueError(
                    f"bbox {axis}min ({lo}) must be strictly less than {axis}max ({hi})"
                )


@dataclass
class Frame:
    """Canonical in-memory representation of one trajectory frame.

    Adapters populate this from format-specific bytes; filter_by_bbox operates
    only on positions and atom_lines.  header_lines and extras are passed
    through verbatim to write_frame.
    """

    timestep: Optional[int]
    natoms: int
    positions: list[tuple[float, float, float]]  # parallel to atom_lines
    atom_lines: list[str]                         # raw text, written verbatim for kept atoms
    header_lines: list[str]                       # everything before atom data; passed through
    extras: dict = field(default_factory=dict)    # e.g. atoms_header_str, column indices


def filter_by_bbox(frame: Frame, bbox: BBox) -> Frame:
    """Return a new Frame containing only atoms whose (x, y, z) lie inside bbox.

    When NumPy is available the bounding-box test is a single vectorized C-level
    operation on the entire position array (no Python loop over atoms):

        xyz  = np.array(frame.positions)                 # (N, 3)
        lo   = np.array([xmin, ymin, zmin])
        hi   = np.array([xmax, ymax, zmax])
        mask = np.all((xyz >= lo) & (xyz <= hi), axis=1) # (N,) bool

    Atom lines are written verbatim for kept atoms — coordinates, ids, types,
    and all other columns are preserved at original precision and column layout.
    header_lines, timestep, and extras are copied unchanged.
    natoms is updated to reflect the kept count.
    """
    if _HAS_NUMPY and frame.natoms > 0:
        xyz = np.array(frame.positions, dtype=np.float64)        # (N, 3)
        lo  = np.array([bbox.xmin, bbox.ymin, bbox.zmin], dtype=np.float64)
        hi  = np.array([bbox.xmax, bbox.ymax, bbox.zmax], dtype=np.float64)
        mask = np.all((xyz >= lo) & (xyz <= hi), axis=1)          # (N,) bool

        # NumPy object-array indexing for string selection — no Python loop
        atom_arr   = np.array(frame.atom_lines, dtype=object)
        kept_lines: list[str] = atom_arr[mask].tolist()
        kept_positions: list[tuple[float, float, float]] = [
            tuple(row) for row in xyz[mask].tolist()              # type: ignore[misc]
        ]
    else:
        # Pure-Python fallback (no NumPy, or empty frame)
        kept_positions = []
        kept_lines = []
        for pos, line in zip(frame.positions, frame.atom_lines):
            x, y, z = pos
            if (
                bbox.xmin <= x <= bbox.xmax
                and bbox.ymin <= y <= bbox.ymax
                and bbox.zmin <= z <= bbox.zmax
            ):
                kept_positions.append(pos)
                kept_lines.append(line)

    return Frame(
        timestep=frame.timestep,
        natoms=len(kept_lines),
        positions=kept_positions,
        atom_lines=kept_lines,
        header_lines=frame.header_lines,
        extras=frame.extras,
    )
