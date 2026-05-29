"""Format-agnostic bbox filtering for trajectory frames.

Assumptions (enforced by callers, not checked here):
  - natoms is constant across all source frames (NVT/NVE; grand-canonical not supported)
  - Every bbox-filtered frame contains at least one atom
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional


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

    Coordinates are taken from frame.positions (parallel to frame.atom_lines).
    Original atom line text is preserved verbatim for kept atoms.
    header_lines, timestep, and extras are copied unchanged.
    natoms is updated to reflect the kept count.
    """
    kept_positions: list[tuple[float, float, float]] = []
    kept_lines: list[str] = []

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
