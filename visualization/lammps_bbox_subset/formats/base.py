"""Adapter protocol and shared FrameSpec dataclass."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Protocol, TextIO

from frame_filter import Frame


@dataclass
class FrameSpec:
    """Structural constants probed from the first frame of a SOURCE trajectory.

    Valid only for source files (natoms constant).  Do NOT reuse for part
    files, which have variable filtered natoms per frame.
    """

    natoms: int
    lines_per_frame: int   # header_lines + natoms — constant for source
    header_lines: int      # lines before atom data in each frame
    x_col: int             # 0-indexed token position in each atom line
    y_col: int
    z_col: int
    timestep_line: int     # 0-indexed line within a frame that holds the timestep value
    extra: dict = field(default_factory=dict)  # e.g. atoms_header_str, column name list


class TrajectoryAdapter(Protocol):
    """Interface every format adapter must satisfy."""

    def detect(self, path: str) -> bool:
        """Return True if the file at path matches this format."""
        ...

    def probe(self, path: str) -> FrameSpec:
        """Read the first frame; return FrameSpec valid for all source frames."""
        ...

    def read_frame(self, f: TextIO, spec: FrameSpec) -> Optional[Frame]:
        """Read and return one Frame from f, or None at EOF."""
        ...

    def read_frame_header(self, f: TextIO) -> Optional[tuple[int, Optional[int]]]:
        """Read frame header only; skip atom lines.

        Returns (natoms, timestep) or None at EOF.
        Used by merge validation (V3/V4) — never parses coordinates.
        """
        ...

    def skip_frames(self, f: TextIO, spec: FrameSpec, n: int) -> None:
        """Advance f by exactly n * spec.lines_per_frame lines.

        Uses bare readline() with no float parsing — a fast line-skip.
        For source files only (spec.lines_per_frame is constant).
        """
        ...

    def count_frames(self, path: str, spec: FrameSpec) -> int:
        """Return total frame count using line-count arithmetic.

        total_lines(path) // spec.lines_per_frame.  Source files only.
        """
        ...

    def write_frame(self, f: TextIO, frame: Frame, spec: FrameSpec) -> None:
        """Serialize one filtered frame to f."""
        ...
