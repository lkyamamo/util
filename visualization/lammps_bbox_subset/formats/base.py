"""Adapter protocol and FrameSpec dataclass.

Each format adapter module must implement:
  detect(path)                                   -> bool
  probe(path)                                    -> FrameSpec
  build_index(path, spec)                        -> list[int]
  count_frames(path, spec)                       -> int
  skip_frames(f, spec, n)                        -> None
  read_frame_header(f)                           -> (natoms, timestep) | None
  stream_filter_frame(f_in, f_out, spec, bbox, batch_size) -> (n_kept, timestep) | None
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class FrameSpec:
    """Structural constants probed from the first frame of a SOURCE trajectory.

    Valid only for source files where natoms is constant across all frames.
    Do NOT reuse for part files, which have variable filtered natoms per frame.
    """

    natoms: int
    lines_per_frame: int   # header_lines + natoms — constant for source
    header_lines: int      # lines before atom data in each frame
    x_col: int             # 0-indexed token position in each atom line
    y_col: int
    z_col: int
    timestep_line: int     # 0-indexed line within a frame holding the timestep value
    extra: dict = field(default_factory=dict)
