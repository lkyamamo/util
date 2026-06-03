"""Adapter protocol, FrameSpec dataclass, and shared index builder.

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

import numpy as np


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


def build_index_from_line_count(path: str, lines_per_frame: int) -> list[int]:
    """Return the byte offset of each frame start using pure numpy newline counting.

    Because natoms is constant, lines_per_frame = header_lines + natoms is fixed
    for every frame.  The start of frame i is simply the byte that follows the
    (i * lines_per_frame)-th newline in the file — no regex, no subprocess.

    Algorithm (single pass, reads at memory-bandwidth speed):
      • Scan the file in 8 MB binary chunks.
      • Use numpy to locate newline bytes (\\n = 0x0A) in each chunk.
      • When line_in_frame reaches 0 (start of a new frame), record the current
        byte position as a frame offset.
      • When line_in_frame reaches lines_per_frame, reset it to 0 and advance
        frame_start to the byte after the last newline of the completed frame.
    """
    CHUNK = 1 << 23   # 8 MB — large enough to amortize Python overhead

    offsets: list[int] = []
    byte_pos: int = 0
    line_in_frame: int = 0
    frame_start: int = 0

    with open(path, "rb") as f:
        while True:
            chunk = f.read(CHUNK)
            if not chunk:
                break

            arr = np.frombuffer(chunk, dtype=np.uint8)
            nl_positions = np.where(arr == 0x0A)[0]

            for nl_pos in nl_positions:
                if line_in_frame == 0:
                    offsets.append(frame_start)
                line_in_frame += 1
                if line_in_frame == lines_per_frame:
                    frame_start = byte_pos + int(nl_pos) + 1
                    line_in_frame = 0

            byte_pos += len(chunk)

    return offsets
