"""Adapter for the dielectric-concatenated XYZ format (lammps_xyz).

Frame structure (2 header lines + natoms atom lines):
    {natoms}
    Atoms. Timestep: {timestep}
    {type} {x} {y} {z}

Auto-detect: line 1 = integer natoms, line 2 starts with 'Atoms. Timestep:'.

NOTE — no box metadata:
    This format carries no simulation cell information.  The output file will
    not display a bounding box in OVITO.  Use lammps_custom if you need cell
    context for visualization.

Requires: numpy
"""

from __future__ import annotations

from typing import Optional, TextIO

import numpy as np

from formats.base import FrameSpec

HEADER_LINES = 2
X_COL, Y_COL, Z_COL = 1, 2, 3

DEFAULT_BATCH_SIZE = 100_000


def detect(path: str) -> bool:
    with open(path, "r") as f:
        line1 = f.readline().strip()
        line2 = f.readline().strip()
    try:
        int(line1)
    except ValueError:
        return False
    return line2.startswith("Atoms. Timestep:")


def probe(path: str) -> FrameSpec:
    with open(path, "r") as f:
        line1 = f.readline().strip()
        line2 = f.readline().strip()

    try:
        natoms = int(line1)
    except ValueError:
        raise ValueError(f"lammps_xyz probe: expected integer natoms, got {line1!r}")
    if not line2.startswith("Atoms. Timestep:"):
        raise ValueError(f"lammps_xyz probe: expected 'Atoms. Timestep: N', got {line2!r}")

    return FrameSpec(
        natoms=natoms,
        lines_per_frame=HEADER_LINES + natoms,
        header_lines=HEADER_LINES,
        x_col=X_COL, y_col=Y_COL, z_col=Z_COL,
        timestep_line=1,
        extra={},
    )


def build_index(path: str, spec: FrameSpec) -> list[int]:
    """Return byte offset of each frame start.

    Delegates to build_index_from_line_count: because natoms is constant,
    lines_per_frame is fixed, and every lines_per_frame-th newline marks a
    frame boundary.  No subprocess or regex needed.
    """
    from formats.base import build_index_from_line_count
    return build_index_from_line_count(path, spec.lines_per_frame)


def count_frames(path: str, spec: FrameSpec) -> int:
    return len(build_index(path, spec))


def skip_frames(f: TextIO, spec: FrameSpec, n: int) -> None:
    for _ in range(n * spec.lines_per_frame):
        if f.readline() == "":
            raise ValueError(f"lammps_xyz skip_frames: unexpected EOF skipping {n} frames")


def read_frame_header(f: TextIO) -> Optional[tuple[int, Optional[int]]]:
    natoms_line = f.readline()
    if natoms_line == "":
        return None
    try:
        natoms = int(natoms_line.strip())
    except ValueError:
        raise ValueError(f"lammps_xyz read_frame_header: expected integer natoms, got {natoms_line!r}")

    comment = f.readline()
    if comment == "":
        raise ValueError("lammps_xyz read_frame_header: EOF after natoms line")

    try:
        timestep: Optional[int] = int(comment.split(":")[1].strip())
    except (IndexError, ValueError):
        timestep = None

    for i in range(natoms):
        if f.readline() == "":
            raise ValueError(f"lammps_xyz read_frame_header: EOF skipping atom {i+1}/{natoms}")

    return natoms, timestep


def stream_filter_frame(
    f_in: TextIO,
    f_out: Optional[TextIO],
    spec: FrameSpec,
    bbox,
    batch_size: int = DEFAULT_BATCH_SIZE,
) -> Optional[tuple[int, Optional[int]]]:
    """Read one frame, filter atoms in batches, write kept atoms.

    Returns (n_kept, timestep) or None at EOF.
    Pass f_out=None to skip (stride path).
    See lammps_custom.stream_filter_frame for memory/seek-back details.
    """
    natoms_line = f_in.readline()
    if natoms_line == "":
        return None
    try:
        natoms = int(natoms_line.strip())
    except ValueError:
        raise ValueError(f"lammps_xyz: expected integer natoms, got {natoms_line!r}")
    if natoms != spec.natoms:
        raise ValueError(
            f"lammps_xyz: natoms mismatch — expected {spec.natoms}, got {natoms}"
        )

    comment_line = f_in.readline()
    if comment_line == "":
        raise ValueError("lammps_xyz: EOF after natoms line")
    try:
        timestep: Optional[int] = int(comment_line.split(":")[1].strip())
    except (IndexError, ValueError):
        timestep = None

    if f_out is None:
        for _ in range(natoms):
            f_in.readline()
        return 0, timestep

    natoms_width = len(natoms_line.rstrip())
    count_pos = f_out.tell()
    f_out.write(natoms_line)
    f_out.write(comment_line)

    lo = np.array([bbox.xmin, bbox.ymin, bbox.zmin], dtype=np.float64)
    hi = np.array([bbox.xmax, bbox.ymax, bbox.zmax], dtype=np.float64)

    n_kept = 0
    remaining = natoms

    while remaining > 0:
        n_batch = min(batch_size, remaining)

        batch_lines: list[str] = []
        for _ in range(n_batch):
            line = f_in.readline()
            if line == "":
                raise ValueError("lammps_xyz: unexpected EOF in atom data")
            batch_lines.append(line)

        xyz = np.array(
            [ln.split() for ln in batch_lines], dtype=object
        )[:, [X_COL, Y_COL, Z_COL]].astype(np.float64)

        mask = np.all((xyz >= lo) & (xyz <= hi), axis=1)
        n_batch_kept = int(mask.sum())

        if n_batch_kept > 0:
            f_out.writelines(np.array(batch_lines, dtype=object)[mask].tolist())

        n_kept += n_batch_kept
        remaining -= n_batch

    end_pos = f_out.tell()
    f_out.seek(count_pos)
    f_out.write(f"{n_kept:>{natoms_width}d}\n")
    f_out.seek(end_pos)

    return n_kept, timestep
