"""Adapter for the dielectric-concatenated XYZ format (lammps_xyz).

Frame structure (2 header lines + natoms atom lines):
    {natoms}                        <- line 0 of frame
    Atoms. Timestep: {timestep}     <- line 1 (timestep_line = 1)
    {type} {x} {y} {z}             <- lines 2 .. 2+natoms-1

Auto-detect: line 1 = integer natoms, line 2 starts with 'Atoms. Timestep:'.

NOTE — no box metadata:
  This format carries no simulation cell information.  The output file will
  not display a bounding box in OVITO or other tools.  Only the filtered atoms
  at their original coordinates are written.  Use lammps_custom if you need
  cell context for visualization.

Reference: old_to_new_mpi.py read_xyz_frame
"""

from __future__ import annotations

import subprocess
from typing import Optional, TextIO

from formats.base import FrameSpec
from frame_filter import Frame

HEADER_LINES = 2      # natoms line + comment line
TIMESTEP_LINE = 1     # 0-indexed line within frame containing the timestep

# Atom columns: {type} {x} {y} {z}
X_COL = 1
Y_COL = 2
Z_COL = 3


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
    """Read the first frame of path and return a FrameSpec for all source frames."""
    with open(path, "r") as f:
        line1 = f.readline().strip()
        line2 = f.readline().strip()

    try:
        natoms = int(line1)
    except ValueError:
        raise ValueError(f"lammps_xyz probe: expected integer natoms on line 1, got {line1!r}")

    if not line2.startswith("Atoms. Timestep:"):
        raise ValueError(
            f"lammps_xyz probe: expected 'Atoms. Timestep: N' on line 2, got {line2!r}"
        )

    return FrameSpec(
        natoms=natoms,
        lines_per_frame=HEADER_LINES + natoms,
        header_lines=HEADER_LINES,
        x_col=X_COL,
        y_col=Y_COL,
        z_col=Z_COL,
        timestep_line=TIMESTEP_LINE,
        extra={},
    )


def _parse_timestep(comment_line: str) -> Optional[int]:
    """Extract integer timestep from 'Atoms. Timestep: N'."""
    try:
        return int(comment_line.split(":")[1].strip())
    except (IndexError, ValueError):
        return None


def read_frame(f: TextIO, spec: FrameSpec) -> Optional[Frame]:
    """Stream one frame from f.  Returns None at EOF."""
    natoms_line = f.readline()
    if natoms_line == "":
        return None

    try:
        natoms = int(natoms_line.strip())
    except ValueError:
        raise ValueError(
            f"lammps_xyz read_frame: expected integer natoms, got {natoms_line!r}"
        )

    if natoms != spec.natoms:
        raise ValueError(
            f"lammps_xyz read_frame: natoms mismatch — expected {spec.natoms}, got {natoms}. "
            "Source trajectory must have constant atom count."
        )

    comment_line = f.readline()
    if comment_line == "":
        raise ValueError("lammps_xyz read_frame: EOF after natoms line")

    timestep = _parse_timestep(comment_line.strip())
    if timestep is None:
        raise ValueError(
            f"lammps_xyz read_frame: could not parse timestep from {comment_line!r}"
        )

    header_lines = [natoms_line, comment_line]

    positions: list[tuple[float, float, float]] = []
    atom_lines: list[str] = []
    for i in range(natoms):
        atom_line = f.readline()
        if atom_line == "":
            raise ValueError(
                f"lammps_xyz read_frame: EOF while reading atom {i + 1} of {natoms}"
            )
        parts = atom_line.split()
        positions.append((float(parts[X_COL]), float(parts[Y_COL]), float(parts[Z_COL])))
        atom_lines.append(atom_line)

    return Frame(
        timestep=timestep,
        natoms=natoms,
        positions=positions,
        atom_lines=atom_lines,
        header_lines=header_lines,
        extras={},
    )


def read_frame_header(f: TextIO) -> Optional[tuple[int, Optional[int]]]:
    """Read frame header and skip atom lines.  Returns (natoms, timestep) or None at EOF.

    Never parses atom coordinates.  Used by merge validation (V3/V4).
    """
    natoms_line = f.readline()
    if natoms_line == "":
        return None

    try:
        natoms = int(natoms_line.strip())
    except ValueError:
        raise ValueError(
            f"lammps_xyz read_frame_header: expected integer natoms, got {natoms_line!r}"
        )

    comment_line = f.readline()
    if comment_line == "":
        raise ValueError("lammps_xyz read_frame_header: EOF after natoms line")

    timestep = _parse_timestep(comment_line.strip())

    for i in range(natoms):
        if f.readline() == "":
            raise ValueError(
                f"lammps_xyz read_frame_header: EOF while skipping atom {i + 1} of {natoms}"
            )

    return natoms, timestep


def skip_frames(f: TextIO, spec: FrameSpec, n: int) -> None:
    """Skip n complete frames using fast readline() — no float parsing."""
    total_lines = n * spec.lines_per_frame
    for _ in range(total_lines):
        if f.readline() == "":
            raise ValueError(
                f"lammps_xyz skip_frames: unexpected EOF while skipping {n} frames "
                f"({total_lines} lines)"
            )


def count_frames(path: str, spec: FrameSpec) -> int:
    """Return frame count via line-count arithmetic.  Source files only."""
    try:
        result = subprocess.run(
            ["wc", "-l", path],
            capture_output=True, text=True, check=True,
        )
        total_lines = int(result.stdout.split()[0])
    except Exception:
        with open(path, "r") as f:
            total_lines = sum(1 for _ in f)
    return total_lines // spec.lines_per_frame


def write_frame(f: TextIO, frame: Frame, spec: FrameSpec) -> None:
    """Write one filtered frame.

    Comment line is written verbatim (timestep preserved).
    Only the natoms count line is updated to the filtered count.
    """
    f.write(f"{frame.natoms}\n")
    # comment line (header_lines[1]) passed through unchanged
    f.write(frame.header_lines[1])
    for atom_line in frame.atom_lines:
        f.write(atom_line)
