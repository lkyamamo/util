"""Adapter for LAMMPS custom dump format.

Frame structure (9 header lines + natoms atom lines = lines_per_frame):
    ITEM: TIMESTEP          <- line 0 of frame
    {timestep}              <- line 1 (timestep_line = 1)
    ITEM: NUMBER OF ATOMS   <- line 2
    {natoms}                <- line 3
    ITEM: BOX BOUNDS ...    <- line 4
    {xlo xhi}               <- line 5
    {ylo yhi}               <- line 6
    {zlo zhi}               <- line 7
    ITEM: ATOMS {col ...}   <- line 8
    {atom lines ...}        <- lines 9 .. 9+natoms-1

Auto-detect: first non-empty line starts with 'ITEM: TIMESTEP'.

Filtering semantics for lammps_custom:
  - Original ITEM: BOX BOUNDS header + 3 box lines are ALWAYS written unchanged.
  - Atom coordinates are NEVER shifted or rescaled.
  - Only ITEM: NUMBER OF ATOMS count and the atom lines change.
  - Visual result in OVITO: full simulation cell with a floating sliver of
    atoms at their original lab-frame positions.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path
from typing import Optional, TextIO

from formats.base import FrameSpec
from frame_filter import Frame

HEADER_LINES = 9       # lines before atom data
TIMESTEP_LINE = 1      # 0-indexed line within frame containing the timestep value

try:
    import numpy as np
    _HAS_NUMPY = True
except ImportError:
    _HAS_NUMPY = False


def detect(path: str) -> bool:
    with open(path, "r") as f:
        for line in f:
            stripped = line.strip()
            if stripped:
                return stripped.startswith("ITEM: TIMESTEP")
    return False


def _parse_atoms_header(header: str) -> tuple[int, int, int, list[str]]:
    """Parse 'ITEM: ATOMS col1 col2 ...' and return (x_col, y_col, z_col, columns)."""
    if not header.startswith("ITEM: ATOMS"):
        raise ValueError(f"Expected 'ITEM: ATOMS ...', got {header!r}")
    columns = header.split()[2:]
    for name in ("x", "y", "z"):
        if name not in columns:
            raise ValueError(
                f"lammps_custom: 'ITEM: ATOMS' header missing required column '{name}'. "
                f"Got columns: {columns}"
            )
    x_col = columns.index("x")
    y_col = columns.index("y")
    z_col = columns.index("z")
    return x_col, y_col, z_col, columns


def probe(path: str) -> FrameSpec:
    """Read the first frame of path and return a FrameSpec for all source frames."""
    with open(path, "r") as f:
        # line 0
        line = f.readline()
        if not line or not line.strip().startswith("ITEM: TIMESTEP"):
            raise ValueError(f"lammps_custom probe: expected 'ITEM: TIMESTEP', got {line!r}")
        # line 1 — timestep
        f.readline()
        # line 2
        if f.readline().strip() != "ITEM: NUMBER OF ATOMS":
            raise ValueError("lammps_custom probe: expected 'ITEM: NUMBER OF ATOMS'")
        # line 3 — natoms
        natoms = int(f.readline().strip())
        # line 4 — ITEM: BOX BOUNDS
        bounds_hdr = f.readline().strip()
        if not bounds_hdr.startswith("ITEM: BOX BOUNDS"):
            raise ValueError(f"lammps_custom probe: expected 'ITEM: BOX BOUNDS', got {bounds_hdr!r}")
        # lines 5-7 — box lines
        f.readline(); f.readline(); f.readline()
        # line 8 — ITEM: ATOMS ...
        atoms_hdr = f.readline().strip()
        x_col, y_col, z_col, columns = _parse_atoms_header(atoms_hdr)

    return FrameSpec(
        natoms=natoms,
        lines_per_frame=HEADER_LINES + natoms,
        header_lines=HEADER_LINES,
        x_col=x_col,
        y_col=y_col,
        z_col=z_col,
        timestep_line=TIMESTEP_LINE,
        extra={"atoms_header": atoms_hdr, "columns": columns},
    )


def build_index(path: str, spec: FrameSpec) -> list[int]:
    """Return byte offset of each frame start in path.

    Uses 'grep -b ^ITEM: TIMESTEP' (C-level SIMD search, ~1-2 GB/s) when
    available.  Falls back to a Python binary readline scan if grep is absent
    or returns a non-zero exit code.

    The index is built once and shared by all Slurm workers, replacing the
    O(i * lines_per_frame) readline skip each worker previously had to perform
    before reaching its start frame.
    """
    try:
        result = subprocess.run(
            ["grep", "-b", "^ITEM: TIMESTEP", path],
            capture_output=True, text=True,
        )
        if result.returncode in (0, 1):  # 1 = no matches (empty file)
            offsets = []
            for line in result.stdout.splitlines():
                if line:
                    offsets.append(int(line.split(":", 1)[0]))
            return offsets
    except FileNotFoundError:
        pass  # grep not available; fall through to Python scan

    # Python fallback: binary readline scan, O(total_lines)
    offsets = []
    with open(path, "rb") as f:
        while True:
            pos = f.tell()
            first = f.readline()
            if not first:
                break
            offsets.append(pos)
            for _ in range(spec.lines_per_frame - 1):
                if f.readline() == b"":
                    return offsets
    return offsets


def read_frame(f: TextIO, spec: FrameSpec) -> Optional[Frame]:
    """Stream one frame from f.  Returns None at EOF."""
    line0 = f.readline()
    if line0 == "":
        return None
    if not line0.strip().startswith("ITEM: TIMESTEP"):
        raise ValueError(f"lammps_custom read_frame: expected 'ITEM: TIMESTEP', got {line0!r}")

    timestep = int(f.readline().strip())

    line2 = f.readline().strip()
    if line2 != "ITEM: NUMBER OF ATOMS":
        raise ValueError(f"lammps_custom read_frame: expected 'ITEM: NUMBER OF ATOMS', got {line2!r}")

    natoms_line = f.readline()
    natoms = int(natoms_line.strip())
    if natoms != spec.natoms:
        raise ValueError(
            f"lammps_custom read_frame: natoms mismatch — expected {spec.natoms}, got {natoms}. "
            "Source trajectory must have constant atom count."
        )

    bounds_hdr = f.readline()
    if not bounds_hdr.strip().startswith("ITEM: BOX BOUNDS"):
        raise ValueError(f"lammps_custom read_frame: expected 'ITEM: BOX BOUNDS', got {bounds_hdr!r}")

    box_line0 = f.readline()
    box_line1 = f.readline()
    box_line2 = f.readline()
    if "" in (box_line0, box_line1, box_line2):
        raise ValueError("lammps_custom read_frame: unexpected EOF in box bounds")

    atoms_hdr_line = f.readline()
    if not atoms_hdr_line.strip().startswith("ITEM: ATOMS"):
        raise ValueError(f"lammps_custom read_frame: expected 'ITEM: ATOMS', got {atoms_hdr_line!r}")

    # header_lines stores all 9 preamble lines verbatim for write_frame
    header_lines = [
        line0,
        f"{timestep}\n",
        "ITEM: NUMBER OF ATOMS\n",
        natoms_line,
        bounds_hdr,
        box_line0,
        box_line1,
        box_line2,
        atoms_hdr_line,
    ]

    # Read all atom lines upfront
    atom_lines: list[str] = []
    for _ in range(natoms):
        atom_line = f.readline()
        if atom_line == "":
            raise ValueError(
                f"lammps_custom read_frame: unexpected EOF while reading atoms "
                f"(expected {natoms})"
            )
        atom_lines.append(atom_line)

    # Extract positions — numpy batch parse is ~5x faster than per-line float()
    # for large natoms; falls back to pure Python if numpy is unavailable.
    x_col, y_col, z_col = spec.x_col, spec.y_col, spec.z_col
    if _HAS_NUMPY:
        tokens = [line.split() for line in atom_lines]
        arr = np.array(tokens, dtype=object)[:, [x_col, y_col, z_col]].astype(np.float64)
        positions = list(map(tuple, arr.tolist()))
    else:
        positions = []
        for line in atom_lines:
            parts = line.split()
            positions.append((float(parts[x_col]), float(parts[y_col]), float(parts[z_col])))

    return Frame(
        timestep=timestep,
        natoms=natoms,
        positions=positions,
        atom_lines=atom_lines,
        header_lines=header_lines,
        extras=spec.extra,
    )


def read_frame_header(f: TextIO) -> Optional[tuple[int, Optional[int]]]:
    """Read frame header and skip atom lines.  Returns (natoms, timestep) or None at EOF.

    Never parses atom coordinates.  Used by merge validation (V3/V4).
    """
    line0 = f.readline()
    if line0 == "":
        return None
    if not line0.strip().startswith("ITEM: TIMESTEP"):
        raise ValueError(f"lammps_custom read_frame_header: expected 'ITEM: TIMESTEP', got {line0!r}")

    timestep = int(f.readline().strip())

    if f.readline().strip() != "ITEM: NUMBER OF ATOMS":
        raise ValueError("lammps_custom read_frame_header: expected 'ITEM: NUMBER OF ATOMS'")

    natoms = int(f.readline().strip())

    # skip ITEM: BOX BOUNDS + 3 box lines + ITEM: ATOMS header
    for _ in range(5):
        f.readline()

    # skip atom lines without parsing
    for _ in range(natoms):
        if f.readline() == "":
            raise ValueError(
                f"lammps_custom read_frame_header: EOF while skipping {natoms} atom lines"
            )

    return natoms, timestep


def skip_frames(f: TextIO, spec: FrameSpec, n: int) -> None:
    """Skip n complete frames using fast readline() — no float parsing.

    Prefer the byte-offset index + f.seek() approach when available
    (see lammps_bbox_subset.py cmd_subset).  This fallback is retained
    for single-process runs without a pre-built index.
    """
    total_lines = n * spec.lines_per_frame
    for _ in range(total_lines):
        if f.readline() == "":
            raise ValueError(
                f"lammps_custom skip_frames: unexpected EOF while skipping {n} frames "
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
        # fallback: count in Python
        with open(path, "r") as f:
            total_lines = sum(1 for _ in f)
    return total_lines // spec.lines_per_frame


def write_frame(f: TextIO, frame: Frame, spec: FrameSpec) -> None:
    """Write one filtered frame.

    Box bounds and all other header lines are passed through verbatim.
    Only ITEM: NUMBER OF ATOMS count (header_lines[3]) is rewritten.
    """
    # Rewrite the natoms line (index 3) to the filtered count
    for i, hline in enumerate(frame.header_lines):
        if i == 3:
            f.write(f"{frame.natoms}\n")
        else:
            f.write(hline)
    for atom_line in frame.atom_lines:
        f.write(atom_line)
