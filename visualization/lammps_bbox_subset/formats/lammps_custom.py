"""Adapter for LAMMPS custom dump format.

Frame structure (9 header lines + natoms atom lines):
    ITEM: TIMESTEP
    {timestep}
    ITEM: NUMBER OF ATOMS
    {natoms}
    ITEM: BOX BOUNDS pp pp pp
    {xlo xhi}
    {ylo yhi}
    {zlo zhi}
    ITEM: ATOMS {col ...}
    {atom lines ...}

Auto-detect: first non-empty line starts with 'ITEM: TIMESTEP'.

Filtering semantics (floating-sliver):
    ITEM: BOX BOUNDS and all box lines are always written verbatim.
    Atom coordinates are never shifted or rescaled.
    Only ITEM: NUMBER OF ATOMS count and the atom lines change.
    Result in OVITO: full simulation cell with atoms at original lab positions.

Requires: numpy
"""

from __future__ import annotations

import subprocess
from typing import Optional, TextIO

import numpy as np

from formats.base import FrameSpec

HEADER_LINES = 9

DEFAULT_BATCH_SIZE = 100_000


def detect(path: str) -> bool:
    with open(path, "r") as f:
        for line in f:
            s = line.strip()
            if s:
                return s.startswith("ITEM: TIMESTEP")
    return False


def _parse_atoms_header(header: str) -> tuple[int, int, int, list[str]]:
    if not header.startswith("ITEM: ATOMS"):
        raise ValueError(f"Expected 'ITEM: ATOMS ...', got {header!r}")
    columns = header.split()[2:]
    for name in ("x", "y", "z"):
        if name not in columns:
            raise ValueError(
                f"lammps_custom: missing required column '{name}'. Got: {columns}"
            )
    return columns.index("x"), columns.index("y"), columns.index("z"), columns


def probe(path: str) -> FrameSpec:
    with open(path, "r") as f:
        line = f.readline()
        if not line.strip().startswith("ITEM: TIMESTEP"):
            raise ValueError(f"lammps_custom probe: expected 'ITEM: TIMESTEP', got {line!r}")
        f.readline()                              # timestep value
        if f.readline().strip() != "ITEM: NUMBER OF ATOMS":
            raise ValueError("lammps_custom probe: expected 'ITEM: NUMBER OF ATOMS'")
        natoms = int(f.readline().strip())
        bounds_hdr = f.readline().strip()
        if not bounds_hdr.startswith("ITEM: BOX BOUNDS"):
            raise ValueError(f"lammps_custom probe: expected 'ITEM: BOX BOUNDS', got {bounds_hdr!r}")
        f.readline(); f.readline(); f.readline()  # three box lines
        atoms_hdr = f.readline().strip()
        x_col, y_col, z_col, columns = _parse_atoms_header(atoms_hdr)

    return FrameSpec(
        natoms=natoms,
        lines_per_frame=HEADER_LINES + natoms,
        header_lines=HEADER_LINES,
        x_col=x_col, y_col=y_col, z_col=z_col,
        timestep_line=1,
        extra={"atoms_header": atoms_hdr, "columns": columns},
    )


def build_index(path: str, spec: FrameSpec) -> list[int]:
    """Return byte offset of each frame start.

    Uses grep -b (C-level SIMD, ~1-2 GB/s) when available.
    Falls back to a Python binary readline scan.
    """
    try:
        result = subprocess.run(
            ["grep", "-b", "^ITEM: TIMESTEP", path],
            capture_output=True, text=True,
        )
        if result.returncode in (0, 1):
            return [
                int(line.split(":", 1)[0])
                for line in result.stdout.splitlines()
                if line
            ]
    except FileNotFoundError:
        pass

    offsets: list[int] = []
    with open(path, "rb") as f:
        while True:
            pos = f.tell()
            if not f.readline():
                break
            offsets.append(pos)
            for _ in range(spec.lines_per_frame - 1):
                if f.readline() == b"":
                    return offsets
    return offsets


def count_frames(path: str, spec: FrameSpec) -> int:
    try:
        result = subprocess.run(
            ["wc", "-l", path], capture_output=True, text=True, check=True,
        )
        total_lines = int(result.stdout.split()[0])
    except Exception:
        with open(path, "r") as f:
            total_lines = sum(1 for _ in f)
    return total_lines // spec.lines_per_frame


def skip_frames(f: TextIO, spec: FrameSpec, n: int) -> None:
    """Skip n source frames with bare readline() — no parsing.

    Only used when no frame index is available.
    """
    for _ in range(n * spec.lines_per_frame):
        if f.readline() == "":
            raise ValueError(f"lammps_custom skip_frames: unexpected EOF skipping {n} frames")


def read_frame_header(f: TextIO) -> Optional[tuple[int, Optional[int]]]:
    """Read one frame header, skip atom lines.  Returns (natoms, timestep) or None at EOF.

    No coordinate parsing — used only by merge validation (V3/V4).
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

    for _ in range(5):          # BOX BOUNDS header + 3 box lines + ITEM: ATOMS header
        f.readline()
    for _ in range(natoms):
        if f.readline() == "":
            raise ValueError("lammps_custom read_frame_header: EOF while skipping atom lines")

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
    Pass f_out=None to skip the frame (stride path) — reads and discards
    atom lines with bare readline(), no numpy.

    Memory per batch (batch_size=100_000):
        batch_lines list  ~  4 MB
        tokens list[list] ~ 38 MB  (transient; freed after each batch)
        xyz numpy array   ~  2 MB
        ─────────────────────────
        Peak per batch    ~ 44 MB  regardless of total natoms per frame.

    Atom count placeholder:
        ITEM: NUMBER OF ATOMS is written with the source natoms as a
        placeholder.  After all batches complete, f_out is seeked back to
        that position and the actual kept count is written, right-aligned
        to the same character width so all subsequent byte offsets are
        unchanged.  Safe on Linux text-mode files with ASCII content.
    """
    # ── read header ──────────────────────────────────────────────────────────
    line0 = f_in.readline()
    if line0 == "":
        return None
    if not line0.strip().startswith("ITEM: TIMESTEP"):
        raise ValueError(f"lammps_custom: expected 'ITEM: TIMESTEP', got {line0!r}")

    timestep_line = f_in.readline()
    timestep = int(timestep_line.strip())
    item_natoms_line = f_in.readline()   # "ITEM: NUMBER OF ATOMS\n"
    natoms_line = f_in.readline()        # "{natoms}\n"
    natoms = int(natoms_line.strip())
    if natoms != spec.natoms:
        raise ValueError(
            f"lammps_custom: natoms mismatch — expected {spec.natoms}, got {natoms}"
        )
    bounds_hdr = f_in.readline()
    box0 = f_in.readline()
    box1 = f_in.readline()
    box2 = f_in.readline()
    atoms_hdr = f_in.readline()

    # ── stride skip ──────────────────────────────────────────────────────────
    if f_out is None:
        for _ in range(natoms):
            f_in.readline()
        return 0, timestep

    # ── write header with placeholder natoms ─────────────────────────────────
    f_out.write(line0)
    f_out.write(timestep_line)
    f_out.write(item_natoms_line)
    natoms_width = len(natoms_line.rstrip())   # digit count; used for same-width overwrite
    count_pos = f_out.tell()
    f_out.write(natoms_line)                   # placeholder = source natoms
    f_out.write(bounds_hdr)
    f_out.write(box0)
    f_out.write(box1)
    f_out.write(box2)
    f_out.write(atoms_hdr)

    # ── filter in batches ────────────────────────────────────────────────────
    lo = np.array([bbox.xmin, bbox.ymin, bbox.zmin], dtype=np.float64)
    hi = np.array([bbox.xmax, bbox.ymax, bbox.zmax], dtype=np.float64)
    cols = [spec.x_col, spec.y_col, spec.z_col]

    n_kept = 0
    remaining = natoms

    while remaining > 0:
        n_batch = min(batch_size, remaining)

        batch_lines: list[str] = []
        for _ in range(n_batch):
            line = f_in.readline()
            if line == "":
                raise ValueError("lammps_custom: unexpected EOF in atom data")
            batch_lines.append(line)

        xyz = np.array(
            [ln.split() for ln in batch_lines], dtype=object
        )[:, cols].astype(np.float64)

        mask = np.all((xyz >= lo) & (xyz <= hi), axis=1)
        n_batch_kept = int(mask.sum())

        if n_batch_kept > 0:
            f_out.writelines(np.array(batch_lines, dtype=object)[mask].tolist())

        n_kept += n_batch_kept
        remaining -= n_batch

    # ── overwrite placeholder with actual atom count ──────────────────────────
    end_pos = f_out.tell()
    f_out.seek(count_pos)
    f_out.write(f"{n_kept:>{natoms_width}d}\n")
    f_out.seek(end_pos)

    return n_kept, timestep
