#!/usr/bin/env python3
"""Count surface silanols and export analysis artifacts with OVITO.

Usage:
    python silanol.py (--local-path PATH | --remote-path PATH) [options]

Primary options:
    --local-path PATH            Analyze a local trajectory file.
    --remote-path PATH           Analyze a remote trajectory file via SSH/SFTP.
    --frames START END           Optional frame window [START, END).
    --surface-axis {x,y,z}       Axis used to detect the Si-gap interfaces.
    --surface-thickness FLOAT    Symmetric +/- window (Angstrom) around each
                                 Si-gap interface for surface filtering.
    --no-surface-filter          Disable interface filtering and count all
                                 silanol candidates.
    --type-si/--type-o/--type-h  Optional explicit particle type IDs.
    --require-h / --no-require-h Require at least one bonded H on silanol O.

Output behavior:
    - Creates output directory next to trajectory source:
      local:  <local parent>/output
      remote: <remote parent>/output
    - Writes:
      boundary_si_ids.csv
      silanol_ids_per_frame.jsonl
      silanol_expression_selection.txt
      silanol_counts_per_frame.csv
      silanol_statistics.txt
"""

from __future__ import annotations

import argparse
import math
import tempfile
from dataclasses import dataclass, field
from pathlib import Path, PurePosixPath
from typing import Dict, List, Optional, Tuple

import paramiko
from ovito.io import import_file
from ovito.modifiers import CreateBondsModifier


# ---------------------------------------------------------------------------
# Statistics container
# ---------------------------------------------------------------------------

@dataclass
class SilanoilStatistics:
    """Summary statistics over all analysed frames."""

    counts_per_frame: List[int] = field(default_factory=list)

    # Populated by finalize().
    n_frames: int = 0
    mean: float = 0.0
    median: float = 0.0
    std: float = 0.0
    sem: float = 0.0          # standard error of the mean
    variance: float = 0.0
    minimum: int = 0
    maximum: int = 0
    count_range: int = 0      # max - min
    q1: float = 0.0
    q3: float = 0.0
    iqr: float = 0.0

    def finalize(self) -> None:
        """Compute all derived statistics from counts_per_frame."""
        n = len(self.counts_per_frame)
        if n == 0:
            return

        self.n_frames = n
        self.minimum = min(self.counts_per_frame)
        self.maximum = max(self.counts_per_frame)
        self.count_range = self.maximum - self.minimum

        self.mean = sum(self.counts_per_frame) / n
        self.variance = sum((x - self.mean) ** 2 for x in self.counts_per_frame) / n
        self.std = math.sqrt(self.variance)
        self.sem = self.std / math.sqrt(n) if n > 1 else 0.0

        sorted_counts = sorted(self.counts_per_frame)
        self.median = _percentile(sorted_counts, 50)
        self.q1 = _percentile(sorted_counts, 25)
        self.q3 = _percentile(sorted_counts, 75)
        self.iqr = self.q3 - self.q1

    def as_text(self) -> str:
        lines = [
            "Silanol count statistics",
            "========================",
            f"  Frames analysed : {self.n_frames}",
            f"  Mean            : {self.mean:.4f}",
            f"  Median          : {self.median:.4f}",
            f"  Std dev         : {self.std:.4f}",
            f"  Variance        : {self.variance:.4f}",
            f"  Std error       : {self.sem:.4f}",
            f"  Minimum         : {self.minimum}",
            f"  Maximum         : {self.maximum}",
            f"  Range           : {self.count_range}",
            f"  Q1 (25th pct)   : {self.q1:.4f}",
            f"  Q3 (75th pct)   : {self.q3:.4f}",
            f"  IQR             : {self.iqr:.4f}",
        ]
        return "\n".join(lines) + "\n"

    def as_dict(self) -> Dict[str, float]:
        return {
            "n_frames": self.n_frames,
            "mean": self.mean,
            "median": self.median,
            "std": self.std,
            "variance": self.variance,
            "sem": self.sem,
            "minimum": self.minimum,
            "maximum": self.maximum,
            "range": self.count_range,
            "q1": self.q1,
            "q3": self.q3,
            "iqr": self.iqr,
        }


def _percentile(sorted_data: List[int], pct: float) -> float:
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
    """Ensure the pipeline has a CreateBondsModifier with desired cutoffs."""
    cb = next((m for m in pipeline.modifiers if isinstance(m, CreateBondsModifier)), None)
    if cb is None:
        cb = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)
        pipeline.modifiers.append(cb)

    if type_si is None or type_o is None or type_h is None:
        data0 = pipeline.compute(0)
        available = [(int(t.id), str(t.name)) for t in data0.particles.particle_types.types]

        def _norm(s: str) -> str:
            return "".join(ch for ch in s.lower() if ch.isalpha())

        norm_to_id: dict[str, int] = {}
        for tid, tname in available:
            n = _norm(tname)
            if n:
                norm_to_id.setdefault(n, tid)

        if type_si is None:
            type_si = norm_to_id.get("si")
        if type_o is None:
            type_o = norm_to_id.get("o")
        if type_h is None:
            type_h = norm_to_id.get("h")

        if type_si is None or type_o is None or type_h is None:
            raise ValueError(
                "Could not infer particle type IDs for Si/O/H from the file. "
                "Pass them explicitly with --type-si/--type-o/--type-h. "
                f"Available types: {available}"
            )

    t_si = int(type_si)
    t_o = int(type_o)
    t_h = int(type_h)
    cb.set_pairwise_cutoff(t_si, t_o, si_o_cutoff)
    cb.set_pairwise_cutoff(t_o, t_h, o_h_cutoff)
    return t_si, t_o, t_h


# ---------------------------------------------------------------------------
# PBC helpers
# ---------------------------------------------------------------------------

def _mic_dist(a: float, b: float, box_length: float) -> float:
    """Minimum-image distance between two coordinates along a single periodic axis."""
    d = abs(a - b)
    return min(d, box_length - d)


def _si_gap_boundaries(
    si_axis_values: List[float],
    si_ids: List[int],
    *,
    axis: str,
    box_length: float,
) -> Tuple[int, int, float, float]:
    """Return boundary Si IDs and coordinates for the largest Si-Si gap along axis.

    Considers both interior gaps and the periodic wraparound gap so that slabs
    whose vacuum region straddles the box boundary are handled correctly.
    """
    if len(si_axis_values) < 2:
        raise ValueError(
            f"Need at least 2 Si atoms to detect a gap along axis {axis!r}; "
            f"found {int(len(si_axis_values))}."
        )

    paired = sorted(zip(si_axis_values, si_ids), key=lambda t: t[0])

    max_gap = -1.0
    lower_coord, lower_id = paired[0]
    upper_coord, upper_id = paired[1]

    for i in range(len(paired) - 1):
        gap = paired[i + 1][0] - paired[i][0]
        if gap > max_gap:
            max_gap = gap
            lower_coord, lower_id = paired[i]
            upper_coord, upper_id = paired[i + 1]

    # Wraparound gap: vacuum straddles the periodic boundary.
    wrap_gap = (paired[0][0] + box_length) - paired[-1][0]
    if wrap_gap > max_gap:
        lower_coord, lower_id = paired[-1]
        upper_coord, upper_id = paired[0]

    return int(lower_id), int(upper_id), float(lower_coord), float(upper_coord)


# ---------------------------------------------------------------------------
# Core analysis
# ---------------------------------------------------------------------------

def compute_silanols(
    pipeline,
    *,
    frames: Optional[Tuple[int, int]] = None,
    si_o_cutoff: float = 1.9,
    o_h_cutoff: float = 1.2,
    allow_other_neighbors: bool = True,
    return_identifiers_if_available: bool = True,
    type_si: Optional[int] = None,
    type_o: Optional[int] = None,
    type_h: Optional[int] = None,
    surface_axis: str = "x",
    apply_surface_filter: bool = True,
    surface_thickness_angstrom: float = 5.0,
    require_at_least_one_h: bool = True,
) -> Tuple[SilanoilStatistics, List[List[int]], List[Tuple[int, int]]]:
    """Return (statistics, ids_per_frame, boundary_si_ids_per_frame).

    ids are OVITO particle identifiers if present, otherwise particle indices (0-based).
    """
    t_Si, t_O, t_H = _ensure_bonds(
        pipeline,
        si_o_cutoff=si_o_cutoff,
        o_h_cutoff=o_h_cutoff,
        type_si=type_si,
        type_o=type_o,
        type_h=type_h,
    )

    stats = SilanoilStatistics()
    ids: List[List[int]] = []
    boundary_si_ids_per_frame: List[Tuple[int, int]] = []

    start_frame = frames[0] if frames else 0
    end_frame = frames[1] if frames else len(pipeline.frames)

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

        box_length = float(data.cell[ax, ax])

        si_indices = [i for i in range(data.particles.count) if ptype[i] == t_Si]
        if len(si_indices) == 0:
            raise ValueError("No Si particles found; cannot define surface by Si-gap interfaces.")

        si_axis_values = [float(pos[i, ax]) for i in si_indices]
        has_pid = "Particle Identifier" in data.particles
        if has_pid:
            pid_arr = data.particles["Particle Identifier"].array
            si_ids = [int(pid_arr[i]) for i in si_indices]
        else:
            si_ids = [int(i) for i in si_indices]

        boundary_lo_id, boundary_hi_id, si_lo, si_hi = _si_gap_boundaries(
            si_axis_values,
            si_ids,
            axis=axis,
            box_length=box_length,
        )
        boundary_si_ids_per_frame.append((boundary_lo_id, boundary_hi_id))

        neighbors: List[List[int]] = [[] for _ in range(data.particles.count)]
        for i, j in bonds.topology:
            neighbors[i].append(j)
            neighbors[j].append(i)

        matched_indices: List[int] = []
        for i in range(data.particles.count):
            if ptype[i] != t_O:
                continue

            nb = neighbors[i]
            n_si = sum(ptype[j] == t_Si for j in nb)
            if n_si != 1:
                continue

            if require_at_least_one_h:
                n_h = sum(ptype[j] == t_H for j in nb)
                if n_h < 1:
                    continue

            if not allow_other_neighbors:
                if any((ptype[j] != t_Si and ptype[j] != t_H) for j in nb):
                    continue

            matched_indices.append(i)

        if apply_surface_filter:
            if surface_thickness_angstrom <= 0:
                matched_indices = []
            else:
                thickness = float(surface_thickness_angstrom)
                matched_indices = [
                    i
                    for i in matched_indices
                    if (
                        _mic_dist(float(pos[i, ax]), si_lo, box_length) <= thickness
                        or _mic_dist(float(pos[i, ax]), si_hi, box_length) <= thickness
                    )
                ]

        if return_identifiers_if_available and has_pid:
            matched_ids = [int(pid_arr[i]) for i in matched_indices]
        else:
            matched_ids = matched_indices

        # Diagnostics on the first frame only.
        if frame == start_frame:
            import numpy as np
            print(f"\n=== DIAGNOSTICS for frame {frame} ===")
            print(f"  Particle Identifier property present: {has_pid}")
            if has_pid:
                unique_pids = len(set(int(x) for x in pid_arr))
                print(f"  Total particles: {data.particles.count}, Unique PIDs: {unique_pids}")
                if unique_pids < data.particles.count:
                    print("  WARNING: Duplicate Particle Identifiers detected!")
            print(f"  Type IDs -> Si={t_Si}, O={t_O}, H={t_H}")
            print(f"  Box length along {axis!r}: {box_length:.4f} Å")
            print(f"  Boundary Si coords: si_lo={si_lo:.4f}, si_hi={si_hi:.4f}")
            print(f"  Selected {len(matched_indices)} silanols (indices) -> {len(matched_ids)} IDs")
            for idx in matched_indices[:20]:
                nb = neighbors[idx]
                n_si_nb = sum(ptype[j] == t_Si for j in nb)
                n_h_nb  = sum(ptype[j] == t_H  for j in nb)
                n_o_nb  = sum(ptype[j] == t_O  for j in nb)
                pid_val = int(pid_arr[idx]) if has_pid else idx
                si_dists = [
                    float(np.linalg.norm(
                        data.particles.positions[j] - data.particles.positions[idx]
                    ))
                    for j in nb if ptype[j] == t_Si
                ]
                print(
                    f"  idx={idx} PID={pid_val} type={int(ptype[idx])} "
                    f"pos=({pos[idx,0]:.2f},{pos[idx,1]:.2f},{pos[idx,2]:.2f}) "
                    f"n_Si={n_si_nb} n_H={n_h_nb} n_O={n_o_nb} Si_dists={si_dists}"
                )
            print("=== END DIAGNOSTICS ===\n")

        stats.counts_per_frame.append(len(matched_ids))
        ids.append(matched_ids)

    stats.finalize()
    return stats, ids, boundary_si_ids_per_frame


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


def _counts_per_frame_csv(stats: SilanoilStatistics) -> str:
    lines = ["frame,silanol_count"]
    lines.extend(
        f"{frame},{count}"
        for frame, count in enumerate(stats.counts_per_frame)
    )
    return "\n".join(lines) + "\n"


def _silanol_ids_jsonl(ids_per_frame: List[List[int]]) -> str:
    lines = [
        f'{{"frame": {int(frame)}, "silanol_ids": [{", ".join(str(int(pid)) for pid in frame_ids)}]}}'
        for frame, frame_ids in enumerate(ids_per_frame)
    ]
    return "\n".join(lines) + "\n"


def _ovito_expression_for_particle_identifiers(frame_ids: List[int]) -> str:
    if not frame_ids:
        return "(ParticleIdentifier == -1) && (ParticleIdentifier == -2)"
    return " || ".join(f"(ParticleIdentifier == {int(pid)})" for pid in frame_ids)


def _silanol_expression_selection_text(ids_per_frame: List[List[int]]) -> str:
    lines: List[str] = []
    for frame, frame_ids in enumerate(ids_per_frame):
        lines.append(f"# frame {frame}")
        lines.append(_ovito_expression_for_particle_identifiers(frame_ids))
        lines.append("")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def _parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute silanol count and RDF (OVITO).")
    p.add_argument("--base-dir", type=Path, default=None,
                   help="Deprecated: no longer used for output location.")
    p.add_argument("--host", default="endeavour.usc.edu", help="SSH hostname.")
    p.add_argument("--username", default="lkyamamo", help="SSH username.")
    p.add_argument("--type-si", type=int, default=None, help="Particle type ID for Si.")
    p.add_argument("--type-o",  type=int, default=None, help="Particle type ID for O.")
    p.add_argument("--type-h",  type=int, default=None, help="Particle type ID for H.")
    p.add_argument("--require-h", action=argparse.BooleanOptionalAction, default=True,
                   help="Require at least one H bonded to the silanol O (Si–O–H).")
    p.add_argument("--surface-axis", default="x", choices=["x", "y", "z"],
                   help="Surface normal axis used to identify top/bottom surfaces.")
    p.add_argument("--surface-thickness", type=float, default=5.0,
                   help="Symmetric half-width in Å around each Si-gap interface.")
    p.add_argument("--no-surface-filter", action="store_true",
                   help="Count all silanols (do not restrict to surface slabs).")

    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument("--remote-path",
                        help="Remote trajectory path; creates '<remote parent>/output'.")
    source.add_argument("--local-path", type=Path,
                        help="Local trajectory path; creates '<local parent>/output'.")
    p.add_argument("--frames", nargs=2, type=int, default=None,
                   metavar=("START", "END"),
                   help="Analyse only frames [START, END) of a LAMMPS dump.")
    return p.parse_args(argv)


# ---------------------------------------------------------------------------
# Remote I/O helpers
# ---------------------------------------------------------------------------

def _download_remote_file_to_temp(*, host: str, username: str, remote_path: str) -> Path:
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(host, username=username)
    sftp = ssh.open_sftp()
    try:
        with tempfile.NamedTemporaryFile(
            prefix="silanol_input_", suffix=".data", delete=False
        ) as tmp:
            local_path = Path(tmp.name)
            with sftp.open(remote_path, "rb") as src:
                while True:
                    chunk = src.read(8 * 1024 * 1024)
                    if not chunk:
                        break
                    tmp.write(chunk)
    finally:
        sftp.close()
        ssh.close()
    return local_path


def _download_lammps_dump_frames_to_temp(
    *,
    host: str,
    username: str,
    remote_path: str,
    frames: Tuple[int, int],
) -> Path:
    start, end = frames
    if start < 0 or end <= start:
        raise ValueError(f"Invalid frames range: {frames}. Expected START >= 0 and END > START.")

    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(host, username=username)
    sftp = ssh.open_sftp()
    try:
        with tempfile.NamedTemporaryFile(
            prefix="silanol_input_", suffix=".dump", delete=False, mode="w"
        ) as tmp:
            local_path = Path(tmp.name)

            with sftp.open(remote_path, "r") as src:
                first = src.readline()
                if not first.startswith("ITEM: TIMESTEP"):
                    raise ValueError(
                        "Expected LAMMPS dump format beginning with 'ITEM: TIMESTEP'. "
                        f"First line was: {first[:80]!r}"
                    )
                lines_per_frame = 1
                while True:
                    line = src.readline()
                    if line == "":
                        raise ValueError(
                            "Reached EOF before finding the second 'ITEM: TIMESTEP' marker."
                        )
                    if line.startswith("ITEM: TIMESTEP"):
                        break
                    lines_per_frame += 1

            skip_lines = start * lines_per_frame
            copy_lines = (end - start) * lines_per_frame

            with sftp.open(remote_path, "r") as src:
                for _ in range(skip_lines):
                    if src.readline() == "":
                        raise ValueError(
                            f"Reached EOF while skipping to START frame {start}. "
                            f"(lines_per_frame={lines_per_frame})"
                        )
                for _ in range(copy_lines):
                    line = src.readline()
                    if line == "":
                        raise ValueError(
                            f"Reached EOF while copying frame range {frames}. "
                            f"(lines_per_frame={lines_per_frame})"
                        )
                    tmp.write(line)
    finally:
        sftp.close()
        ssh.close()
    return local_path


def _ensure_remote_output_dir(*, host: str, username: str, remote_path: str) -> str:
    rpath = PurePosixPath(remote_path)
    parent = str(rpath.parent) if str(rpath.parent) else "."
    remote_output = str(PurePosixPath(parent) / "output")

    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(host, username=username)
    sftp = ssh.open_sftp()
    try:
        current = ""
        for part in PurePosixPath(remote_output).parts:
            if part == "/":
                current = "/"
                continue
            current = f"{current}/{part}" if current not in ("", "/") else f"{current}{part}"
            try:
                sftp.stat(current)
            except OSError:
                sftp.mkdir(current)
    finally:
        sftp.close()
        ssh.close()
    return remote_output


def _write_remote_text_file(
    *, host: str, username: str, remote_file_path: str, contents: str
) -> None:
    ssh = paramiko.SSHClient()
    ssh.set_missing_host_key_policy(paramiko.AutoAddPolicy())
    ssh.connect(host, username=username)
    sftp = ssh.open_sftp()
    try:
        with sftp.open(remote_file_path, "w") as f:
            f.write(contents)
    finally:
        sftp.close()
        ssh.close()


# ---------------------------------------------------------------------------
# Output writing (local / remote)
# ---------------------------------------------------------------------------

def _write_outputs_local(
    output_dir: Path,
    *,
    stats: SilanoilStatistics,
    ids: List[List[int]],
    boundary_si_ids_per_frame: List[Tuple[int, int]],
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "boundary_si_ids.csv").write_text(_boundary_si_csv(boundary_si_ids_per_frame))
    (output_dir / "silanol_counts_per_frame.csv").write_text(_counts_per_frame_csv(stats))
    (output_dir / "silanol_ids_per_frame.jsonl").write_text(_silanol_ids_jsonl(ids))
    (output_dir / "silanol_expression_selection.txt").write_text(
        _silanol_expression_selection_text(ids)
    )
    (output_dir / "silanol_statistics.txt").write_text(stats.as_text())
    print(f"output_dir={output_dir}")


def _write_outputs_remote(
    remote_output_dir: str,
    *,
    host: str,
    username: str,
    stats: SilanoilStatistics,
    ids: List[List[int]],
    boundary_si_ids_per_frame: List[Tuple[int, int]],
) -> None:
    base = PurePosixPath(remote_output_dir)
    files = {
        "boundary_si_ids.csv":              _boundary_si_csv(boundary_si_ids_per_frame),
        "silanol_counts_per_frame.csv":     _counts_per_frame_csv(stats),
        "silanol_ids_per_frame.jsonl":      _silanol_ids_jsonl(ids),
        "silanol_expression_selection.txt": _silanol_expression_selection_text(ids),
        "silanol_statistics.txt":           stats.as_text(),
    }
    for filename, contents in files.items():
        _write_remote_text_file(
            host=host,
            username=username,
            remote_file_path=str(base / filename),
            contents=contents,
        )
    print(f"remote_output_dir={remote_output_dir}")


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(argv: Optional[List[str]] = None) -> None:
    args = _parse_args(argv)

    cleanup_local_input = False
    analysis_frames: Optional[Tuple[int, int]] = None
    output_dir: Optional[Path] = None
    remote_output_dir: Optional[str] = None

    if args.local_path is not None:
        local_input_path = args.local_path.expanduser().resolve()
        output_dir = local_input_path.parent / "output"
        analysis_frames = (
            (int(args.frames[0]), int(args.frames[1])) if args.frames is not None else None
        )
    else:
        cleanup_local_input = True
        remote_output_dir = _ensure_remote_output_dir(
            host=args.host,
            username=args.username,
            remote_path=args.remote_path,
        )
        if args.frames is None:
            local_input_path = _download_remote_file_to_temp(
                host=args.host,
                username=args.username,
                remote_path=args.remote_path,
            )
            analysis_frames = None
        else:
            requested = (int(args.frames[0]), int(args.frames[1]))
            local_input_path = _download_lammps_dump_frames_to_temp(
                host=args.host,
                username=args.username,
                remote_path=args.remote_path,
                frames=requested,
            )
            analysis_frames = (0, requested[1] - requested[0])

    try:
        pipeline = import_file(str(local_input_path))

        stats, ids, boundary_si_ids_per_frame = compute_silanols(
            pipeline,
            frames=analysis_frames,
            type_si=args.type_si,
            type_o=args.type_o,
            type_h=args.type_h,
            surface_axis=args.surface_axis,
            apply_surface_filter=not args.no_surface_filter,
            surface_thickness_angstrom=args.surface_thickness,
            require_at_least_one_h=args.require_h,
        )

        print(stats.as_text())
        print(f"boundary_si_ids[:20]={boundary_si_ids_per_frame[:20]}")

        if output_dir is not None:
            _write_outputs_local(
                output_dir,
                stats=stats,
                ids=ids,
                boundary_si_ids_per_frame=boundary_si_ids_per_frame,
            )
        elif remote_output_dir is not None:
            _write_outputs_remote(
                remote_output_dir,
                host=args.host,
                username=args.username,
                stats=stats,
                ids=ids,
                boundary_si_ids_per_frame=boundary_si_ids_per_frame,
            )
    finally:
        if cleanup_local_input:
            try:
                local_input_path.unlink(missing_ok=True)
            except Exception:
                pass


if __name__ == "__main__":
    main()