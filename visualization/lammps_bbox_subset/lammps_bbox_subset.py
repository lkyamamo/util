#!/usr/bin/env python3
"""LAMMPS trajectory bounding-box subset tool.

Subcommands
-----------
subset        Filter atoms by Cartesian bbox; write new trajectory file.
count-frames  Print total frame count for a source trajectory (integer).
build-index   Scan source and write a byte-offset frame index (.idx).
merge         Concatenate ordered part files with V1-V7 validation.

Performance notes
-----------------
Without a frame index, each Slurm worker skips its start offset via
readline() — O(start_frame × lines_per_frame) calls.  For N workers the
total skip work scales as O(N² × CHUNK × lines_per_frame).

With a frame index the Slurm script builds it once (O(file_size), using
grep at ~1-2 GB/s for lammps_custom), then every worker seeks in O(1).

Typical workflow:
    python3 lammps_bbox_subset.py build-index INPUT --output INPUT.idx
    python3 lammps_bbox_subset.py subset INPUT --frame-index INPUT.idx ...

Assumptions (must hold for source trajectories):
  - natoms is constant across all frames (NVT/NVE; grand-canonical not supported).
  - Every bbox-filtered frame contains at least one atom.
  - lammps_custom: ITEM: BOX BOUNDS is always written unchanged (floating-sliver
    semantics — full simulation cell preserved, only atoms outside bbox deleted).
"""

from __future__ import annotations

import argparse
import glob
import os
import re
import struct
import sys
import time
from pathlib import Path
from typing import Optional

# Allow running from the lammps_bbox_subset/ directory directly.
sys.path.insert(0, str(Path(__file__).parent))

import formats as fmt_registry
from frame_filter import BBox, filter_by_bbox


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _resolve(path: str) -> Path:
    return Path(path).resolve()


def _numeric_sort_key(p: str) -> int:
    m = re.search(r"(\d+)", Path(p).stem)
    return int(m.group(1)) if m else 0


def _sorted_parts(pattern_or_list) -> list[str]:
    """Glob a pattern or sort an explicit list, both by numeric stem."""
    if isinstance(pattern_or_list, str):
        paths = glob.glob(pattern_or_list)
    else:
        paths = list(pattern_or_list)
    return sorted(paths, key=_numeric_sort_key)


def _default_output(input_path: str) -> Path:
    p = Path(input_path)
    return p.parent / f"{p.stem}_bbox{p.suffix}"


def _abort(msg: str, check_id: str = "") -> None:
    prefix = f"merge validation failed: {check_id} — " if check_id else ""
    print(f"ERROR: {prefix}{msg}", file=sys.stderr)
    sys.exit(1)


# ---------------------------------------------------------------------------
# Frame index I/O  (unsigned 64-bit little-endian byte offsets)
# ---------------------------------------------------------------------------

def _save_index(offsets: list[int], path: str) -> None:
    """Write frame byte offsets as a compact binary file."""
    with open(path, "wb") as f:
        f.write(struct.pack(f"<{len(offsets)}Q", *offsets))


def _load_index(path: str) -> list[int]:
    """Load frame byte offsets from a binary index file."""
    with open(path, "rb") as f:
        data = f.read()
    n = len(data) // 8
    return list(struct.unpack(f"<{n}Q", data))


# ---------------------------------------------------------------------------
# Subcommand: count-frames
# ---------------------------------------------------------------------------

def cmd_count_frames(args: argparse.Namespace) -> None:
    adapter = fmt_registry.resolve_adapter(args.format, args.input)
    spec = adapter.probe(args.input)  # type: ignore[attr-defined]
    n = adapter.count_frames(args.input, spec)  # type: ignore[attr-defined]
    print(n)


# ---------------------------------------------------------------------------
# Subcommand: build-index
# ---------------------------------------------------------------------------

def cmd_build_index(args: argparse.Namespace) -> None:
    """Scan the source trajectory and write a binary frame byte-offset index.

    The index allows Slurm workers to seek directly to their start frame
    (O(1)) instead of skipping via readline() (O(start * lines_per_frame)).
    For lammps_custom the scan uses grep at ~1-2 GB/s.  For lammps_xyz it
    falls back to a Python binary readline scan.

    Index format: packed little-endian uint64 values, one per frame.
    """
    adapter = fmt_registry.resolve_adapter(args.format, args.input)
    spec = adapter.probe(args.input)  # type: ignore[attr-defined]

    if not hasattr(adapter, "build_index"):
        _abort(f"format adapter for {args.input!r} does not support build-index")

    offsets = adapter.build_index(args.input, spec)  # type: ignore[attr-defined]

    out_path = args.output if args.output else args.input + ".idx"
    _save_index(offsets, out_path)
    print(f"build-index: {len(offsets)} frames → {out_path}")


# ---------------------------------------------------------------------------
# Subcommand: subset
# ---------------------------------------------------------------------------

def cmd_subset(args: argparse.Namespace) -> None:
    bbox = BBox(
        xmin=args.xmin, xmax=args.xmax,
        ymin=args.ymin, ymax=args.ymax,
        zmin=args.zmin, zmax=args.zmax,
    )
    try:
        bbox.validate()
    except ValueError as e:
        _abort(str(e))

    input_path = _resolve(args.input)
    output_path = _resolve(args.output) if args.output else _resolve(str(_default_output(args.input)))

    if input_path == output_path:
        _abort(f"input and output paths are the same file: {input_path}")

    adapter = fmt_registry.resolve_adapter(args.format, str(input_path))
    spec = adapter.probe(str(input_path))  # type: ignore[attr-defined]

    start_frame: int = args.frames[0] if args.frames else 0
    end_frame: Optional[int] = args.frames[1] if args.frames else None
    total = adapter.count_frames(str(input_path), spec)  # type: ignore[attr-defined]
    if end_frame is None:
        end_frame = total

    if start_frame >= end_frame:
        _abort(f"--frames START ({start_frame}) must be less than END ({end_frame})")

    # Load byte-offset index for O(1) seeking, if provided
    frame_offsets: Optional[list[int]] = None
    if args.frame_index:
        frame_offsets = _load_index(args.frame_index)

    stride: int = args.stride
    verbose: bool = args.verbose

    frames_read = 0
    frames_written = 0
    atoms_kept_total = 0
    atoms_total = 0
    t0 = time.perf_counter()

    with open(input_path, "r") as infile, open(output_path, "w") as outfile:
        if start_frame > 0:
            if frame_offsets is not None:
                # O(1) seek — the entire point of the index
                infile.seek(frame_offsets[start_frame])
            else:
                # Fallback: O(start * lines_per_frame) readline skip
                adapter.skip_frames(infile, spec, start_frame)  # type: ignore[attr-defined]

        for frame_idx in range(start_frame, end_frame):
            frame = adapter.read_frame(infile, spec)  # type: ignore[attr-defined]
            if frame is None:
                break

            frames_read += 1
            atoms_total += frame.natoms

            # Stride: skip frames that don't fall on the stride boundary
            if (frame_idx - start_frame) % stride != 0:
                continue

            filtered = filter_by_bbox(frame, bbox)
            adapter.write_frame(outfile, filtered, spec)  # type: ignore[attr-defined]
            atoms_kept_total += filtered.natoms
            frames_written += 1

            if verbose:
                print(
                    f"  frame {frame_idx:>6d}  ts={frame.timestep!s:>10s}  "
                    f"kept {filtered.natoms:>8,d} / {frame.natoms:,d}",
                    flush=True,
                )

    elapsed = time.perf_counter() - t0
    fps = frames_written / elapsed if elapsed > 0 else 0
    print(
        f"subset: {frames_written} frames written to {output_path}  "
        f"({atoms_kept_total:,}/{atoms_total:,} atom-slots kept)  "
        f"[{elapsed:.1f}s  {fps:.1f} frames/s  stride={stride}]"
    )


# ---------------------------------------------------------------------------
# Merge validation helpers (V1-V7)
# ---------------------------------------------------------------------------

def _stream_frame_headers(path: str, adapter) -> list[tuple[int, Optional[int]]]:
    """Stream all frame headers from path; return list of (natoms, timestep).
    Uses read_frame_header — never parses atom coordinates.
    """
    results = []
    with open(path, "r") as f:
        while True:
            hdr = adapter.read_frame_header(f)  # type: ignore[attr-defined]
            if hdr is None:
                break
            results.append(hdr)
    return results


def _validate_parts_pre(
    sorted_part_paths: list[str],
    expected_parts: Optional[int],
    output_path: Path,
    adapter,
) -> dict[str, list[tuple[int, Optional[int]]]]:
    """Run V1–V5.  Returns per-part header lists (cached from V3 pass)."""

    # V5: output path safety
    for p in sorted_part_paths:
        if _resolve(p) == output_path:
            _abort(f"output path {output_path} matches part file {p}", "V5")

    # V1: part file count
    if expected_parts is not None and len(sorted_part_paths) != expected_parts:
        found_ids = [Path(p).stem for p in sorted_part_paths]
        _abort(
            f"expected {expected_parts} parts, found {len(sorted_part_paths)}. "
            f"Found stems: {found_ids}",
            "V1",
        )

    # V2: non-empty parts
    for p in sorted_part_paths:
        if not os.path.exists(p) or os.path.getsize(p) == 0:
            _abort(f"part file is missing or zero-byte: {p}", "V2")

    # V3 + V4 combined: stream frame headers per part
    part_headers: dict[str, list[tuple[int, Optional[int]]]] = {}
    for i, part in enumerate(sorted_part_paths):
        try:
            headers = _stream_frame_headers(part, adapter)
        except Exception as e:
            _abort(f"part {part} failed frame-header validation: {e}", "V3")
        if not headers:
            _abort(f"part {part} contains no frames", "V3")
        part_headers[part] = headers

    # V4: boundary timestep sanity — last(part_i) <= first(part_i+1)
    for i in range(len(sorted_part_paths) - 1):
        part_a = sorted_part_paths[i]
        part_b = sorted_part_paths[i + 1]
        headers_a = part_headers[part_a]
        headers_b = part_headers[part_b]
        last_ts = headers_a[-1][1]
        first_ts = headers_b[0][1]
        if last_ts is None or first_ts is None:
            # lammps_xyz can have None if comment line was unparseable; skip
            continue
        if last_ts > first_ts:
            _abort(
                f"boundary violation between {Path(part_a).name} "
                f"(last timestep {last_ts}) and {Path(part_b).name} "
                f"(first timestep {first_ts})",
                "V4",
            )

    return part_headers


def _validate_merged_post(
    output_path: Path,
    part_headers: dict[str, list[tuple[int, Optional[int]]]],
    sorted_part_paths: list[str],
    adapter,
) -> None:
    """Run V6 and V7 in a single header-stream pass over the merged output."""
    expected_frame_count = sum(len(part_headers[p]) for p in sorted_part_paths)

    merged_headers = _stream_frame_headers(str(output_path), adapter)

    # V6: frame count accounting
    if len(merged_headers) != expected_frame_count:
        _abort(
            f"merged output has {len(merged_headers)} frames; "
            f"expected {expected_frame_count} (sum of part frame counts)",
            "V6",
        )

    # V7: strict global monotonic timesteps
    prev_ts: Optional[int] = None
    for frame_idx, (_, ts) in enumerate(merged_headers):
        if ts is None:
            continue  # format doesn't carry timestep; skip ordering check
        if prev_ts is not None and ts <= prev_ts:
            _abort(
                f"timestep ordering violation at merged frame {frame_idx}: "
                f"timestep {ts} follows {prev_ts} (must be strictly increasing)",
                "V7",
            )
        prev_ts = ts


# ---------------------------------------------------------------------------
# Subcommand: merge
# ---------------------------------------------------------------------------

def cmd_merge(args: argparse.Namespace) -> None:
    # Resolve and numerically sort parts
    raw = args.parts if args.parts else []
    sorted_part_paths = _sorted_parts(raw)

    if not sorted_part_paths:
        _abort("no part files found — check --parts argument")

    output_path = _resolve(args.output)
    temp_path = output_path.parent / (output_path.name + ".merging")

    # Resolve adapter from first part file
    adapter = fmt_registry.resolve_adapter(args.format, sorted_part_paths[0])

    # --- Pre-merge validation (V1–V5) ---
    part_headers = _validate_parts_pre(
        sorted_part_paths,
        args.expected_parts,
        output_path,
        adapter,
    )

    # --- Concatenate to temp file ---
    try:
        with open(temp_path, "w") as outfile:
            for part in sorted_part_paths:
                with open(part, "r") as infile:
                    for chunk in iter(lambda: infile.read(1 << 20), ""):
                        outfile.write(chunk)
    except Exception as e:
        if temp_path.exists():
            temp_path.unlink()
        _abort(f"concatenation failed: {e}")

    # Atomically replace output
    temp_path.replace(output_path)

    # --- Post-merge validation (V6, V7) ---
    _validate_merged_post(output_path, part_headers, sorted_part_paths, adapter)

    total_frames = sum(len(part_headers[p]) for p in sorted_part_paths)
    print(f"merge: {total_frames} frames → {output_path}")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def _add_bbox_args(p: argparse.ArgumentParser) -> None:
    for axis in ("x", "y", "z"):
        p.add_argument(f"--{axis}min", type=float, required=True,
                       help=f"Minimum {axis} coordinate (Å)")
        p.add_argument(f"--{axis}max", type=float, required=True,
                       help=f"Maximum {axis} coordinate (Å)")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="lammps_bbox_subset",
        description="Extract a spatial subset of a LAMMPS trajectory by bounding box.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # subset
    p_sub = sub.add_parser("subset", help="Filter atoms by bbox and write new trajectory.")
    p_sub.add_argument("input", help="Source trajectory file (read-only).")
    p_sub.add_argument("--format", default="auto",
                       choices=["auto", "lammps_custom", "lammps_xyz"],
                       help="File format (default: auto-detect).")
    p_sub.add_argument("--output", default=None,
                       help="Output path (default: {stem}_bbox{ext} beside input).")
    p_sub.add_argument("--frames", nargs=2, type=int, metavar=("START", "END"),
                       help="Process only frames [START, END). Used by Slurm workers.")
    p_sub.add_argument("--frame-index", default=None, dest="frame_index", metavar="PATH",
                       help="Binary frame-offset index (.idx) produced by build-index. "
                            "Enables O(1) seek to start frame instead of O(n) readline skip.")
    p_sub.add_argument("--stride", type=int, default=1,
                       help="Write every Nth frame (default: 1 = all frames). "
                            "Frames are still read sequentially; only writing is skipped.")
    p_sub.add_argument("--verbose", action="store_true",
                       help="Print per-frame atom counts and timing.")
    _add_bbox_args(p_sub)

    # count-frames
    p_count = sub.add_parser("count-frames",
                              help="Print total frame count to stdout.")
    p_count.add_argument("input", help="Source trajectory file.")
    p_count.add_argument("--format", default="auto",
                         choices=["auto", "lammps_custom", "lammps_xyz"])

    # build-index
    p_idx = sub.add_parser(
        "build-index",
        help="Scan source trajectory and write a binary byte-offset frame index (.idx). "
             "Pass the resulting file to subset via --frame-index for O(1) seeking.",
    )
    p_idx.add_argument("input", help="Source trajectory file.")
    p_idx.add_argument("--format", default="auto",
                       choices=["auto", "lammps_custom", "lammps_xyz"])
    p_idx.add_argument("--output", default=None,
                       help="Index output path (default: INPUT.idx).")

    # merge
    p_merge = sub.add_parser("merge",
                              help="Concatenate part files with V1-V7 validation.")
    p_merge.add_argument("--parts", nargs="+", required=True,
                         help="Part files to merge (numerically sorted).")
    p_merge.add_argument("--output", required=True,
                         help="Output merged trajectory path.")
    p_merge.add_argument("--format", default="auto",
                         choices=["auto", "lammps_custom", "lammps_xyz"])
    p_merge.add_argument("--expected-parts", type=int, default=None,
                         dest="expected_parts",
                         help="Expected number of part files (V1 check). "
                              "Set to $SLURM_NTASKS in Slurm scripts.")

    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    if args.command == "subset":
        cmd_subset(args)
    elif args.command == "count-frames":
        cmd_count_frames(args)
    elif args.command == "build-index":
        cmd_build_index(args)
    elif args.command == "merge":
        cmd_merge(args)


if __name__ == "__main__":
    main()
