#!/usr/bin/env python3
"""LAMMPS trajectory bounding-box subset tool.

Subcommands
-----------
subset        Filter atoms by Cartesian bbox.

              Single-process:
                  --output PATH  → all frames written to one file.

              Parallel worker (round-robin, one file per frame):
                  --worker W --num-workers N --output-dir DIR --frame-index IDX
                  Worker W processes frames W, W+N, W+2N, …
                  Each frame written to DIR/frame_{i:07d}.{ext}.
                  Requires --frame-index for O(1) per-frame seeking.

count-frames  Print total source frame count to stdout.
build-index   Scan source and write a binary byte-offset frame index (.idx).
merge         Sort frame files in --parts-dir by index, validate, concatenate.

Assumptions (must hold for source trajectories):
  - natoms is constant across all frames (NVT/NVE).
  - Every bbox-filtered frame contains at least one atom.
  - lammps_custom: ITEM: BOX BOUNDS written unchanged (floating-sliver).

Requires: numpy
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

sys.path.insert(0, str(Path(__file__).parent))

import formats as fmt_registry
from frame_filter import BBox


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _resolve(path: str) -> Path:
    return Path(path).resolve()


def _abort(msg: str, check_id: str = "") -> None:
    prefix = f"merge validation failed: {check_id} — " if check_id else ""
    print(f"ERROR: {prefix}{msg}", file=sys.stderr)
    sys.exit(1)


def _default_output(input_path: str) -> Path:
    p = Path(input_path)
    return p.parent / f"{p.stem}_bbox{p.suffix}"


def _input_ext(input_path: str) -> str:
    return Path(input_path).suffix.lstrip(".")


# ---------------------------------------------------------------------------
# Frame index I/O
# ---------------------------------------------------------------------------

def _save_index(offsets: list[int], path: str) -> None:
    with open(path, "wb") as f:
        f.write(struct.pack(f"<{len(offsets)}Q", *offsets))


def _load_index(path: str) -> list[int]:
    with open(path, "rb") as f:
        data = f.read()
    n = len(data) // 8
    return list(struct.unpack(f"<{n}Q", data))


# ---------------------------------------------------------------------------
# Subcommand: count-frames
# ---------------------------------------------------------------------------

def cmd_count_frames(args: argparse.Namespace) -> None:
    # If the index already exists alongside the source, read its length directly
    # (O(1) stat + small binary read) rather than re-scanning the source file.
    default_idx = args.input + ".idx"
    if Path(default_idx).exists():
        print(len(_load_index(default_idx)))
        return
    adapter = fmt_registry.resolve_adapter(args.format, args.input)
    spec = adapter.probe(args.input)
    print(adapter.count_frames(args.input, spec))


# ---------------------------------------------------------------------------
# Subcommand: build-index
# ---------------------------------------------------------------------------

def cmd_build_index(args: argparse.Namespace) -> None:
    adapter = fmt_registry.resolve_adapter(args.format, args.input)
    spec = adapter.probe(args.input)
    offsets = adapter.build_index(args.input, spec)
    out_path = args.output or (args.input + ".idx")
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
    bbox.validate()

    input_path = _resolve(args.input)
    adapter = fmt_registry.resolve_adapter(args.format, str(input_path))
    spec = adapter.probe(str(input_path))
    total = adapter.count_frames(str(input_path), spec)

    batch_size: int = args.batch_size
    verbose: bool = args.verbose
    t0 = time.perf_counter()

    if args.worker is not None:
        # ── Parallel round-robin mode ────────────────────────────────────────
        # Worker W processes frames W, W+N, W+2N, …
        # Each frame → output_dir/frame_{i:07d}.{ext}
        worker: int = args.worker
        num_workers: int = args.num_workers
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        ext = _input_ext(str(input_path))

        if not args.frame_index:
            _abort("--frame-index is required in worker mode (needed for O(1) per-frame seek)")
        frame_offsets = _load_index(args.frame_index)

        frames_written = 0
        atoms_kept_total = 0

        with open(input_path, "r") as infile:
            for frame_idx in range(worker, total, num_workers):
                out_path = output_dir / f"frame_{frame_idx:07d}.{ext}"
                infile.seek(frame_offsets[frame_idx])
                with open(out_path, "w") as outfile:
                    result = adapter.stream_filter_frame(
                        infile, outfile, spec, bbox, batch_size
                    )
                if result is None:
                    break
                n_kept, timestep = result
                frames_written += 1
                atoms_kept_total += n_kept
                if verbose:
                    print(
                        f"  frame {frame_idx:>7d}  ts={timestep!s:>10s}  "
                        f"kept {n_kept:>10,d} / {spec.natoms:,d}",
                        flush=True,
                    )

        elapsed = time.perf_counter() - t0
        fps = frames_written / elapsed if elapsed > 0 else 0
        print(
            f"subset worker {worker}/{num_workers}: {frames_written} frames written  "
            f"({atoms_kept_total:,} atoms kept)  [{elapsed:.1f}s  {fps:.1f} frames/s]"
        )

    else:
        # ── Single-process mode ──────────────────────────────────────────────
        output_path = _resolve(args.output) if args.output else _resolve(str(_default_output(args.input)))
        if input_path == output_path:
            _abort(f"input and output resolve to the same file: {input_path}")

        frame_offsets: Optional[list[int]] = None
        if args.frame_index:
            frame_offsets = _load_index(args.frame_index)

        stride: int = args.stride
        frames_written = 0
        atoms_kept_total = 0

        with open(input_path, "r") as infile, open(output_path, "w") as outfile:
            for frame_idx in range(total):
                should_write = frame_idx % stride == 0
                f_out = outfile if should_write else None

                if frame_offsets is not None:
                    infile.seek(frame_offsets[frame_idx])
                    result = adapter.stream_filter_frame(
                        infile, f_out, spec, bbox, batch_size
                    )
                else:
                    result = adapter.stream_filter_frame(
                        infile, f_out, spec, bbox, batch_size
                    )

                if result is None:
                    break
                n_kept, timestep = result

                if should_write:
                    frames_written += 1
                    atoms_kept_total += n_kept
                    if verbose:
                        print(
                            f"  frame {frame_idx:>7d}  ts={timestep!s:>10s}  "
                            f"kept {n_kept:>10,d} / {spec.natoms:,d}",
                            flush=True,
                        )

        elapsed = time.perf_counter() - t0
        fps = frames_written / elapsed if elapsed > 0 else 0
        print(
            f"subset: {frames_written} frames → {output_path}  "
            f"({atoms_kept_total:,} atoms kept)  "
            f"[{elapsed:.1f}s  {fps:.1f} frames/s  stride={stride}]"
        )


# ---------------------------------------------------------------------------
# Subcommand: merge
# ---------------------------------------------------------------------------

def cmd_merge(args: argparse.Namespace) -> None:
    parts_dir = Path(args.parts_dir)
    output_path = _resolve(args.output)
    temp_path = output_path.parent / (output_path.name + ".merging")

    # Detect extension from output path, find frame files
    ext = output_path.suffix.lstrip(".")
    pattern = str(parts_dir / f"frame_*.{ext}")
    frame_files = sorted(
        glob.glob(pattern),
        key=lambda p: int(re.search(r"frame_(\d+)", Path(p).name).group(1)),
    )

    if not frame_files:
        _abort(f"no frame files found matching {pattern}")

    # V5: output not inside parts_dir under a colliding name
    if _resolve(str(parts_dir)) == output_path.parent and output_path.name.startswith("frame_"):
        _abort(f"output {output_path} would collide with frame files in {parts_dir}", "V5")

    # V1: expected frame count
    if args.expected_frames is not None and len(frame_files) != args.expected_frames:
        _abort(
            f"expected {args.expected_frames} frames, found {len(frame_files)}",
            "V1",
        )

    # V2: non-empty files
    for p in frame_files:
        if os.path.getsize(p) == 0:
            _abort(f"frame file is zero-byte: {p}", "V2")

    # Concatenate to temp, then atomically replace
    try:
        with open(temp_path, "w") as out:
            for p in frame_files:
                with open(p) as f:
                    for chunk in iter(lambda: f.read(1 << 20), ""):
                        out.write(chunk)
    except Exception as e:
        if temp_path.exists():
            temp_path.unlink()
        _abort(f"concatenation failed: {e}")

    temp_path.replace(output_path)

    # V7: strict global timestep monotonicity (single pass over merged output)
    adapter = fmt_registry.resolve_adapter(args.format, str(output_path))
    prev_ts: Optional[int] = None
    frame_count = 0
    with open(output_path) as f:
        while True:
            hdr = adapter.read_frame_header(f)
            if hdr is None:
                break
            _, ts = hdr
            frame_count += 1
            if ts is not None:
                if prev_ts is not None and ts <= prev_ts:
                    _abort(
                        f"timestep disorder at frame {frame_count}: "
                        f"{ts} follows {prev_ts} (must be strictly increasing)",
                        "V7",
                    )
                prev_ts = ts

    # V6: merged frame count matches file count
    if frame_count != len(frame_files):
        _abort(
            f"merged output has {frame_count} frames; expected {len(frame_files)}",
            "V6",
        )

    print(f"merge: {frame_count} frames → {output_path}")


# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------

def _add_bbox_args(p: argparse.ArgumentParser) -> None:
    for axis in ("x", "y", "z"):
        p.add_argument(f"--{axis}min", type=float, required=True)
        p.add_argument(f"--{axis}max", type=float, required=True)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="lammps_bbox_subset",
        description="Extract a spatial subset of a LAMMPS trajectory by bounding box.",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    # ── subset ────────────────────────────────────────────────────────────────
    p_sub = sub.add_parser("subset",
                            help="Filter atoms by bbox.")
    p_sub.add_argument("input")
    p_sub.add_argument("--format", default="auto",
                       choices=["auto", "lammps_custom", "lammps_xyz"])
    p_sub.add_argument("--frame-index", default=None, dest="frame_index",
                       metavar="PATH",
                       help="Binary frame-offset index (.idx) from build-index.")
    p_sub.add_argument("--batch-size", type=int, default=100_000, dest="batch_size",
                       metavar="N",
                       help="Atom lines per batch (default 100000). "
                            "Peak RAM ≈ 27 MB × N/100000 per worker.")

    # Single-process options
    p_sub.add_argument("--output", default=None,
                       help="Single-process: output file path.")
    p_sub.add_argument("--stride", type=int, default=1,
                       help="Single-process: write every Nth frame (default 1).")
    p_sub.add_argument("--verbose", action="store_true")

    # Parallel worker options
    p_sub.add_argument("--worker", type=int, default=None, metavar="W",
                       help="Parallel: worker index (0-based).")
    p_sub.add_argument("--num-workers", type=int, default=None, metavar="N",
                       help="Parallel: total number of workers.")
    p_sub.add_argument("--output-dir", default=None, dest="output_dir",
                       metavar="DIR",
                       help="Parallel: directory for per-frame output files "
                            "(frame_{i:07d}.{ext}).")
    _add_bbox_args(p_sub)

    # ── count-frames ──────────────────────────────────────────────────────────
    p_count = sub.add_parser("count-frames")
    p_count.add_argument("input")
    p_count.add_argument("--format", default="auto",
                         choices=["auto", "lammps_custom", "lammps_xyz"])

    # ── build-index ───────────────────────────────────────────────────────────
    p_idx = sub.add_parser("build-index")
    p_idx.add_argument("input")
    p_idx.add_argument("--format", default="auto",
                       choices=["auto", "lammps_custom", "lammps_xyz"])
    p_idx.add_argument("--output", default=None,
                       help="Index path (default: INPUT.idx).")

    # ── merge ─────────────────────────────────────────────────────────────────
    p_merge = sub.add_parser("merge",
                              help="Concatenate per-frame files → single trajectory.")
    p_merge.add_argument("--parts-dir", required=True, dest="parts_dir",
                         metavar="DIR",
                         help="Directory containing frame_*.{ext} files.")
    p_merge.add_argument("--output", required=True)
    p_merge.add_argument("--format", default="auto",
                         choices=["auto", "lammps_custom", "lammps_xyz"])
    p_merge.add_argument("--expected-frames", type=int, default=None,
                         dest="expected_frames",
                         help="Expected number of frame files (= total source frames). "
                              "Pass $TOTAL from Slurm script.")

    return parser


def main() -> None:
    args = build_parser().parse_args()
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
