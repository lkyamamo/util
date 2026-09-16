"""
Render reaction-path LAMMPS data frames into one GIF per site (OVITO + Tachyon).

Walks a run tree, finds every directory holding numerically-named LAMMPS data
files (e.g. `0264/setup/atom1/frames/{1..7}.data`), and renders each such frame
set to an animated GIF.

The camera is placed once per site and held fixed for the whole sequence, so the
only thing that moves in the GIF is the structure. By default the view is cut
down to a sphere around the reacting atoms (found from `manifest.csv`, or from
whichever atoms actually move between the first and last frame), which is the
part anyone wants to look at — the full slab renders as an undifferentiated blob.

Typical use, from a run tree such as .../runs/reaction_energy:

    python reaction_path_gifs.py .                 # -> ./gifs/0264_atom1.gif, ...
    python reaction_path_gifs.py . --radius 16 --fps 3 --pingpong
    python reaction_path_gifs.py 0266 --out-dir /tmp/g --jobs 4

To animate the minimized outputs instead of the setup path:

    python reaction_path_gifs.py . --pattern '*_min.data'
"""

from __future__ import annotations

import argparse
import csv
import math
import multiprocessing
import re
import sys
import tempfile
from collections.abc import Iterable, Mapping, Sequence
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np

# Sphere radii (Å) keyed by element symbol. Kept well under the covalent radii so
# the bonds stay visible — these are ball-and-stick, not space-filling.
ATOM_DISPLAY_RADII: dict[str, float] = {"Si": 0.50, "O": 0.36, "H": 0.20}

# Bonds are drawn when a pair is closer than this (Å). Pairs absent here get none.
BOND_PAIR_CUTOFFS: dict[tuple[str, str], float] = {
    ("Si", "O"): 2.10,
    ("O", "H"): 1.25,
}
BOND_VISUAL_WIDTH = 0.22

# Colors for the atoms being tracked, so the incoming molecule stands out
# against the identically-typed atoms of the slab.
TRACKED_COLORS: dict[str, tuple[float, float, float]] = {
    "O": (0.15, 0.45, 1.00),
    "H": (0.55, 0.85, 1.00),
    "*": (0.20, 0.90, 0.35),
}
TRACKED_RADIUS_SCALE = 1.25

# Element symbol by atomic mass (amu), for data files whose types are unnamed.
_MASS_TO_SYMBOL: tuple[tuple[float, str], ...] = (
    (1.008, "H"),
    (12.011, "C"),
    (14.007, "N"),
    (15.999, "O"),
    (22.990, "Na"),
    (26.982, "Al"),
    (28.085, "Si"),
    (30.974, "P"),
    (32.06, "S"),
    (35.45, "Cl"),
    (39.098, "K"),
    (40.078, "Ca"),
)

_ELEMENT_NAME = re.compile(r"^[A-Z][a-z]?$")


# --------------------------------------------------------------------------
# Discovery
# --------------------------------------------------------------------------


@dataclass(frozen=True)
class FrameSet:
    """One directory of numbered data files that becomes one GIF."""

    directory: Path
    files: tuple[Path, ...]
    name: str

    @property
    def manifest(self) -> Path:
        return self.directory / "manifest.csv"


def _frame_index(path: Path) -> int | None:
    """Leading integer of a frame filename (`12.data`, `12_min.data`) or None."""
    match = re.match(r"^(\d+)", path.stem)
    return int(match.group(1)) if match else None


def _frame_set_name(directory: Path, root: Path) -> str:
    """Flatten a path into a GIF basename: 0264/setup/atom1/frames -> 0264_atom1."""
    try:
        parts = list(directory.relative_to(root).parts)
    except ValueError:
        parts = [directory.name]
    parts = [p for p in parts if p not in {"frames", "setup", "run", "."}]
    if not parts:
        parts = [root.resolve().name or directory.name]
    return "_".join(parts)


def discover_frame_sets(root: Path, pattern: str, min_frames: int) -> list[FrameSet]:
    """Every directory under `root` holding at least `min_frames` numbered files."""
    found: list[FrameSet] = []
    for directory in sorted({p.parent for p in root.rglob(pattern)}):
        numbered = [(i, p) for p in directory.glob(pattern) if (i := _frame_index(p)) is not None]
        if len(numbered) < min_frames:
            continue
        numbered.sort(key=lambda pair: (pair[0], pair[1].name))
        found.append(
            FrameSet(
                directory=directory,
                files=tuple(p for _, p in numbered),
                name=_frame_set_name(directory, root),
            )
        )
    return found


# --------------------------------------------------------------------------
# Where to point the camera, and what to highlight
# --------------------------------------------------------------------------


@dataclass
class SiteInfo:
    """The region of interest for one frame set."""

    center: np.ndarray  # (3,) point the camera looks at
    path_direction: np.ndarray | None  # unit vector the tracked atoms travel along
    tracked_ids: frozenset[int]  # particle identifiers to highlight
    labels: tuple[str, ...]  # per-frame caption, or empty


def _read_manifest(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def _manifest_site(manifest_rows: Sequence[Mapping[str, str]], n_frames: int) -> SiteInfo | None:
    """Read the moving oxygen's track and the water ids straight out of manifest.csv."""
    if len(manifest_rows) != n_frames:
        return None
    try:
        positions = np.array(
            [[float(row["o_x"]), float(row["o_y"]), float(row["o_z"])] for row in manifest_rows]
        )
    except (KeyError, ValueError):
        return None

    ids: set[int] = set()
    for key in ("o_id", "h1_id", "h2_id"):
        for row in manifest_rows:
            try:
                ids.add(int(row[key]))
            except (KeyError, ValueError):
                pass

    labels = []
    for i, row in enumerate(manifest_rows, start=1):
        caption = f"frame {i}/{n_frames}"
        distance = row.get("o_si_distance")
        if distance:
            try:
                caption += f"    d = {float(distance):.2f} A"
            except ValueError:
                pass
        labels.append(caption)

    return SiteInfo(
        center=positions.mean(axis=0),
        path_direction=_unit(positions[-1] - positions[0]),
        tracked_ids=frozenset(ids),
        labels=tuple(labels),
    )


def _motion_site(pipeline, n_frames: int, min_displacement: float) -> SiteInfo:
    """Fallback: find the atoms that move between the first and last frame."""
    first = pipeline.compute(0)
    last = pipeline.compute(n_frames - 1)
    cell = np.asarray(first.cell[...])
    labels = tuple(f"frame {i}/{n_frames}" for i in range(1, n_frames + 1))
    whole_cell = SiteInfo(
        center=_cell_center(cell),
        path_direction=None,
        tracked_ids=frozenset(),
        labels=labels,
    )

    ids_first = np.asarray(first.particles["Particle Identifier"][...])
    ids_last = np.asarray(last.particles["Particle Identifier"][...])
    order_first = np.argsort(ids_first)
    order_last = np.argsort(ids_last)
    if not np.array_equal(ids_first[order_first], ids_last[order_last]):
        # Atoms were added or removed along the way; there is no correspondence
        # to difference, so fall back to viewing the whole cell.
        return whole_cell

    start = np.asarray(first.particles.positions[...])[order_first]
    end = np.asarray(last.particles.positions[...])[order_last]
    delta = _minimum_image(end - start, cell[:, :3], tuple(bool(f) for f in first.cell.pbc))

    moving = np.linalg.norm(delta, axis=1) > min_displacement
    if not moving.any():
        return whole_cell

    midpoints = start[moving] + 0.5 * delta[moving]
    return SiteInfo(
        center=midpoints.mean(axis=0),
        path_direction=_unit(delta[moving].mean(axis=0)),
        tracked_ids=frozenset(int(i) for i in ids_first[order_first][moving]),
        labels=labels,
    )


def _cell_center(cell: np.ndarray) -> np.ndarray:
    """Center of an OVITO 3x4 cell matrix (column 3 is the origin)."""
    return cell[:, :3].sum(axis=0) / 2.0 + cell[:, 3]


def _minimum_image(delta: np.ndarray, cell: np.ndarray, pbc: Sequence[bool]) -> np.ndarray:
    """Wrap displacements into the primary image along each periodic axis."""
    if not any(pbc):
        return delta
    try:
        inverse = np.linalg.inv(cell)
    except np.linalg.LinAlgError:
        return delta
    fractional = delta @ inverse
    for axis, periodic in enumerate(pbc):
        if periodic:
            fractional[:, axis] -= np.round(fractional[:, axis])
    return fractional @ cell


def _unit(vector: np.ndarray) -> np.ndarray | None:
    norm = float(np.linalg.norm(vector))
    return None if norm < 1e-9 else np.asarray(vector, dtype=float) / norm


# --------------------------------------------------------------------------
# Pipeline
# --------------------------------------------------------------------------


def _symbol_for_type(particle_type) -> str:
    """Element symbol for an OVITO particle type, from its name or its mass."""
    name = str(particle_type.name).strip()
    if _ELEMENT_NAME.match(name):
        return name
    mass = float(getattr(particle_type, "mass", 0.0) or 0.0)
    if mass > 0.0:
        symbol, error = min(
            ((sym, abs(mass - ref)) for ref, sym in _MASS_TO_SYMBOL), key=lambda pair: pair[1]
        )
        if error < 0.5:
            return symbol
    return name


def _make_local_cut(site: SiteInfo, radius: float):
    """Keep only atoms within `radius` of the site, unwrapped around it."""
    center = np.asarray(site.center, dtype=float)

    def modify(frame: int, data):  # noqa: ARG001 — OVITO modifier signature
        cell_matrix = np.asarray(data.cell[...])[:, :3]
        pbc = tuple(bool(f) for f in data.cell.pbc)
        delta = _minimum_image(np.asarray(data.particles.positions[...]) - center, cell_matrix, pbc)

        keep = np.linalg.norm(delta, axis=1) <= radius
        if not keep.any():
            return
        data.particles_.delete_indices(np.flatnonzero(~keep))
        # Positions are now contiguous around the site, so bonds and the camera
        # do not have to reason about periodic images.
        data.particles_.positions_[...] = center + delta[keep]
        data.cell_.pbc = (False, False, False)

    return modify


def _make_style(site: SiteInfo, radii: Mapping[str, float], highlight: bool):
    """Assign per-type display radii and recolor the tracked atoms."""
    tracked = site.tracked_ids

    def modify(frame: int, data):  # noqa: ARG001 — OVITO modifier signature
        types = data.particles_.particle_types_
        symbols: dict[int, str] = {}
        for particle_type in types.types_:
            symbol = _symbol_for_type(particle_type)
            symbols[int(particle_type.id)] = symbol
            if symbol in radii:
                particle_type.radius = float(radii[symbol])

        if not (highlight and tracked) or data.particles.count == 0:
            return
        identifiers = data.particles["Particle Identifier"]
        if identifiers is None:
            return
        selected = np.isin(np.asarray(identifiers[...]), np.fromiter(tracked, dtype=np.int64))
        if not selected.any():
            return

        colors = data.particles_.create_property("Color")
        sizes = data.particles_.create_property("Radius")
        type_ids = np.asarray(data.particles["Particle Type"][...])
        base_colors = {int(t.id): np.asarray(t.color, dtype=float) for t in types.types_}
        base_radii = {int(t.id): float(t.radius) for t in types.types_}
        with colors, sizes:
            for index in np.flatnonzero(selected):
                type_id = int(type_ids[index])
                symbol = symbols.get(type_id, "")
                colors[index] = TRACKED_COLORS.get(symbol, TRACKED_COLORS["*"])
                sizes[index] = base_radii.get(type_id, 0.5) * TRACKED_RADIUS_SCALE
            for index in np.flatnonzero(~selected):
                type_id = int(type_ids[index])
                colors[index] = base_colors.get(type_id, np.array([0.7, 0.7, 0.7]))
                sizes[index] = base_radii.get(type_id, 0.5)

    return modify


def _make_labeler(labels: Sequence[str]):
    def modify(frame: int, data):
        data.attributes["FrameLabel"] = labels[frame] if frame < len(labels) else f"frame {frame + 1}"

    return modify


def _add_bonds(pipeline, cutoffs: Mapping[tuple[str, str], float], width: float) -> None:
    """Pairwise CreateBonds, skipped when no cutoff matches the types present."""
    from ovito.modifiers import CreateBondsModifier

    data = pipeline.compute(0)
    available = {_symbol_for_type(t) for t in data.particles["Particle Type"].types}

    modifier = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)
    n_set = 0
    for (a, b), cutoff in cutoffs.items():
        if cutoff <= 0.0 or a not in available or b not in available:
            continue
        modifier.set_pairwise_cutoff(a, b, float(cutoff))
        n_set += 1
    if n_set == 0:
        return
    # OVITO renamed BondsVis.radius to .width (diameter) in 3.8.
    if hasattr(modifier.vis, "width"):
        modifier.vis.width = float(width)
    else:
        modifier.vis.radius = float(width) / 2.0
    pipeline.modifiers.append(modifier)


# --------------------------------------------------------------------------
# Camera
# --------------------------------------------------------------------------


def _camera_direction(site: SiteInfo, override: np.ndarray | None) -> np.ndarray:
    """
    Look at the site from a direction perpendicular to the reaction path, so the
    atoms travel across the image rather than toward the viewer.
    """
    if override is not None:
        return override
    path = site.path_direction
    if path is None:
        return _unit(np.array([-1.0, -0.45, -0.30]))

    for reference in (np.array([0.0, 0.0, 1.0]), np.array([0.0, 1.0, 0.0])):
        perpendicular = np.cross(path, reference)
        direction = _unit(perpendicular)
        if direction is not None:
            # Tilt slightly off-axis so the structure reads as three-dimensional.
            return _unit(direction + 0.25 * np.cross(direction, path))
    return _unit(np.array([-1.0, -0.45, -0.30]))


def _place_camera(viewport, site: SiteInfo, radius: float, direction: np.ndarray, fov_deg: float):
    fov = math.radians(fov_deg)
    distance = radius / math.tan(fov / 2.0) * 1.15
    viewport.camera_dir = tuple(float(v) for v in direction)
    viewport.camera_pos = tuple(float(v) for v in (site.center - direction * distance))
    viewport.fov = fov


# --------------------------------------------------------------------------
# Rendering
# --------------------------------------------------------------------------


def _write_gif(png_paths: Sequence[Path], out_path: Path, fps: float, pingpong: bool, hold: int):
    """
    Assemble the rendered PNGs into a looping GIF.

    All frames are quantized against one palette built from the whole sequence,
    otherwise per-frame adaptive palettes make the slab shimmer between frames.
    """
    from PIL import Image

    frames = [Image.open(path).convert("RGB") for path in png_paths]
    if not frames:
        raise ValueError("no frames to write")

    stacked = np.concatenate([np.asarray(frame) for frame in frames], axis=0)
    palette = Image.fromarray(stacked).quantize(colors=255, method=Image.Quantize.MEDIANCUT)
    # No dithering: these are smooth shaded spheres on a flat background, and
    # dither noise is very visible on the few small, brightly colored atoms.
    quantized = [frame.quantize(palette=palette, dither=Image.Dither.NONE) for frame in frames]

    order = list(range(len(frames)))
    if pingpong and len(frames) > 2:
        order += list(reversed(range(1, len(frames) - 1)))

    milliseconds = max(20, int(round(1000.0 / fps)))
    durations = [milliseconds] * len(order)
    if hold > 1:
        durations[len(frames) - 1] *= hold  # linger on the end of the forward pass
        if pingpong:
            durations[0] *= hold  # ...and on the start, at the loop seam

    out_path.parent.mkdir(parents=True, exist_ok=True)
    quantized[order[0]].save(
        out_path,
        save_all=True,
        append_images=[quantized[i] for i in order[1:]],
        duration=durations,
        loop=0,
        optimize=True,
    )


def render_frame_set(frame_set: FrameSet, options: argparse.Namespace) -> Path:
    """Build the pipeline for one directory of frames and write its GIF."""
    from ovito.io import import_file
    from ovito.vis import TachyonRenderer, TextLabelOverlay, Viewport

    pipeline = import_file([str(p) for p in frame_set.files])
    n_frames = pipeline.num_frames

    site = _manifest_site(_read_manifest(frame_set.manifest), n_frames)
    if site is None:
        site = _motion_site(pipeline, n_frames, options.min_displacement)

    radius = options.radius
    if radius <= 0.0:
        cell = np.asarray(pipeline.compute(0).cell[...])
        site.center = cell[:, :3].sum(axis=0) / 2.0 + cell[:, 3]
        radius = float(np.linalg.norm(cell[:, :3].sum(axis=0))) / 2.0
    else:
        pipeline.modifiers.append(_make_local_cut(site, radius))

    pipeline.modifiers.append(_make_style(site, ATOM_DISPLAY_RADII, not options.no_highlight))
    if not options.no_bonds:
        _add_bonds(pipeline, BOND_PAIR_CUTOFFS, BOND_VISUAL_WIDTH)
    if site.labels and not options.no_label:
        pipeline.modifiers.append(_make_labeler(site.labels))

    source = pipeline.source.data
    if source.cell is not None:
        source.cell.vis.enabled = options.show_cell

    pipeline.add_to_scene()
    try:
        viewport = Viewport(type=Viewport.Type.Perspective)
        _place_camera(viewport, site, radius, _camera_direction(site, options.camera_dir), options.fov)

        if site.labels and not options.no_label:
            viewport.overlays.append(
                TextLabelOverlay(
                    text="[FrameLabel]",
                    source_pipeline=pipeline,
                    alignment=0x01 | 0x20,  # Qt.AlignLeft | Qt.AlignTop
                    offset_x=0.02,
                    offset_y=0.02,
                    font_size=0.035,
                    text_color=(0.1, 0.1, 0.1),
                )
            )

        renderer = TachyonRenderer(ambient_occlusion=not options.no_ao, shadows=not options.no_ao)
        out_path = options.out_dir / f"{frame_set.name}.gif"
        with tempfile.TemporaryDirectory(prefix="reaction_gif_") as tmp:
            png_paths = []
            for frame in range(n_frames):
                png = Path(tmp) / f"{frame:05d}.png"
                viewport.render_image(
                    size=options.size,
                    filename=str(png),
                    frame=frame,
                    renderer=renderer,
                    background=(1.0, 1.0, 1.0),
                    crop=False,
                )
                png_paths.append(png)
                if options.keep_pngs:
                    keep_dir = options.out_dir / frame_set.name
                    keep_dir.mkdir(parents=True, exist_ok=True)
                    keep_dir.joinpath(png.name).write_bytes(png.read_bytes())
            _write_gif(png_paths, out_path, options.fps, options.pingpong, options.hold)
    finally:
        pipeline.remove_from_scene()

    return out_path


def _render_worker(args: tuple[FrameSet, argparse.Namespace]) -> tuple[str, str | None]:
    frame_set, options = args
    try:
        return str(render_frame_set(frame_set, options)), None
    except Exception as exc:  # noqa: BLE001 — one bad frame set must not sink the batch
        return frame_set.name, f"{type(exc).__name__}: {exc}"


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------


def _parse_vector(text: str) -> np.ndarray:
    parts = [float(v) for v in re.split(r"[,\s]+", text.strip()) if v]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("expected three comma-separated numbers, e.g. -1,0,-0.3")
    vector = _unit(np.array(parts))
    if vector is None:
        raise argparse.ArgumentTypeError("direction must be non-zero")
    return vector


def _parse_size(text: str) -> tuple[int, int]:
    match = re.fullmatch(r"(\d+)\s*[xX×]\s*(\d+)", text.strip())
    if not match:
        raise argparse.ArgumentTypeError("expected WIDTHxHEIGHT, e.g. 800x600")
    return int(match.group(1)), int(match.group(2))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("root", nargs="?", default=".", type=Path, help="tree to search (default: .)")
    parser.add_argument("--pattern", default="*.data", help="frame filename glob (default: *.data)")
    parser.add_argument("--out-dir", type=Path, help="GIF output directory (default: ROOT/gifs)")
    parser.add_argument("--min-frames", type=int, default=2, help="skip directories with fewer frames")
    parser.add_argument(
        "--radius",
        type=float,
        default=12.0,
        help="Å around the reaction site to render; 0 renders the whole cell (default: 12)",
    )
    parser.add_argument("--size", type=_parse_size, default=(800, 600), help="render size WxH")
    parser.add_argument("--fps", type=float, default=2.5, help="GIF frame rate (default: 2.5)")
    parser.add_argument("--hold", type=int, default=3, help="repeat the last frame N times")
    parser.add_argument("--pingpong", action="store_true", help="play forward then back")
    parser.add_argument("--fov", type=float, default=35.0, help="vertical field of view, degrees")
    parser.add_argument("--camera-dir", type=_parse_vector, help="viewing direction, e.g. -1,0,-0.3")
    parser.add_argument(
        "--min-displacement",
        type=float,
        default=0.5,
        help="Å of motion that marks an atom as reacting, when there is no manifest.csv",
    )
    parser.add_argument("--no-bonds", action="store_true", help="draw atoms only")
    parser.add_argument("--no-highlight", action="store_true", help="do not recolor the moving atoms")
    parser.add_argument("--no-label", action="store_true", help="omit the frame/distance caption")
    parser.add_argument("--no-ao", action="store_true", help="disable shadows and ambient occlusion")
    parser.add_argument("--show-cell", action="store_true", help="draw the simulation cell box")
    parser.add_argument("--keep-pngs", action="store_true", help="also keep the rendered PNG frames")
    parser.add_argument("--overwrite", action="store_true", help="re-render GIFs that already exist")
    parser.add_argument("--jobs", type=int, default=1, help="frame sets to render in parallel")
    parser.add_argument("--dry-run", action="store_true", help="list what would be rendered and exit")
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    options = build_parser().parse_args(argv)

    root = options.root.resolve()
    if not root.is_dir():
        print(f"error: {root} is not a directory", file=sys.stderr)
        return 2
    options.out_dir = (options.out_dir or root / "gifs").resolve()

    frame_sets = discover_frame_sets(root, options.pattern, options.min_frames)
    if not frame_sets:
        print(f"No directories under {root} contain >= {options.min_frames} '{options.pattern}' frames.")
        return 1

    if not options.overwrite:
        skipped = [fs for fs in frame_sets if (options.out_dir / f"{fs.name}.gif").exists()]
        for frame_set in skipped:
            print(f"skip {frame_set.name} (exists; --overwrite to redo)")
        frame_sets = [fs for fs in frame_sets if fs not in skipped]

    if options.dry_run:
        for frame_set in frame_sets:
            print(
                f"{options.out_dir / (frame_set.name + '.gif')}"
                f"  <- {len(frame_set.files)} frames from {frame_set.directory}"
            )
        return 0
    if not frame_sets:
        return 0

    options.out_dir.mkdir(parents=True, exist_ok=True)
    work = [(fs, options) for fs in frame_sets]

    failures = 0
    if options.jobs > 1 and len(work) > 1:
        context = multiprocessing.get_context("spawn")
        with ProcessPoolExecutor(max_workers=options.jobs, mp_context=context) as pool:
            results: Iterable[tuple[str, str | None]] = pool.map(_render_worker, work)
            for label, error in results:
                failures += _report(label, error)
    else:
        for item in work:
            label, error = _render_worker(item)
            failures += _report(label, error)

    return 1 if failures else 0


def _report(label: str, error: str | None) -> int:
    if error:
        print(f"FAILED {label}: {error}", file=sys.stderr)
        return 1
    print(f"wrote {label}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
