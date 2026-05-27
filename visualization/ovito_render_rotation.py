"""Orbit camera and export GIFs (see `ovito_rotation_config` for radii and bonds)."""

from __future__ import annotations

import sys
from pathlib import Path

_viz_dir = Path(__file__).resolve().parent
if str(_viz_dir) not in sys.path:
    sys.path.insert(0, str(_viz_dir))

from collections.abc import Iterable
import math
import tempfile

from ovito.vis import TachyonRenderer, Viewport

from ovito_rotation_config import (
    CAMERA_LOOK_AT,
    CAMERA_ORBIT_POLAR_ANGLE_DEG,
    CAMERA_ORBIT_RADIUS,
    CAMERA_ORBIT_START_ANGLE_DEG,
    DUMP_PATH,
    TRAJECTORY_FRAMES,
    LAMMPS_TYPE_ID_TO_SYMBOL,
    ATOM_DISPLAY_RADII,
    BOND_PAIR_CUTOFFS,
    BOND_VISUAL_RADIUS,
    ROTATION_SUBFRAMES,
    GIF_BASENAME,
    GIF_FPS,
    KEEP_PNG_FRAMES,
    RENDER_SIZE,
)
from ovito_visual_pipeline import load_styled_pipeline


def normalize_trajectory_frames(spec: int | Iterable[int]) -> tuple[int, ...]:
    """Return a sorted tuple of unique trajectory indices."""
    if isinstance(spec, int):
        return (spec,)
    if isinstance(spec, Iterable) and not isinstance(spec, (str, bytes)):
        return tuple(sorted(set(int(i) for i in spec)))
    raise TypeError("TRAJECTORY_FRAMES must be an int or iterable of ints")


def gif_path_for_traj(traj_frame: int, multi_traj: bool) -> Path:
    if multi_traj:
        return Path(f"{GIF_BASENAME}_t{traj_frame:05d}.gif")
    return Path(f"{GIF_BASENAME}.gif")


def write_gif_from_png_paths(paths: list[Path], out_path: Path, *, fps: float) -> None:
    try:
        import imageio.v2 as imageio
    except ImportError as e:
        raise ImportError(
            "Writing GIF requires the `imageio` package. Install with: pip install imageio"
        ) from e

    dt = 1.0 / fps
    with imageio.get_writer(str(out_path), mode="I", duration=dt, loop=0) as writer:
        for p in paths:
            writer.append_data(imageio.imread(p))


def orbit_camera_ray(
    look_at: tuple[float, float, float],
    azimuth: float,
    orbit_radius: float,
    polar_rad: float,
) -> tuple[tuple[float, float, float], tuple[float, float, float]]:
    """
    Spherical orbit around look_at: polar_rad is colatitude from +Z; azimuth rotates in the ring.
    polar_rad = π/2 recovers horizontal motion in the z = lz plane.
    """
    lx, ly, lz = look_at
    s_horiz = orbit_radius * math.sin(polar_rad)
    px = lx + s_horiz * math.sin(azimuth)
    py = ly - s_horiz * math.cos(azimuth)
    pz = lz + orbit_radius * math.cos(polar_rad)

    dx = lx - px
    dy = ly - py
    dz = lz - pz
    n = math.hypot(math.hypot(dx, dy), dz)
    if n <= 1e-12:
        return (px, py, pz), (-1.0, 0.0, 0.0)
    inv = 1.0 / n
    return (px, py, pz), (dx * inv, dy * inv, dz * inv)


def render_rotation() -> None:
    pipeline = load_styled_pipeline(
        DUMP_PATH,
        atom_display_radii=ATOM_DISPLAY_RADII,
        pairwise_bond_cutoffs=BOND_PAIR_CUTOFFS,
        bond_visual_radius=BOND_VISUAL_RADIUS,
        lammps_type_id_to_symbol=LAMMPS_TYPE_ID_TO_SYMBOL,
    )
    pipeline.add_to_scene()

    vp = Viewport()
    vp.type = Viewport.Type.Perspective

    traj_frames = normalize_trajectory_frames(TRAJECTORY_FRAMES)
    multi_traj = len(traj_frames) > 1
    cwd = Path.cwd()

    for traj_frame in traj_frames:
        if not (0 <= traj_frame < pipeline.num_frames):
            raise ValueError(
                f"trajectory frame {traj_frame} out of range [0, {pipeline.num_frames - 1}]"
            )

        gif_out = cwd / gif_path_for_traj(traj_frame, multi_traj)

        def orbit_one_trajectory(tmpdir: Path | None) -> list[Path]:
            phase0 = math.radians(CAMERA_ORBIT_START_ANGLE_DEG)
            polar = math.radians(CAMERA_ORBIT_POLAR_ANGLE_DEG)
            paths: list[Path] = []
            for i_rot in range(ROTATION_SUBFRAMES):
                azimuth = phase0 + 2 * math.pi * i_rot / ROTATION_SUBFRAMES
                pos, direction = orbit_camera_ray(
                    CAMERA_LOOK_AT,
                    azimuth,
                    CAMERA_ORBIT_RADIUS,
                    polar,
                )
                vp.camera_pos = pos
                vp.camera_dir = direction

                if tmpdir is not None:
                    name = tmpdir / f"{i_rot:05d}.png"
                elif multi_traj:
                    name = cwd / f"frame_t{traj_frame:05d}_r{i_rot:04d}.png"
                else:
                    name = cwd / f"frame_{i_rot:04d}.png"

                vp.render_image(
                    size=RENDER_SIZE,
                    filename=str(name),
                    frame=traj_frame,
                    renderer=TachyonRenderer(),
                )
                paths.append(name)
            return paths

        if KEEP_PNG_FRAMES:
            ordered_pngs = orbit_one_trajectory(None)
            write_gif_from_png_paths(ordered_pngs, gif_out, fps=GIF_FPS)
        else:
            with tempfile.TemporaryDirectory(prefix="ovito_rot_") as tmp:
                ordered_pngs = orbit_one_trajectory(Path(tmp))
                write_gif_from_png_paths(ordered_pngs, gif_out, fps=GIF_FPS)

        print(f"Wrote {gif_out}")


if __name__ == "__main__":
    render_rotation()
