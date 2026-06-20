import argparse
import glob
import re
import sys

import numpy as np
import h5py
from scipy.optimize import least_squares


def fit_sphere(points):
    """Least-squares sphere fit to (N, 3) points. Returns (cx, cy, cz, R)."""
    A = np.column_stack([-2.0 * points, np.ones(len(points))])
    b = -(points ** 2).sum(axis=1)
    result, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
    cx, cy, cz = result[:3]
    R = float(np.sqrt(max(cx**2 + cy**2 + cz**2 - result[3], 0.0)))
    return cx, cy, cz, R


def crater_analysis(out, initial_cutoff, secondary_cutoff):
    si_surface = out['si_surface'][-1][:]    # (ny, nz), last frame
    voxel_size = float(out.attrs['voxel_size'])
    ylo        = float(out.attrs['ylo'])
    zlo        = float(out.attrs['zlo'])
    ny, nz     = si_surface.shape

    valid = ~np.isnan(si_surface)
    if not valid.any():
        print("WARNING: si_surface is entirely NaN — skipping crater analysis", file=sys.stderr)
        return

    x_vals = si_surface[valid]

    # reference plane: mean x of undisturbed surface (above median - initial_cutoff)
    x_median  = float(np.median(x_vals))
    undisturbed_mask = valid & (si_surface > x_median - initial_cutoff)
    x_ref = float(si_surface[undisturbed_mask].mean())

    # initial crater candidates: bins that dip below the reference plane
    crater_mask = valid & (si_surface < x_ref - initial_cutoff)
    if not crater_mask.any():
        print("WARNING: no crater detected with initial_cutoff — skipping crater analysis", file=sys.stderr)
        return

    iy_idx, iz_idx = np.where(crater_mask)
    y_coords = (iy_idx + 0.5) * voxel_size + ylo
    z_coords = (iz_idx + 0.5) * voxel_size + zlo
    x_coords = si_surface[crater_mask]
    points   = np.column_stack([x_coords, y_coords, z_coords])

    # first sphere fit
    cx, cy, cz, R = fit_sphere(points)

    # refine: keep points within secondary_cutoff of sphere surface
    dist_to_sphere = np.abs(np.linalg.norm(points - np.array([cx, cy, cz]), axis=1) - R)
    refined = dist_to_sphere < secondary_cutoff

    if refined.sum() >= 4:
        cx, cy, cz, R = fit_sphere(points[refined])
    else:
        print("WARNING: fewer than 4 points after secondary_cutoff — using initial sphere fit", file=sys.stderr)

    # depth: how far the deepest point of the sphere sits below the reference plane
    depth = float(x_ref - (cx - R))

    # crater_radius: radius of the circle where the sphere intersects the reference plane (x = x_ref)
    d_center_to_ref = x_ref - cx
    crater_radius = float(np.sqrt(max(R**2 - d_center_to_ref**2, 0.0))) if abs(d_center_to_ref) < R else 0.0

    # 2D surface and depth map
    depth_map = np.where(valid, x_ref - si_surface, np.nan).astype(np.float32)

    grp = out.require_group('crater')
    grp.create_dataset('surface_2d',    data=si_surface.astype(np.float32))
    grp.create_dataset('depth_map',     data=depth_map)
    grp.create_dataset('sphere_center', data=np.array([cx, cy, cz], dtype=np.float32))
    grp.create_dataset('sphere_radius', data=np.float32(R))
    grp.create_dataset('depth',         data=np.float32(depth))
    grp.create_dataset('crater_radius', data=np.float32(crater_radius))
    grp.attrs['reference_x']       = x_ref
    grp.attrs['initial_cutoff']    = initial_cutoff
    grp.attrs['secondary_cutoff']  = secondary_cutoff

    print(f"Crater: depth={depth:.2f} Å  sphere_radius={R:.2f} Å  crater_radius={crater_radius:.2f} Å  center=({cx:.1f}, {cy:.1f}, {cz:.1f})")


def merge(output_dir, final_h5, initial_cutoff, secondary_cutoff):
    files = sorted(glob.glob(f"{output_dir}/output_*.h5"),
                   key=lambda f: int(re.search(r'output_(\d+)\.h5$', f).group(1)))
    if not files:
        raise FileNotFoundError(f"No output_*.h5 files found in {output_dir}")

    T = len(files)
    print(f"Merging {T} files -> {final_h5}")

    with h5py.File(final_h5, 'w') as out:
        for i, path in enumerate(files):
            with h5py.File(path, 'r') as src:
                for name in src:
                    data = src[name][:]
                    if i == 0:
                        shape = (T,) + data.shape
                        out.create_dataset(name, shape=shape, dtype=data.dtype)
                        for key, val in src.attrs.items():
                            out.attrs[key] = val
                    out[name][i] = data
            if (i + 1) % 50 == 0 or (i + 1) == T:
                print(f"  {i + 1}/{T}")

        crater_analysis(out, initial_cutoff, secondary_cutoff)

    print("Done.")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("output_dir")
    p.add_argument("final_h5")
    p.add_argument("--initial-cutoff",   type=float, required=True,
                   help="Å below reference plane to classify a bin as crater (first pass)")
    p.add_argument("--secondary-cutoff", type=float, required=True,
                   help="Å from fitted sphere surface to keep a point (second pass)")
    args = p.parse_args()

    merge(args.output_dir, args.final_h5, args.initial_cutoff, args.secondary_cutoff)
