import argparse
import glob
import re
import sys

import numpy as np
import h5py


def fit_sphere_constrained(points, cy, cz):
    """Least-squares sphere fit with fixed y/z center. Returns (cx, R).

    For each point (xi, yi, zi):
        (xi - cx)^2 + di^2 = R^2   where di^2 = (yi-cy)^2 + (zi-cz)^2
    Rearranging:  -2*cx*xi + (cx^2 - R^2) = -(xi^2 + di^2)
    Linear system: A @ [cx, cx^2-R^2]^T = rhs, solved via lstsq.
    """
    x = points[:, 0]
    di2 = (points[:, 1] - cy) ** 2 + (points[:, 2] - cz) ** 2
    A = np.column_stack([-2.0 * x, np.ones(len(x))])
    rhs = -(x ** 2 + di2)
    result, _, _, _ = np.linalg.lstsq(A, rhs, rcond=None)
    cx = result[0]
    R = float(np.sqrt(max(cx ** 2 - result[1], 0.0)))
    return cx, R


def crater_analysis(out, initial_cutoff, secondary_cutoff, sphere_cy, sphere_cz):
    si_surface    = out['si_surface']       # (T, ny, nz)
    si_surf_mask  = out['si_surface_mask']  # (T, nx, ny, nz)
    T, ny, nz     = si_surface.shape
    nx            = si_surf_mask.shape[1]
    voxel_size    = float(out.attrs['voxel_size'])
    ylo           = float(out.attrs['ylo'])
    zlo           = float(out.attrs['zlo'])

    ref_x_arr       = np.full(T,          np.nan, dtype=np.float32)
    depth_map_arr   = np.full((T, ny, nz), np.nan, dtype=np.float32)
    crater_mask_arr = np.zeros((T, nx, ny, nz),    dtype=np.uint8)
    sc_arr          = np.full((T, 3),     np.nan, dtype=np.float32)
    sr_arr          = np.full(T,          np.nan, dtype=np.float32)
    depth_arr       = np.full(T,          np.nan, dtype=np.float32)
    cr_arr          = np.full(T,          np.nan, dtype=np.float32)

    for t in range(T):
        ss = si_surface[t][:]   # (ny, nz)
        valid = ~np.isnan(ss)
        if not valid.any():
            continue

        x_vals   = ss[valid]
        x_median = float(np.median(x_vals))
        undist   = valid & (ss > x_median - initial_cutoff)
        if not undist.any():
            continue
        x_ref = float(ss[undist].mean())
        ref_x_arr[t] = x_ref

        depth_map_arr[t] = np.where(valid, x_ref - ss, np.nan).astype(np.float32)

        crater_bins = valid & (ss < x_ref - initial_cutoff)
        if not crater_bins.any():
            continue

        # mark surface Si voxels that fall in crater bins
        sm = si_surf_mask[t][:]  # (nx, ny, nz)
        cm = np.zeros((nx, ny, nz), dtype=np.uint8)
        iy_idx, iz_idx = np.where(crater_bins)
        for iy, iz in zip(iy_idx, iz_idx):
            col = sm[:, iy, iz]
            hits = np.where(col == 1)[0]
            if hits.size:
                cm[hits[0], iy, iz] = 1
        crater_mask_arr[t] = cm

        # constrained sphere fit (fixed cy=sphere_cy, cz=sphere_cz)
        y_coords = (iy_idx + 0.5) * voxel_size + ylo
        z_coords = (iz_idx + 0.5) * voxel_size + zlo
        x_coords = ss[crater_bins]
        points   = np.column_stack([x_coords, y_coords, z_coords])

        if len(points) < 2:
            print(f"  t={t}: fewer than 2 crater points — skipping sphere fit", file=sys.stderr)
            continue

        cx, R = fit_sphere_constrained(points, sphere_cy, sphere_cz)

        dist = np.abs(np.sqrt((points[:, 0] - cx) ** 2
                              + (points[:, 1] - sphere_cy) ** 2
                              + (points[:, 2] - sphere_cz) ** 2) - R)
        refined = dist < secondary_cutoff
        if refined.sum() >= 2:
            cx, R = fit_sphere_constrained(points[refined], sphere_cy, sphere_cz)
        else:
            print(f"  t={t}: fewer than 2 points after secondary_cutoff — using initial fit",
                  file=sys.stderr)

        depth         = float(x_ref - (cx - R))
        d_ctr         = x_ref - cx
        crater_radius = float(np.sqrt(max(R ** 2 - d_ctr ** 2, 0.0))) if abs(d_ctr) < R else 0.0

        sc_arr[t]    = [cx, sphere_cy, sphere_cz]
        sr_arr[t]    = R
        depth_arr[t] = depth
        cr_arr[t]    = crater_radius

        print(f"  t={t}: depth={depth:.2f} Å  R={R:.2f} Å  crater_r={crater_radius:.2f} Å  cx={cx:.1f}")

    grp = out.require_group('crater')
    grp.create_dataset('reference_x',    data=ref_x_arr)
    grp.create_dataset('depth_map',      data=depth_map_arr)
    grp.create_dataset('si_crater_mask', data=crater_mask_arr)
    grp.create_dataset('sphere_center',  data=sc_arr)
    grp.create_dataset('sphere_radius',  data=sr_arr)
    grp.create_dataset('depth',          data=depth_arr)
    grp.create_dataset('crater_radius',  data=cr_arr)
    grp.attrs['initial_cutoff']   = initial_cutoff
    grp.attrs['secondary_cutoff'] = secondary_cutoff
    grp.attrs['sphere_cy']        = sphere_cy
    grp.attrs['sphere_cz']        = sphere_cz

    print(f"Crater analysis complete for {T} frames.")


def merge(output_dir, final_h5, initial_cutoff, secondary_cutoff, sphere_cy, sphere_cz):
    files = sorted(glob.glob(f"{output_dir}/output_*.h5"),
                   key=lambda f: int(re.search(r'output_(\d+)\.h5$', f).group(1)))
    if not files:
        raise FileNotFoundError(f"No output_*.h5 files found in {output_dir}")

    T = len(files)
    print(f"Merging {T} files -> {final_h5}")

    with h5py.File(final_h5, 'w') as out:
        out.create_dataset('timestep', shape=(T,), dtype=np.int64)

        for i, path in enumerate(files):
            with h5py.File(path, 'r') as src:
                for name in src:
                    data = src[name][()]
                    if i == 0:
                        shape = (T,) + data.shape
                        out.create_dataset(name, shape=shape, dtype=data.dtype)
                        for key, val in src.attrs.items():
                            out.attrs[key] = val
                    out[name][i] = data
                out['timestep'][i] = src.attrs['timestep']
            if (i + 1) % 50 == 0 or (i + 1) == T:
                print(f"  {i + 1}/{T}")

        crater_analysis(out, initial_cutoff, secondary_cutoff, sphere_cy, sphere_cz)

    print("Done.")


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("output_dir")
    p.add_argument("final_h5")
    p.add_argument("--initial-cutoff",   type=float, required=True,
                   help="Å below reference plane to classify a bin as crater (first pass)")
    p.add_argument("--secondary-cutoff", type=float, required=True,
                   help="Å from fitted sphere surface to keep a point (second pass)")
    p.add_argument("--sphere-y-center",  type=float, required=True,
                   help="Fixed y-center for constrained sphere fit (Å)")
    p.add_argument("--sphere-z-center",  type=float, required=True,
                   help="Fixed z-center for constrained sphere fit (Å)")
    args = p.parse_args()

    merge(args.output_dir, args.final_h5,
          args.initial_cutoff, args.secondary_cutoff,
          args.sphere_y_center, args.sphere_z_center)
