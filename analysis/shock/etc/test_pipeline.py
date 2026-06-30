"""
End-to-end pipeline test with two synthetic frames on a 5x6x6 voxel system.

Box: 50 x 60 x 60 Å  (VOXEL_SIZE=10 Å)
Hydronium rod: y_center=30, z_center=30

Run:
    python test_pipeline.py

Produces test_output/ with per-frame h5 files and trajectory.h5,
then prints PASS/FAIL for each property.

Strategy: compute raw per-voxel values independently, apply the same
NaN-aware uniform_filter smoothing the pipeline uses, then compare
against the h5 output.  This lets the test verify both the raw
calculation and the smoothing in a single pass.
"""

import math
import os
import subprocess
import sys
import shutil

import numpy as np
import h5py
from scipy.ndimage import uniform_filter

TEST_DIR   = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.dirname(TEST_DIR)
OUT_DIR    = os.path.join(TEST_DIR, 'test_output')
DUMP_DIR   = os.path.join(OUT_DIR, 'dumps')
H5_DIR     = os.path.join(OUT_DIR, 'h5')
TRAJ_H5    = os.path.join(OUT_DIR, 'trajectory.h5')

Y_CENTER = 30.0
Z_CENTER = 30.0

NX, NY, NZ   = 5, 6, 6
VOXEL_SIZE   = 10.0
JET_THRESHOLD = 30.0

# unit constants (from voxel_analysis.py)
AMU_A3_TO_G_CM3     = 1.6605
AMU_ANGS2_PS2_TO_EV = 1.0364e-4
EV_A3_TO_GPA        = 160.2176
K_B                 = 8.617333e-5

O_MASS  = 15.999
H_MASS  = 1.008
SI_MASS = 28.085


# -----------------------------------------------------------------------
# dump generation

def write_dump(path, timestep, atoms):
    """atoms: list of (id, type, x, y, z, vx, vy, vz, fx, fy, fz), sorted by x"""
    with open(path, 'w') as f:
        f.write("ITEM: TIMESTEP\n")
        f.write(f"{timestep}\n")
        f.write("ITEM: NUMBER OF ATOMS\n")
        f.write(f"{len(atoms)}\n")
        f.write("ITEM: BOX BOUNDS pp pp pp\n")
        f.write("0.0 50.0\n")
        f.write("0.0 60.0\n")
        f.write("0.0 60.0\n")
        f.write("ITEM: ATOMS id type x y z vx vy vz fx fy fz\n")
        for a in atoms:
            f.write(" ".join(str(v) for v in a) + "\n")


# -----------------------------------------------------------------------
# atom layouts (sorted by x so streaming loop sees atoms in ix order)
#
# Frame 1 voxel assignments:
#   ix=0, iy=0, iz=0 : 1 Si                → silica, si_surface
#   ix=1, iy=2, iz=2 : O + 3H (1.2 Å)     → hydronium
#   ix=2, iy=2, iz=2 : 2 O, ±v, ±f        → temperature / pressure
#   ix=3, iy=3, iz=3 : 1 O, v=0            → placeholder (speed < 30)
#   ix=4, iy=3, iz=3 : O+H, v=40           → jet tip
#
# Frame 2 changes:
#   Si moved to x=8  → si_surface[0,0] = 8
#   ix=3 O speed=35  → jet tip moves to ix=3
#   ix=4 atoms v=0   → no longer jet tip

FRAME1_ATOMS = [
    # id  type   x      y      z     vx    vy  vz  fx   fy  fz
    (  1,   1,  5.0,  5.0,  5.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # Si  ix=0
    (  2,   2, 15.0, 25.0, 25.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # O   ix=1 hydronium
    (  3,   3, 15.0, 25.8, 25.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # H
    (  4,   3, 15.0, 25.0, 25.8,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # H
    (  5,   3, 15.8, 25.0, 25.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # H
    (  6,   2, 24.0, 25.0, 25.0,  10.0, 0.0, 0.0, 1.0, 0.0, 0.0),  # O   ix=2 temp/pressure
    (  7,   2, 26.0, 25.0, 25.0, -10.0, 0.0, 0.0,-1.0, 0.0, 0.0),  # O   ix=2
    (  8,   2, 35.0, 35.0, 35.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # O   ix=3 placeholder
    (  9,   2, 45.0, 35.0, 35.0,  40.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # O   ix=4 jet tip
    ( 10,   3, 45.9, 35.0, 35.0,  40.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # H   ix=4
]

FRAME2_ATOMS = [
    (  1,   1,  8.0,  5.0,  5.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # Si  ix=0 (x=8)
    (  2,   2, 15.0, 25.0, 25.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # O   ix=1 hydronium
    (  3,   3, 15.0, 25.8, 25.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    (  4,   3, 15.0, 25.0, 25.8,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    (  5,   3, 15.8, 25.0, 25.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
    (  6,   2, 24.0, 25.0, 25.0,  10.0, 0.0, 0.0, 1.0, 0.0, 0.0),  # O   ix=2
    (  7,   2, 26.0, 25.0, 25.0, -10.0, 0.0, 0.0,-1.0, 0.0, 0.0),
    (  8,   2, 35.0, 35.0, 35.0,  35.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # O   ix=3 speed=35 → jet tip
    (  9,   2, 45.0, 35.0, 35.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # O   ix=4 v=0
    ( 10,   3, 45.9, 35.0, 35.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # H   ix=4 v=0
]


# -----------------------------------------------------------------------
# raw voxel calculations (mirrors process_voxel exactly)

def _vcom(masses, velocities):
    total = sum(masses)
    return [sum(m * v[i] for m, v in zip(masses, velocities)) / total for i in range(3)]

def raw_density(masses, V=1000.0):
    return sum(masses) / V * AMU_A3_TO_G_CM3

def raw_temperature(masses, vel_th):
    N = len(masses)
    if N <= 1:
        return np.nan
    KE = sum(m * sum(vi**2 for vi in v) for m, v in zip(masses, vel_th))
    KE *= AMU_ANGS2_PS2_TO_EV
    return KE / ((3 * N - 3) * K_B)

def raw_pressure(masses, velocities, positions, forces, V=1000.0):
    vc = _vcom(masses, velocities)
    total = sum(masses)
    rc = [sum(m * p[i] for m, p in zip(masses, positions)) / total for i in range(3)]
    virial  = sum((p[i] - rc[i]) * f[i]
                  for p, f in zip(positions, forces) for i in range(3))
    KE = sum(m * sum(v[i]**2 for i in range(3)) for m, v in zip(masses, velocities))
    KE *= AMU_ANGS2_PS2_TO_EV
    return (KE + virial) / (3 * V) * EV_A3_TO_GPA, virial / (3 * V) * EV_A3_TO_GPA

def raw_avg_speed(velocities):
    return float(np.mean([math.sqrt(sum(v**2 for v in vs)) for vs in velocities]))

def raw_number_density(N, V=1000.0):
    return N / V


# -----------------------------------------------------------------------
# pre-computed raw values for the 5 populated voxels in frame 1
# (frame 2 is identical except (0,0,0) Si position and (3,3,3)/(4,3,3) speeds)

def frame_raw_values(si_x, v33_speed, v43_speed):
    """Return dict of (ix,iy,iz) -> dict of property -> raw value."""

    # voxel (0,0,0): 1 Si
    d = {}
    d[(0,0,0)] = dict(
        density         = raw_density([SI_MASS]),
        pressure        = 0.0,
        virial_pressure = 0.0,
        temperature     = np.nan,
        avg_speed       = 0.0,
        avg_O_speed     = np.nan,
        number_density  = raw_number_density(1),
        voxel_type      = 2,
        si_surface_x    = si_x,      # max-x Si in bin (iy=0, iz=0)
    )

    # voxel (1,2,2): O + 3H hydronium, all v=0
    masses = [O_MASS, H_MASS, H_MASS, H_MASS]
    vels   = [[0,0,0]]*4
    vc     = [0,0,0]
    d[(1,2,2)] = dict(
        density         = raw_density(masses),
        pressure        = 0.0,
        virial_pressure = 0.0,
        temperature     = raw_temperature(masses, vels),   # 0 K (all v=0)
        avg_speed       = 0.0,
        avg_O_speed     = 0.0,
        number_density  = raw_number_density(4),
        voxel_type      = 1,
    )

    # voxel (2,2,2): 2 O with opposite velocities and forces
    masses2 = [O_MASS, O_MASS]
    vels2   = [[10,0,0], [-10,0,0]]
    pos2    = [[24,25,25], [26,25,25]]
    frc2    = [[1,0,0], [-1,0,0]]
    vc2     = _vcom(masses2, vels2)
    vth2    = [[v[i]-vc2[i] for i in range(3)] for v in vels2]
    P, VP   = raw_pressure(masses2, vels2, pos2, frc2)
    d[(2,2,2)] = dict(
        density         = raw_density(masses2),
        pressure        = P,
        virial_pressure = VP,
        temperature     = raw_temperature(masses2, vth2),
        avg_speed       = raw_avg_speed(vels2),
        avg_O_speed     = raw_avg_speed(vels2),
        number_density  = raw_number_density(2),
        voxel_type      = 1,
    )

    # voxel (3,3,3): 1 O (N=1 → temperature=NaN; pressure uses raw velocity)
    v33_pos = [[35, 35, 35]]
    v33_vel = [[v33_speed, 0, 0]]
    v33_frc = [[0, 0, 0]]
    P33, VP33 = raw_pressure([O_MASS], v33_vel, v33_pos, v33_frc)
    d[(3,3,3)] = dict(
        density         = raw_density([O_MASS]),
        pressure        = P33,
        virial_pressure = VP33,
        temperature     = np.nan,
        avg_speed       = float(v33_speed),
        avg_O_speed     = float(v33_speed),
        number_density  = raw_number_density(1),
        voxel_type      = 1,
    )

    # voxel (4,3,3): O + H at same speed (v_COM = v43_speed → thermal v = 0 → T=0)
    masses4  = [O_MASS, H_MASS]
    v43_pos  = [[45, 35, 35], [45.9, 35, 35]]
    vels4    = [[v43_speed, 0, 0], [v43_speed, 0, 0]]
    frc4     = [[0, 0, 0], [0, 0, 0]]
    vc4      = _vcom(masses4, vels4)
    vth4     = [[v[i]-vc4[i] for i in range(3)] for v in vels4]
    P43, VP43 = raw_pressure(masses4, vels4, v43_pos, frc4)
    d[(4,3,3)] = dict(
        density         = raw_density(masses4),
        pressure        = P43,
        virial_pressure = VP43,
        temperature     = raw_temperature(masses4, vth4),
        avg_speed       = float(v43_speed),
        avg_O_speed     = float(v43_speed),
        number_density  = raw_number_density(2),
        voxel_type      = 1,
    )

    return d


# -----------------------------------------------------------------------
# smoothing (mirrors voxel_analysis._smooth_3d)

def smooth_3d(arr):
    nan_mask = np.isnan(arr)
    filled   = np.where(nan_mask, 0.0, arr)
    weights  = np.where(nan_mask, 0.0, 1.0)
    s = uniform_filter(filled,  size=7, mode='constant', cval=0.0)
    c = uniform_filter(weights, size=7, mode='constant', cval=0.0)
    return np.where(c > 0, s / c, np.nan).astype(np.float32)


def build_smoothed_grids(raw_vals):
    """Build (NX, NY, NZ) smoothed grids from a raw_vals dict."""
    PROPS = ['density', 'pressure', 'virial_pressure', 'temperature',
             'avg_speed', 'avg_O_speed', 'number_density']
    grids = {p: np.full((NX, NY, NZ), np.nan) for p in PROPS}

    for (ix, iy, iz), vals in raw_vals.items():
        for p in PROPS:
            if p in vals:
                grids[p][ix, iy, iz] = vals[p]

    return {p: smooth_3d(grids[p]) for p in PROPS}


# -----------------------------------------------------------------------
# pipeline execution

def run_pipeline():
    shutil.rmtree(OUT_DIR, ignore_errors=True)
    os.makedirs(DUMP_DIR)
    os.makedirs(H5_DIR)

    dump1 = os.path.join(DUMP_DIR, 'dump.1000')
    dump2 = os.path.join(DUMP_DIR, 'dump.2000')
    h5_1  = os.path.join(H5_DIR, 'output_1000.h5')
    h5_2  = os.path.join(H5_DIR, 'output_2000.h5')

    write_dump(dump1, 1000, FRAME1_ATOMS)
    write_dump(dump2, 2000, FRAME2_ATOMS)

    voxel_script = os.path.join(SCRIPT_DIR, 'voxel_analysis.py')
    merge_script = os.path.join(SCRIPT_DIR, 'merge_h5.py')

    python = '/Users/loganyamamoto/util/analysis/.venv/bin/python'

    for dump, h5 in [(dump1, h5_1), (dump2, h5_2)]:
        r = subprocess.run(
            [python, voxel_script, dump, h5, str(Y_CENTER), str(Z_CENTER)],
            capture_output=True, text=True
        )
        if r.returncode != 0:
            print(f"ERROR running voxel_analysis on {dump}:\n{r.stderr}")
            sys.exit(1)

    r = subprocess.run(
        [python, merge_script, H5_DIR, TRAJ_H5,
         '--initial-cutoff', '2.0', '--secondary-cutoff', '1.0',
         '--sphere-y-center', str(Y_CENTER), '--sphere-z-center', str(Z_CENTER)],
        capture_output=True, text=True
    )
    if r.returncode != 0:
        print(f"ERROR running merge_h5:\n{r.stderr}")
        sys.exit(1)
    print(r.stdout.strip())


# -----------------------------------------------------------------------
# verification

def chk(label, got, expected, tol=1e-3):
    if np.isnan(expected):
        ok = np.isnan(got)
        got_str = 'NaN' if np.isnan(got) else f'{got:.6f}'
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}: got={got_str}  expected=NaN")
    else:
        ok = abs(float(got) - float(expected)) <= tol
        print(f"  [{'PASS' if ok else 'FAIL'}] {label}: got={float(got):.6f}  expected={float(expected):.6f}")
    return ok


def verify():
    passes = fails = 0

    raw1 = frame_raw_values(si_x=5.0,  v33_speed=0.0,  v43_speed=40.0)
    raw2 = frame_raw_values(si_x=8.0,  v33_speed=35.0, v43_speed=0.0)
    sm1  = build_smoothed_grids(raw1)
    sm2  = build_smoothed_grids(raw2)

    SMOOTHED_PROPS = list(sm1.keys())
    VOXELS_TO_CHECK = [(0,0,0), (1,2,2), (2,2,2), (3,3,3), (4,3,3)]

    with h5py.File(TRAJ_H5, 'r') as f:

        for frame_idx, (label, sm, raw) in enumerate([
                ('Frame 0 (timestep=1000)', sm1, raw1),
                ('Frame 1 (timestep=2000)', sm2, raw2),
        ]):
            print(f"\n=== {label} — smoothed voxel properties ===")
            for vox in VOXELS_TO_CHECK:
                ix, iy, iz = vox
                if vox not in raw:
                    continue
                print(f"\n  Voxel {vox}:")

                # voxel_type (not smoothed)
                got = int(f['voxel_type'][frame_idx, ix, iy, iz])
                exp = raw[vox]['voxel_type']
                ok  = got == exp
                print(f"  [{'PASS' if ok else 'FAIL'}] voxel_type: got={got}  expected={exp}")
                passes += ok; fails += not ok

                # smoothed scalar properties
                for prop in SMOOTHED_PROPS:
                    if prop not in raw[vox]:
                        continue
                    got = float(f[prop][frame_idx, ix, iy, iz])
                    exp = float(sm[prop][ix, iy, iz])
                    ok  = chk(prop, got, exp)
                    passes += ok; fails += not ok

        # --- frame-level scalars ---
        print("\n=== Frame 0 scalars ===")
        ok = chk('jet_tip_x', float(f['jet_tip_x'][0]), 45.0)
        passes += ok; fails += not ok
        ok = chk('si_surface[0,0]', float(f['si_surface'][0, 0, 0]), 5.0)
        passes += ok; fails += not ok
        hc = int(f['hydronium_count'][0, 1])
        ok = hc == 1
        print(f"  [{'PASS' if ok else 'FAIL'}] hydronium_count[1]: got={hc}  expected=1")
        passes += ok; fails += not ok

        print("\n=== Frame 1 scalars ===")
        ok = chk('jet_tip_x', float(f['jet_tip_x'][1]), 35.0)
        passes += ok; fails += not ok
        ok = chk('si_surface[0,0]', float(f['si_surface'][1, 0, 0]), 8.0)
        passes += ok; fails += not ok
        hc = int(f['hydronium_count'][1, 1])
        ok = hc == 1
        print(f"  [{'PASS' if ok else 'FAIL'}] hydronium_count[1]: got={hc}  expected=1")
        passes += ok; fails += not ok

        # --- si_surface_mask ---
        print("\n=== si_surface_mask ===")
        # Frame 0: Si at x=5.0, y=5.0, z=5.0 → ix=0, iy=0, iz=0
        # Frame 1: Si at x=8.0, y=5.0, z=5.0 → ix=0, iy=0, iz=0
        for fi, si_x in [(0, 5.0), (1, 8.0)]:
            mask = f['si_surface_mask'][fi]  # (nx, ny, nz)
            ok = int(mask[0, 0, 0]) == 1
            print(f"  [{'PASS' if ok else 'FAIL'}] Frame {fi}: si_surface_mask[0,0,0]=={int(mask[0,0,0])} expected 1")
            passes += ok; fails += not ok
            total_set = int(mask.sum())
            ok = total_set == 1
            print(f"  [{'PASS' if ok else 'FAIL'}] Frame {fi}: total set voxels=={total_set} expected 1")
            passes += ok; fails += not ok

        # --- crater datasets ---
        print("\n=== crater datasets ===")
        ok = 'crater' in f
        print(f"  [{'PASS' if ok else 'FAIL'}] crater group exists")
        passes += ok; fails += not ok

        if 'crater' in f:
            cg = f['crater']

            ok = 'reference_x' in cg and cg['reference_x'].shape == (2,)
            print(f"  [{'PASS' if ok else 'FAIL'}] crater/reference_x shape (2,)")
            passes += ok; fails += not ok

            ok = 'depth_map' in cg and cg['depth_map'].shape == (2, NY, NZ)
            print(f"  [{'PASS' if ok else 'FAIL'}] crater/depth_map shape (2,{NY},{NZ})")
            passes += ok; fails += not ok

            ok = 'sphere_center' in cg and cg['sphere_center'].shape == (2, 3)
            print(f"  [{'PASS' if ok else 'FAIL'}] crater/sphere_center shape (2,3)")
            passes += ok; fails += not ok

            ok = 'sphere_radius' in cg and cg['sphere_radius'].shape == (2,)
            print(f"  [{'PASS' if ok else 'FAIL'}] crater/sphere_radius shape (2,)")
            passes += ok; fails += not ok

            # with only 1 Si bin and cutoff=2.0, no crater detected → reference_x valid, no sphere
            ref_x = cg['reference_x'][()]
            ok = not np.isnan(ref_x[0]) and not np.isnan(ref_x[1])
            print(f"  [{'PASS' if ok else 'FAIL'}] reference_x not NaN: {ref_x}")
            passes += ok; fails += not ok

            ok = float(ref_x[0]) == pytest_approx(5.0, 1e-3) if False else abs(float(ref_x[0]) - 5.0) < 0.1
            print(f"  [{'PASS' if ok else 'FAIL'}] reference_x[0]≈5.0: got {float(ref_x[0]):.3f}")
            passes += ok; fails += not ok

            ok = abs(float(ref_x[1]) - 8.0) < 0.1
            print(f"  [{'PASS' if ok else 'FAIL'}] reference_x[1]≈8.0: got {float(ref_x[1]):.3f}")
            passes += ok; fails += not ok

            # sphere fit should be NaN (no crater bins with depth > 2.0)
            sr = cg['sphere_radius'][()]
            ok = np.isnan(float(sr[0])) and np.isnan(float(sr[1]))
            print(f"  [{'PASS' if ok else 'FAIL'}] sphere_radius NaN for both frames (no crater): {sr}")
            passes += ok; fails += not ok

    print(f"\n{'='*40}")
    print(f"Results: {passes} PASS  {fails} FAIL")
    return fails == 0


if __name__ == '__main__':
    print("Running pipeline...")
    run_pipeline()
    print(f"Output written to {OUT_DIR}")
    ok = verify()
    sys.exit(0 if ok else 1)
