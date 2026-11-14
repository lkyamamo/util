"""
End-to-end pipeline test with two synthetic frames on a 5x6x6 voxel system.

Box: 50 x 60 x 60 Å  (VOXEL_SIZE=10 Å)
Hydronium rod: y_center=30, z_center=30

Run:
    pytest shock/etc/test_pipeline.py   (from analysis/)

Writes the dumps, per-frame h5 files and trajectory.h5 to a temporary
directory, runs voxel_analysis.py and merge_h5.py on them, and checks every
property.

Strategy: compute raw per-voxel values independently, apply the same
NaN-aware uniform_filter smoothing the pipeline uses, then compare
against the h5 output.  This lets the test verify both the raw
calculation and the smoothing in a single pass.
"""

import math
import os

import numpy as np
import h5py
import pytest
from scipy.ndimage import uniform_filter

from testkit import load_script, run_script, write_dump as write_lammps_dump

TEST_DIR   = os.path.dirname(os.path.abspath(__file__))
SCRIPT_DIR = os.path.dirname(TEST_DIR)
VOXEL_SCRIPT = os.path.join(SCRIPT_DIR, 'voxel_analysis.py')
MERGE_SCRIPT = os.path.join(SCRIPT_DIR, 'merge_h5.py')

Y_CENTER = 30.0
Z_CENTER = 30.0
BUBBLE_RADIUS = 10.0   # jet tip rod half-width; both ix=3,4 atoms at y,z=35,35
                        # (distance 5 from center) are well within this rod

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

COLUMNS = ("id", "type", "x", "y", "z", "vx", "vy", "vz", "fx", "fy", "fz")


def write_dump(path, timestep, atoms):
    """atoms: list of (id, type, x, y, z, vx, vy, vz, fx, fy, fz), sorted by x"""
    table = list(zip(*atoms))
    frame = {c: table[i] for i, c in enumerate(COLUMNS)}
    write_lammps_dump(path, [frame], ((0.0, 0.0, 0.0), (50.0, 60.0, 60.0)), COLUMNS,
                      timesteps=[timestep])


# -----------------------------------------------------------------------
# atom layouts (sorted by x so streaming loop sees atoms in ix order)
#
# Frame 1 voxel assignments:
#   ix=0, iy=0, iz=0 : 1 Si                → silica, si_surface
#   ix=1, iy=2, iz=2 : O + 3H (1.2 Å)     → hydronium
#   ix=2, iy=2, iz=2 : 2 O, ±v, ±f        → temperature / pressure
#   ix=3, iy=3, iz=3 : 1 O, v=35           → low density, not a jet tip candidate
#   ix=4, iy=3, iz=3 : O+H, v=40           → low density, not a jet tip candidate
# detect_jet_tip is now density-based (JET_TIP_DENSITY_THRESHOLD=0.15 g/cm^3),
# and runs on SMOOTHED density. This test grid is sparse (5 populated voxels
# in a 5x6x6 grid) with a 7-wide smoothing kernel that spans the whole grid,
# so every voxel smooths to roughly the same low density (~0.037 g/cm^3) —
# well below threshold. Both frames should report jet_tip_x = NaN, which
# confirms isolated low-density atoms no longer register as a false jet tip.
# The rod-restriction / smallest-x / density-threshold logic itself is
# covered by a focused unit test on detect_jet_tip below.
#
# Frame 2 changes:
#   Si moved to x=8  → si_surface[0,0] = 8
#   ix=3 O speed=35  → unchanged (no longer relevant — density based now)
#   ix=4 atoms v=0   → unchanged

FRAME1_ATOMS = [
    # id  type   x      y      z     vx    vy  vz  fx   fy  fz
    (  1,   1,  5.0,  5.0,  5.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # Si  ix=0
    (  2,   2, 15.0, 25.0, 25.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # O   ix=1 hydronium
    (  3,   3, 15.0, 25.8, 25.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # H
    (  4,   3, 15.0, 25.0, 25.8,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # H
    (  5,   3, 15.8, 25.0, 25.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # H
    (  6,   2, 24.0, 25.0, 25.0,  10.0, 0.0, 0.0, 1.0, 0.0, 0.0),  # O   ix=2 temp/pressure
    (  7,   2, 26.0, 25.0, 25.0, -10.0, 0.0, 0.0,-1.0, 0.0, 0.0),  # O   ix=2
    (  8,   2, 35.0, 35.0, 35.0,  35.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # O   ix=3 jet tip candidate
    (  9,   2, 45.0, 35.0, 35.0,  40.0, 0.0, 0.0, 0.0, 0.0, 0.0),  # O   ix=4 jet tip candidate
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

@pytest.fixture(scope="module")
def traj(tmp_path_factory):
    """Run voxel_analysis on both frames and merge_h5 on the results; yield trajectory.h5."""
    out = tmp_path_factory.mktemp("shock")
    h5_dir = out / "h5"
    h5_dir.mkdir()
    for step, atoms in [(1000, FRAME1_ATOMS), (2000, FRAME2_ATOMS)]:
        dump = out / f"dump.{step}"
        write_dump(dump, step, atoms)
        run_script(VOXEL_SCRIPT, dump, h5_dir / f"output_{step}.h5",
                   Y_CENTER, Z_CENTER, BUBBLE_RADIUS)
    traj_h5 = out / "trajectory.h5"
    run_script(MERGE_SCRIPT, h5_dir, traj_h5,
               "--initial-cutoff", "2.0", "--secondary-cutoff", "1.0",
               "--sphere-y-center", Y_CENTER, "--sphere-z-center", Z_CENTER)
    with h5py.File(traj_h5, "r") as f:
        yield f


# -----------------------------------------------------------------------
# tests

FRAMES = [
    (frame_raw_values(si_x=5.0, v33_speed=35.0, v43_speed=40.0)),
    (frame_raw_values(si_x=8.0, v33_speed=35.0, v43_speed=0.0)),
]
VOXELS_TO_CHECK = [(0,0,0), (1,2,2), (2,2,2), (3,3,3), (4,3,3)]


@pytest.mark.parametrize("frame_idx", [0, 1])
def test_smoothed_voxel_properties(traj, checker, frame_idx):
    raw = FRAMES[frame_idx]
    sm = build_smoothed_grids(raw)
    for vox in VOXELS_TO_CHECK:
        ix, iy, iz = vox
        checker.chk(f"{vox} voxel_type", int(traj['voxel_type'][frame_idx, ix, iy, iz]),
                    raw[vox]['voxel_type'])
        for prop in sm:
            if prop in raw[vox]:
                checker.chk(f"{vox} {prop}", float(traj[prop][frame_idx, ix, iy, iz]),
                            float(sm[prop][ix, iy, iz]))


@pytest.mark.parametrize("frame_idx,si_x", [(0, 5.0), (1, 8.0)])
def test_frame_scalars(traj, checker, frame_idx, si_x):
    checker.chk("jet_tip_x", float(traj['jet_tip_x'][frame_idx]), float('nan'))
    checker.chk("si_surface[0,0]", float(traj['si_surface'][frame_idx, 0, 0]), si_x)
    checker.chk("hydronium_count[1]", int(traj['hydronium_count'][frame_idx, 1]), 1)
    # Si at (si_x, 5, 5) is in voxel (0, 0, 0), the only surface voxel.
    mask = traj['si_surface_mask'][frame_idx]
    checker.chk("si_surface_mask[0,0,0]", int(mask[0, 0, 0]), 1)
    checker.chk("si_surface_mask total", int(mask.sum()), 1)


def test_crater_datasets(traj, checker):
    assert 'crater' in traj
    cg = traj['crater']
    checker.chk("reference_x shape", cg['reference_x'].shape, (2,))
    checker.chk("depth_map shape", cg['depth_map'].shape, (2, NY, NZ))
    checker.chk("sphere_center shape", cg['sphere_center'].shape, (2, 3))
    checker.chk("sphere_radius shape", cg['sphere_radius'].shape, (2,))
    # One Si bin and cutoff 2.0: no crater, so reference_x is valid and there is no sphere.
    ref_x = cg['reference_x'][()]
    checker.chk("reference_x[0]", float(ref_x[0]), 5.0, tol=0.1)
    checker.chk("reference_x[1]", float(ref_x[1]), 8.0, tol=0.1)
    sr = cg['sphere_radius'][()]
    checker.chk("sphere_radius[0]", float(sr[0]), float('nan'))
    checker.chk("sphere_radius[1]", float(sr[1]), float('nan'))


def test_jet_tip_detection():
    """Rod restriction, density threshold and smallest-x selection on detect_jet_tip directly.

      ix=0, y,z=(5,5)   outside rod:  density 1.0  -> excluded by the rod
      ix=1, y,z=(35,35) inside rod:   density 0.05 -> excluded by the density threshold
      ix=2, y,z=(35,35) inside rod:   density 1.0  -> valid, smallest qualifying
      ix=3, y,z=(35,35) inside rod:   density 1.0  -> valid, but larger ix
    Expected tip: ix=2, at (2 + 0.5) * 10 = 25 A.
    """
    va = load_script(VOXEL_SCRIPT)
    nx, ny, nz, vs = 5, 6, 6, 10.0
    attrs = {'nx': nx, 'ny': ny, 'nz': nz, 'voxel_size': vs,
             'xlo': 0.0, 'ylo': 0.0, 'zlo': 0.0}
    density = np.zeros((nx, ny, nz), dtype=np.float32)
    voxel_type = np.zeros((nx, ny, nz), dtype=np.uint8)
    voxel_type[0, 0, 0] = 1; density[0, 0, 0] = 1.0
    voxel_type[1, 3, 3] = 1; density[1, 3, 3] = 0.05
    voxel_type[2, 3, 3] = 1; density[2, 3, 3] = 1.0
    voxel_type[3, 3, 3] = 1; density[3, 3, 3] = 1.0
    with h5py.File('jet_tip_unit_test.h5', 'w', driver='core', backing_store=False) as h5:
        h5.create_dataset('density', data=density)
        h5.create_dataset('voxel_type', data=voxel_type)
        tip_x = va.detect_jet_tip(h5, attrs, y_center=30.0, z_center=30.0, bubble_radius=10.0)
    assert float(tip_x) == pytest.approx((2 + 0.5) * vs, abs=1e-3)
