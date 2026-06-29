import math
import os
import sys
from collections import defaultdict

import numpy as np
import h5py
from scipy.spatial import cKDTree
from scipy.ndimage import uniform_filter

# --- Configuration ---
# DUMP_FILE, OUTPUT_FILE, HYDRONIUM_Y_CENTER, HYDRONIUM_Z_CENTER are CLI args
# Usage: python voxel_analysis.py <dump_file> <output_file> <y_center> <z_center>
VOXEL_SIZE  = 10.0      # Å (10 Å = 1 nm)

# Atom type IDs (1-based, must match LAMMPS dump)
SI_TYPE = 1
O_TYPE  = 2
H_TYPE  = 3

# Mass lookup array indexed by type int (index 0 unused)
TYPE_TO_MASS = np.array([0.0, 28.085, 15.999, 1.008])  # amu: Si, O, H

# Unit conversion constants (LAMMPS metal units)
# positions: Å, velocities: Å/ps, forces: eV/Å, masses: amu
K_B                 = 8.617333e-5   # eV/K
AMU_ANGS2_PS2_TO_EV = 1.0364e-4    # amu*(Å/ps)² -> eV (for kinetic energy)
EV_A3_TO_GPA        = 160.2176     # eV/Å³ -> GPa (for pressure)
AMU_A3_TO_G_CM3     = 1.6605       # amu/Å³ -> g/cm³ (for density)

# --- Hydronium detection ---
HYDRONIUM_VOXEL_X = 10.0   # Å, voxel depth along x
HYDRONIUM_VOXEL_Y = 50.0   # Å, voxel cross-section in y
HYDRONIUM_VOXEL_Z = 50.0   # Å, voxel cross-section in z
HYDRONIUM_PADDING =  1.3   # Å, padding around rod for H atom collection
OH_CUTOFF         =  1.2   # Å, O-H bond distance cutoff

# --- Jet tip detection ---
JET_TIP_SPEED_THRESHOLD = 30.0  # Å/ps, minimum avg_speed to count as jet front


def parse_header(file):
    file.readline()                        # skip ITEM: TIMESTEP
    timestep = int(file.readline())        # read timestep
    file.readline()                        # skip ITEM: NUMBER OF ATOMS
    N = int(file.readline())               # read number of atoms
    file.readline()                        # skip ITEM: BOX BOUNDS
    xlo, xhi = map(float, file.readline().split())
    ylo, yhi = map(float, file.readline().split())
    zlo, zhi = map(float, file.readline().split())
    file.readline()                        # skip ITEM: ATOMS header
    return timestep, N, xlo, xhi, ylo, yhi, zlo, zhi

def preallocate_hdf5(path, nx, ny, nz, attrs):

    h5file = h5py.File(path, 'w')

    scalar_datasets = ['density', 'pressure', 'virial_pressure', 'temperature', 'avg_speed', 'avg_O_speed', 'number_density']
    for name in scalar_datasets:
        h5file.create_dataset(
            name,
            shape=(nx, ny, nz),
            dtype=np.float32,
            chunks=(1, ny, nz),
            fillvalue=np.nan,
        )

    h5file.create_dataset(
        'voxel_type',
        shape=(nx, ny, nz),
        dtype=np.uint8,
        chunks=(1, ny, nz),
        fillvalue=0,
    )

    h5file.create_dataset(
        'v_COM',
        shape=(nx, ny, nz, 3),
        dtype=np.float32,
        chunks=(1, ny, nz, 3),
        fillvalue=np.nan,
    )

    h5file.create_dataset('si_surface',      shape=(ny, nz),                    dtype=np.float32, fillvalue=np.nan)
    h5file.create_dataset('si_surface_mask', shape=(nx, ny, nz),                dtype=np.uint8,   fillvalue=0)
    h5file.create_dataset('hydronium_count', shape=(attrs['n_hydronium_x'],),   dtype=np.int32,   fillvalue=0)
    h5file.create_dataset('jet_tip_x',       shape=(),                          dtype=np.float32)

    for key, val in attrs.items():
        h5file.attrs[key] = val

    return h5file


def process_voxel(arr, masses, V):
    """
    Parameters
    ----------
    arr : np.ndarray, shape (N, 11), float64
        Columns: id, type, x, y, z, vx, vy, vz, fx, fy, fz
    masses : np.ndarray, shape (N,), float64
        Per-atom masses looked up from TYPE_TO_MASS
    V : float
        Voxel volume in Å³

    Returns
    -------
    density        : float  (g/cm³)
    pressure       : float  (GPa, kinetic + virial)
    virial_pressure: float  (GPa, virial only)
    temperature    : float  (K)
    avg_speed      : float  (Å/ps)
    avg_O_speed    : float  (Å/ps, nan if no O atoms)
    voxel_type     : int    (1 = water, 2 = silica)
    v_COM          : np.ndarray shape (3,) float64  (Å/ps)
    """
    N = len(arr)
    types = arr[:, 1].astype(int)

    total_mass = masses.sum()
    density = (total_mass / V) * AMU_A3_TO_G_CM3

    v_COM = np.average(arr[:, 5:8], axis=0, weights=masses)

    v_th = arr[:, 5:8] - v_COM
    speeds_sq = np.sum(v_th**2, axis=1)
    thermal_KE = np.dot(masses, speeds_sq) * AMU_ANGS2_PS2_TO_EV
    # dof = 3N - 3: subtract 3 translational DOF removed by v_COM, matching LAMMPS compute temp/com
    temperature = thermal_KE / ((3 * N - 3) * K_B) if N > 1 else np.nan

    r_COM = np.average(arr[:, 2:5], axis=0, weights=masses)
    r_prime = arr[:, 2:5] - r_COM
    virial_sum = np.sum(r_prime * arr[:, 8:11])
    kinetic_term = np.dot(masses, np.sum(arr[:, 5:8]**2, axis=1)) * AMU_ANGS2_PS2_TO_EV
    pressure        = ((kinetic_term + virial_sum) / (3 * V)) * EV_A3_TO_GPA
    virial_pressure = (virial_sum / (3 * V)) * EV_A3_TO_GPA

    all_speeds = np.sqrt(np.sum(arr[:, 5:8]**2, axis=1))
    avg_speed = all_speeds.mean()

    o_mask = types == O_TYPE
    avg_O_speed = all_speeds[o_mask].mean() if o_mask.any() else np.nan

    voxel_type = 2 if SI_TYPE in types else 1
    number_density = N / V   # atoms/Å³, unnormalized — normalization done at visualization time

    return density, pressure, virial_pressure, temperature, avg_speed, avg_O_speed, voxel_type, v_COM, number_density


def flush_layer(layer_buf, ix, h5file, attrs):
    ny = attrs['ny']
    nz = attrs['nz']
    vs = attrs['voxel_size']

    box_x = attrs['xhi'] - attrs['xlo']
    box_y = attrs['yhi'] - attrs['ylo']
    box_z = attrs['zhi'] - attrs['zlo']

    # for edge voxels, use the actual distance to the edge of the box
    dx = min(vs, box_x - ix * vs)

    # accumulate results into layer-sized arrays, then write once per dataset
    out_density          = np.full((ny, nz), np.nan,  dtype=np.float32)
    out_pressure         = np.full((ny, nz), np.nan,  dtype=np.float32)
    out_virial_pressure  = np.full((ny, nz), np.nan,  dtype=np.float32)
    out_temperature      = np.full((ny, nz), np.nan,  dtype=np.float32)
    out_avg_speed        = np.full((ny, nz), np.nan,  dtype=np.float32)
    out_avg_O_speed      = np.full((ny, nz), np.nan,  dtype=np.float32)
    out_number_density   = np.full((ny, nz), np.nan,  dtype=np.float32)
    out_voxel_type       = np.zeros((ny, nz),          dtype=np.uint8)
    out_v_COM            = np.full((ny, nz, 3), np.nan, dtype=np.float32)

    for iy in range(ny):
        # for edge voxels, use the actual distance to the edge of the box
        dy = min(vs, box_y - iy * vs)

        for iz in range(nz):
            atoms = layer_buf.get((iy, iz))
            if not atoms:
                continue

            # for edge voxels, use the actual distance to the edge of the box
            dz = min(vs, box_z - iz * vs)
            V  = dx * dy * dz

            arr    = np.array(atoms, dtype=np.float64)
            masses = TYPE_TO_MASS[arr[:, 1].astype(int)]

            density, pressure, virial_pressure, temperature, avg_speed, avg_O_speed, voxel_type, v_COM, number_density = \
                process_voxel(arr, masses, V)

            out_density         [iy, iz]    = density
            out_pressure        [iy, iz]    = pressure
            out_virial_pressure [iy, iz]    = virial_pressure
            out_temperature     [iy, iz]    = temperature
            out_avg_speed       [iy, iz]    = avg_speed
            out_avg_O_speed     [iy, iz]    = avg_O_speed
            out_number_density  [iy, iz]    = number_density
            out_voxel_type      [iy, iz]    = voxel_type
            out_v_COM           [iy, iz, :] = v_COM

    # single write per dataset — one chunk I/O instead of ny*nz
    h5file['density'        ][ix] = out_density
    h5file['pressure'       ][ix] = out_pressure
    h5file['virial_pressure'][ix] = out_virial_pressure
    h5file['temperature'    ][ix] = out_temperature
    h5file['avg_speed'      ][ix] = out_avg_speed
    h5file['avg_O_speed'    ][ix] = out_avg_O_speed
    h5file['number_density' ][ix] = out_number_density
    h5file['voxel_type'     ][ix] = out_voxel_type
    h5file['v_COM'          ][ix] = out_v_COM


def streaming_loop(file, attrs, h5file, y_center, z_center):
    vs = attrs['voxel_size']
    xlo = attrs['xlo']
    ylo = attrs['ylo']
    zlo = attrs['zlo']
    nx = attrs['nx']
    ny = attrs['ny']
    nz = attrs['nz']

    # hydronium rod bounds (no PBC — absolute coordinates)
    rod_y_lo = y_center - HYDRONIUM_VOXEL_Y / 2.0
    rod_y_hi = y_center + HYDRONIUM_VOXEL_Y / 2.0
    rod_z_lo = z_center - HYDRONIUM_VOXEL_Z / 2.0
    rod_z_hi = z_center + HYDRONIUM_VOXEL_Z / 2.0
    pad_y_lo = rod_y_lo - HYDRONIUM_PADDING
    pad_y_hi = rod_y_hi + HYDRONIUM_PADDING
    pad_z_lo = rod_z_lo - HYDRONIUM_PADDING
    pad_z_hi = rod_z_hi + HYDRONIUM_PADDING

    si_surface = np.full((ny, nz), np.nan, dtype=np.float64)
    rod_o = []   # (x, y, z) of O atoms inside the rod cross-section
    rod_h = []   # (x, y, z) of H atoms inside the padded cross-section

    layer_buf = defaultdict(list)
    current_ix = None

    def _process_atom(fields):
        atype = int(fields[1])
        x, y, z = float(fields[2]), float(fields[3]), float(fields[4])

        ix = min(int((x - xlo) / vs), nx - 1)
        iy = min(int((y - ylo) / vs), ny - 1)
        iz = min(int((z - zlo) / vs), nz - 1)

        # si surface: track max-x Si atom per (y,z) bin
        if atype == SI_TYPE:
            if np.isnan(si_surface[iy, iz]) or x > si_surface[iy, iz]:
                si_surface[iy, iz] = x

        # hydronium rod collection
        if atype == O_TYPE and rod_y_lo <= y <= rod_y_hi and rod_z_lo <= z <= rod_z_hi:
            rod_o.append((x, y, z))
        elif atype == H_TYPE and pad_y_lo <= y <= pad_y_hi and pad_z_lo <= z <= pad_z_hi:
            rod_h.append((x, y, z))

        return ix, iy, iz

    fields = file.readline().split()
    ix, iy, iz = _process_atom(fields)
    current_ix = ix
    layer_buf[(iy, iz)].append([float(v) for v in fields])

    for line in file:
        fields = line.split()
        ix, iy, iz = _process_atom(fields)

        if ix != current_ix:
            flush_layer(layer_buf, current_ix, h5file, attrs)
            layer_buf = defaultdict(list)
            current_ix = ix

        layer_buf[(iy, iz)].append([float(v) for v in fields])

    if layer_buf:
        flush_layer(layer_buf, current_ix, h5file, attrs)

    return si_surface, rod_o, rod_h


def detect_hydronium(rod_o, rod_h, attrs):
    n_hx = attrs['n_hydronium_x']
    counts = np.zeros(n_hx, dtype=np.int32)

    if not rod_o or len(rod_h) < 3:
        return counts

    o_arr = np.array(rod_o, dtype=np.float64)   # (N_o, 3)
    h_arr = np.array(rod_h, dtype=np.float64)   # (N_h, 3)

    k = min(3, len(h_arr))
    tree = cKDTree(h_arr)
    dists, _ = tree.query(o_arr, k=k, distance_upper_bound=OH_CUTOFF)

    if k < 3:
        return counts

    is_hydronium = ~np.isinf(dists).any(axis=1)
    hydronium_x = o_arr[is_hydronium, 0]

    xlo = attrs['xlo']
    for ox in hydronium_x:
        ix = min(int((ox - xlo) / HYDRONIUM_VOXEL_X), n_hx - 1)
        counts[ix] += 1

    return counts


def detect_jet_tip(h5file, attrs):
    avg_speed  = h5file['avg_speed'][:]    # (nx, ny, nz)
    voxel_type = h5file['voxel_type'][:]   # (nx, ny, nz)

    water_above_threshold = (voxel_type == 1) & (avg_speed > JET_TIP_SPEED_THRESHOLD)
    x_layers_hit = np.any(water_above_threshold, axis=(1, 2))   # (nx,)

    if not x_layers_hit.any():
        return np.float32(np.nan)

    tip_ix = int(np.where(x_layers_hit)[0].max())
    return np.float32((tip_ix + 0.5) * attrs['voxel_size'])


SMOOTH_DATASETS = [
    'density', 'pressure', 'virial_pressure', 'temperature',
    'avg_speed', 'avg_O_speed', 'number_density',
]

SMOOTH_KERNEL_SIZE = 7  # box width per axis: 3 neighbors each side + itself

def _smooth_3d(arr):
    """NaN-aware uniform_filter smoothing. Empty voxels receive the weighted
    average of their neighbors; non-empty voxels include themselves."""
    nan_mask = np.isnan(arr)
    filled   = np.where(nan_mask, 0.0, arr)
    weights  = np.where(nan_mask, 0.0, 1.0)
    smoothed = uniform_filter(filled,   size=SMOOTH_KERNEL_SIZE, mode='constant', cval=0.0)
    counts   = uniform_filter(weights,  size=SMOOTH_KERNEL_SIZE, mode='constant', cval=0.0)
    return np.where(counts > 0, smoothed / counts, np.nan).astype(arr.dtype)


def smooth_voxel_data(h5file):
    for name in SMOOTH_DATASETS:
        arr = h5file[name][:]
        h5file[name][:] = _smooth_3d(arr)

    v_com = h5file['v_COM'][:]   # (nx, ny, nz, 3)
    for i in range(3):
        v_com[..., i] = _smooth_3d(v_com[..., i])
    h5file['v_COM'][:] = v_com


def main():
    dump_file   = sys.argv[1]
    output_file = sys.argv[2]
    y_center    = float(sys.argv[3])
    z_center    = float(sys.argv[4])

    with open(dump_file, "r") as f:
        timestep, N, xlo, xhi, ylo, yhi, zlo, zhi = parse_header(f)

        nx = math.ceil((xhi - xlo) / VOXEL_SIZE)
        ny = math.ceil((yhi - ylo) / VOXEL_SIZE)
        nz = math.ceil((zhi - zlo) / VOXEL_SIZE)
        n_hydronium_x = math.ceil((xhi - xlo) / HYDRONIUM_VOXEL_X)

        attrs = {
            'timestep': timestep,
            'N': N,
            'xlo': xlo, 'xhi': xhi,
            'ylo': ylo, 'yhi': yhi,
            'zlo': zlo, 'zhi': zhi,
            'voxel_size': VOXEL_SIZE,
            'nx': nx, 'ny': ny, 'nz': nz,
            'si_type': SI_TYPE,
            'o_type': O_TYPE,
            'h_type': H_TYPE,
            'n_hydronium_x': n_hydronium_x,
            'hydronium_y_center': y_center,
            'hydronium_z_center': z_center,
        }

        tmp_file = output_file + ".tmp"
        h5file = preallocate_hdf5(tmp_file, nx, ny, nz, attrs)

        si_surface, rod_o, rod_h = streaming_loop(f, attrs, h5file, y_center, z_center)

        h5file['si_surface'][:]      = si_surface.astype(np.float32)

        # mark the surface Si voxel for each (iy, iz) bin
        si_surface_mask = np.zeros((nx, ny, nz), dtype=np.uint8)
        vs = VOXEL_SIZE
        for iy in range(ny):
            for iz in range(nz):
                sx = si_surface[iy, iz]
                if not np.isnan(sx):
                    six = min(int((sx - attrs['xlo']) / vs), nx - 1)
                    si_surface_mask[six, iy, iz] = 1
        h5file['si_surface_mask'][:] = si_surface_mask

        h5file['hydronium_count'][:] = detect_hydronium(rod_o, rod_h, attrs)
        h5file['jet_tip_x'][()]      = detect_jet_tip(h5file, attrs)

        smooth_voxel_data(h5file)

        h5file.close()
        os.rename(tmp_file, output_file)


if __name__ == "__main__":
    main()
