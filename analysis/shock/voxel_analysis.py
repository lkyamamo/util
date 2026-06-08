import math
from collections import defaultdict

import numpy as np
import h5py

# --- Configuration ---
DUMP_FILE   = "shock_analysis_test.dump"
OUTPUT_FILE = "output.h5"
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

    scalar_datasets = ['density', 'pressure', 'virial_pressure', 'temperature', 'avg_speed', 'avg_O_speed']
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

    return density, pressure, virial_pressure, temperature, avg_speed, avg_O_speed, voxel_type, v_COM


def flush_layer(layer_buf, ix, h5file, attrs):
    ny = attrs['ny']
    nz = attrs['nz']
    vs = attrs['voxel_size']
    V  = vs ** 3

    for iy in range(ny):
        for iz in range(nz):
            atoms = layer_buf.get((iy, iz))
            if not atoms:
                continue

            arr    = np.array(atoms, dtype=np.float64)
            masses = TYPE_TO_MASS[arr[:, 1].astype(int)]

            density, pressure, virial_pressure, temperature, avg_speed, avg_O_speed, voxel_type, v_COM = \
                process_voxel(arr, masses, V)

            h5file['density'         ][ix, iy, iz]    = np.float32(density)
            h5file['pressure'        ][ix, iy, iz]    = np.float32(pressure)
            h5file['virial_pressure' ][ix, iy, iz]    = np.float32(virial_pressure)
            h5file['temperature'     ][ix, iy, iz]    = np.float32(temperature)
            h5file['avg_speed'       ][ix, iy, iz]    = np.float32(avg_speed)
            h5file['avg_O_speed'     ][ix, iy, iz]    = np.float32(avg_O_speed)
            h5file['voxel_type'      ][ix, iy, iz]    = np.uint8(voxel_type)
            h5file['v_COM'           ][ix, iy, iz, :] = v_COM.astype(np.float32)


def streaming_loop(file, attrs, h5file):

    vs = attrs['voxel_size']
    xlo = attrs['xlo']
    ylo = attrs['ylo']
    zlo = attrs['zlo']
    nx = attrs['nx']
    ny = attrs['ny']
    nz = attrs['nz']

    layer_buf = defaultdict(list)
    current_ix = None
    
    fields = file.readline().split()
    x, y, z = float(fields[2]), float(fields[3]), float(fields[4])

    ix = min(int((x - xlo) / vs), nx - 1)
    iy = min(int((y - ylo) / vs), ny - 1)
    iz = min(int((z - zlo) / vs), nz - 1)

    current_ix = ix

    layer_buf[(iy, iz)].append([float(v) for v in fields])

    for line in file:
        fields = line.split()
        x, y, z = float(fields[2]), float(fields[3]), float(fields[4])

        ix = min(int((x - xlo) / vs), nx - 1)
        iy = min(int((y - ylo) / vs), ny - 1)
        iz = min(int((z - zlo) / vs), nz - 1)

        if ix != current_ix:
            flush_layer(layer_buf, current_ix, h5file, attrs)
            layer_buf = defaultdict(list)
            current_ix = ix

        layer_buf[(iy, iz)].append([float(v) for v in fields])

    if layer_buf:
        flush_layer(layer_buf, current_ix, h5file, attrs)



def main():
    with open(DUMP_FILE, "r") as f:
        timestep, N, xlo, xhi, ylo, yhi, zlo, zhi = parse_header(f)

        nx = math.ceil((xhi - xlo) / VOXEL_SIZE)
        ny = math.ceil((yhi - ylo) / VOXEL_SIZE)
        nz = math.ceil((zhi - zlo) / VOXEL_SIZE)

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
            'h_type': H_TYPE
        }

        h5file = preallocate_hdf5(OUTPUT_FILE, nx, ny, nz, attrs)
        streaming_loop(f, attrs, h5file)
        h5file.close()


if __name__ == "__main__":
    main()

