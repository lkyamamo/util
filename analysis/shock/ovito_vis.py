"""
Ovito visualizer: LAMMPS dump atoms + h5 voxel grid in one scene.

Run via:
    ovito --script ovito_vis.py -- <dump_dir> <dump_glob> <trajectory_h5>
Or use the provided run_ovito_vis.sh wrapper.

Dump files are sorted numerically by the timestep embedded in the filename
(e.g. full.4000.lammpstrj → 4000). This order matches the h5 frame index:
h5 frame 0 = lowest-timestep dump file, frame 1 = next, etc.

Voxel grid is hidden by default — toggle visibility in Ovito's pipeline panel.
"""

import sys
import os
import glob
import re

import numpy as np
import h5py

from ovito.io import import_file
from ovito.data import VoxelGrid, SimulationCell, DataCollection


# h5 datasets exposed as voxel grid properties
VOXEL_PROPERTIES = [
    'density', 'pressure', 'virial_pressure', 'temperature',
    'avg_speed', 'avg_O_speed', 'number_density', 'voxel_type',
]


def _parse_args():
    # When run via `ovito --script script.py -- arg1 arg2 arg3`,
    # sys.argv is [script_path, arg1, arg2, arg3]
    if len(sys.argv) < 4:
        print("Usage: ovito --script ovito_vis.py -- <dump_dir> <dump_glob> <trajectory_h5>",
              file=sys.stderr)
        sys.exit(1)
    return sys.argv[1], sys.argv[2], sys.argv[3]


def _sorted_dumps(dump_dir, dump_glob):
    """Return dump files sorted numerically by the timestep in the filename."""
    files = glob.glob(os.path.join(dump_dir, dump_glob))
    if not files:
        print(f"ERROR: No files matching '{dump_glob}' in '{dump_dir}'", file=sys.stderr)
        sys.exit(1)
    return sorted(files, key=lambda f: int(re.search(r'\.(\d+)\.[^.]+$', f).group(1)))


# -----------------------------------------------------------------------
dump_dir, dump_glob, h5_path = _parse_args()
dump_files = _sorted_dumps(dump_dir, dump_glob)
print(f"Found {len(dump_files)} dump files")
print(f"Trajectory h5: {h5_path}")


# -----------------------------------------------------------------------
# Atom pipeline — Ovito reads LAMMPS dump files natively in sorted order.
# Frame index in Ovito matches h5 frame index since both use the same sort.
pipeline = import_file(dump_files, input_format='lammps/dump')


# -----------------------------------------------------------------------
# Voxel modifier — appends a VoxelGrid to the pipeline's DataCollection
# at each frame by reading the corresponding h5 frame.

def _add_voxel_grid(frame, data: DataCollection):
    with h5py.File(h5_path, 'r') as f:
        n_frames = f['density'].shape[0]
        t = min(frame, n_frames - 1)

        voxel_type = f['voxel_type'][t][:]   # (nx, ny, nz)
        nx, ny, nz = voxel_type.shape

        vs   = float(f.attrs.get('voxel_size', 10.0))
        xlo  = float(f.attrs.get('xlo', 0.0))
        ylo  = float(f.attrs.get('ylo', 0.0))
        zlo  = float(f.attrs.get('zlo', 0.0))

        # Build voxel grid
        grid = VoxelGrid()
        grid.shape = (nx, ny, nz)

        # Spatial domain: align with the simulation box stored in the h5.
        # The dump file's SimulationCell is already in `data` from the atom
        # pipeline, but the voxel grid needs its own explicit domain so that
        # it occupies the correct region in world space.
        grid_cell = SimulationCell()
        grid_cell.pbc = (False, False, False)
        grid_cell.matrix = np.array([
            [nx * vs,      0,      0, xlo],
            [     0, ny * vs,      0, ylo],
            [     0,      0, nz * vs, zlo],
        ], dtype=float)
        grid.domain_ = grid_cell

        # Attach each h5 property as a flattened grid property.
        # NaN (empty voxels) is replaced with 0 since Ovito's color coding
        # does not handle NaN. The voxel_type==0 sentinel still identifies
        # empty voxels for display purposes.
        for name in VOXEL_PROPERTIES:
            if name not in f:
                continue
            arr = f[name][t][:].astype(np.float32).flatten()
            arr = np.where(np.isnan(arr), 0.0, arr)
            grid.create_property(name, data=arr)

        data.objects.append(grid)


# Append the voxel modifier to the existing atom pipeline so both share
# the same frame counter — no separate alignment needed.
pipeline.modifiers.append(_add_voxel_grid)
pipeline.add_to_scene()


# -----------------------------------------------------------------------
# Hide the voxel grid by default. Compute frame 0 to get the VoxelGrid
# vis object, then disable it. The setting persists across frames.
output = pipeline.compute(0)
for obj in output.objects:
    if isinstance(obj, VoxelGrid):
        obj.vis.enabled = False
        print("Voxel grid loaded and hidden — toggle in pipeline panel to show.")
        break
