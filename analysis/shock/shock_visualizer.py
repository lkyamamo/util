import sys
import tkinter as tk
from tkinter import simpledialog

import numpy as np
import h5py
import pyvista as pv
from scipy.ndimage import uniform_filter

# --- Configuration ---
H5_FILE    = sys.argv[1] if len(sys.argv) > 1 else "trajectory.h5"
QUANTITIES = ['density', 'pressure', 'virial_pressure', 'temperature',
              'avg_speed', 'avg_O_speed', 'voxel_type']
CMAPS = {
    'density':          'viridis',
    'pressure':         'RdBu_r',
    'virial_pressure':  'RdBu_r',
    'temperature':      'plasma',
    'avg_speed':        'viridis',
    'avg_O_speed':      'viridis',
    'voxel_type':       'Set1',
}


# ---------------------------------------------------------------------------
# State
# ---------------------------------------------------------------------------
class ViewerState:
    def __init__(self):
        self.timestep_idx = 0
        self.quantity_idx = 0
        self.smooth_kernel = 1
        self.clip_active   = False
        self.clip_normal   = np.array([1.0, 0.0, 0.0])
        self.clip_origin   = np.array([0.0, 0.0, 0.0])

    @property
    def quantity(self):
        return QUANTITIES[self.quantity_idx]


# ---------------------------------------------------------------------------
# Data helpers
# ---------------------------------------------------------------------------
def load_grid(h5file, state, meta):
    nx, ny, nz     = meta['nx'], meta['ny'], meta['nz']
    voxel_size     = meta['voxel_size']
    xlo, ylo, zlo  = meta['xlo'], meta['ylo'], meta['zlo']

    raw = h5file[state.quantity][state.timestep_idx]   # (nx, ny, nz)
    if state.smooth_kernel > 1:
        raw = uniform_filter(raw.astype(np.float64), size=state.smooth_kernel)

    grid = pv.ImageData()
    grid.dimensions = (nx + 1, ny + 1, nz + 1)        # cell-centred needs dims = shape + 1
    grid.spacing    = (voxel_size, voxel_size, voxel_size)
    grid.origin     = (xlo, ylo, zlo)
    grid.cell_data[state.quantity] = raw.flatten(order='F')
    return grid


def read_meta(h5file):
    """Read scalar metadata from the first timestep's HDF5 attributes."""
    attrs = dict(h5file.attrs) if h5file.attrs else {}
    # fall back to inferring from dataset shape
    ds = h5file[QUANTITIES[0]]
    T, nx, ny, nz = ds.shape
    return {
        'T':          T,
        'nx':         nx,
        'ny':         ny,
        'nz':         nz,
        'voxel_size': float(attrs.get('voxel_size', 10.0)),
        'xlo':        float(attrs.get('xlo', 0.0)),
        'ylo':        float(attrs.get('ylo', 0.0)),
        'zlo':        float(attrs.get('zlo', 0.0)),
    }


# ---------------------------------------------------------------------------
# Clip plane dialog  (ax + by + cz = d)
# ---------------------------------------------------------------------------
def ask_clip_plane(state):
    root = tk.Tk()
    root.withdraw()
    raw = simpledialog.askstring(
        "Clip plane",
        "Enter plane equation coefficients:\n  a  b  c  d   (for ax+by+cz=d)",
        parent=root,
    )
    root.destroy()
    if raw is None:
        return
    try:
        a, b, c, d = map(float, raw.split())
    except ValueError:
        print("Invalid input — expected four floats: a b c d")
        return
    normal = np.array([a, b, c])
    norm_sq = np.dot(normal, normal)
    if norm_sq == 0:
        print("Normal vector (a,b,c) cannot be zero.")
        return
    state.clip_normal  = normal
    state.clip_origin  = normal * (d / norm_sq)
    state.clip_active  = True
    print(f"Clip plane set: {a}x + {b}y + {c}z = {d}")


# ---------------------------------------------------------------------------
# Main viewer
# ---------------------------------------------------------------------------
def main():
    h5file = h5py.File(H5_FILE, 'r')
    meta   = read_meta(h5file)
    state  = ViewerState()
    T      = meta['T']

    plotter = pv.Plotter(title="Shock Voxel Visualizer")
    plotter.set_background("black")

    # ---- render helper ----
    def refresh():
        plotter.clear()

        grid = load_grid(h5file, state, meta)
        cmap = CMAPS.get(state.quantity, 'viridis')

        if state.clip_active:
            clipped = grid.clip(normal=state.clip_normal, origin=state.clip_origin)
            plotter.add_mesh(clipped, scalars=state.quantity, cmap=cmap,
                             show_scalar_bar=True)
        else:
            plotter.add_volume(grid, scalars=state.quantity, cmap=cmap,
                               opacity='linear', shade=False)

        ts_val = ""
        if 'timestep' in (h5file.attrs or {}):
            ts_val = f"  (step {int(h5file.attrs['timestep'])})"
        plotter.add_text(
            f"Quantity : {state.quantity}\n"
            f"Timestep : {state.timestep_idx}/{T - 1}{ts_val}\n"
            f"Smoothing: {state.smooth_kernel}\n"
            f"Clip     : {'ON' if state.clip_active else 'OFF'}\n\n"
            f"[n/N] next/prev timestep\n"
            f"[q]   cycle quantity\n"
            f"[s/S] smoothing +/-\n"
            f"[p]   set clip plane  (ax+by+cz=d)\n"
            f"[c]   toggle clip\n"
            f"[r]   reset clip",
            position="upper_left", font_size=9, color="white",
        )
        plotter.render()

    # ---- key bindings ----
    def next_timestep():
        state.timestep_idx = min(state.timestep_idx + 1, T - 1)
        refresh()

    def prev_timestep():
        state.timestep_idx = max(state.timestep_idx - 1, 0)
        refresh()

    def cycle_quantity():
        state.quantity_idx = (state.quantity_idx + 1) % len(QUANTITIES)
        refresh()

    def smooth_up():
        state.smooth_kernel = min(state.smooth_kernel + 2, 15)
        refresh()

    def smooth_down():
        state.smooth_kernel = max(state.smooth_kernel - 2, 1)
        refresh()

    def set_clip_plane():
        ask_clip_plane(state)
        refresh()

    def toggle_clip():
        state.clip_active = not state.clip_active
        refresh()

    def reset_clip():
        state.clip_active = False
        refresh()

    plotter.add_key_event('n', next_timestep)
    plotter.add_key_event('N', prev_timestep)
    plotter.add_key_event('q', cycle_quantity)
    plotter.add_key_event('s', smooth_up)
    plotter.add_key_event('S', smooth_down)
    plotter.add_key_event('p', set_clip_plane)
    plotter.add_key_event('c', toggle_clip)
    plotter.add_key_event('r', reset_clip)

    refresh()
    plotter.show()
    h5file.close()


if __name__ == "__main__":
    main()
