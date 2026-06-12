import sys
import threading
import time

import h5py
import numpy as np
import pyvista as pv
from scipy.ndimage import uniform_filter
from vtkmodules.vtkCommonDataModel import vtkPiecewiseFunction
from vtkmodules.vtkRenderingCore import vtkColorTransferFunction, vtkVolume, vtkVolumeProperty
from vtkmodules.vtkRenderingVolume import vtkSmartVolumeMapper
from vtkmodules.util import numpy_support

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
H5_FILE = sys.argv[1] if len(sys.argv) > 1 else "trajectory.h5"

QUANTITIES = [
    'density', 'pressure', 'virial_pressure',
    'temperature', 'avg_speed', 'avg_O_speed', 'voxel_type',
]

CMAPS = {
    'density':         'viridis',
    'pressure':        'RdBu_r',
    'virial_pressure': 'RdBu_r',
    'temperature':     'plasma',
    'avg_speed':       'viridis',
    'avg_O_speed':     'viridis',
    'voxel_type':      'tab10',
}

SMOOTH_KERNEL  = 3
VIEW_TYPES     = ['all', 'water', 'silica']
PLAY_FPS_STEPS = [1, 2, 5, 10, 15, 24, 30, 60]


# ---------------------------------------------------------------------------
# Data layer
# ---------------------------------------------------------------------------
class TrajectoryData:
    def __init__(self, h5_path):
        print(f"Loading {h5_path} ...", flush=True)
        self.data = {}

        with h5py.File(h5_path, 'r') as f:
            attrs = dict(f.attrs)

            for qty in QUANTITIES:
                if qty not in f:
                    continue
                self.data[qty] = f[qty][:].astype(np.float32)

            self.voxel_type = f['voxel_type'][:] if 'voxel_type' in f else None

            # opacity from number_density, normalized to [0, 1] globally
            if 'number_density' in f:
                nd = f['number_density'][:].astype(np.float32)
                np.nan_to_num(nd, nan=0.0, copy=False)
                nd_max = nd.max()
                self.opacity = nd / nd_max if nd_max > 0 else nd
            else:
                # fallback: binary mask from voxel_type
                self.opacity = (self.voxel_type != 0).astype(np.float32) \
                    if self.voxel_type is not None \
                    else np.ones_like(self.data[QUANTITIES[0]], dtype=np.float32)

            if 'timesteps' in f:
                self.timestep_labels = f['timesteps'][:]
            elif 'timestep_labels' in f:
                self.timestep_labels = f['timestep_labels'][:]
            else:
                self.timestep_labels = None

        self.available = [q for q in QUANTITIES if q in self.data]

        self.type_ids = {
            'water':  int(attrs.get('h_type',  3)),
            'silica': int(attrs.get('si_type', 1)),
        }

        # data range over filled cells only
        filled_mask = self.opacity > 0
        self.data_range = {}
        for q in self.available:
            vals = self.data[q][filled_mask]
            if vals.size > 0:
                self.data_range[q] = (float(vals.min()), float(vals.max()))
            else:
                self.data_range[q] = (0.0, 1.0)

        # grid metadata
        first = self.data[self.available[0]]
        T, nx, ny, nz = first.shape
        voxel_size = float(attrs.get('voxel_size', 10.0))
        xlo = float(attrs.get('xlo', 0.0))
        ylo = float(attrs.get('ylo', 0.0))
        zlo = float(attrs.get('zlo', 0.0))
        self.meta = {
            'T':   T,
            'nx':  nx,  'ny':  ny,  'nz':  nz,
            'xlo': xlo, 'ylo': ylo, 'zlo': zlo,
            'xhi': float(attrs.get('xhi', xlo + nx * voxel_size)),
            'yhi': float(attrs.get('yhi', ylo + ny * voxel_size)),
            'zhi': float(attrs.get('zhi', zlo + nz * voxel_size)),
        }
        print(
            f"Loaded: {T} timesteps, "
            f"grid {nx}×{ny}×{nz}, "
            f"quantities: {self.available}",
            flush=True,
        )

    def get(self, quantity, t):
        return self.data[quantity][t]          # (nx, ny, nz)

    def get_voxel_type(self, t):
        return self.voxel_type[t] if self.voxel_type is not None else None

    def get_opacity(self, t):
        return self.opacity[t]                 # (nx, ny, nz) float [0, 1]

    def timestep_label(self, t):
        if self.timestep_labels is not None and t < len(self.timestep_labels):
            return f"  (step {int(self.timestep_labels[t])})"
        return ""

    @property
    def T(self):
        return self.meta['T']


# ---------------------------------------------------------------------------
# Parameter layer
# ---------------------------------------------------------------------------
class ViewParams:
    def __init__(self, first_quantity):
        self.quantity      = first_quantity
        self.view_type     = 'all'
        self.clip_axis     = None
        self.clip_position = 0.0
        self.clip_side     = 'below'
        self.smoothing     = False

    def copy(self):
        p = ViewParams(self.quantity)
        p.view_type     = self.view_type
        p.clip_axis     = self.clip_axis
        p.clip_position = self.clip_position
        p.clip_side     = self.clip_side
        p.smoothing     = self.smoothing
        return p


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _build_clip_mask(meta, clip_axis, clip_position, clip_side):
    nx, ny, nz = meta['nx'], meta['ny'], meta['nz']
    xlo, ylo, zlo = meta['xlo'], meta['ylo'], meta['zlo']
    xhi, yhi, zhi = meta['xhi'], meta['yhi'], meta['zhi']
    dx = (xhi - xlo) / nx
    dy = (yhi - ylo) / ny
    dz = (zhi - zlo) / nz

    if clip_axis == 'x':
        centers = xlo + (np.arange(nx) + 0.5) * dx
        mask_1d = (centers <= clip_position) if clip_side == 'below' else (centers >= clip_position)
        return np.broadcast_to(mask_1d[:, None, None], (nx, ny, nz)).copy()
    elif clip_axis == 'y':
        centers = ylo + (np.arange(ny) + 0.5) * dy
        mask_1d = (centers <= clip_position) if clip_side == 'below' else (centers >= clip_position)
        return np.broadcast_to(mask_1d[None, :, None], (nx, ny, nz)).copy()
    elif clip_axis == 'z':
        centers = zlo + (np.arange(nz) + 0.5) * dz
        mask_1d = (centers <= clip_position) if clip_side == 'below' else (centers >= clip_position)
        return np.broadcast_to(mask_1d[None, None, :], (nx, ny, nz)).copy()
    else:
        return np.ones((nx, ny, nz), dtype=bool)


def _build_ctf(quantity, trajectory):
    dmin, dmax   = trajectory.data_range[quantity]
    cmap_name    = CMAPS.get(quantity, 'viridis')
    lut          = pv.LookupTable(cmap=cmap_name, n_colors=256)
    lut.SetRange(dmin, dmax)
    ctf = vtkColorTransferFunction()
    for i in range(256):
        t   = dmin + i * (dmax - dmin) / 255.0 if dmax > dmin else dmin
        rgb = lut.map_value(t)[:3]
        ctf.AddRGBPoint(t, rgb[0] / 255.0, rgb[1] / 255.0, rgb[2] / 255.0)
    return ctf, dmin, dmax


# ---------------------------------------------------------------------------
# Cache layer
# ---------------------------------------------------------------------------
class FrameCache:
    def __init__(self, T, N):
        self.frames   = np.zeros((T, N, 2), dtype=np.float32)
        self.progress = set()
        self.ready    = False
        self._thread  = None
        self._cancel  = threading.Event()

    def build(self, trajectory, params, start_t=0):
        self._cancel.set()
        if self._thread is not None:
            self._thread.join()
        self._cancel.clear()
        self.progress = set()
        self.ready    = False

        self._thread = threading.Thread(
            target=self._worker,
            args=(trajectory, params.copy(), start_t),
            daemon=True,
        )
        self._thread.start()

    def _worker(self, trajectory, params, start_t):
        meta = trajectory.meta
        T    = trajectory.T

        clip_mask = _build_clip_mask(
            meta, params.clip_axis, params.clip_position, params.clip_side
        )

        for t in list(range(start_t, T)) + list(range(0, start_t)):
            if self._cancel.is_set():
                return

            data      = trajectory.get(params.quantity, t)    # (nx, ny, nz)
            opacity_t = trajectory.get_opacity(t)             # (nx, ny, nz) [0, 1]

            # 1. type filter
            vt = trajectory.get_voxel_type(t)
            if params.view_type == 'all' or vt is None:
                type_mask = np.ones(data.shape, dtype=bool)
            else:
                type_mask = (vt == trajectory.type_ids[params.view_type])

            # 2. clip
            visible = type_mask & clip_mask

            # 3. filled = has atoms AND passes visibility
            filled = (opacity_t > 0) & visible

            # component 0: color
            if params.smoothing:
                num      = uniform_filter(data,                          size=SMOOTH_KERNEL, mode='constant', cval=0.0)
                den      = uniform_filter((opacity_t > 0).astype(np.float32), size=SMOOTH_KERNEL, mode='constant', cval=0.0)
                smoothed = np.where(den > 0, num / den, 0.0)
                color    = np.where(filled, smoothed, 0.0)
            else:
                color    = np.where(filled, data, 0.0)

            # component 1: opacity — zero out cells that don't pass visibility
            opacity = np.where(filled, opacity_t, 0.0)

            self.frames[t, :, 0] = color.flatten(order='F')
            self.frames[t, :, 1] = opacity.flatten(order='F')
            self.progress.add(t)

        self.ready = True


# ---------------------------------------------------------------------------
# Render layer
# ---------------------------------------------------------------------------
def _make_grid(meta):
    nx, ny, nz = meta['nx'], meta['ny'], meta['nz']
    grid = pv.ImageData()
    grid.dimensions = (nx + 1, ny + 1, nz + 1)
    grid.spacing    = (
        (meta['xhi'] - meta['xlo']) / nx,
        (meta['yhi'] - meta['ylo']) / ny,
        (meta['zhi'] - meta['zlo']) / nz,
    )
    grid.origin = (meta['xlo'], meta['ylo'], meta['zlo'])
    return grid


def main():
    trajectory = TrajectoryData(H5_FILE)
    meta       = trajectory.meta
    T          = trajectory.T
    N          = meta['nx'] * meta['ny'] * meta['nz']

    params = ViewParams(trajectory.available[0])
    cache  = FrameCache(T, N)
    cache.build(trajectory, params, start_t=0)

    state = {
        't':        0,
        'playing':  False,
        'fps_idx':  PLAY_FPS_STEPS.index(10),
        'timer_id': None,
        'qty_idx':  0,
        'vtype_idx': 0,
    }

    # plotter + grid — built once, never rebuilt
    plotter = pv.Plotter(title="Shock Voxel Visualizer")
    plotter.set_background("#1a1a2e")

    grid = _make_grid(meta)

    # 2-component VTK cell array: col 0 = color, col 1 = opacity
    from vtkmodules.vtkCommonCore import vtkFloatArray
    vtk_arr = vtkFloatArray()
    vtk_arr.SetNumberOfComponents(2)
    vtk_arr.SetNumberOfTuples(N)
    vtk_arr.SetName('display')
    vtk_arr.Fill(0.0)
    grid.GetCellData().AddArray(vtk_arr)
    grid.GetCellData().SetActiveScalars('display')

    vtk_np_view = numpy_support.vtk_to_numpy(vtk_arr)   # (N, 2) view into VTK memory

    # volume setup
    mapper = vtkSmartVolumeMapper()
    mapper.SetInputData(grid)
    mapper.SetBlendModeToComposite()

    prop = vtkVolumeProperty()
    prop.IndependentComponentsOn()

    # component 1: opacity only — defined once, never changes
    ctf1 = vtkColorTransferFunction()
    ctf1.AddRGBPoint(0.0, 0.0, 0.0, 0.0)
    ctf1.AddRGBPoint(1.0, 0.0, 0.0, 0.0)
    otf1 = vtkPiecewiseFunction()
    otf1.AddPoint(0.0, 0.0)
    otf1.AddPoint(1.0, 1.0)
    prop.SetColor(1, ctf1)
    prop.SetScalarOpacity(1, otf1)
    prop.SetComponentWeight(1, 0.0)

    def _apply_quantity_to_prop(quantity):
        ctf0, dmin, dmax = _build_ctf(quantity, trajectory)
        otf0 = vtkPiecewiseFunction()
        otf0.AddPoint(dmin, 1.0)
        otf0.AddPoint(dmax, 1.0)
        prop.SetColor(0, ctf0)
        prop.SetScalarOpacity(0, otf0)
        prop.SetComponentWeight(0, 1.0)

    _apply_quantity_to_prop(params.quantity)

    vol = vtkVolume()
    vol.SetMapper(mapper)
    vol.SetProperty(prop)
    plotter.renderer.AddVolume(vol)

    plotter.add_axes(xlabel='X', ylabel='Y', zlabel='Z', color='white', line_width=3)
    plotter.show_bounds(
        bounds=(meta['xlo'], meta['xhi'], meta['ylo'], meta['yhi'], meta['zlo'], meta['zhi']),
        grid='back', location='outer', ticks='outside',
        font_size=10, color='white', fmt='%.0f',
    )

    # HUD
    def _hud_text():
        fps  = PLAY_FPS_STEPS[state['fps_idx']]
        t    = state['t']
        clip = f"{params.clip_axis} @ {params.clip_position:.1f} ({params.clip_side})" \
               if params.clip_axis else 'OFF'
        prog = f"{len(cache.progress)}/{T}" if not cache.ready else "ready"
        return (
            f"Quantity  : {params.quantity}\n"
            f"Timestep  : {t}/{T - 1}{trajectory.timestep_label(t)}\n"
            f"View type : {params.view_type}\n"
            f"Clip      : {clip}\n"
            f"Smoothing : {'ON' if params.smoothing else 'OFF'}\n"
            f"Play      : {'▶' if state['playing'] else '⏸'}  {fps} fps\n"
            f"Cache     : {prog}\n\n"
            f"[space] play/pause    [←/→] step\n"
            f"[[ / ]] speed -/+    [v] quantity\n"
            f"[w] view type        [x/y/z] clip axis (again=off)\n"
            f"[+/-] clip pos       [f] smoothing"
        )

    hud = plotter.add_text(
        _hud_text(), position='upper_left', font_size=9, color='#e0e0e0'
    )

    def _fast_update(t):
        if t not in cache.progress:
            return
        vtk_np_view[:] = cache.frames[t]      # (N, 2) in-place write
        grid.GetCellData().Modified()
        hud.SetText(2, _hud_text())
        plotter.render_window.Render()

    def _on_param_change():
        _apply_quantity_to_prop(params.quantity)
        cache.build(trajectory, params, start_t=state['t'])
        _fast_update(state['t'])

    # timer
    def _advance_frame(step=None):
        if not state['playing']:
            return
        state['t'] = (state['t'] + 1) % T
        _fast_update(state['t'])

    def _restart_timer():
        if state['timer_id'] is not None:
            plotter.remove_timer_event(state['timer_id'])
        fps = PLAY_FPS_STEPS[state['fps_idx']]
        state['timer_id'] = plotter.add_timer_event(
            max_steps=10_000_000,
            duration=int(1000 / fps),
            callback=_advance_frame,
        )

    # key bindings
    def toggle_play():
        state['playing'] = not state['playing']
        if state['playing']:
            _restart_timer()
        else:
            if state['timer_id'] is not None:
                plotter.remove_timer_event(state['timer_id'])
                state['timer_id'] = None
        hud.SetText(2, _hud_text())
        plotter.render_window.Render()

    def step_forward():
        state['playing'] = False
        state['t'] = min(state['t'] + 1, T - 1)
        _fast_update(state['t'])

    def step_back():
        state['playing'] = False
        state['t'] = max(state['t'] - 1, 0)
        _fast_update(state['t'])

    def speed_up():
        state['fps_idx'] = min(state['fps_idx'] + 1, len(PLAY_FPS_STEPS) - 1)
        if state['playing']:
            _restart_timer()
        hud.SetText(2, _hud_text())
        plotter.render_window.Render()

    def slow_down():
        state['fps_idx'] = max(state['fps_idx'] - 1, 0)
        if state['playing']:
            _restart_timer()
        hud.SetText(2, _hud_text())
        plotter.render_window.Render()

    def cycle_quantity():
        state['qty_idx'] = (state['qty_idx'] + 1) % len(trajectory.available)
        params.quantity  = trajectory.available[state['qty_idx']]
        _on_param_change()

    def cycle_view_type():
        state['vtype_idx'] = (state['vtype_idx'] + 1) % len(VIEW_TYPES)
        params.view_type   = VIEW_TYPES[state['vtype_idx']]
        _on_param_change()

    def _set_clip_axis(axis):
        if params.clip_axis == axis:
            params.clip_axis = None
        else:
            params.clip_axis     = axis
            params.clip_position = (meta[axis + 'lo'] + meta[axis + 'hi']) / 2.0
        _on_param_change()

    def clip_pos_increase():
        if params.clip_axis is None:
            return
        axis = params.clip_axis
        step = (meta[axis + 'hi'] - meta[axis + 'lo']) / meta['n' + axis]
        params.clip_position = min(params.clip_position + step, meta[axis + 'hi'])
        _on_param_change()

    def clip_pos_decrease():
        if params.clip_axis is None:
            return
        axis = params.clip_axis
        step = (meta[axis + 'hi'] - meta[axis + 'lo']) / meta['n' + axis]
        params.clip_position = max(params.clip_position - step, meta[axis + 'lo'])
        _on_param_change()

    def toggle_smoothing():
        params.smoothing = not params.smoothing
        _on_param_change()

    plotter.add_key_event('space',        toggle_play)
    plotter.add_key_event('Right',        step_forward)
    plotter.add_key_event('Left',         step_back)
    plotter.add_key_event('bracketright', speed_up)
    plotter.add_key_event('bracketleft',  slow_down)
    plotter.add_key_event('v',            cycle_quantity)
    plotter.add_key_event('w',            cycle_view_type)
    plotter.add_key_event('x',            lambda: _set_clip_axis('x'))
    plotter.add_key_event('y',            lambda: _set_clip_axis('y'))
    plotter.add_key_event('z',            lambda: _set_clip_axis('z'))
    plotter.add_key_event('equal',        clip_pos_increase)
    plotter.add_key_event('minus',        clip_pos_decrease)
    plotter.add_key_event('f',            toggle_smoothing)

    # wait until t=0 is ready then show
    while 0 not in cache.progress:
        time.sleep(0.01)
    _fast_update(0)

    plotter.show()


if __name__ == "__main__":
    main()
