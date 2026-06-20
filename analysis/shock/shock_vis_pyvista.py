import numpy as np
import pyvista as pv
import h5py

VOXEL_SIZE = 1.0
VOXEL_GAP  = 0.05   # fractional gap between voxels
#FILEPATH   = "/tmp/shock_vis_verify.h5"
FILEPATH  = '/Users/loganyamamoto/Desktop/Research/grants/geo_sciences/finalized-bubble-collapse/analysis/small_system/voxel_r0137/trajectory.h5'

PROPERTIES = ['density', 'pressure', 'virial_pressure', 'temperature',
              'avg_speed', 'avg_O_speed', 'voxel_type', 'v_COM']

ZERO_BASED_RANGE_PROPERTIES = {'density', 'temperature', 'avg_speed', 'avg_O_speed'}

# Percentile clip for auto display ranges: e.g. 1.0 uses 1st-99th percentile.
# Set to 0.0 to use true min/max (may bunch up if outliers exist).
PERCENTILE_CLIP = 1.0

# Optional per-property display range in real units: (min, max) or None to use data range
PROPERTY_DISPLAY_RANGES = {
    'density':         None,
    'pressure':        None,
    'virial_pressure': None,
    'temperature':     (0,1000),
    'avg_speed':       None,
    'avg_O_speed':     None,
    'voxel_type':      None,
    'v_COM':           None,
}

# Playback speed options in ms per frame (lower = faster)
PLAY_SPEEDS_MS = [500, 250, 100, 50, 20]

# cube corner offsets (unit cube, scaled by half_size later)
_CORNERS = np.array([
    [-1,-1,-1],[+1,-1,-1],[+1,+1,-1],[-1,+1,-1],
    [-1,-1,+1],[+1,-1,+1],[+1,+1,+1],[-1,+1,+1],
], dtype=np.float32) * 0.5

# 6 faces (quads), each referencing 4 of the 8 corners
_FACES = np.array([
    [0,3,2,1],[4,5,6,7],  # -z, +z
    [0,1,5,4],[3,7,6,2],  # -y, +y
    [0,4,7,3],[1,2,6,5],  # -x, +x
], dtype=np.int32)

KEY_HELP = (
    "--- Navigation ---\n"
    "Right / Left   : next / prev frame\n"
    "Up / Down      : next / prev property\n"
    "Space          : play / pause\n"
    "f / g          : faster / slower playback\n"
    "0-9 + Return   : jump to timestep\n"
    "Escape         : cancel jump input\n"
    "--- Voxel Filter ---\n"
    "t              : toggle voxel type filter\n"
    "] / [          : increase / decrease voxel type\n"
    "--- Slicing ---\n"
    "s              : toggle slice\n"
    "x / y / z      : slice axis (yz / xz / xy)\n"
    ". / ,          : move slice forward / backward\n"
    "= / -          : thicker / thinner slice\n"
    "--- Camera ---\n"
    "h              : toggle this help text"
)


class VoxelGrid:
    def __init__(self, data):
        self.data = data
        self.ntimesteps, self.nx, self.ny, self.nz = data['density'].shape
        self.current_property  = 'density'
        self.current_timestep  = 0
        self.voxel_type_filter = None

        # slicing
        self.slice_enabled   = False
        self.slice_axis      = 'xy'
        self.slice_center    = 0
        self.slice_thickness = 1

        # per-property data ranges
        self.ranges = {}
        for prop in PROPERTIES:
            arr = data[prop]
            if prop == 'v_COM':
                arr = np.linalg.norm(arr, axis=-1)
            valid = arr[~np.isnan(arr)]
            if PERCENTILE_CLIP > 0.0:
                lo = 0.0 if prop in ZERO_BASED_RANGE_PROPERTIES else float(np.percentile(valid, PERCENTILE_CLIP))
                hi = float(np.percentile(valid, 100.0 - PERCENTILE_CLIP))
            else:
                lo = 0.0 if prop in ZERO_BASED_RANGE_PROPERTIES else float(valid.min())
                hi = float(valid.max())
            self.ranges[prop] = (lo, hi)

        self.display_ranges = {}
        for prop in PROPERTIES:
            override = PROPERTY_DISPLAY_RANGES.get(prop)
            self.display_ranges[prop] = override if override is not None else self.ranges[prop]

        self.nan_mask = np.isnan(data['density'])  # (ntimesteps, nx, ny, nz)

    # ------------------------------------------------------------------
    def _compute_mask(self):
        t    = self.current_timestep
        mask = self.nan_mask[t].copy()
        if self.voxel_type_filter is not None:
            mask |= (self.data['voxel_type'][t] != self.voxel_type_filter)
        if self.slice_enabled:
            lo = self.slice_center - self.slice_thickness
            hi = self.slice_center + self.slice_thickness
            if self.slice_axis == 'xy':
                iz = np.arange(self.nz)
                mask |= ((iz < lo) | (iz > hi))[np.newaxis, np.newaxis, :]
            elif self.slice_axis == 'xz':
                iy = np.arange(self.ny)
                mask |= ((iy < lo) | (iy > hi))[np.newaxis, :, np.newaxis]
            elif self.slice_axis == 'yz':
                ix = np.arange(self.nx)
                mask |= ((ix < lo) | (ix > hi))[:, np.newaxis, np.newaxis]
        return mask

    @property
    def slice_axis_size(self):
        return {'xy': self.nz, 'xz': self.ny, 'yz': self.nx}[self.slice_axis]

    def _to_scalar(self, prop, arr):
        """Return raw scalar values, computing magnitude for v_COM."""
        if prop == 'v_COM':
            return np.linalg.norm(arr, axis=-1)
        return arr

    def _normalize_density(self, arr):
        lo, hi = self.display_ranges['density']
        denom  = (hi - lo) if hi != lo else 1.0
        return np.clip((arr - lo) / denom, 0.0, 1.0)

    # ------------------------------------------------------------------
    def build_mesh(self):
        """Build a PyVista PolyData for currently visible voxels."""
        t    = self.current_timestep
        prop = self.current_property

        color_arr   = self._to_scalar(prop, self.data[prop][t].astype(np.float32))
        density_arr = self._normalize_density(self.data['density'][t].astype(np.float32))

        mask = self._compute_mask() | np.isnan(color_arr) | np.isnan(density_arr)
        visible = ~mask

        vox_idx = np.argwhere(visible)  # (N, 3)
        N = len(vox_idx)

        if N == 0:
            lo, hi = self.display_ranges[prop]
            return pv.PolyData(), lo, hi

        half    = VOXEL_SIZE * (1.0 - VOXEL_GAP) / 2.0
        centers = (vox_idx + 0.5) * VOXEL_SIZE  # (N, 3)

        # points: (N*8, 3)
        points = (centers[:, np.newaxis, :] + _CORNERS[np.newaxis, :, :] * (half * 2)).reshape(-1, 3)

        # faces in PyVista flat format: [4, v0,v1,v2,v3, 4, ...]
        base  = np.arange(N, dtype=np.int32) * 8
        quads = (base[:, np.newaxis, np.newaxis] + _FACES[np.newaxis, :, :]).reshape(-1, 4)
        faces_flat = np.empty(N * 6 * 5, dtype=np.int32)
        faces_flat[0::5] = 4
        faces_flat[1::5] = quads[:, 0]
        faces_flat[2::5] = quads[:, 1]
        faces_flat[3::5] = quads[:, 2]
        faces_flat[4::5] = quads[:, 3]

        mesh = pv.PolyData(points.astype(np.float32), faces_flat)

        ci, cj, ck      = vox_idx[:, 0], vox_idx[:, 1], vox_idx[:, 2]
        color_per_vox   = color_arr[ci, cj, ck]    # (N,) real values
        density_per_vox = density_arr[ci, cj, ck]  # (N,) normalized [0,1]

        mesh.cell_data['scalars'] = np.repeat(color_per_vox,   6)
        # alpha: 0.1 at min density, 1.0 at max
        mesh.cell_data['alpha']   = np.clip(0.1 + 0.9 * np.repeat(density_per_vox, 6), 0.0, 1.0)

        lo, hi = self.display_ranges[prop]
        return mesh, lo, hi


# -----------------------------------------------------------------------
# load data
with h5py.File(FILEPATH, 'r') as f:
    full_data = {prop: np.array(f[prop]) for prop in PROPERTIES}

grid = VoxelGrid(full_data)

# -----------------------------------------------------------------------
# plotter setup
pl = pv.Plotter(window_size=(1000, 800))
pl.set_background('white')
pl.enable_parallel_projection()

mesh, lo, hi = grid.build_mesh()
actor = pl.add_mesh(
    mesh,
    scalars='scalars',
    clim=[lo, hi],
    cmap='coolwarm',
    show_scalar_bar=False,
    opacity='alpha',
)

sbar = pl.add_scalar_bar(
    title=grid.current_property,
    n_labels=5,
    vertical=True,
    position_x=0.88,
    position_y=0.1,
    width=0.05,
    height=0.8,
    label_font_size=14,
    title_font_size=16,
    color='black',
)

pl.show_bounds(color='black', xlabel='X', ylabel='Y', zlabel='Z')

# face the XZ plane: camera sits at -Y, looks in +Y, Z is up → +X goes right
_cx = grid.nx * VOXEL_SIZE / 2.0
_cy = grid.ny * VOXEL_SIZE / 2.0
_cz = grid.nz * VOXEL_SIZE / 2.0
pl.camera.focal_point = (_cx, _cy, _cz)
pl.camera.position    = (_cx, _cy - max(grid.nx, grid.nz) * VOXEL_SIZE * 3.0, _cz)
pl.camera.up          = (0.0, 0.0, 1.0)

# -----------------------------------------------------------------------
# state
current_timestep  = [0]
current_prop_idx  = [0]
voxel_type_filter = [None]
active_voxel_type = [1]
playing           = [False]
play_speed_idx    = [0]   # index into PLAY_SPEEDS_MS
jump_buffer       = ['']  # digits typed for jump-to-frame
timer_id          = [None]
help_visible      = [True]


# -----------------------------------------------------------------------
# text actors
def _info_text():
    speed_ms = PLAY_SPEEDS_MS[play_speed_idx[0]]
    play_str = f"[PLAYING {speed_ms}ms/frame]" if playing[0] else "[paused]"
    jump_str = f"  Jump: {jump_buffer[0]}_" if jump_buffer[0] else ''
    return (f"Property : {grid.current_property}\n"
            f"Timestep : {grid.current_timestep} / {grid.ntimesteps - 1}  {play_str}{jump_str}")

info_actor = pl.add_text(_info_text(), position='upper_left', font_size=12, color='black')
help_actor = pl.add_text(KEY_HELP, position='lower_left', font_size=10, color='black')


# -----------------------------------------------------------------------
def refresh():
    global actor
    mesh, lo, hi = grid.build_mesh()
    pl.remove_actor(actor)
    actor = pl.add_mesh(
        mesh,
        scalars='scalars',
        clim=[lo, hi],
        cmap='coolwarm',
        show_scalar_bar=False,
        opacity='alpha',
    )
    sbar.SetTitle(grid.current_property)
    sbar.GetLookupTable().SetTableRange(lo, hi)
    info_actor.SetText(2, _info_text())
    pl.render()
    pl.update()


def _set_timestep(t):
    current_timestep[0] = max(0, min(t, grid.ntimesteps - 1))
    grid.current_timestep = current_timestep[0]


# -----------------------------------------------------------------------
# playback timer
iren = pl.iren

def _on_timer(obj, event):
    if playing[0]:
        if current_timestep[0] >= grid.ntimesteps - 1:
            _stop_playback()
            return
        _set_timestep(current_timestep[0] + 1)
        refresh()

iren._iren.AddObserver('TimerEvent', _on_timer)


def _start_playback():
    playing[0] = True
    timer_id[0] = iren._iren.CreateRepeatingTimer(PLAY_SPEEDS_MS[play_speed_idx[0]])

def _stop_playback():
    playing[0] = False
    if timer_id[0] is not None:
        iren._iren.DestroyTimer(timer_id[0])
        timer_id[0] = None
    info_actor.SetText(2, _info_text())
    pl.render()
    pl.update()

def _restart_timer_if_playing():
    if playing[0]:
        _stop_playback()
        _start_playback()


# -----------------------------------------------------------------------
# key callbacks

# timestep
def next_frame():
    _stop_playback()
    _set_timestep(current_timestep[0] + 1)
    refresh()

def prev_frame():
    _stop_playback()
    _set_timestep(current_timestep[0] - 1)
    refresh()

def toggle_play():
    if playing[0]:
        _stop_playback()
    else:
        _start_playback()
    info_actor.SetText(2, _info_text())
    pl.render()
    pl.update()

def faster():
    play_speed_idx[0] = min(play_speed_idx[0] + 1, len(PLAY_SPEEDS_MS) - 1)
    print(f"Playback speed: {PLAY_SPEEDS_MS[play_speed_idx[0]]} ms/frame")
    _restart_timer_if_playing()
    info_actor.SetText(2, _info_text())
    pl.render()
    pl.update()

def slower():
    play_speed_idx[0] = max(play_speed_idx[0] - 1, 0)
    print(f"Playback speed: {PLAY_SPEEDS_MS[play_speed_idx[0]]} ms/frame")
    _restart_timer_if_playing()
    info_actor.SetText(2, _info_text())
    pl.render()
    pl.update()

# jump-to-frame digit accumulation
def _digit(d):
    def _cb():
        jump_buffer[0] += str(d)
        info_actor.SetText(2, _info_text())
        pl.render()
        pl.update()
    return _cb

def confirm_jump():
    if jump_buffer[0]:
        t = int(jump_buffer[0])
        jump_buffer[0] = ''
        _stop_playback()
        _set_timestep(t)
        refresh()
    else:
        info_actor.SetText(2, _info_text())
        pl.render()
        pl.update()

def cancel_jump():
    jump_buffer[0] = ''
    info_actor.SetText(2, _info_text())
    pl.render()
    pl.update()

# property
def next_prop():
    current_prop_idx[0] = (current_prop_idx[0] + 1) % len(PROPERTIES)
    grid.current_property = PROPERTIES[current_prop_idx[0]]
    refresh()

def prev_prop():
    current_prop_idx[0] = (current_prop_idx[0] - 1) % len(PROPERTIES)
    grid.current_property = PROPERTIES[current_prop_idx[0]]
    refresh()

# voxel type filter
def toggle_voxel_filter():
    if voxel_type_filter[0] is None:
        voxel_type_filter[0] = active_voxel_type[0]
    else:
        voxel_type_filter[0] = None
    grid.voxel_type_filter = voxel_type_filter[0]
    refresh()

def inc_voxel_type():
    active_voxel_type[0] += 1
    if voxel_type_filter[0] is not None:
        voxel_type_filter[0] = active_voxel_type[0]
        grid.voxel_type_filter = voxel_type_filter[0]
        refresh()

def dec_voxel_type():
    active_voxel_type[0] -= 1
    if voxel_type_filter[0] is not None:
        voxel_type_filter[0] = active_voxel_type[0]
        grid.voxel_type_filter = voxel_type_filter[0]
        refresh()

# slicing
def toggle_slice():
    grid.slice_enabled = not grid.slice_enabled
    if grid.slice_enabled:
        grid.slice_center = grid.slice_axis_size // 2
    refresh()

def slice_axis_yz():
    grid.slice_axis = 'yz'; grid.slice_center = grid.nx // 2; refresh()

def slice_axis_xz():
    grid.slice_axis = 'xz'; grid.slice_center = grid.ny // 2; refresh()

def slice_axis_xy():
    grid.slice_axis = 'xy'; grid.slice_center = grid.nz // 2; refresh()

def slice_forward():
    grid.slice_center = min(grid.slice_center + 1, grid.slice_axis_size - 1); refresh()

def slice_backward():
    grid.slice_center = max(grid.slice_center - 1, 0); refresh()

def slice_thicker():
    grid.slice_thickness += 1; refresh()

def slice_thinner():
    grid.slice_thickness = max(0, grid.slice_thickness - 1); refresh()

# help toggle
def toggle_help():
    help_visible[0] = not help_visible[0]
    help_actor.SetVisibility(help_visible[0])
    pl.render()
    pl.update()


# -----------------------------------------------------------------------
# key bindings
pl.add_key_event('Right',        next_frame)
pl.add_key_event('Left',         prev_frame)
pl.add_key_event('Up',           next_prop)
pl.add_key_event('Down',         prev_prop)
pl.add_key_event('space',        toggle_play)
pl.add_key_event('f',            faster)
pl.add_key_event('g',            slower)
pl.add_key_event('Return',       confirm_jump)
pl.add_key_event('Escape',       cancel_jump)
for _d in range(10):
    pl.add_key_event(str(_d), _digit(_d))
pl.add_key_event('t',            toggle_voxel_filter)
pl.add_key_event('bracketright', inc_voxel_type)
pl.add_key_event('bracketleft',  dec_voxel_type)
pl.add_key_event('s',            toggle_slice)
pl.add_key_event('w',            lambda: None)
pl.add_key_event('x',            slice_axis_yz)
pl.add_key_event('y',            slice_axis_xz)
pl.add_key_event('z',            slice_axis_xy)
pl.add_key_event('period',       slice_forward)
pl.add_key_event('comma',        slice_backward)
pl.add_key_event('equal',        slice_thicker)
pl.add_key_event('minus',        slice_thinner)
pl.add_key_event('h',            toggle_help)

pl.show()
