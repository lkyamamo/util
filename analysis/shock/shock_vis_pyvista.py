import numpy as np
import pyvista as pv
import h5py

VOXEL_SIZE = 1.0
VOXEL_GAP  = 0.05   # fractional gap between voxels
#FILEPATH   = "/tmp/shock_vis_verify.h5"
FILEPATH  = '/scratch1/lkyamamo/finalized-bubble-collapse/analysis/small_interface/voxel_0147/trajectory.h5'

PROPERTIES = ['density', 'pressure', 'virial_pressure', 'temperature',
              'avg_speed', 'avg_O_speed', 'voxel_type', 'v_COM']

# Additional properties shown only in the voxel info panel (not cycled
# through with Up/Down or used for the colored mesh display)
INFO_EXTRA_PROPERTIES = ['number_density']

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

# Velocity arrow length range in display units.
# Largest-magnitude voxel → ARROW_LENGTH_MAX; smallest → ARROW_LENGTH_MIN.
ARROW_LENGTH_MIN = 0.5
ARROW_LENGTH_MAX = 3.0

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
    "k              : toggle slice\n"
    "x / y / z      : slice axis (yz / xz / xy)\n"
    ". / ,          : move slice forward / backward\n"
    "= / -          : thicker / thinner slice\n"
    "--- Overlays ---\n"
    "a              : toggle density opacity\n"
    "m              : toggle voxel mesh\n"
    "p              : toggle Si surface plane\n"
    "s              : toggle crater sphere cap\n"
    "c              : toggle crater/surface coloring (Si mode)\n"
    "v              : toggle velocity vectors\n"
    "< / >          : scale velocity arrows down / up\n"
    "--- Voxel Inspection ---\n"
    "left-click     : select voxel, show its properties\n"
    "i              : clear voxel selection\n"
    "--- Camera ---\n"
    "F1/F2          : face +X / -X\n"
    "F3/F4          : face +Y / -Y\n"
    "F5/F6          : face +Z / -Z\n"
    "h              : toggle this help text"
)


class VoxelGrid:
    def __init__(self, data, opt):
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

        # optional overlay data
        self.si_surface          = opt.get('si_surface')           # (T, ny, nz) Å
        self.si_surface_mask     = opt.get('si_surface_mask')      # (T, nx, ny, nz) uint8
        self.si_crater_mask      = opt.get('crater/si_crater_mask')# (T, nx, ny, nz) uint8
        self.crater_reference_x  = opt.get('crater/reference_x')  # (T,) Å
        self.crater_sphere_center= opt.get('crater/sphere_center') # (T, 3) Å
        self.crater_sphere_radius= opt.get('crater/sphere_radius') # (T,) Å

        # coordinate conversion: Å → display (display = voxel index * VOXEL_SIZE)
        self.xlo            = opt.get('xlo', 0.0)
        self.ylo            = opt.get('ylo', 0.0)
        self.zlo            = opt.get('zlo', 0.0)
        self.voxel_size_real= opt.get('voxel_size_real', 10.0)

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

        self.nan_mask = np.isnan(data['density']) | (data['voxel_type'] == 0)  # (T, nx, ny, nz)

    def angstrom_to_display(self, x, y, z):
        """Convert real Å coordinates to display voxel coordinates."""
        dx = (x - self.xlo) / self.voxel_size_real * VOXEL_SIZE
        dy = (y - self.ylo) / self.voxel_size_real * VOXEL_SIZE
        dz = (z - self.zlo) / self.voxel_size_real * VOXEL_SIZE
        return dx, dy, dz

    def x_angstrom_to_display(self, x):
        return (x - self.xlo) / self.voxel_size_real * VOXEL_SIZE

    def radius_to_display(self, r):
        return r / self.voxel_size_real * VOXEL_SIZE

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
        if prop == 'v_COM':
            return np.linalg.norm(arr, axis=-1)
        return arr

    def _normalize_density(self, arr):
        lo, hi = self.display_ranges['density']
        denom  = (hi - lo) if hi != lo else 1.0
        return np.clip((arr - lo) / denom, 0.0, 1.0)

    # ------------------------------------------------------------------
    def build_mesh(self, crater_colors=False):
        """Build a PyVista PolyData for currently visible voxels.

        Returns (mesh, lo, hi, is_categorical).
        is_categorical=True when crater coloring overrides the property scalar.
        """
        t    = self.current_timestep
        prop = self.current_property

        use_categorical = (
            crater_colors
            and self.voxel_type_filter == 2
            and self.si_surface_mask is not None
        )

        if use_categorical:
            surf_mask = self.si_surface_mask[t]  # (nx, ny, nz)
            crat_mask = self.si_crater_mask[t] if self.si_crater_mask is not None else None
            color_arr = np.zeros((self.nx, self.ny, self.nz), dtype=np.float32)
            color_arr[surf_mask == 1] = 1.0
            if crat_mask is not None:
                color_arr[(surf_mask == 1) & (crat_mask == 1)] = 2.0
            nan_mask_color = np.zeros((self.nx, self.ny, self.nz), dtype=bool)
            lo, hi = 0.0, 2.0
        else:
            color_arr      = self._to_scalar(prop, self.data[prop][t].astype(np.float32))
            nan_mask_color = np.isnan(color_arr)
            lo, hi         = self.display_ranges[prop]

        density_arr = self._normalize_density(self.data['density'][t].astype(np.float32))
        mask    = self._compute_mask() | nan_mask_color | np.isnan(density_arr)
        visible = ~mask

        vox_idx = np.argwhere(visible)  # (N, 3)
        N = len(vox_idx)

        if N == 0:
            return pv.PolyData(), lo, hi, use_categorical

        half    = VOXEL_SIZE * (1.0 - VOXEL_GAP) / 2.0
        centers = (vox_idx + 0.5) * VOXEL_SIZE  # (N, 3)

        points = (centers[:, np.newaxis, :] + _CORNERS[np.newaxis, :, :] * (half * 2)).reshape(-1, 3)

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
        color_per_vox   = color_arr[ci, cj, ck]
        density_per_vox = density_arr[ci, cj, ck]

        mesh.cell_data['scalars'] = np.repeat(color_per_vox,   6)
        mesh.cell_data['alpha']   = np.clip(0.1 + 0.9 * np.repeat(density_per_vox, 6), 0.0, 1.0)

        return mesh, lo, hi, use_categorical


# -----------------------------------------------------------------------
# load data
with h5py.File(FILEPATH, 'r') as f:
    full_data = {prop: np.array(f[prop]) for prop in PROPERTIES}
    for prop in INFO_EXTRA_PROPERTIES:
        if prop in f:
            full_data[prop] = np.array(f[prop])

    opt = {}
    opt['xlo']            = float(f.attrs.get('xlo', 0.0))
    opt['ylo']            = float(f.attrs.get('ylo', 0.0))
    opt['zlo']            = float(f.attrs.get('zlo', 0.0))
    opt['voxel_size_real']= float(f.attrs.get('voxel_size', 10.0))

    for key in ('si_surface', 'si_surface_mask'):
        if key in f:
            opt[key] = np.array(f[key])

    if 'crater' in f:
        cg = f['crater']
        for key in ('reference_x', 'si_crater_mask', 'sphere_center', 'sphere_radius'):
            if key in cg:
                opt[f'crater/{key}'] = np.array(cg[key])

grid = VoxelGrid(full_data, opt)

# -----------------------------------------------------------------------
# plotter setup
pl = pv.Plotter(window_size=(1000, 800))
pl.set_background('white')
pl.enable_parallel_projection()

mesh, lo, hi, is_cat = grid.build_mesh()
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

pl.show_bounds(color='black', xtitle='X', ytitle='Y', ztitle='Z')

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
play_speed_idx    = [0]
jump_buffer       = ['']
timer_id          = [None]
help_visible      = [True]

# overlay toggles
density_opacity    = [True]   # True = opacity mapped from density, False = fully opaque
show_mesh          = [True]
show_surface_plane = [False]
show_sphere_cap    = [False]
crater_color_mode  = [False]
show_velocity_vecs = [False]
velocity_scale     = [1.0]   # global multiplier on top of [ARROW_LENGTH_MIN, ARROW_LENGTH_MAX] normalization

# voxel selection / inspection
selected_voxel = [None]   # (ix, iy, iz) or None

# overlay actors
surface_plane_actor = [None]
sphere_cap_actor    = [None]
velocity_actor      = [None]


# -----------------------------------------------------------------------
# text actors
def _info_text():
    speed_ms = PLAY_SPEEDS_MS[play_speed_idx[0]]
    play_str = f"[PLAYING {speed_ms}ms/frame]" if playing[0] else "[paused]"
    jump_str = f"  Jump: {jump_buffer[0]}_" if jump_buffer[0] else ''
    return (f"Property : {grid.current_property}\n"
            f"Timestep : {grid.current_timestep} / {grid.ntimesteps - 1}  {play_str}{jump_str}")

def _voxel_info_text():
    if selected_voxel[0] is None:
        return ""
    ix, iy, iz = selected_voxel[0]
    t = grid.current_timestep
    if not (0 <= ix < grid.nx and 0 <= iy < grid.ny and 0 <= iz < grid.nz):
        return ""
    lines = [f"Voxel ({ix}, {iy}, {iz})"]
    for prop in PROPERTIES + INFO_EXTRA_PROPERTIES:
        if prop not in grid.data:
            continue
        val = grid.data[prop][t, ix, iy, iz]
        if prop == 'v_COM':
            vx, vy, vz = (float(v) for v in val)
            mag = float(np.linalg.norm(val))
            lines.append(f"v_COM: ({vx:.2f}, {vy:.2f}, {vz:.2f})  |v|={mag:.2f}")
        elif prop == 'voxel_type':
            lines.append(f"voxel_type: {int(val)}")
        else:
            lines.append(f"{prop}: {float(val):.4g}")
    return "\n".join(lines)


info_actor = pl.add_text(_info_text(), position='upper_left', font_size=12, color='black')
help_actor = pl.add_text(KEY_HELP, position='lower_left', font_size=10, color='black')
voxel_info_actor = pl.add_text(_voxel_info_text(), position='upper_right', font_size=11, color='black')


# -----------------------------------------------------------------------
# overlay helpers

def _update_surface_plane():
    if surface_plane_actor[0] is not None:
        pl.remove_actor(surface_plane_actor[0])
        surface_plane_actor[0] = None

    if not show_surface_plane[0]:
        return

    t = grid.current_timestep
    ref_x_real = None
    if grid.crater_reference_x is not None:
        val = float(grid.crater_reference_x[t])
        if not np.isnan(val):
            ref_x_real = val
    if ref_x_real is None and grid.si_surface is not None:
        val = float(np.nanmedian(grid.si_surface[t]))
        if not np.isnan(val):
            ref_x_real = val
    if ref_x_real is None:
        return

    ref_x_d = grid.x_angstrom_to_display(ref_x_real)
    cy_d = grid.ny * VOXEL_SIZE / 2.0
    cz_d = grid.nz * VOXEL_SIZE / 2.0

    plane = pv.Plane(
        center=(ref_x_d, cy_d, cz_d),
        direction=(1, 0, 0),
        i_size=grid.ny * VOXEL_SIZE * 1.5,
        j_size=grid.nz * VOXEL_SIZE * 1.5,
    )
    surface_plane_actor[0] = pl.add_mesh(
        plane, color='lightgray', opacity=0.4, show_scalar_bar=False, pickable=False,
    )


def _update_sphere_cap():
    if sphere_cap_actor[0] is not None:
        pl.remove_actor(sphere_cap_actor[0])
        sphere_cap_actor[0] = None

    if not show_sphere_cap[0]:
        return

    if grid.crater_sphere_center is None or grid.crater_sphere_radius is None:
        return

    t = grid.current_timestep
    center = grid.crater_sphere_center[t]
    radius = float(grid.crater_sphere_radius[t])

    if np.any(np.isnan(center)) or np.isnan(radius):
        return

    ref_x_real = None
    if grid.crater_reference_x is not None:
        val = float(grid.crater_reference_x[t])
        if not np.isnan(val):
            ref_x_real = val
    if ref_x_real is None:
        return

    cx_d, cy_d, cz_d = grid.angstrom_to_display(float(center[0]), float(center[1]), float(center[2]))
    r_d    = grid.radius_to_display(radius)
    ref_x_d = grid.x_angstrom_to_display(ref_x_real)

    sphere = pv.Sphere(radius=r_d, center=(cx_d, cy_d, cz_d),
                       theta_resolution=60, phi_resolution=60)
    cap = sphere.clip(normal=(1, 0, 0), origin=(ref_x_d, cy_d, cz_d), invert=True)

    if cap.n_points > 0:
        sphere_cap_actor[0] = pl.add_mesh(
            cap, color='coral', opacity=0.5, show_scalar_bar=False, pickable=False,
        )


def _update_velocity_vectors():
    if velocity_actor[0] is not None:
        pl.remove_actor(velocity_actor[0])
        velocity_actor[0] = None

    if not show_velocity_vecs[0]:
        return

    t = grid.current_timestep
    visible = ~grid._compute_mask()

    vox_idx = np.argwhere(visible)
    if len(vox_idx) == 0:
        return

    ci, cj, ck = vox_idx[:, 0], vox_idx[:, 1], vox_idx[:, 2]
    v_com = grid.data['v_COM'][t][ci, cj, ck]  # (N, 3)

    valid_v = ~np.any(np.isnan(v_com), axis=1)
    if not valid_v.any():
        return

    centers = ((vox_idx[valid_v] + 0.5) * VOXEL_SIZE).astype(np.float32)
    v_com   = v_com[valid_v].astype(np.float32)
    magnitudes = np.linalg.norm(v_com, axis=1).astype(np.float32)

    # map magnitudes to arrow lengths in [ARROW_LENGTH_MIN, ARROW_LENGTH_MAX]
    mag_min, mag_max = magnitudes.min(), magnitudes.max()
    if mag_max > mag_min:
        t_norm = (magnitudes - mag_min) / (mag_max - mag_min)
    else:
        t_norm = np.full_like(magnitudes, 0.5)
    arrow_lengths = (ARROW_LENGTH_MIN + t_norm * (ARROW_LENGTH_MAX - ARROW_LENGTH_MIN)).astype(np.float32)

    points = pv.PolyData(centers)
    points['vectors']      = v_com
    points['arrow_length'] = arrow_lengths
    points['magnitude']    = magnitudes

    arrows = points.glyph(orient='vectors', scale='arrow_length', factor=velocity_scale[0])
    velocity_actor[0] = pl.add_mesh(
        arrows, scalars='magnitude', cmap='coolwarm', show_scalar_bar=False, pickable=False,
    )


# -----------------------------------------------------------------------
def refresh():
    global actor
    mesh, lo, hi, is_cat = grid.build_mesh(crater_colors=crater_color_mode[0])
    pl.remove_actor(actor)

    opacity = 'alpha' if density_opacity[0] else 1.0

    if is_cat:
        actor = pl.add_mesh(
            mesh,
            scalars='scalars',
            clim=[0, 2],
            cmap=['steelblue', 'orange', 'red'],
            n_colors=3,
            show_scalar_bar=False,
            opacity=opacity,
        )
    else:
        actor = pl.add_mesh(
            mesh,
            scalars='scalars',
            clim=[lo, hi],
            cmap='coolwarm',
            show_scalar_bar=False,
            opacity=opacity,
        )
        sbar.SetTitle(grid.current_property)
        sbar.GetLookupTable().SetTableRange(lo, hi)

    actor.SetVisibility(show_mesh[0])

    _update_surface_plane()
    _update_sphere_cap()
    _update_velocity_vectors()

    info_actor.SetText(2, _info_text())
    voxel_info_actor.SetText(3, _voxel_info_text())
    pl.render()
    pl.update()


def _set_timestep(t):
    current_timestep[0] = max(0, min(t, grid.ntimesteps - 1))
    grid.current_timestep = current_timestep[0]


# -----------------------------------------------------------------------
# playback timer
iren = pl.iren
vtk_iren = iren.interactor

def _on_timer(obj, event):
    if playing[0]:
        if current_timestep[0] >= grid.ntimesteps - 1:
            _stop_playback()
            return
        _set_timestep(current_timestep[0] + 1)
        refresh()

vtk_iren.AddObserver('TimerEvent', _on_timer)


def _start_playback():
    playing[0] = True
    timer_id[0] = vtk_iren.CreateRepeatingTimer(PLAY_SPEEDS_MS[play_speed_idx[0]])

def _stop_playback():
    playing[0] = False
    if timer_id[0] is not None:
        vtk_iren.DestroyTimer(timer_id[0])
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

def next_prop():
    current_prop_idx[0] = (current_prop_idx[0] + 1) % len(PROPERTIES)
    grid.current_property = PROPERTIES[current_prop_idx[0]]
    refresh()

def prev_prop():
    current_prop_idx[0] = (current_prop_idx[0] - 1) % len(PROPERTIES)
    grid.current_property = PROPERTIES[current_prop_idx[0]]
    refresh()

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

# slicing (moved from s to k)
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

# overlay toggles
def toggle_density_opacity():
    density_opacity[0] = not density_opacity[0]
    refresh()

def toggle_mesh():
    show_mesh[0] = not show_mesh[0]
    actor.SetVisibility(show_mesh[0])
    pl.render()
    pl.update()

def toggle_surface_plane():
    show_surface_plane[0] = not show_surface_plane[0]
    _update_surface_plane()
    pl.render()
    pl.update()

def toggle_sphere_cap():
    show_sphere_cap[0] = not show_sphere_cap[0]
    _update_sphere_cap()
    pl.render()
    pl.update()

def toggle_crater_colors():
    crater_color_mode[0] = not crater_color_mode[0]
    refresh()

def toggle_velocity_vectors():
    show_velocity_vecs[0] = not show_velocity_vecs[0]
    _update_velocity_vectors()
    pl.render()
    pl.update()

def scale_velocity_down():
    velocity_scale[0] *= 0.5
    print(f"Velocity arrow scale: {velocity_scale[0]:.3f}")
    if show_velocity_vecs[0]:
        _update_velocity_vectors()
        pl.render()
        pl.update()

def scale_velocity_up():
    velocity_scale[0] *= 2.0
    print(f"Velocity arrow scale: {velocity_scale[0]:.3f}")
    if show_velocity_vecs[0]:
        _update_velocity_vectors()
        pl.render()
        pl.update()

def toggle_help():
    help_visible[0] = not help_visible[0]
    help_actor.SetVisibility(help_visible[0])
    pl.render()
    pl.update()

# camera view presets
def _set_camera_view(axis):
    cx = grid.nx * VOXEL_SIZE / 2.0
    cy = grid.ny * VOXEL_SIZE / 2.0
    cz = grid.nz * VOXEL_SIZE / 2.0
    dist = max(grid.nx, grid.ny, grid.nz) * VOXEL_SIZE * 3.0
    offsets = {
        '+x': (-dist, 0, 0), '-x': (dist, 0, 0),
        '+y': (0, -dist, 0), '-y': (0, dist, 0),
        '+z': (0, 0, -dist), '-z': (0, 0, dist),
    }
    ups = {
        '+x': (0, 0, 1), '-x': (0, 0, 1),
        '+y': (0, 0, 1), '-y': (0, 0, 1),
        '+z': (0, 1, 0), '-z': (0, 1, 0),
    }
    ox, oy, oz = offsets[axis]
    pl.camera.focal_point = (cx, cy, cz)
    pl.camera.position    = (cx + ox, cy + oy, cz + oz)
    pl.camera.up          = ups[axis]
    pl.render()
    pl.update()


# voxel selection / inspection
def _on_pick_point(point):
    if point is None:
        return
    ix = int(point[0] // VOXEL_SIZE)
    iy = int(point[1] // VOXEL_SIZE)
    iz = int(point[2] // VOXEL_SIZE)
    if not (0 <= ix < grid.nx and 0 <= iy < grid.ny and 0 <= iz < grid.nz):
        return
    selected_voxel[0] = (ix, iy, iz)
    voxel_info_actor.SetText(3, _voxel_info_text())
    pl.render()
    pl.update()

def clear_voxel_selection():
    selected_voxel[0] = None
    voxel_info_actor.SetText(3, "")
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
pl.add_key_event('k',            toggle_slice)
pl.add_key_event('w',            lambda: None)
pl.add_key_event('x',            slice_axis_yz)
pl.add_key_event('y',            slice_axis_xz)
pl.add_key_event('z',            slice_axis_xy)
pl.add_key_event('period',       slice_forward)
pl.add_key_event('comma',        slice_backward)
pl.add_key_event('equal',        slice_thicker)
pl.add_key_event('minus',        slice_thinner)
pl.add_key_event('a',            toggle_density_opacity)
pl.add_key_event('m',            toggle_mesh)
pl.add_key_event('p',            toggle_surface_plane)
pl.add_key_event('s',            toggle_sphere_cap)
pl.add_key_event('c',            toggle_crater_colors)
pl.add_key_event('v',            toggle_velocity_vectors)
pl.add_key_event('less',         scale_velocity_down)
pl.add_key_event('greater',      scale_velocity_up)
pl.add_key_event('h',            toggle_help)
pl.add_key_event('i',            clear_voxel_selection)
pl.add_key_event('F1',           lambda: _set_camera_view('+x'))
pl.add_key_event('F2',           lambda: _set_camera_view('-x'))
pl.add_key_event('F3',           lambda: _set_camera_view('+y'))
pl.add_key_event('F4',           lambda: _set_camera_view('-y'))
pl.add_key_event('F5',           lambda: _set_camera_view('+z'))
pl.add_key_event('F6',           lambda: _set_camera_view('-z'))

# use_picker=False avoids passing (point, picker) to callback (would raise TypeError)
pl.enable_point_picking(callback=_on_pick_point, use_picker=False,
                         show_message=False, left_clicking=True)

pl.show()
