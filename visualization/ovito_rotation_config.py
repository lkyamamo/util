"""
User-tunable knobs for OVITO trajectory rendering.

Edit paths, camera, GIF output, particle display radii, and pairwise bond
distance cutoffs (neighbors with separation below the cutoff get a bond).
"""

from pathlib import Path

# --- Camera (trajectory coords, usually Å) ---
# Look-at center. Orbit: camera lies on a circle (radius CAMERA_ORBIT_RADIUS) around CAMERA_LOOK_AT.
# Polar angle φ is colatitude from +Z; azimuth θ advances each subframe:
# θ = radians(CAMERA_ORBIT_START_ANGLE_DEG) + 2π * i_rot / ROTATION_SUBFRAMES.
# Eye offset from look-at: (−horizontal along sin/cos of θ plus vertical component from φ).
CAMERA_LOOK_AT = (0.0, 0.0, 0.0)
CAMERA_ORBIT_RADIUS = 50.0
# Polar / colatitude from +Z (degrees). 90 = horizontal ring; 0 = directly on +Z above look_at.
CAMERA_ORBIT_POLAR_ANGLE_DEG = 15.0
# Azimuth at frame 0 (degrees). 0 matches the previous default ring position in XY.
CAMERA_ORBIT_START_ANGLE_DEG = 0.0

# --- Data source ---
DUMP_PATH = Path(
    "/Users/loganyamamoto/Desktop/Research/grants/geo_sciences/bubble_collapse/data/"
    "systems/sioh_cov/a-SiO/2_surface/N6192/z-plane/20260418/2.npt_nvt/"
    "2.int_npt_nvt.dump"
)

# Simulation snapshot indices (0 .. pipeline.num_frames - 1).
TRAJECTORY_FRAMES = range(15000, 20000, 1000)

# --- LAMMPS type column "type" → chemical symbol (ParticleType names in OVITO match dump strings
# before remap, e.g. "1"; this map renames them for ATOM_DISPLAY_RADII / BOND_PAIR_CUTOFFS).
# First-frame counts for this DUMP_PATH (6192 atoms): type 1 = 864, type 2 = 2928, type 3 = 2400.
# Confirm correspondence with Masses/pair_coeff in your LAMMPS data — swap entries if chemistry is wrong.
LAMMPS_TYPE_ID_TO_SYMBOL = {
    1: "Si",
    2: "O",
    3: "H",
}

# --- Particle styles (sphere radii per ParticleType.name) ---
# Symbols must match chemical labels after LAMMPS_TYPE_ID_TO_SYMBOL (Å).
ATOM_DISPLAY_RADII = {
    "Si": 1.1,
    "O": 0.75,
    "H": 0.42,
}

# --- Bonds: maximum center–center distance per type pair (Å) ---
# Only one ordering per pair is required; the pipeline mirrors (a,b) → (b,a) by default.
# Omit pairs that should never bond — the pipeline skips cutoffs <= 0.
BOND_PAIR_CUTOFFS = {
    ("Si", "O"): 1.9,
    ("O", "H"): 1.2,
}

# Cylinder radius for drawn bonds (Å). Use None for OVITO default width.
BOND_VISUAL_RADIUS = None

ROTATION_SUBFRAMES = 5000
GIF_FPS = 24
GIF_BASENAME = "rotation"
KEEP_PNG_FRAMES = False

RENDER_SIZE = (800, 600)
