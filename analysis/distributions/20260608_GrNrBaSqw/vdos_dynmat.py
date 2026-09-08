"""
vdos_dynmat.py — Vibrational Density of States from the LAMMPS dynamical matrix.

An independent, harmonic route to the same quantity vdos.py produces from the
velocity autocorrelation of an MD trajectory. Where vdos.py measures what the
atoms actually did at the simulation temperature, this measures the curvature of
the potential-energy surface at a single minimum: LAMMPS minimizes the final MD
configuration into its nearest inherent structure, builds the mass-weighted
force-constant matrix by finite displacement (dynamical_matrix command), and this
script diagonalizes it and histograms the eigenvalues.

Run both on the same trajectory and overlay the CSVs (identical column layout by
design). Peak positions should agree; this curve is sharper — no thermal
broadening, no correlation-length resolution limit — and usually slightly
blue-shifted, since it has no anharmonic softening.

QUICK START
-----------
1. Run the pipeline with RUN_VDOS_DYNMAT=1, which sets LAMMPS's RUN_DYNMAT=1 and
   produces dynmat.dat in the run directory.
2. Set VDOS_DYNMAT_MAX_FREQUENCY (the one required variable).
3. Run:  python vdos_dynmat.py

CONFIGURATION
-------------
Each setting comes from exactly one environment variable — no fallback names,
and no silent defaults for anything that changes the result.

REQUIRED (unset or empty is a hard error; the script will not guess):
  VDOS_DYNMAT_MAX_FREQUENCY  upper limit of the DOS grid, in VDOS_DYNMAT_XUNIT.
                             Required for the same reason vdos.py's
                             VDOS_MAX_FREQUENCY_EV is: a default that is too low
                             silently truncates the spectrum instead of failing

Optional:
  DYNMAT_FILE                matrix written by LAMMPS      (default dynmat.dat)
  DYNMAT_BINARY              yes | no — must match what    (default no)
                             the LAMMPS run used
  TRAJ                       dump read for the element     (default dump.lammpstrj)
                             column; first frame only
  VDOS_DYNMAT_XUNIT          meV | THz | cm-1 | eV         (default meV)
  VDOS_DYNMAT_BINS           frequency grid points         (default 500)
  VDOS_DYNMAT_SMEARING       Gaussian FWHM in XUNIT;       (default 0)
                             0 = plain histogram
  VDOS_DYNMAT_MATRIX_STYLE   regular | eskm — must match   (default regular)
                             dynamical_matrix's style
  VDOS_DYNMAT_NORMALIZATION  phonon | unit_area            (default phonon)
  VDOS_DYNMAT_PARTIAL        yes | no — per-element curves (default yes)
  VDOS_DYNMAT_ASR            none | simple                 (default none)
  VDOS_DYNMAT_THREADS        BLAS threads, 0/unset = leave (default 0)
                             OMP_NUM_THREADS alone
  VDOS_DYNMAT_OUTPUT         output basename               (default vdos_dynmat)

Mode-character analysis (see MODE CHARACTER below), all optional:
  VDOS_DYNMAT_CHARACTER          yes | no                  (default no)
  DYNMAT_REF_TRAJ                minimized coordinates     (default dynmat_ref.lammpstrj)
  VDOS_DYNMAT_BRIDGE_ELEMENT     the bridging atom         (default O)
  VDOS_DYNMAT_NEIGHBOR_ELEMENT   its two neighbours        (default Si)
  VDOS_DYNMAT_BOND_CUTOFF        bridge-neighbour max, A   (default 2.2)

BINS only changes how smooth the curve looks, so it has a default; MAX_FREQUENCY
changes which physics is in the plot at all, so it does not.

WHERE THE ELEMENT LABELS COME FROM
----------------------------------
dynamical_matrix iterates atoms by global ID, so matrix row 3i+alpha belongs to
atom ID i+1. The pipeline's dump.lammpstrj is written with `dump_modify sort id`,
so its first frame lists atoms in that same order. Since the dynamical matrix is
built on the final configuration of the run that wrote that dump — same atoms,
same IDs, no replicate in between — the dump's element column indexes the matrix
rows directly, and no extra output from LAMMPS is needed. The 9*N^2 element count
of the matrix is checked against the dump's atom count, which is what catches this
assumption being broken later (a group other than `all`, a replicate after the
dump, a mismatched pair of files).

METHOD
------
1. Read the matrix as a flat stream of 9*N^2 floats and reshape to (3N, 3N):
   LAMMPS writes it three numbers per line, ordered j fastest, then alpha, then i.
2. Symmetrize, D = (D + D.T)/2. Finite differences leave a small asymmetry; the
   ratio max|D - D.T| / max|D| is reported as a quality check.
3. Optionally impose the acoustic sum rule (ASR='simple'): un-mass-weight, zero
   each atom's force-constant row sum onto its own diagonal block, re-weight.
   This pushes the three acoustic modes to exactly zero. Off by default because
   it modifies the matrix LAMMPS produced.
4. Diagonalize. eigh when partial DOS is wanted, eigvalsh when it is not (which
   halves both peak memory and runtime).
5. Convert eigenvalues to frequencies. D is mass-weighted, so in `units metal`:
       regular : lambda in eV/(A^2*amu), nu[THz] = sqrt(lambda) * 15.6333042
       eskm    : lambda in 1/ps^2,       nu[THz] = sqrt(lambda) / (2*pi)
   The two agree — LAMMPS's eskm conversion factor is 9648.5, and
   sqrt(9648.5)/(2*pi) = 15.6333. Negative eigenvalues (imaginary modes) are
   reported as negative frequencies and excluded from the spectrum; see
   DIAGNOSTICS. Modes within one bin width of zero — the acoustic translations,
   scattered across zero by finite-difference noise — are snapped to exactly
   zero first, so they all land in the first bin instead of half of them falling
   out of range and costing the spectrum weight it should have.
6. Bin the real modes over [0, MAX_FREQUENCY] with BINS bins — unweighted for the
   total, weighted by each mode's per-element participation for the partials.
   With SMEARING > 0, sum Gaussians of that FWHM straight onto the bin centers
   instead, which avoids the artifacts of binning and then smoothing.

Partial weights are w[e,m] = sum over atoms i of element e and directions alpha of
|v[3i+alpha, m]|^2. D is mass-weighted, so its eigenvectors are orthonormal and
these weights sum to 1 for every mode.

WEIGHTING
---------
VDOS_DYNMAT_WEIGHTING sets how much each element contributes to a total, and is
a different axis from VDOS_DYNMAT_NORMALIZATION below: weighting decides the
relative species contributions, normalization decides the sum rule. It takes a
SEMICOLON-separated list and emits one DoS(Total_<weighting>) column per entry:

  unity       k_el = 1                                    (the plain total)
  coherent    k_el = sigma_coh_el / m_el
  incoherent  k_el = sigma_inc_el / m_el
  total       k_el = (sigma_coh_el + sigma_inc_el) / m_el

k_el multiplies each element's per-mode participation before binning, so the
weight enters per mode rather than by rescaling finished curves, and is
normalized so sum_el k_el*c_el = 1 — which keeps the 'phonon' 3-per-atom sum
rule intact. Everything but 'unity' needs VDOS_DYNMAT_PARTIAL=yes, since without
eigenvectors there are no per-element participations to weight.

Careful when comparing with vdos.py: its partials are per-element *shapes* and
its weight carries the concentration, w_el = c_el*sigma/m. Here the partials
already carry it (each integrates to 3*c_el under 'phonon'), so k_el must not
include another factor of c_el. Both scripts end up giving each element the same
share of the total, c_el(sigma/m)_el / sum(c*sigma/m), which is what the startup
table prints and what makes the two directly comparable.

Why sigma/m and not a scattering length: what inelastic neutron scattering
measures is the generalized DOS, in which the one-phonon incoherent cross
section carries a factor sigma/m per species. That is a different quantity from
the coherent scattering length b that weights diffraction, and sigma_inc cannot
be derived from b_coh at all. The mass comes from ELEMENT_MASSES — the same
masses the dynamical matrix was mass-weighted with, so the weighting cannot
drift from the physics the eigenvectors describe.

Not applied: the Debye-Waller factor exp(-2W), which is Q-dependent while this
DOS is not Q-resolved.

Hydrogen is refused: with H or D present, any weighting other than 'unity' exits
with an explanation, for the same reason as vdos.py — the isotope cannot be read
off a dump that labels both 'H', and the two differ by ~21x in sigma/m.

NORMALIZATION
-------------
'phonon' (default) — every curve is divided by the atom count, so the total
   integrates to exactly 3 per atom and each partial to 3*N_el/N_total. The
   partials therefore sum to the total exactly. This holds with smearing on as
   well, because each mode's Gaussian is normalized over the plotted range
   rather than analytically — see bin_spectrum. It does *not* hold if modes lie
   above MAX_FREQUENCY, which the script warns about explicitly.
'unit_area' — every curve (each partial and the total) is independently rescaled
   so integral(curve) dnu = 1. Same convention as vdos.py's option of that name.

Note that vdos.py's 'phonon' normalization reaches ~3 per atom by a different
route (msd.cpp's 6/pi mole-fraction-weighted sum over C(0)=1 VACFs), so absolute
y-values from the two scripts are close but not identical. Peak positions and
curve shapes are what compare directly; use 'unit_area' on both if you want the
overlay to line up in height as well.

DIAGNOSTICS
-----------
Three zero modes are expected — the acoustic translations. Any imaginary mode
surviving the zero-snap in METHOD step 5 means the minimization did not reach a
true minimum, so the script prints the six lowest signed frequencies, the number
snapped to zero, and the number still negative, rather than folding them silently
into the spectrum. If imaginary modes appear, tighten DYNMAT_MIN_FTOL or reduce
DYNMAT_DISPLACEMENT in the LAMMPS stage before trusting anything downstream.

PARALLELIZATION AND MEMORY
--------------------------
The expensive half of this method is on the LAMMPS side and is already MPI
parallel: dynamical_matrix's 6N force evaluations use the same domain
decomposition as any other LAMMPS force computation, so they scale across all the
ranks the trajectory job was given. Its one caveat is memory rather than speed —
LAMMPS gathers 9*N doubles on *every* rank, so per-rank memory does not fall as
ranks are added.

This script is single-node and thread parallel. np.linalg.eigh is a LAPACK
dsyevd call that the underlying BLAS (OpenBLAS or MKL) threads across cores;
distribution_submit.slurm's OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK is what drives
it, and VDOS_DYNMAT_THREADS overrides that for hand runs. There is no distributed
diagonalization — that would need ScaLAPACK or ELPA — so a cell whose matrix does
not fit one node is out of reach by this route. Everything after the
diagonalization (weights, binning) is negligible by comparison.

The matrix is (3N)^2 float64: 1.9 GB at N=5184 (replicate 6 6 6). Symmetrizing
and diagonalizing peak at roughly three times that with partials on, or two times
with VDOS_DYNMAT_PARTIAL=no. Sized for a --mem=0 --exclusive node.

OUTPUT
------
- <date>_vdos_dynmat.csv — freq_meV, freq_THz, freq_cm-1, freq_eV, then
                           DoS(<element>) per element and DoS(Total). Same column
                           layout as vdos.py's CSV so the two can be overlaid.
- <date>_vdos_dynmat.png — all curves on one axes, x-axis in XUNIT

DEPENDENCIES
------------
  pip install numpy matplotlib
"""

import os

# Thread count must be set BEFORE numpy is imported — OpenBLAS and MKL both read
# their thread environment once, at load time, and ignore later changes. This is
# the whole parallelization story on the Python side: np.linalg.eigh is a LAPACK
# dsyevd call, which the underlying BLAS threads across cores on one node. There
# is no MPI here (that would need ScaLAPACK/ELPA), so the diagonalization is
# bounded by one node's core count.
#
# distribution_submit.slurm already exports OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK,
# so the normal pipeline path is threaded without setting anything; this variable
# exists for running the script by hand, and to make the choice visible in the log.
_THREADS = os.environ.get("VDOS_DYNMAT_THREADS", "")
if _THREADS not in ("", "0"):
    for _var in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS",
                 "VECLIB_MAXIMUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[_var] = _THREADS

from datetime import date

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# =============================================================================
# CONFIGURATION — edit these variables between runs
# =============================================================================

# Matrix written by LAMMPS's dynamical_matrix command, and the dump whose first
# frame supplies the element label for each matrix row. Both are symlinked into
# this directory by submit_pipeline.sh.
DYNMAT_FILE = os.environ.get("DYNMAT_FILE", "dynmat.dat")
DUMP_FILE   = os.environ.get("TRAJ", "dump.lammpstrj")

# Required settings have no default: an unset or empty value is a hard error,
# never a silently substituted number. Missing ones are collected so a single run
# reports every one of them at once instead of failing one at a time. Same
# helpers, and the same reasoning, as vdos.py and msd.py.
_MISSING = []

def _require(name, description):
    """Value of `name`, or None after recording it as missing."""
    value = os.environ.get(name, "")
    if value == "":
        _MISSING.append(f"  {name:<26} {description}")
        return None
    return value

def _check_required(script):
    if _MISSING:
        raise SystemExit(
            f"{script}: required environment variable(s) not set:\n"
            + "\n".join(_MISSING)
            + "\n\nThese determine the numbers this script produces, so it will not "
              "guess them.\nSet them in submit_pipeline.conf (driven by "
              "submit_pipeline.sh / submit_pipeline_local.sh)\nor in the Analysis "
              "parameters block of distribution_run.sh / distribution_submit.slurm."
        )

def _env(name, default):
    """Optional setting: value of `name` if non-empty, else `default`. Used only
    where the default cannot silently distort the result — a documented algorithm
    choice, or a grid density that changes smoothness but not physics."""
    value = os.environ.get(name, "")
    return value if value != "" else default


# x-axis unit. Set before MAX_FREQUENCY because that value is expressed in it.
XUNIT = _env("VDOS_DYNMAT_XUNIT", "meV")
_VALID_XUNITS = ('meV', 'THz', 'cm-1', 'eV')
if XUNIT not in _VALID_XUNITS:
    raise ValueError(f"Unknown VDOS_DYNMAT_XUNIT={XUNIT!r}; use one of {_VALID_XUNITS}.")

# --- Required. Upper limit of the DOS grid, in XUNIT. A default here would
# silently cut the spectrum short — the same failure mode vdos.py's
# VDOS_MAX_FREQUENCY_EV was made required to avoid.
_MAX_FREQUENCY_ENV = _require("VDOS_DYNMAT_MAX_FREQUENCY", f"upper limit of the DOS grid, in {XUNIT}")
_check_required("vdos_dynmat.py")

MAX_FREQUENCY = float(_MAX_FREQUENCY_ENV)
if MAX_FREQUENCY <= 0:
    raise ValueError(f"VDOS_DYNMAT_MAX_FREQUENCY must be positive, got {MAX_FREQUENCY}.")

# Number of frequency bins. Only changes how smooth the curve looks.
BINS = int(_env("VDOS_DYNMAT_BINS", "500"))
if BINS < 2:
    raise ValueError(f"VDOS_DYNMAT_BINS must be at least 2, got {BINS}.")

# Gaussian FWHM in XUNIT applied to the discrete mode spectrum; 0 = plain
# histogram. A finite cell has finitely many modes, so some smearing is usually
# what makes the curve readable.
SMEARING = float(_env("VDOS_DYNMAT_SMEARING", "0"))
if SMEARING < 0:
    raise ValueError(f"VDOS_DYNMAT_SMEARING must be >= 0, got {SMEARING}.")

# Which dynamical_matrix style wrote the file — sets the eigenvalue-to-frequency
# conversion, so a mismatch rescales the whole spectrum.
MATRIX_STYLE = _env("VDOS_DYNMAT_MATRIX_STYLE", "regular")
if MATRIX_STYLE not in ("regular", "eskm"):
    raise ValueError(f"Unknown VDOS_DYNMAT_MATRIX_STYLE={MATRIX_STYLE!r}; use 'regular' or 'eskm'.")

# Whether the file is raw float64 (binary yes) or text (binary no) — must match
# the `binary` keyword the LAMMPS run used.
BINARY = _env("DYNMAT_BINARY", "no")
if BINARY not in ("yes", "no"):
    raise ValueError(f"Unknown DYNMAT_BINARY={BINARY!r}; use 'yes' or 'no'.")

# How curves are scaled — see NORMALIZATION in the module docstring.
NORMALIZATION = _env("VDOS_DYNMAT_NORMALIZATION", "phonon")
if NORMALIZATION not in ("phonon", "unit_area"):
    raise ValueError(f"Unknown VDOS_DYNMAT_NORMALIZATION={NORMALIZATION!r}; use 'phonon' or 'unit_area'.")

# How much each element contributes to a total — see WEIGHTING in the module
# docstring. SEMICOLON-separated; one total column per entry. A different axis
# from VDOS_DYNMAT_NORMALIZATION above, which is a sum rule, not a weighting.
WEIGHTING = _env("VDOS_DYNMAT_WEIGHTING", "unity")

# Per-element partial DOS needs eigenvectors; turning it off lets the script use
# eigvalsh instead of eigh, halving peak memory and runtime on large cells.
PARTIAL = _env("VDOS_DYNMAT_PARTIAL", "yes")
if PARTIAL not in ("yes", "no"):
    raise ValueError(f"Unknown VDOS_DYNMAT_PARTIAL={PARTIAL!r}; use 'yes' or 'no'.")

# Acoustic sum rule — see METHOD step 3. Off by default: it modifies the matrix
# LAMMPS produced, which is a choice worth making explicitly.
ASR = _env("VDOS_DYNMAT_ASR", "none")
if ASR not in ("none", "simple"):
    raise ValueError(f"Unknown VDOS_DYNMAT_ASR={ASR!r}; use 'none' or 'simple'.")

# Mode-character analysis — see MODE CHARACTER in the module docstring. Off by
# default: it needs positions at the minimum, which only exist if the LAMMPS
# stage wrote them, and it costs an extra pass over the eigenvectors.
CHARACTER = _env("VDOS_DYNMAT_CHARACTER", "no")
if CHARACTER not in ("yes", "no"):
    raise ValueError(f"Unknown VDOS_DYNMAT_CHARACTER={CHARACTER!r}; use 'yes' or 'no'.")

# Coordinates at the minimized geometry, written by the LAMMPS stage's
# write_dump immediately after `minimize`. Only read when CHARACTER='yes' — this
# has to be the post-minimization geometry, since the whole analysis is a
# projection onto bond directions and the relaxation rotates them.
DYNMAT_REF_TRAJ = os.environ.get("DYNMAT_REF_TRAJ", "dynmat_ref.lammpstrj")

# The local frame is built at each BRIDGE atom from its two NEIGHBOR-element
# neighbours: for a-SiO2 that is an oxygen bridging two silicons, the Si-O-Si
# unit whose three orthogonal displacement directions are the standard band
# assignment for silica glass.
BRIDGE_ELEMENT   = _env("VDOS_DYNMAT_BRIDGE_ELEMENT", "O")
NEIGHBOR_ELEMENT = _env("VDOS_DYNMAT_NEIGHBOR_ELEMENT", "Si")
BOND_CUTOFF      = float(_env("VDOS_DYNMAT_BOND_CUTOFF", "2.2"))   # Angstrom
if BOND_CUTOFF <= 0:
    raise ValueError(f"VDOS_DYNMAT_BOND_CUTOFF must be positive, got {BOND_CUTOFF}.")

# Output basename; .csv and .png are appended (set either to None to skip).
_OUTPUT_BASE = _env("VDOS_DYNMAT_OUTPUT", "vdos_dynmat")
OUTPUT_CSV   = f"{_OUTPUT_BASE}.csv"
OUTPUT_PLOT  = f"{_OUTPUT_BASE}.png"
# Written only when CHARACTER='yes'.
OUTPUT_MODES     = f"{_OUTPUT_BASE}_modes.csv"
OUTPUT_CHARACTER = f"{_OUTPUT_BASE}_character.png"

# Plot appearance
PLOT_DPI = 150

# =============================================================================
# END CONFIGURATION
# =============================================================================

# Prepend today's date (YYYYMMDD_) to every output filename.
def _dated(filename):
    return None if filename is None else f"{date.today():%Y%m%d}_{filename}"

OUTPUT_CSV       = _dated(OUTPUT_CSV)
OUTPUT_PLOT      = _dated(OUTPUT_PLOT)
OUTPUT_MODES     = _dated(OUTPUT_MODES)
OUTPUT_CHARACTER = _dated(OUTPUT_CHARACTER)

# The three orthogonal directions the bridging atom can move in, in the order
# they are reported. Names follow the silica-glass literature.
CHARACTER_NAMES = ('stretch', 'bend', 'rock')

THZ_PER_CM1 = 0.0299792458   # 1 cm^-1 = 0.0299792458 THz
EV_PER_THZ  = 0.0041356677   # 1 THz = h * 1e12 Hz = 4.1356677e-3 eV (E = h*nu)

# nu[THz] per sqrt(eigenvalue), by dynamical_matrix style, for `units metal`.
#
# 'regular' leaves the mass-weighted force constants in eV/(A^2*amu). One of
# those is 1 eV / (1 A^2 * 1 amu) = 1.602176634e-19 / (1e-20 * 1.66053907e-27)
# = 9.648533e27 s^-2, so omega = 9.822695e13 rad/s and nu = omega/(2*pi) =
# 15.6333042 THz.
#
# 'eskm' means LAMMPS already applied that same conversion itself (its metal
# force->mvv2e = 1.0364269e-4, so 1/mvv2e = 9648.53 ps^-2), leaving lambda in
# 1/ps^2 — from which sqrt is omega in rad/ps and nu is just omega/(2*pi).
#
# The two rows are therefore the same physical conversion against different
# input units, and agree to ~1e-6: LAMMPS's mvv2e is itself a rounded constant,
# so the eskm path inherits that rounding rather than the CODATA value used here.
# Far below any frequency this method resolves, but it is why the two styles do
# not agree to the last digit.
THZ_PER_SQRT_EIGENVALUE = {
    'regular': 15.6333042,
    'eskm':    1.0 / (2.0 * np.pi),
}

# XUNIT value of 1 THz, i.e. multiply a THz frequency by this to get XUNIT.
UNIT_PER_THZ = {
    'THz':  1.0,
    'meV':  EV_PER_THZ * 1000.0,
    'eV':   EV_PER_THZ,
    'cm-1': 1.0 / THZ_PER_CM1,
}

FREQ_UNIT_LABELS = {'meV': 'E (meV)', 'THz': 'ν (THz)', 'cm-1': 'ν (cm⁻¹)', 'eV': 'E (eV)'}

# Only needed for ASR='simple', which has to undo the mass weighting. These are
# the masses OH-therm.input/b-SiO-therm.input set.
ELEMENT_MASSES = {'O': 15.9994, 'H': 1.00784, 'Si': 28.0855}

# Coherent and incoherent neutron scattering cross-sections in barn, from NIST
# (https://www.ncnr.nist.gov/resources/n-lengths/). Used only by
# VDOS_DYNMAT_WEIGHTING; see WEIGHTING in the module docstring.
#
# The mass in the sigma/m weight deliberately comes from ELEMENT_MASSES above
# rather than from a second table here: that is the mass the dynamical matrix
# was actually mass-weighted with, so the weighting cannot drift away from the
# physics the eigenvectors describe. Adding an element means adding it to both.
NEUTRON_CROSS_SECTIONS = {
    'H':  (1.7568,  80.26),      # protium
    'D':  (5.592,    2.05),
    'C':  (5.551,    0.001),
    'N':  (11.01,    0.5),
    'O':  (4.232,    0.0008),
    'Na': (1.66,     1.62),
    'Mg': (3.631,    0.08),
    'Al': (1.495,    0.0082),
    'Si': (2.163,    0.004),
    'P':  (3.307,    0.005),
    'S':  (1.0186,   0.007),
    'Cl': (11.5257,  5.3),
    'K':  (1.69,     0.27),
    'Ca': (2.78,     0.05),
    'Fe': (11.22,    0.4),
    'Ni': (13.3,     5.2),
    'Zr': (6.44,     0.02),
    'Ba': (3.23,     0.15),
}

WEIGHTING_EQUATIONS = {
    'unity':      'k_el = 1',
    'coherent':   'k_el = sigma_coh_el / m_el',
    'incoherent': 'k_el = sigma_inc_el / m_el',
    'total':      'k_el = (sigma_coh_el + sigma_inc_el) / m_el',
}


def parse_weightings(value):
    """
    Parse the semicolon-separated VDOS_DYNMAT_WEIGHTING list.

    Semicolon and not comma because submit_pipeline.sh passes settings through
    `sbatch --export`, which is comma-delimited and silently truncates a value at
    the first embedded comma.
    """
    if ',' in value:
        raise ValueError(
            f"VDOS_DYNMAT_WEIGHTING={value!r} uses ',' but the separator is ';' — a comma "
            f"would be truncated by `sbatch --export` in submit_pipeline.sh. Write it as "
            f"{value.replace(',', ';')!r}."
        )
    keys = [k.strip() for k in value.split(';') if k.strip()]
    if not keys:
        raise ValueError(f"VDOS_DYNMAT_WEIGHTING is empty; choose from {list(WEIGHTING_EQUATIONS)}.")
    unknown = [k for k in keys if k not in WEIGHTING_EQUATIONS]
    if unknown:
        raise ValueError(
            f"Unknown VDOS_DYNMAT_WEIGHTING {unknown}; choose from {list(WEIGHTING_EQUATIONS)}.")
    return keys


def check_weighting_prerequisites(unique_elements, weightings, partial):
    """
    Refuse the cases where a neutron-weighted DOS would be meaningless here.

    Hydrogen: protium and deuterium differ by ~21x in sigma/m and a LAMMPS dump
    labels both 'H', while under 'incoherent' weighting H carries >99.9% of the
    weight — the result would be set almost entirely by the species whose
    isotope is unknown. Same guard as vdos.py.

    PARTIAL='no': without eigenvectors there are no per-element participations,
    so there is nothing to weight — only the plain total exists.
    """
    neutron = [w for w in weightings if w != 'unity']
    if not neutron:
        return

    if partial != 'yes':
        raise SystemExit(
            f"vdos_dynmat.py: VDOS_DYNMAT_WEIGHTING={neutron} needs per-element mode\n"
            f"  participations, which require eigenvectors. Set VDOS_DYNMAT_PARTIAL=yes,\n"
            f"  or use VDOS_DYNMAT_WEIGHTING=unity."
        )

    present = sorted({'H', 'D'} & set(unique_elements))
    if present:
        raise SystemExit(
            f"\nvdos_dynmat.py: refusing to neutron-weight a hydrogen-bearing system.\n"
            f"  elements present:     {present}\n"
            f"  weightings requested: {neutron}\n\n"
            f"  Protium and deuterium differ by ~21x in sigma/m (81.4 vs 3.8 barn/amu)\n"
            f"  and a LAMMPS dump labels both 'H', so the isotope cannot be determined\n"
            f"  from the trajectory. Under 'incoherent' weighting H would carry >99.9%\n"
            f"  of the weight, so the result would be dominated by exactly the species\n"
            f"  whose treatment is undecided. This needs a careful implementation that\n"
            f"  is deliberately not attempted here.\n\n"
            f"  Use VDOS_DYNMAT_WEIGHTING=unity for this system.\n"
        )


def species_weights(weighting, elements):
    """
    Per-element factors k_el to apply to the partial mode participations.

    NOTE the difference from vdos.py, which is easy to get wrong: vdos.py's
    partials are per-element *shapes* and its weight carries the concentration,
    w_el = c_el*sigma/m. Here the partials already carry it — partial_weights
    sums |v|^2 over an element's rows, so under 'phonon' each integrates to
    3*c_el — and multiplying by another c_el would count concentration twice.

    So k_el is sigma/m normalized so that sum_el k_el*c_el = 1, which keeps the
    total integrating to 3 per atom. The resulting share of the total carried by
    each element is k_el*c_el = c_el(sigma/m)_el / sum(c*sigma/m), the same
    quantity vdos.py prints — which is what makes the two directly comparable.

    Returns (k, shares) with shares[el] = k_el * c_el summing to 1.
    """
    counts = {el: int((elements == el).sum()) for el in sorted(set(elements.tolist()))}
    n_total = sum(counts.values())
    concentration = {el: n / n_total for el, n in counts.items()}

    if weighting == 'unity':
        k = {el: 1.0 for el in counts}
        return k, dict(concentration)

    no_sigma = [el for el in counts if el not in NEUTRON_CROSS_SECTIONS]
    no_mass = [el for el in counts if el not in ELEMENT_MASSES]
    if no_sigma or no_mass:
        raise SystemExit(
            f"vdos_dynmat.py: VDOS_DYNMAT_WEIGHTING={weighting!r} needs both a cross-section "
            f"and a mass for every element.\n"
            f"  missing from NEUTRON_CROSS_SECTIONS: {no_sigma or 'none'}\n"
            f"  missing from ELEMENT_MASSES:         {no_mass or 'none'}\n"
            f"  Add them at the top of this file (masses matching the LAMMPS input's "
            f"`mass` lines), or use VDOS_DYNMAT_WEIGHTING=unity."
        )

    raw = {}
    for el in counts:
        sigma_coh, sigma_inc = NEUTRON_CROSS_SECTIONS[el]
        sigma = {'coherent': sigma_coh,
                 'incoherent': sigma_inc,
                 'total': sigma_coh + sigma_inc}[weighting]
        raw[el] = sigma / ELEMENT_MASSES[el]

    denominator = sum(concentration[el] * raw[el] for el in counts)
    if denominator <= 0:
        raise SystemExit(
            f"vdos_dynmat.py: VDOS_DYNMAT_WEIGHTING={weighting!r} gives zero total weight — "
            f"every cross-section involved is zero."
        )
    k = {el: raw[el] / denominator for el in counts}
    return k, {el: k[el] * concentration[el] for el in counts}


def print_weighting_table(shares_by_key, normalization):
    """Print each total's defining equation and the share of it each element carries."""
    sum_rule = ('integral = 3 per atom' if normalization == 'phonon' else 'integral = 1')
    print("\nWeighted totals (check these shares against the curves):")
    for key, shares in shares_by_key.items():
        listed = ', '.join(f'{el}={s:.4f}' for el, s in sorted(shares.items()))
        print(f"  DoS(Total_{key})")
        print(f"    {WEIGHTING_EQUATIONS[key]}, normalized so sum_el k_el*c_el = 1")
        print(f"    share of total: {listed}   sum={sum(shares.values()):.6f}   {sum_rule}")
    print("  Debye-Waller factor exp(-2W) is NOT applied: it is Q-dependent, while this DOS\n"
          "  is not Q-resolved. Shares are normalized, so they say nothing about absolute\n"
          "  signal: for light elements sigma_inc is tiny and a real measurement is\n"
          "  coherent-dominated.")


def read_elements(filename):
    """
    Read only the element column of the dump's FIRST frame and stop — this file
    is the trajectory, and everything after frame 0 is irrelevant here.

    Returns a (n_atoms,) array of element labels in atom-ID order, which is the
    order dynamical_matrix uses for the matrix rows. Requires
    `dump_modify ... sort id`, already used throughout this pipeline.
    """
    with open(filename) as f:
        line = f.readline()
        if not line:
            raise ValueError(f"{filename} is empty.")

        # TIMESTEP value
        f.readline()

        # NUMBER OF ATOMS
        f.readline()
        n_atoms = int(f.readline().strip())

        # BOX BOUNDS (header + 3 dims) — not needed here
        for _ in range(4):
            f.readline()

        # ATOMS header — parse column positions dynamically
        header = f.readline().split()   # ['ITEM:', 'ATOMS', 'id', 'element', ...]
        cols = header[2:]
        if 'element' not in cols:
            raise ValueError(
                f"{filename} has no 'element' column — vdos_dynmat.py takes the "
                f"per-atom element from this dump (dump ... id element x y z ...) "
                f"so the LAMMPS stage does not have to write a second reference file."
            )
        col_element = cols.index('element')
        col_id = cols.index('id') if 'id' in cols else None

        rows = [f.readline().split() for _ in range(n_atoms)]
        elements = np.array([r[col_element] for r in rows])

        # The row-to-atom mapping this whole script rests on assumes the dump is
        # in ascending atom-ID order, which is what dynamical_matrix iterates in.
        # Cheap to check, and silent garbage if it is ever wrong.
        if col_id is not None:
            ids = np.array([int(r[col_id]) for r in rows])
            if not np.array_equal(ids, np.arange(1, n_atoms + 1)):
                raise ValueError(
                    f"{filename} frame 0 is not in ascending atom-ID order 1..{n_atoms}. "
                    f"dynamical_matrix orders its rows by global atom ID, so the element "
                    f"labels would be assigned to the wrong matrix rows.\n"
                    f"Add 'dump_modify <id> sort id' to the LAMMPS input."
                )

    return elements


def read_frame_with_positions(filename):
    """
    Read elements, positions and box from the first (only) frame of the
    reference dump written after `minimize`.

    Returns (elements, positions, box_lengths, tilt) with positions in the
    dump's own order, which `read_elements`' ID check has already established is
    ascending atom ID — the matrix row order.
    """
    with open(filename) as f:
        if not f.readline():
            raise ValueError(f"{filename} is empty.")
        f.readline()                                   # timestep
        f.readline()
        n_atoms = int(f.readline().strip())

        bounds_header = f.readline().split()           # ITEM: BOX BOUNDS [xy xz yz] pp pp pp
        lo, hi, tilt = np.empty(3), np.empty(3), np.zeros(3)
        for axis in range(3):
            parts = f.readline().split()
            lo[axis], hi[axis] = float(parts[0]), float(parts[1])
            if len(parts) > 2:                         # triclinic: xy, xz, yz
                tilt[axis] = float(parts[2])

        header = f.readline().split()
        cols = header[2:]
        missing = [c for c in ('element', 'x', 'y', 'z') if c not in cols]
        if missing:
            raise ValueError(
                f"{filename} is missing column(s) {missing}. The mode-character "
                f"analysis needs 'id element x y z'; see the write_dump line in the "
                f"RUN_DYNMAT block of OH-therm.input / b-SiO-therm.input."
            )
        ce, cx, cy, cz = (cols.index(c) for c in ('element', 'x', 'y', 'z'))

        rows = [f.readline().split() for _ in range(n_atoms)]

    elements  = np.array([r[ce] for r in rows])
    positions = np.array([[float(r[cx]), float(r[cy]), float(r[cz])] for r in rows])

    # LAMMPS triclinic bounds are the bounding box of the tilted cell, not the
    # cell lengths; correcting for that is only worth doing if it comes up.
    if np.any(tilt != 0.0):
        raise ValueError(
            f"{filename} describes a triclinic box (tilt factors {tilt}). The "
            f"mode-character analysis assumes an orthogonal cell for its minimum-"
            f"image bond vectors, so it would silently mis-assign bonds across the "
            f"periodic boundary. Set VDOS_DYNMAT_CHARACTER=no, or extend "
            f"find_bridges() to pass the tilt through to freud."
        )

    return elements, positions, hi - lo


def find_bridges(elements, positions, box_lengths, bridge_element,
                 neighbor_element, cutoff):
    """
    Find every bridge_element atom bonded to exactly two neighbor_element atoms,
    and build the local orthonormal frame at each.

    With r1, r2 the unit vectors from the bridge to its two neighbours, the three
    directions are

        stretch = norm(r1 - r2)   along the neighbour-neighbour line; one bond
                                  lengthens as the other shortens
        bend    = norm(r1 + r2)   along the bisector; opens and closes the
                                  neighbour-bridge-neighbour angle
        rock    = norm(r1 x r2)   perpendicular to the plane

    r1-r2 and r1+r2 are orthogonal whenever |r1| = |r2|, and the cross product is
    perpendicular to both, so the three form a complete orthonormal basis for the
    bridge atom's motion — this is a decomposition, not a heuristic split.

    Returns (bridge_indices, frames, census) where frames is
    (n_bridge, 3, 3) indexed [atom, direction, xyz] in CHARACTER_NAMES order.
    """
    import freud   # only needed for CHARACTER='yes'; keeps the base script numpy-only

    box = freud.box.Box(Lx=box_lengths[0], Ly=box_lengths[1], Lz=box_lengths[2])
    is_bridge   = elements == bridge_element
    is_neighbor = elements == neighbor_element
    if not is_bridge.any() or not is_neighbor.any():
        raise ValueError(
            f"No {bridge_element} or no {neighbor_element} atoms in the reference "
            f"dump (elements present: {sorted(set(elements.tolist()))}). Set "
            f"VDOS_DYNMAT_BRIDGE_ELEMENT / VDOS_DYNMAT_NEIGHBOR_ELEMENT for this system."
        )

    bridge_idx   = np.flatnonzero(is_bridge)
    neighbor_idx = np.flatnonzero(is_neighbor)
    wrapped = box.wrap(positions)

    # Query direction, and the use of nl.vectors, both follow bad_freud.py's
    # _query(): building over the neighbours and probing with the bridge atoms
    # makes freud return the bonds already segmented by bridge atom, and its
    # vectors are minimum-image corrected so no manual wrapping is needed.
    nl = freud.locality.AABBQuery(box, wrapped[neighbor_idx]).query(
        wrapped[bridge_idx], {'r_max': cutoff, 'r_min': 1e-6, 'exclude_ii': False}
    ).toNeighborList()

    counts   = np.asarray(nl.neighbor_counts).astype(np.intp)
    segments = np.asarray(nl.segments).astype(np.intp)
    vectors  = np.asarray(nl.vectors) / np.asarray(nl.distances)[:, None]

    census = {
        'bridge_total':   int(len(bridge_idx)),
        'coordination':   {int(c): int(n) for c, n in
                           zip(*np.unique(counts, return_counts=True))},
        'neighbor_total': int(len(neighbor_idx)),
    }

    keep = counts == 2
    kept_bridge = bridge_idx[keep]
    if not keep.any():
        raise ValueError(
            f"No {bridge_element} atom has exactly two {neighbor_element} neighbours "
            f"within {cutoff} A, so no local frame can be built. Check "
            f"VDOS_DYNMAT_BOND_CUTOFF against the first peak of the {bridge_element}-"
            f"{neighbor_element} RDF."
        )

    # The two bond directions of each two-coordinated bridge atom sit at
    # segments[i] and segments[i]+1 in the flat bond list.
    first = segments[keep]
    r1, r2 = vectors[first], vectors[first + 1]

    def _unit(v):
        return v / np.linalg.norm(v, axis=1, keepdims=True)

    angles = np.degrees(np.arccos(np.clip((r1 * r2).sum(axis=1), -1.0, 1.0)))

    # Built as stretch -> rock -> bend rather than the three literal expressions,
    # because rock x stretch is proportional to (r1 + r2) exactly (expand the
    # triple product) while staying orthonormal by construction instead of by
    # coincidence of three separate normalizations.
    stretch = _unit(r1 - r2)
    cross = np.cross(r1, r2)
    cross_norm = np.linalg.norm(cross, axis=1)          # = |sin(angle)|

    # A linear bridge has no plane: at 180 degrees both r1 + r2 and r1 x r2
    # vanish, and bend and rock become physically degenerate — any two
    # perpendicular directions across the axis are equivalent by symmetry. Ideal
    # beta-cristobalite is exactly this case, so it has to be handled rather than
    # left to divide by zero.
    #
    # The threshold is about the precision of the input coordinates, not machine
    # epsilon. The reference dump stores positions to ~6 decimals, so a bond
    # direction carries ~1e-6 of error; the cross product's *direction* is then
    # only meaningful while |sin(angle)| stays well above that. At 1e-3 (within
    # 0.06 degrees of linear) the direction is good to ~0.1%, and below it the
    # cross is noise and the arbitrary perpendicular is strictly better.
    degenerate = cross_norm < 1e-3
    if degenerate.any():
        fallback = np.zeros((int(degenerate.sum()), 3))
        axis = np.argmin(np.abs(stretch[degenerate]), axis=1)     # least-aligned axis
        fallback[np.arange(len(fallback)), axis] = 1.0
        cross[degenerate] = np.cross(stretch[degenerate], fallback)

    # r1 x r2 is perpendicular to r1 - r2 analytically, but that cancellation is
    # ill-conditioned exactly where the cross is small — normalizing then amplifies
    # the residual. Projecting the component along stretch back out costs nothing
    # and makes the frame orthonormal to machine precision at every angle.
    cross -= (cross * stretch).sum(axis=1, keepdims=True) * stretch
    rock = _unit(cross)
    bend = np.cross(rock, stretch)                       # unit and orthogonal by construction

    frames = np.stack([stretch, bend, rock], axis=1)

    census['angles_deg'] = angles
    census['degenerate'] = int(degenerate.sum())
    # Not singular but ill-conditioned: near 180 degrees the plane is poorly
    # defined, so the bend/rock split is noisy even though their sum is not.
    census['near_linear'] = int((angles > 175.0).sum())
    return kept_bridge, frames, census


def mode_character(eigenvectors, bridge_indices, frames, masses, chunk=256):
    """
    Fraction of each mode's bridge-atom motion along stretch, bend and rock.

    The eigenvectors of a mass-weighted matrix are not displacements: with
    D = Phi/sqrt(m_i m_j), the physical displacement is u_i = e_i/sqrt(m_i). The
    projections below are geometric, so they must use u, not e — skipping this
    would systematically misweight the lighter species, which for silica is the
    oxygen that carries all the band character.

    Chunked over bridge atoms so the (n_bridge, 3, n_modes) displacement block
    never has to exist in full.

    Returns {name: (n_modes,) fraction}, each in [0, 1] and summing to 1 across
    the three names (up to floating point).
    """
    n_modes = eigenvectors.shape[1]
    totals = {name: np.zeros(n_modes) for name in CHARACTER_NAMES}
    denominator = np.zeros(n_modes)

    inv_sqrt_m = 1.0 / np.sqrt(masses[bridge_indices])        # (n_bridge,)

    for start in range(0, len(bridge_indices), chunk):
        stop = min(start + chunk, len(bridge_indices))
        rows = (3 * bridge_indices[start:stop, None] + np.arange(3)).ravel()
        u = eigenvectors[rows].reshape(stop - start, 3, n_modes)
        u = u * inv_sqrt_m[start:stop, None, None]            # e -> displacement

        denominator += np.einsum('oam,oam->m', u, u)
        for axis, name in enumerate(CHARACTER_NAMES):
            projection = np.einsum('oam,oa->om', u, frames[start:stop, axis, :])
            totals[name] += np.einsum('om,om->m', projection, projection)

    safe = np.where(denominator > 0, denominator, 1.0)
    return {name: value / safe for name, value in totals.items()}


def participation_ratio(eigenvectors, n_atoms, chunk=512):
    """
    PR(m) = 1 / (N * sum_i |e_i|^4), with |e_i|^2 the atom's share of mode m.

    Ranges from 1/N (all the motion on one atom) to 1 (every atom moving
    equally), so it separates extended modes from localized ones — the standard
    way to pick out the boson-peak region and the localized high-frequency modes
    in a glass.

    Defined on the mass-weighted eigenvector, which is the usual convention: e is
    orthonormal, so |e_i|^2 is atom i's fraction of the mode's kinetic energy and
    the sum over atoms is exactly 1.
    """
    n_modes = eigenvectors.shape[1]
    quartic = np.zeros(n_modes)
    for start in range(0, n_atoms, chunk):
        stop = min(start + chunk, n_atoms)
        block = eigenvectors[3 * start:3 * stop].reshape(stop - start, 3, n_modes)
        share = np.einsum('oam,oam->om', block, block)        # |e_i|^2
        quartic += np.einsum('om,om->m', share, share)
    return 1.0 / (n_atoms * quartic)


def reduced_dos(total_curve, grid):
    """
    g(nu)/nu^2 — the conventional way to display the boson peak, the excess over
    the Debye prediction that g(nu) ~ nu^2 flattens into a constant. The nu = 0
    bin is NaN rather than infinity so it simply does not plot.
    """
    with np.errstate(divide='ignore', invalid='ignore'):
        reduced = np.where(grid > 0, total_curve / np.square(grid), np.nan)
    return reduced


def read_matrix(filename, n_atoms, binary):
    """
    Read the dynamical matrix as a flat stream and reshape to (3N, 3N).

    LAMMPS prints three numbers per line, ordered j fastest, then alpha, then i,
    so the flat stream maps straight onto row 3i+alpha, column 3j+beta with no
    reordering. np.fromfile is used rather than np.loadtxt because this file
    reaches gigabytes on realistic cells.

    The length check is the load-bearing assertion of this script: it is what
    fails loudly if the matrix and the dump describe different systems.
    """
    if binary == 'yes':
        values = np.fromfile(filename, dtype=np.float64)
    else:
        values = np.fromfile(filename, sep=' ')

    expected = 9 * n_atoms * n_atoms
    if values.size != expected:
        raise ValueError(
            f"{filename} holds {values.size} values but {n_atoms} atoms in "
            f"{DUMP_FILE} require {expected} (= 9*N^2).\n"
            f"Likely causes: the matrix and the dump come from different runs; "
            f"dynamical_matrix was given a group other than 'all'; or "
            f"DYNMAT_BINARY={binary!r} does not match the `binary` keyword the "
            f"LAMMPS run used."
        )

    return values.reshape(3 * n_atoms, 3 * n_atoms)


def apply_asr(matrix, elements):
    """
    Impose the acoustic sum rule in place: sum_j Phi_ij = 0 for every atom i.

    The matrix LAMMPS writes is mass-weighted (D = Phi/sqrt(m_i*m_j)), so this
    un-weights it, corrects each atom's diagonal 3x3 block by its own row sum,
    and re-weights. The mass scaling is done with broadcasting rather than an
    outer product to avoid allocating a second full (3N, 3N) array.
    """
    unknown = sorted(set(elements.tolist()) - set(ELEMENT_MASSES))
    if unknown:
        raise ValueError(
            f"VDOS_DYNMAT_ASR='simple' needs the mass of every element to undo the "
            f"mass weighting, but {unknown} are not in ELEMENT_MASSES. Add them at "
            f"the top of this file (use the same values as the LAMMPS input's "
            f"`mass` lines), or set VDOS_DYNMAT_ASR=none."
        )

    n_atoms = len(elements)
    sqrt_m = np.repeat(np.sqrt([ELEMENT_MASSES[el] for el in elements]), 3)   # (3N,)

    matrix *= sqrt_m[:, None]
    matrix *= sqrt_m[None, :]                                                 # now Phi

    blocks = matrix.reshape(n_atoms, 3, n_atoms, 3)
    row_sums = blocks.sum(axis=2)                                             # (N, 3, 3)
    idx = np.arange(n_atoms)
    blocks[idx, :, idx, :] -= row_sums

    matrix /= sqrt_m[:, None]
    matrix /= sqrt_m[None, :]                                                 # back to D
    return matrix


def eigenvalues_to_frequencies(eigenvalues, matrix_style, xunit):
    """
    Signed frequencies in `xunit`: nu = sign(lambda) * sqrt(|lambda|) * factor.

    Imaginary modes (lambda < 0) come back negative, which is the usual plotting
    convention and keeps them visible in the diagnostics rather than turning into
    NaN or silently folding onto the real axis.
    """
    freq_THz = np.sign(eigenvalues) * np.sqrt(np.abs(eigenvalues)) * THZ_PER_SQRT_EIGENVALUE[matrix_style]
    return freq_THz * UNIT_PER_THZ[xunit]


def partial_weights(eigenvectors, elements):
    """
    Per-mode participation of each element: w[e, m] = sum over that element's
    rows of |v[row, m]|^2. The eigenvectors of a mass-weighted matrix are
    orthonormal, so these sum to 1 across elements for every mode.

    Computed one element at a time so the squared array is never materialized
    for the whole matrix at once.
    """
    weights = {}
    for el in sorted(set(elements.tolist())):
        rows = np.repeat(elements == el, 3)
        weights[el] = np.einsum('rm,rm->m', eigenvectors[rows], eigenvectors[rows])
    return weights


def bin_spectrum(freqs, weights, max_frequency, bins, smearing):
    """
    Turn a discrete set of mode frequencies into a curve on a common grid.

    smearing == 0 gives a plain histogram; smearing > 0 sums a Gaussian of that
    FWHM per mode straight onto the bin centers instead. Doing it that way rather
    than histogramming and then convolving avoids the double discretization —
    a mode's contribution is placed at its actual frequency, not at its bin's.

    Only modes in [0, max_frequency] contribute, so imaginary modes cannot leak
    into the real spectrum. Returns (grid, {label: curve}), both unnormalized.
    """
    edges = np.linspace(0.0, max_frequency, bins + 1)
    grid = 0.5 * (edges[:-1] + edges[1:])
    in_range = (freqs >= 0.0) & (freqs <= max_frequency)

    if smearing == 0:
        curves = {label: np.histogram(freqs[in_range], bins=edges, weights=w[in_range])[0]
                  for label, w in weights.items()}
        # Histogram bins hold counts; divide by the bin width to make them a
        # density, so both branches of this function return the same units.
        width = edges[1] - edges[0]
        return grid, {label: c / width for label, c in curves.items()}

    sigma = smearing / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    nu = freqs[in_range]
    # (bins, n_modes) kernel, each column a Gaussian centered on one mode. Every
    # curve is then one matrix-vector product against it.
    kernel = np.exp(-0.5 * ((grid[:, None] - nu[None, :]) / sigma) ** 2)
    # Normalize each column over the grid rather than analytically. A mode within
    # a few sigma of an edge — every acoustic mode, sitting at 0 — has part of its
    # Gaussian off the plotted range, and an analytic 1/(sigma*sqrt(2pi)) would
    # silently drop that part: the three zero modes alone would cost the total
    # 1.5 states. Per-column normalization instead says each mode contributes
    # exactly one state spread over the visible range, so the sum rule holds
    # whether or not smearing is on.
    column_area = kernel.sum(axis=0) * (edges[1] - edges[0])
    kernel /= np.where(column_area > 0, column_area, 1.0)
    return grid, {label: kernel @ w[in_range] for label, w in weights.items()}


def normalize(curves, grid, n_atoms, normalization):
    """
    'phonon': divide every curve by the atom count, so the total integrates to
       exactly 3 per atom and each partial to 3*N_el/N_total — the partials sum
       to the total by construction.
    'unit_area': rescale each curve independently to unit integral, matching
       vdos.py's option of the same name.
    """
    if normalization == 'phonon':
        return {label: c / n_atoms for label, c in curves.items()}

    d_nu = grid[1] - grid[0]

    def unit_area(curve):
        area = curve.sum() * d_nu
        return curve / area if area > 0 else curve

    return {label: unit_area(c) for label, c in curves.items()}


def _freq_all_units(grid, xunit):
    """Expand a frequency grid given in `xunit` into all four supported units."""
    freq_THz = grid / UNIT_PER_THZ[xunit]
    return {unit: freq_THz * factor for unit, factor in UNIT_PER_THZ.items()}


def save_csv(results, freq_by_unit, filename, extra=None):
    """Save results dict {element_or_'total': array} to CSV. The first columns are
    identical to vdos.py's so the two methods' outputs can be overlaid; `extra`
    appends further named columns (character-resolved DOS, reduced DOS) after."""
    order = ([el for el in results if not el.startswith('total_')]
             + [el for el in results if el.startswith('total_')])
    header_parts = ['freq_meV', 'freq_THz', 'freq_cm-1', 'freq_eV']
    columns = [freq_by_unit['meV'], freq_by_unit['THz'], freq_by_unit['cm-1'], freq_by_unit['eV']]
    for label in order:
        header_parts.append(f'DoS(Total_{label[len("total_"):]})'
                            if label.startswith('total_') else f'DoS({label})')
        columns.append(results[label])
    for label, values in (extra or {}).items():
        header_parts.append(label)
        columns.append(values)
    header = ','.join(header_parts)
    data = np.column_stack(columns)
    np.savetxt(filename, data, delimiter=',', header=header, comments='', fmt='%.6f')
    print(f"Data table saved to {filename}")


def save_modes_csv(freqs, freq_by_unit_of, pr, character, element_fracs, filename):
    """
    One row per mode: its frequency in all four units, its participation ratio,
    and its character and element fractions. This is what makes 'which modes are
    in this peak' answerable — the DOS curves only show which character dominates
    a region, not which individual modes put it there.
    """
    header_parts = ['freq_meV', 'freq_THz', 'freq_cm-1', 'freq_eV', 'participation_ratio']
    columns = [freq_by_unit_of['meV'], freq_by_unit_of['THz'],
               freq_by_unit_of['cm-1'], freq_by_unit_of['eV'], pr]
    for name in CHARACTER_NAMES:
        header_parts.append(f'frac_{name}')
        columns.append(character[name])
    for element, values in sorted(element_fracs.items()):
        header_parts.append(f'frac_{element}')
        columns.append(values)
    np.savetxt(filename, np.column_stack(columns), delimiter=',',
               header=','.join(header_parts), comments='', fmt='%.6f')
    print(f"Per-mode table saved to {filename}  ({len(freqs)} modes)")


def plot_character(grid, character_curves, total, pr, mode_freqs, filename, xunit=XUNIT):
    """
    Three stacked panels: what kind of motion each band is, how localized its
    modes are, and the reduced DOS that exposes the boson peak.
    """
    fig, axes = plt.subplots(3, 1, figsize=(8, 11), sharex=True)

    ax = axes[0]
    ax.plot(grid, total, color='0.3', linewidth=1.6, label='total')
    for name in CHARACTER_NAMES:
        ax.plot(grid, character_curves[name], linewidth=1.2, label=name)
    ax.set_ylabel('DOS (states / atom / ' + xunit + ')')
    ax.set_title('Character-resolved DOS (bridging-atom motion)')
    ax.legend()

    ax = axes[1]
    # One point per mode: scatter rather than a curve, because the spread of PR
    # at a given frequency is itself the information — a tight low band means
    # every mode there is equally extended, a wide one means they are not.
    ax.plot(mode_freqs, pr, '.', markersize=2, alpha=0.4)
    ax.set_ylabel('participation ratio')
    ax.set_title('Localization (1 = every atom moves, 1/N = one atom moves)')
    ax.set_ylim(0, 1)

    ax = axes[2]
    ax.plot(grid, reduced_dos(total, grid), color='C3', linewidth=1.2)
    ax.set_ylabel(f'g / {xunit}²')
    ax.set_xlabel(FREQ_UNIT_LABELS[xunit])
    ax.set_title('Reduced DOS g(ν)/ν² — a peak here is the boson peak')

    fig.tight_layout()
    fig.savefig(filename, dpi=PLOT_DPI)
    plt.close(fig)
    print(f"Character plot saved to {filename}")


def report_census(census, bridge_element, neighbor_element, cutoff):
    """
    Report the coordination of every bridge atom. In a quenched glass some
    oxygens are non-bridging and some silicons are mis-coordinated; the local
    frame only exists for the two-coordinated ones, so this says how much of the
    structure the decomposition actually covers — and doubles as a glass-quality
    check worth reading on its own.
    """
    total = census['bridge_total']
    print(f"  Coordination census ({bridge_element} by {neighbor_element} neighbours "
          f"within {cutoff} A):")
    for coordination, count in sorted(census['coordination'].items()):
        kind = {0: 'isolated', 1: 'non-bridging', 2: 'bridging'}.get(coordination, 'over-coordinated')
        print(f"    {coordination}-coordinated: {count:6d}  ({100.0*count/total:5.1f}%)  {kind}")
    bridging = census['coordination'].get(2, 0)
    print(f"  Local frames built on {bridging}/{total} {bridge_element} atoms "
          f"({100.0*bridging/total:.1f}% of them); the character fractions describe "
          f"only their motion.")
    if bridging < total:
        print(f"  NOTE: {total - bridging} {bridge_element} atom(s) are not two-coordinated "
              f"and contribute nothing to the character fractions. If that count is large, "
              f"check VDOS_DYNMAT_BOND_CUTOFF against the {bridge_element}-{neighbor_element} "
              f"RDF before reading anything into the decomposition.")
    angles = census['angles_deg']
    if len(angles):
        print(f"  {bridge_element}-{neighbor_element} bridge angle: mean {angles.mean():.1f}"
              f" +/- {angles.std():.1f} deg  (a-SiO2 sits near 144)")
    near_linear = census.get('near_linear', 0)
    exactly_linear = census.get('degenerate', 0)
    if near_linear:
        detail = f" ({exactly_linear} of them exactly linear)" if exactly_linear else ""
        print(f"  NOTE: {near_linear} bridge(s) within 5 deg of linear{detail}. A linear "
              f"bridge has no plane,\n"
              f"        so bend and rock are physically degenerate there and their split is "
              f"arbitrary — only\n"
              f"        their SUM (the transverse motion) is meaningful for those. Stretch is "
              f"unaffected.")


def plot_vdos(results, freq_by_unit, filename, xunit=XUNIT):
    x = freq_by_unit[xunit]

    fig, ax = plt.subplots(figsize=(8, 5))
    for label, curve in results.items():
        ax.plot(x, curve, label=label, linewidth=1.5 if label.startswith('total_') else 1.0)
    ax.set_xlabel(FREQ_UNIT_LABELS[xunit])
    ax.set_ylabel('DOS (states / atom / ' + xunit + ')' if NORMALIZATION == 'phonon'
                  else 'VDOS (unit-area normalized)')
    ax.legend()
    fig.tight_layout()
    fig.savefig(filename, dpi=PLOT_DPI)
    plt.close(fig)
    print(f"Plot saved to {filename}")


def snap_zero_modes(freqs, tolerance):
    """
    Snap modes within `tolerance` of zero to exactly zero, in place.

    The three acoustic translations are zero by symmetry, but finite-difference
    noise scatters them across zero — so without this, one of them lands at, say,
    -0.005 meV, falls outside the [0, MAX_FREQUENCY] range, and quietly costs the
    spectrum 1/N of its total weight. Snapping puts them all in the first bin, so
    the 3-modes-per-atom sum rule comes out exact and a negative frequency that
    survives is a real imaginary mode rather than rounding.

    Tolerance is one bin width: below the plot's own resolution, so this can
    never move a mode the spectrum could have distinguished.
    """
    near_zero = np.abs(freqs) < tolerance
    freqs[near_zero] = 0.0
    return int(near_zero.sum())


def report_diagnostics(freqs, xunit, max_frequency, n_snapped, tolerance):
    """
    Three zero modes are the acoustic translations and are expected. Any other
    imaginary mode means the minimization did not reach a true minimum, so print
    enough to tell those two cases apart instead of quietly dropping them.
    """
    lowest = np.sort(freqs)[:6]
    print("  Lowest 6 signed frequencies ({}): {}".format(
        xunit, ", ".join(f"{v:.4f}" for v in lowest)))
    print(f"  Modes snapped to zero (|nu| < {tolerance:.4g} {xunit}, one bin width): {n_snapped}"
          f"  (3 acoustic translations are expected)")
    if n_snapped != 3:
        print(f"  NOTE: {n_snapped} zero modes rather than 3. More can mean a disconnected "
              f"fragment or a floppy mode; fewer can mean the bin width is too fine to "
              f"absorb the acoustic modes' numerical scatter.")

    n_negative = int((freqs < 0).sum())
    if n_negative:
        print(f"  WARNING: {n_negative} imaginary mode(s), beyond the acoustic ones — the "
              f"minimization did not reach a true minimum, and they are excluded from the "
              f"spectrum.\n"
              f"           Most negative: {freqs.min():.4f} {xunit}. Tighten DYNMAT_MIN_FTOL "
              f"or reduce DYNMAT_DISPLACEMENT in the LAMMPS stage before trusting this.")
    else:
        print("  Imaginary modes: 0")

    n_above = int((freqs > max_frequency).sum())
    if n_above:
        print(f"  WARNING: {n_above} mode(s) lie above VDOS_DYNMAT_MAX_FREQUENCY="
              f"{max_frequency} {xunit} and are excluded — the spectrum is truncated. "
              f"Highest mode: {freqs.max():.4f} {xunit}.")


if __name__ == '__main__':
    import time

    # Validated before anything is read, so a typo costs a second rather than a
    # full diagonalization.
    weightings = parse_weightings(WEIGHTING)

    t0 = time.time()
    positions = box_lengths = None
    if CHARACTER == 'yes':
        # The character analysis is a projection onto bond directions, so it must
        # use the geometry the force constants were built at — the post-minimize
        # dump, not the last MD frame.
        print(f"Reading minimized geometry from: {DYNMAT_REF_TRAJ}")
        try:
            elements, positions, box_lengths = read_frame_with_positions(DYNMAT_REF_TRAJ)
        except FileNotFoundError:
            raise SystemExit(
                f"vdos_dynmat.py: VDOS_DYNMAT_CHARACTER=yes needs {DYNMAT_REF_TRAJ}, the "
                f"coordinates written just after `minimize` in the LAMMPS stage.\n"
                f"It is produced by the write_dump line in the RUN_DYNMAT block of "
                f"OH-therm.input / b-SiO-therm.input, so a matrix from an older run "
                f"predating that line will not have one.\n"
                f"Re-run the LAMMPS stage, or set VDOS_DYNMAT_CHARACTER=no."
            )
        print(f"  Box: {np.round(box_lengths, 4)} A")
    else:
        print(f"Reading elements from: {DUMP_FILE}")
        elements = read_elements(DUMP_FILE)

    n_atoms = len(elements)
    unique_els = sorted(set(elements.tolist()))
    print(f"  Atoms: {n_atoms}, elements: {unique_els}")

    # After element detection, so the guard can name what it actually found.
    print(f"  VDOS_DYNMAT_WEIGHTING={weightings}")
    check_weighting_prerequisites(unique_els, weightings, PARTIAL)

    print(f"Reading dynamical matrix: {DYNMAT_FILE} (binary={BINARY}, style={MATRIX_STYLE})")
    matrix = read_matrix(DYNMAT_FILE, n_atoms, BINARY)
    t1 = time.time()
    print(f"  Matrix: {matrix.shape[0]}x{matrix.shape[1]} ({matrix.nbytes / 1e9:.2f} GB, {t1 - t0:.2f}s)")

    scale = float(np.abs(matrix).max())
    asymmetry = float(np.abs(matrix - matrix.T).max())
    print(f"  Asymmetry max|D - D.T| / max|D|: {asymmetry / scale:.3e}" if scale > 0 else
          "  Matrix is all zeros — check the LAMMPS run.")
    matrix = 0.5 * (matrix + matrix.T)

    if ASR == 'simple':
        print("  Applying acoustic sum rule (ASR='simple')")
        matrix = apply_asr(matrix, elements)

    # Character and participation ratio are both defined on the eigenvectors, so
    # either of them forces the full eigh regardless of VDOS_DYNMAT_PARTIAL.
    need_vectors = PARTIAL == 'yes' or CHARACTER == 'yes'
    if CHARACTER == 'yes' and PARTIAL == 'no':
        print("  Note: VDOS_DYNMAT_CHARACTER=yes needs eigenvectors, so PARTIAL=no "
              "does not save anything here — computing them anyway.")

    print(f"Diagonalizing ({'eigh, with eigenvectors' if need_vectors else 'eigvalsh, eigenvalues only'})...")
    # An accidentally single-threaded LAPACK turns minutes into hours at this
    # size, and it fails silently — so say out loud how many threads it will use.
    print(f"  BLAS threads: {os.environ.get('OMP_NUM_THREADS', 'unset (BLAS default)')}"
          f"  [set VDOS_DYNMAT_THREADS to override]")
    t2 = time.time()
    if need_vectors:
        eigenvalues, eigenvectors = np.linalg.eigh(matrix)
    else:
        eigenvalues, eigenvectors = np.linalg.eigvalsh(matrix), None
    del matrix
    t3 = time.time()
    print(f"  Diagonalization: {t3 - t2:.2f}s ({len(eigenvalues)} modes)")

    freqs = eigenvalues_to_frequencies(eigenvalues, MATRIX_STYLE, XUNIT)
    zero_tol = MAX_FREQUENCY / BINS
    n_snapped = snap_zero_modes(freqs, zero_tol)
    report_diagnostics(freqs, XUNIT, MAX_FREQUENCY, n_snapped, zero_tol)

    # 'unity' is the plain total: every mode counts once, which is already the
    # concentration-weighted sum of the partials. Any other weighting rescales
    # each element's participation by k_el before binning, so the weight enters
    # per mode rather than by rescaling finished curves.
    weights = {'total_unity': np.ones_like(freqs)}
    element_fracs = {}
    shares_by_key = {'unity': species_weights('unity', elements)[1]}
    if eigenvectors is not None and PARTIAL == 'yes':
        element_fracs = partial_weights(eigenvectors, elements)
        weights.update(element_fracs)

        for weighting in (w for w in weightings if w != 'unity'):
            k, shares = species_weights(weighting, elements)
            shares_by_key[weighting] = shares
            weights[f'total_{weighting}'] = sum(
                k[el] * frac for el, frac in element_fracs.items())

    character = pr = None
    if CHARACTER == 'yes':
        print(f"Mode character: local frames at {BRIDGE_ELEMENT} bridging two "
              f"{NEIGHBOR_ELEMENT} within {BOND_CUTOFF} A")
        bridge_indices, frames, census = find_bridges(
            elements, positions, box_lengths, BRIDGE_ELEMENT, NEIGHBOR_ELEMENT, BOND_CUTOFF)
        report_census(census, BRIDGE_ELEMENT, NEIGHBOR_ELEMENT, BOND_CUTOFF)

        masses = np.array([ELEMENT_MASSES[el] for el in elements])
        character = mode_character(eigenvectors, bridge_indices, frames, masses)
        pr = participation_ratio(eigenvectors, n_atoms)
        print(f"  Participation ratio: min {pr.min():.4f}, median {np.median(pr):.4f}, "
              f"max {pr.max():.4f}")
        # Weighting the DOS by each fraction is what turns a per-mode number into
        # a band assignment: the curves say which motion dominates each region.
        weights.update(character)

    if eigenvectors is not None:
        del eigenvectors

    print(f"Binning ({BINS} bins over [0, {MAX_FREQUENCY}] {XUNIT}"
          f"{f', Gaussian FWHM {SMEARING} {XUNIT}' if SMEARING else ', no smearing'})...")
    grid, curves = bin_spectrum(freqs, weights, MAX_FREQUENCY, BINS, SMEARING)
    results = normalize(curves, grid, n_atoms, NORMALIZATION)
    freq_by_unit = _freq_all_units(grid, XUNIT)

    # Order the CSV/legend as vdos.py does: partials first, total last. Character
    # curves are pulled out here so the leading columns stay byte-comparable with
    # vdos.py's CSV; they come back as extra columns after the DoS block.
    character_curves = {name: results.pop(name) for name in CHARACTER_NAMES if name in results}
    totals = {name: results[name] for name in results if name.startswith('total_')}
    results = {el: results[el] for el in unique_els if el in results} | totals

    extra = {f'DoS({name})': curve for name, curve in character_curves.items()}
    # The reduced DOS and the character plot describe the unweighted spectrum.
    extra['g/nu^2'] = reduced_dos(results['total_unity'], grid)

    print(f"  Total: {time.time() - t0:.2f}s")
    print_weighting_table(shares_by_key, NORMALIZATION)

    if OUTPUT_CSV is not None:
        save_csv(results, freq_by_unit, OUTPUT_CSV, extra=extra)
    if OUTPUT_PLOT is not None:
        plot_vdos(results, freq_by_unit, OUTPUT_PLOT)
    if CHARACTER == 'yes':
        if OUTPUT_MODES is not None:
            save_modes_csv(freqs, _freq_all_units(freqs, XUNIT), pr, character,
                           element_fracs, OUTPUT_MODES)
        if OUTPUT_CHARACTER is not None:
            plot_character(grid, character_curves, results['total_unity'], pr, freqs,
                           OUTPUT_CHARACTER)
