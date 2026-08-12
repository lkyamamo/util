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
       regular : lambda in eV/(A^2*amu), nu[THz] = sqrt(lambda) * 15.633302
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

# Output basename; .csv and .png are appended (set either to None to skip).
_OUTPUT_BASE = _env("VDOS_DYNMAT_OUTPUT", "vdos_dynmat")
OUTPUT_CSV   = f"{_OUTPUT_BASE}.csv"
OUTPUT_PLOT  = f"{_OUTPUT_BASE}.png"

# Plot appearance
PLOT_DPI = 150

# =============================================================================
# END CONFIGURATION
# =============================================================================

# Prepend today's date (YYYYMMDD_) to every output filename.
def _dated(filename):
    return None if filename is None else f"{date.today():%Y%m%d}_{filename}"

OUTPUT_CSV  = _dated(OUTPUT_CSV)
OUTPUT_PLOT = _dated(OUTPUT_PLOT)

THZ_PER_CM1 = 0.0299792458   # 1 cm^-1 = 0.0299792458 THz
EV_PER_THZ  = 0.0041356677   # 1 THz = h * 1e12 Hz = 4.1356677e-3 eV (E = h*nu)

# nu[THz] per sqrt(eigenvalue), by dynamical_matrix style, for `units metal`.
# 'regular' leaves the mass-weighted force constants in eV/(A^2*amu); LAMMPS's
# own eskm conversion factor of 9648.5 turns those into 1/ps^2, and
# sqrt(9648.5)/(2*pi) = 15.633302 — so the two rows below are the same number
# expressed against the two different input units.
THZ_PER_SQRT_EIGENVALUE = {
    'regular': 15.633302,
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


def save_csv(results, freq_by_unit, filename):
    """Save results dict {element_or_'total': array} to CSV. Column layout is
    identical to vdos.py's so the two methods' outputs can be overlaid."""
    order = [el for el in results if el != 'total'] + ['total']
    header_parts = ['freq_meV', 'freq_THz', 'freq_cm-1', 'freq_eV']
    columns = [freq_by_unit['meV'], freq_by_unit['THz'], freq_by_unit['cm-1'], freq_by_unit['eV']]
    for label in order:
        header_parts.append('DoS(Total)' if label == 'total' else f'DoS({label})')
        columns.append(results[label])
    header = ','.join(header_parts)
    data = np.column_stack(columns)
    np.savetxt(filename, data, delimiter=',', header=header, comments='', fmt='%.6f')
    print(f"Data table saved to {filename}")


def plot_vdos(results, freq_by_unit, filename, xunit=XUNIT):
    x = freq_by_unit[xunit]

    fig, ax = plt.subplots(figsize=(8, 5))
    for label, curve in results.items():
        ax.plot(x, curve, label=label, linewidth=1.5 if label == 'total' else 1.0)
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

    print(f"Reading elements from: {DUMP_FILE}")
    t0 = time.time()
    elements = read_elements(DUMP_FILE)
    n_atoms = len(elements)
    unique_els = sorted(set(elements.tolist()))
    print(f"  Atoms: {n_atoms}, elements: {unique_els}")

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

    print(f"Diagonalizing ({'eigh, with eigenvectors' if PARTIAL == 'yes' else 'eigvalsh, eigenvalues only'})...")
    # An accidentally single-threaded LAPACK turns minutes into hours at this
    # size, and it fails silently — so say out loud how many threads it will use.
    print(f"  BLAS threads: {os.environ.get('OMP_NUM_THREADS', 'unset (BLAS default)')}"
          f"  [set VDOS_DYNMAT_THREADS to override]")
    t2 = time.time()
    if PARTIAL == 'yes':
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

    weights = {'total': np.ones_like(freqs)}
    if eigenvectors is not None:
        weights.update(partial_weights(eigenvectors, elements))
        del eigenvectors

    print(f"Binning ({BINS} bins over [0, {MAX_FREQUENCY}] {XUNIT}"
          f"{f', Gaussian FWHM {SMEARING} {XUNIT}' if SMEARING else ', no smearing'})...")
    grid, curves = bin_spectrum(freqs, weights, MAX_FREQUENCY, BINS, SMEARING)
    results = normalize(curves, grid, n_atoms, NORMALIZATION)
    freq_by_unit = _freq_all_units(grid, XUNIT)

    # Order the CSV/legend as vdos.py does: partials first, total last.
    results = {el: results[el] for el in unique_els if el in results} | {'total': results['total']}

    print(f"  Total: {time.time() - t0:.2f}s")

    if OUTPUT_CSV is not None:
        save_csv(results, freq_by_unit, OUTPUT_CSV)
    if OUTPUT_PLOT is not None:
        plot_vdos(results, freq_by_unit, OUTPUT_PLOT)
