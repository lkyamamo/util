"""
rdf_freud.py — Radial distribution function, coordination number, and selectable
               neutron correlation functions, using freud

QUICK START
-----------
1. Set DUMP_FILE (or $TRAJ) to your LAMMPS custom dump trajectory.
2. Set R_MAX to less than half the shortest box dimension.
3. Pick the convention you want with $RDF_NORMALIZATION and $RDF_FUNCTIONS (see below).
4. Run:  python rdf_freud.py

OUTPUT
------
- rdfs.csv  — r (Å), partial g_AB(r), one column per requested convention
                                                              (set OUTPUT_CSV=None to skip)
- nrs.csv   — r (Å), n(r) both directions per pair            (set OUTPUT_NR_CSV=None to skip)
- wright.csv — a single T(r) built to overlay directly on a published
                    neutron correlation function        (only when RDF_WRIGHT=yes)

This script writes NO plots. Figures come from the plotting pipeline
(jobs/pipeline/plotting/plot_pipeline.sh), which reads the CSVs above and
writes one PNG per quantity. Run it in this directory, or let the analysis
runner call it via RUN_PLOTS=1.

CONVENTIONS
-----------
The same physical content gets packaged many ways in the literature, and the
symbols collide: Soper's G(r) is a weighted h-sum, the PDF community's G(r) is
4πrρ[g−1], Keen's G(r) is a third function.  Keen, J. Appl. Cryst. 34, 172
(2001) tabulates the conventions against each other.  So nothing here is named
by a bare letter — every combined column is named <function>_<normalization>,
and the script prints each column's defining equation, weight sum, and
asymptotic limits at startup.  Check a printed limit against the curve; never
trust the symbol.

Both keys take SEMICOLON-separated lists ("FZ;absolute"), not commas — see
parse_keys for why.  Both are RDF_-prefixed so the distribution they configure
is explicit; vdos.py has its own unrelated VDOS_NORMALIZATION.

$RDF_NORMALIZATION — how much each partial contributes to the sum.  It sets the
pair weight w_AB, with f = 2 − δ_AB written explicitly rather than folded in:

    unity     w_AB = f c_A c_B                    Σw = 1          dimensionless
    FZ        w_AB = f c_A c_B b_A b_B / <b>²     Σw = 1          dimensionless
    absolute  w_AB = f c_A c_B b_A b_B / 100      Σw = <b>²/100   barn/sr/atom
    formula   w_AB = n f c_A c_B b_A b_B / 100   Σw = n<b>²/100  barn/sr/formula-unit

  unity     Every element scatters identically — b = 1, hence the name, which
            refers to the scattering lengths and not to Σw (FZ also sums to 1).
            The weights are then just mole-fraction products.  This is not a
            measurable quantity: it is the composition-averaged structure,
            useful as a structural summary and as the b-free baseline the
            neutron-weighted curves depart from.  Not called 'number', which
            would collide with Bhatia-Thornton's number-number correlation, a
            different construction.  This is the column previously called
            'total'.

  FZ        Faber-Ziman: divide by <b>², so Σw = 1 and the result approaches 1
            at large r exactly like a partial g_AB does.  Dimensionless, which
            makes samples of different composition superimposable — at the cost
            of dividing by a quantity that nearly cancels for H-rich samples
            (light water: <b>² = 0.0031 barn) and vanishes for a null mixture.

  absolute  No division at all: the weighted sum in the units a measured
            differential cross-section carries, barn/sr/atom (the /100 converts
            fm² to barn).  The scale is physical rather than conventional, so
            the excluded-volume plateau lands at −Σw = −<b>²/100 and is a
            direct check on the composition.  Well conditioned as <b> → 0, so
            this is the one to use for light or null samples.

  formula   'absolute' re-quoted per FORMULA UNIT rather than per atom, which is
            what neutron diffraction papers on compounds usually do. Wright's
            eq. (10) writes the baseline as T0 = 4*pi*r*rho0*(sum_j b_j)^2 with
            rho0 in units/A^3 and sum_j b_j the TOTAL scattering length of one
            unit; since sum_j b_j = n<b> and rho0 = rho_atom/n, that equals
            n*rho_a*<b>^2 — exactly n times the per-atom result, for any
            composition. Needs RDF_ATOMS_PER_FORMULA_UNIT (SiO2 -> 3). This
            factor is the commonest reason a published curve sits a constant
            factor above a per-atom calculation.

COMPARING AGAINST A MEASURED CURVE
----------------------------------
RDF_RESOLUTION_MODE=lorch with RDF_LORCH_QMAX reproduces the modification
function neutron glass diffraction uses,

    M(Q) = sin(dr Q)/(dr Q)   for Q <= Q_max,  0 above

whose cosine transform is the real-space peak function convolved with the
correlation function. dr defaults to pi/Q_max, the standard Lorch choice that
places M's first zero at the truncation.

RDF_RESOLUTION_MODE=modified_lorch with RDF_MODIFIED_LORCH_DELTA is Soper's
eq. (60): instead of a window in Q space, smear h(r) in r space with a uniform
sphere of radius D,

    L'(r, D) = 3/(4 pi D^3)   |r| <= D,   0 above

applied as a true 3D convolution. It is meant for h, the baseline-subtracted
form: h → 0 at large r, so the part of the integral that falls off the end of a
finite grid contributes nothing, whereas a column carrying a baseline loses that
baseline at the edge. Being a top hat in real space it is everywhere
non-negative, so unlike the Q-space Lorch window it cannot push a correlation
function negative between peaks; the trade is that it corresponds to no actual
measurement's truncation, so it smooths rather than reproduces an instrument.
Its resolution is FWHM = sqrt(2) D, so D is 1.41x SMALLER than a quoted
resolution — divide, do not substitute.

Papers quote the RESULTING resolution (the FWHM of that peak function) rather
than dr, and the two differ by ~1.73x: Q_max = 45.2 gives dr = 0.0695 and
FWHM = 0.120 A. Passing the quoted resolution as dr instead doubles the
broadening and makes the kernel double-humped — the script prints the FWHM so
it can be checked against the paper, and warns when the kernel comes out
double-humped.

A Gaussian matches the width but not the shape: the Lorch kernel has negative
side lobes near -5% of the peak, which appear beside strong peaks.

RDF_WRIGHT=yes bundles the three choices such a comparison needs into one
output, so none of them has to be remembered separately:

    RDF_WRIGHT=yes RDF_WRIGHT_QMAX=45.2 RDF_ATOMS_PER_FORMULA_UNIT=3

writes <date>_wright.csv with T(r) Lorch-broadened, per formula unit,
alongside the unbroadened curve and the T0(r) baseline it oscillates about. It
prints Sum(w), the resolution FWHM and the T0 slope 4*pi*rho*Sum(w) — that last
number is the check worth making, since the slope of the paper's average-density
line must equal it and is fixed by composition and density with nothing fitted.

$RDF_FUNCTIONS — what is built from those weights and the partials g_AB:

    g   Σ w g_AB               r→0: 0     r→∞: Σw
    h   Σ w [g_AB − 1]         r→0: −Σw   r→∞: 0
    D   4πrρ Σ w [g_AB − 1]    r→0: 0     r→∞: oscillates about 0
    T   4πrρ Σ w g_AB          r→0: 0     r→∞: 4πrρ Σw

  g   The weighted pair distribution: peak heights are pair density relative to
      the bulk average, and the baseline sits at Σw.
  h   Total correlation function — g with the bulk baseline subtracted off
      first, so peaks sit on zero and the excluded-volume region reads −Σw
      instead of 0.  h × absolute IS Soper's equation (20), the neutron G_n(r).
  D   Differential correlation function.  The 4πr factor offsets the way peak
      amplitude decays with distance, so far-field oscillations stay legible
      rather than flattening into the baseline; oscillates about zero.
      D × FZ is what the PDF community calls G(r).
  T   Total radial distribution function: D but keeping the bulk baseline, so
      it climbs as 4πrρΣw.  The area under a T peak is a coordination number.
      T × FZ reproduces the column this script used to call 't'.

Columns written before this rewrite map as total → g_unity, neutron → g_FZ,
t → T_FZ (T is not in the default $RDF_FUNCTIONS; add it to get that column back).

HYDROGEN IS TREATED AS DEUTERIUM
--------------------------------
NEUTRON_SCATTERING_LENGTHS['H'] holds deuterium's b = 6.671 fm, not protium's
b = −3.7406 fm, because LAMMPS dumps label both 'H'.  Every neutron-weighted
column is therefore for a fully deuterated sample, and the script says so at
startup whenever H is present.  For a light or mixed sample, pass b_override to
pair_weights() with the effective length b_eff = x_H b_H + x_D b_D (null water
is x_H ≈ 0.64).

COLUMN LAYOUT
-------------
Expects LAMMPS custom dump format with at minimum columns: id type element x y z
Column positions are read automatically from the ITEM: ATOMS header line.
Wrapped coordinates (x y z) and a constant atom count are assumed.
"""

import itertools
import os
import time
from datetime import date

import numpy as np
import freud

# freud uses TBB and does NOT read OMP_NUM_THREADS — left alone it takes every
# core on the node, which oversubscribes a shared node and ignores whatever the
# pipeline asked for. 0 is freud's own "all cores" default, so an unset variable
# behaves exactly as before.
#
# Measured on a 134784-atom frame at R_MAX=8 (6 partials), warmed up at each
# thread count: 1 -> 2.01 s, 2 -> 1.68x, 4 -> 2.65x, 6 -> 3.05x, 12 -> 3.30x.
# The speedup is real but strongly sub-linear, and core-seconds climb the whole
# way (2.0 -> 3.0 at 4 threads -> 7.3 at 12). Parallel efficiency is 66% at 4
# threads and 28% at 12, so 2-4 threads is the efficient range and more buys
# little for a lot of allocation. The limit is freud, not this script's serial
# parsing, which is only ~2.6% of the run at that setting.
#
# Note what does NOT help: adding frames. Each frame is a separate compute()
# call, so frames scale the work linearly at unchanged efficiency (3.37x at one
# frame, 3.39x at four). What improves the scaling is more work per call —
# more atoms (25% of them gives only 2.26x) or a larger R_MAX up to ~8.
# Empty is treated as unset, not as an error: distribution_run.sh defaults to
# $(nproc), which does not exist on macOS and exports the variable empty.
freud.parallel.set_num_threads(int(os.environ.get("OMP_NUM_THREADS", "").strip() or 0))
from scipy.ndimage import gaussian_filter1d

# =============================================================================
# CONFIGURATION — edit these variables between runs
# =============================================================================

# Input trajectory file
DUMP_FILE = os.environ.get("TRAJ", "../int_dump.lammpstrj")

# Output data table (CSV); set to None to skip
OUTPUT_CSV = "rdfs.csv"

# RDF parameters
R_MAX = float(os.environ.get("R_MAX", "20.0"))   # maximum r in Angstroms; must be < half shortest box dimension
BINS  = int(os.environ.get("RDF_BINS", "2000"))  # number of bins

# Correlation-function conventions.  Both take SEMICOLON-separated lists (see
# parse_keys) and the cross product is emitted, one column per combination,
# named <function>_<normalization>.  See CONVENTIONS in the module docstring.
#
# Everything here is RDF_-prefixed — env key and Python name alike — because
# 'NORMALIZATION' alone does not say which distribution it belongs to: vdos.py
# has its own VDOS_NORMALIZATION (phonon | unit_area), a completely different
# choice, and the two scripts run in one shared environment.
RDF_NORMALIZATION = os.environ.get("RDF_NORMALIZATION", "FZ;absolute;unity")
RDF_FUNCTIONS     = os.environ.get("RDF_FUNCTIONS", "g;h;D")

# Instrumental resolution matching: Gaussian sigma in Angstroms applied to every
# convention column (Soper used ~0.1 Å to suppress Fourier-termination ripples).
# Written as a *_broadened twin, so the raw curve is never lost.  0 disables.
RDF_RESOLUTION_SIGMA = float(os.environ.get("RDF_RESOLUTION_SIGMA", "0.1"))

# Resolution kernel shape. 'gaussian' is the generic stand-in; 'lorch' reproduces
# the modification function used in neutron glass diffraction (Wright's eq. 10
# and the equations that follow it), where the Q-space window
#     M(Q) = sin(dr Q)/(dr Q)   for Q <= Q_max,  0 above
# is cosine-transformed into a real-space peak function P(r) and convolved with
# the correlation function. A Gaussian matches its WIDTH but not its shape: the
# Lorch kernel has negative side lobes (~-9% of peak) that a Gaussian cannot
# reproduce, which matters between and beside strong peaks.
# 'modified_lorch' is Soper's eq. (60): rather than windowing in Q space, smear
# h(r) in r space with a uniform sphere of radius Δ, L' = 3/(4πΔ³) for |r| <= Δ.
# It is meant for h rather than g because h → 0 at large r, which is what lets
# the convolution be evaluated against a grid that stops at R_MAX — see
# apply_modified_lorch() for the measured difference. It is a 3D top hat, so it
# is everywhere non-negative — it cannot drive a correlation function negative
# between peaks the way the Q-space Lorch window's side lobes can — at the price
# of not corresponding to any real truncation.
RDF_RESOLUTION_MODE = os.environ.get("RDF_RESOLUTION_MODE", "gaussian")
if RDF_RESOLUTION_MODE not in ("gaussian", "lorch", "modified_lorch"):
    raise ValueError(
        f"Unknown RDF_RESOLUTION_MODE={RDF_RESOLUTION_MODE!r}; use 'gaussian', "
        f"'lorch' or 'modified_lorch'.")
RDF_LORCH_QMAX = float(os.environ.get("RDF_LORCH_QMAX", "0") or 0)   # Å⁻¹
# The parameter INSIDE M(Q) = sin(dr Q)/(dr Q). Defaults to the standard Lorch
# choice pi/Q_max, which puts M's first zero exactly at the truncation so the
# window closes smoothly. Papers usually quote the RESULTING real-space
# resolution (the FWHM of P(r)) rather than this parameter, and the two differ
# by about 1.73x — Wright's Q_max = 45.2 with pi/Q_max = 0.0695 yields FWHM
# 0.120 A, which is the number his text quotes. Setting this to the quoted
# resolution instead is the easy mistake, and it doubles the broadening.
RDF_LORCH_DR = float(os.environ.get("RDF_LORCH_DR", "0") or 0)       # Å; 0 = pi/Q_max
# Radius of the smearing sphere for 'modified_lorch', in Angstroms. NOT the
# resolution: the resulting FWHM is sqrt(2)*Δ, so a paper quoting 0.12 Å wants
# Δ = 0.0849. The same trap as RDF_LORCH_DR above, with a different factor.
RDF_MODIFIED_LORCH_DELTA = float(os.environ.get("RDF_MODIFIED_LORCH_DELTA", "0") or 0)
if RDF_RESOLUTION_MODE == "modified_lorch" and RDF_MODIFIED_LORCH_DELTA <= 0:
    raise ValueError(
        "RDF_RESOLUTION_MODE=modified_lorch needs RDF_MODIFIED_LORCH_DELTA (Å), the "
        "radius of the smearing sphere in Soper eq. (60). The resolution it produces "
        "is sqrt(2) x that, so divide a quoted FWHM by 1.4142 to get it.")

if RDF_RESOLUTION_MODE == "lorch":
    if RDF_LORCH_QMAX <= 0:
        raise ValueError(
            "RDF_RESOLUTION_MODE=lorch needs RDF_LORCH_QMAX (Å⁻¹), the truncation of "
            "the Fourier transform in the paper you are comparing against.")
    if RDF_LORCH_DR <= 0:
        RDF_LORCH_DR = np.pi / RDF_LORCH_QMAX

# Atoms per formula unit, needed only by the 'formula' normalization. Neutron
# diffraction papers frequently quote cross-sections per formula unit (per SiO2)
# rather than per atom; that is a factor of n, and it is the single most common
# reason a published curve sits a constant factor above a per-atom calculation.
RDF_ATOMS_PER_FORMULA_UNIT = float(os.environ.get("RDF_ATOMS_PER_FORMULA_UNIT", "0") or 0)



# n(r) output file; set to None to skip
OUTPUT_NR_CSV  = "nrs.csv"

# ---- Wright comparison output -------------------------------------------
# A single curve built to be overlaid directly on a published neutron T(r),
# bundling the three choices such a comparison needs so none has to be
# remembered separately:
#     function      T(r)
#     normalization per FORMULA UNIT   (RDF_ATOMS_PER_FORMULA_UNIT)
#     broadening    Lorch, dr = pi/Q_max   (RDF_WRIGHT_QMAX)
# Emitted to its own CSV/PNG rather than mixed into rdfs.csv, because it is a
# finished comparison artefact with its own units, not another convention column.
RDF_WRIGHT      = os.environ.get("RDF_WRIGHT", "no")
if RDF_WRIGHT not in ("yes", "no"):
    raise ValueError(f"Unknown RDF_WRIGHT={RDF_WRIGHT!r}; use 'yes' or 'no'.")
RDF_WRIGHT_QMAX = float(os.environ.get("RDF_WRIGHT_QMAX", "0") or 0)   # Å⁻¹
if RDF_WRIGHT == "yes":
    if RDF_WRIGHT_QMAX <= 0:
        raise ValueError(
            "RDF_WRIGHT=yes needs RDF_WRIGHT_QMAX (Å⁻¹), the truncation of the Fourier "
            "transform in the paper. Wright's vitreous-silica work used 45.2.")
    if RDF_ATOMS_PER_FORMULA_UNIT <= 0:
        raise ValueError(
            "RDF_WRIGHT=yes needs RDF_ATOMS_PER_FORMULA_UNIT — the number of atoms in "
            "the formula unit the paper quotes its cross-section per (SiO2 -> 3).")
OUTPUT_WRIGHT_CSV  = "wright.csv"

# =============================================================================
# END CONFIGURATION
# =============================================================================

# Prepend today's date (YYYYMMDD_) to every output filename.
def _dated(filename):
    return None if filename is None else f"{date.today():%Y%m%d}_{filename}"

OUTPUT_CSV     = _dated(OUTPUT_CSV)
OUTPUT_NR_CSV  = _dated(OUTPUT_NR_CSV)
OUTPUT_WRIGHT_CSV  = _dated(OUTPUT_WRIGHT_CSV)

# Coherent neutron scattering lengths (fm).  Add elements as needed.
# Values from NIST: https://www.ncnr.nist.gov/resources/n-lengths/
#
# 'H' is DELIBERATELY deuterium: LAMMPS dumps label both isotopes 'H', and this
# pipeline's samples are deuterated.  Protium is b = -3.7406 fm — opposite in
# sign, so every H-containing term would flip and peaks would point the other
# way.  _warn_hydrogen_is_deuterium() prints this at startup; for a light or
# mixed sample use pair_weights(..., b_override={'H': b_eff}).
NEUTRON_SCATTERING_LENGTHS = {
    'H':   6.671,    # deuterium (D), not protium — see note above
    'D':   6.671,
    'C':   6.6460,
    'N':   9.36,
    'O':   5.803,
    'Na':  3.63,
    'Mg':  5.375,
    'Al':  3.449,
    'Si':  4.1491,
    'P':   5.13,
    'S':   2.847,
    'Cl':  9.577,
    'K':   3.67,
    'Ca':  4.70,
    'Fe':  9.45,
    'Ni': 10.3,
    'Zr':  7.16,
    'Ba':  5.07,
}


def read_lammps_dump(filename):
    frames = []

    with open(filename) as f:
        while True:
            line = f.readline()
            if not line:
                break  # EOF

            # TIMESTEP
            timestep = int(f.readline().strip())

            # NUMBER OF ATOMS
            f.readline()
            n_atoms = int(f.readline().strip())

            # BOX BOUNDS — a triclinic dump carries a third tilt column per line,
            # and LAMMPS reports the bounding box of the tilted cell rather than
            # the cell lengths, so freud.box.Box would be silently wrong.  Same
            # guard as vdos_dynmat.py's read_reference_dump().
            f.readline()
            bounds = [f.readline().split() for _ in range(3)]
            tilt   = [float(p[2]) for p in bounds if len(p) > 2]
            if any(t != 0.0 for t in tilt):
                raise ValueError(
                    f"{filename} describes a triclinic box (tilt factors {tilt}). "
                    f"Its ITEM: BOX BOUNDS lines give the bounding box of the tilted "
                    f"cell, not the cell lengths, so the orthogonal box built here "
                    f"would mis-wrap every pair distance. Pass the tilt through to "
                    f"freud.box.Box, or run the analysis on an orthogonal cell."
                )
            (xlo, xhi), (ylo, yhi), (zlo, zhi) = (
                (float(p[0]), float(p[1])) for p in bounds
            )

            # ATOMS header — parse column positions dynamically
            header = f.readline().split()  # ['ITEM:', 'ATOMS', 'id', 'type', 'element', ...]
            cols = header[2:]
            col_element = cols.index('element')
            col_x       = cols.index('x')
            col_y       = cols.index('y')
            col_z       = cols.index('z')

            elements, positions = [], []
            for _ in range(n_atoms):
                parts = f.readline().split()
                elements.append(parts[col_element])
                positions.append([float(parts[col_x]), float(parts[col_y]), float(parts[col_z])])

            box = freud.box.Box(Lx=xhi - xlo, Ly=yhi - ylo, Lz=zhi - zlo)
            positions = np.array(positions)
            center = np.array([(xlo + xhi) / 2, (ylo + yhi) / 2, (zlo + zhi) / 2])
            positions -= center

            frames.append({
                'timestep':       timestep,
                'box':            box,
                'positions':      positions,
                'elements':       np.array(elements),
                'number_density': n_atoms / box.volume,
            })

    return frames


def compute_rdf(frames, get_a, get_b, self_pair=False):
    """
    Compute the frame-averaged RDF and n(r) using freud's built-in accumulation.
    freud accumulates across compute() calls when reset=False, so rdf.rdf and
    rdf.n_r at the end are already the properly normalized frame averages.

    g_AB is symmetric, but n(r) is not: freud counts SYSTEM points around each
    QUERY point, and pos_a is passed as the system, so mean_nr is A around each
    B.  build_coordination() names it accordingly and adds the reverse.

    bin_counts is freud's raw pair-count histogram accumulated over every frame,
    i.e. the number of A-B pairs actually observed in each bin.  That is the
    sampling number M for this partial: the relative statistical error of g(r)
    in a bin is ~1/sqrt(M), so it is what says whether the curve is converged.
    It is a measurement, not a model — nothing here estimates it.

    Returns: r, mean_g, mean_nr, bin_counts
    """
    rdf   = freud.density.RDF(bins=BINS, r_max=R_MAX)
    first = True

    for frame in frames:
        pos_a = get_a(frame)
        pos_b = get_b(frame)

        if len(pos_a) == 0 or len(pos_b) == 0:
            continue

        if self_pair:
            rdf.compute((frame['box'], pos_a), reset=first)
        else:
            rdf.compute((frame['box'], pos_a), query_points=pos_b, reset=first)

        first = False

    if first:
        empty = np.zeros(BINS)
        return rdf.bin_centers, empty, empty, empty

    return rdf.bin_centers, rdf.rdf, rdf.n_r, np.asarray(rdf.bin_counts)


# =============================================================================
# Sampling and cost reporting — see "How much data is enough" in the README
# =============================================================================

# Relative statistical error of a binned average is ~1/sqrt(M), where M is the
# number of independent contributions to that bin.  The ladder:
#     M = 1e2   10%    exploratory only
#     M = 1e3    3%    pass mark: peak positions, coordination numbers
#     M = 1e4    1%    publication / comparison against measured data
#     M = 1e6  0.1%    past the point of diminishing returns
SAMPLING_TARGET = 1e3


def sampling_verdict(m):
    """(relative error, verdict) for a sampling number M."""
    if m <= 0:
        return float('inf'), 'EMPTY'
    error = 1.0 / np.sqrt(m)
    return error, ('ok' if m >= SAMPLING_TARGET else 'LOW')


def report_sampling(rows, note=''):
    """
    Print achieved sampling per row and name the limiting one.

    rows is [(label, where, M)].  Every partial is listed rather than an
    aggregate, because the minority species is normally what limits the result
    and an aggregate hides it entirely.
    """
    print(f"\nSampling achieved (relative error ~ 1/sqrt(M), target M >= {SAMPLING_TARGET:.0e}):")
    width = max((len(r[0]) for r in rows), default=8)
    worst = None
    for label, where, m in rows:
        error, verdict = sampling_verdict(m)
        print(f"  {label.ljust(width)}  {where:<22}  M = {m:9.3g}  {100*error:6.2f}%  {verdict}")
        if worst is None or m < worst[2]:
            worst = (label, where, m)
    if worst is not None:
        error, verdict = sampling_verdict(worst[2])
        print(f"  limiting: {worst[0]} at {100*error:.2f}% ({verdict})")
        if verdict == 'LOW':
            print(f"  -> below the {SAMPLING_TARGET:.0e} pass mark; add frames or atoms"
                  f" (M scales linearly in both)")
    if note:
        print(f"  {note}")


def report_cost(t_serial, t_parallel, threads, m_limiting, scaling=''):
    """
    Print wall time, the threads actually in use, and core-seconds.

    Three numbers because one is not enough: wall is what you wait for, threads
    is what you got (not what you asked for — freud and numba ignore
    OMP_NUM_THREADS), and core-seconds is the portable figure that survives
    comparison across machines and allocations.
    """
    wall = t_serial + t_parallel
    # freud reports 0 for "TBB default" rather than a count, so say that instead
    # of substituting a guess and presenting it as a measurement.
    if threads:
        label, effective = f"{threads} (freud TBB, reported)", threads
    else:
        effective = os.cpu_count() or 1
        label = f"TBB default = all cores (~{effective}, not reported by freud)"
    core_seconds = wall * max(effective, 1)
    serial_fraction = t_serial / wall if wall > 0 else 0.0
    print(f"\nCost: wall {wall:.2f} s   threads {label}   "
          f"core-seconds {core_seconds:.1f}")
    print(f"      parse {t_serial:.2f} s serial ({100*serial_fraction:.1f}%) + "
          f"compute {t_parallel:.2f} s")
    if m_limiting > 0 and wall > 0:
        print(f"      efficiency {m_limiting/wall:.3g} samples/wall-s, "
              f"{m_limiting/core_seconds:.3g} samples/core-s (limiting partial)")
        print(f"      measured: parallel efficiency 66% at 4 threads, 28% at 12; "
              f"frames add work but not scaling")
    if scaling:
        print(f"      {scaling}")


def build_pair_getters(elements):
    """
    Auto-generate {label: (get_a, get_b, self_pair)} for all unique element pairs
    found in the trajectory.  self_pair=True for A-A pairs so freud excludes i=j.
    """
    def make_getter(el):
        return lambda f: f['positions'][f['elements'] == el]

    getters = {el: make_getter(el) for el in elements}

    pairs = {}
    for a, b in itertools.combinations_with_replacement(sorted(elements), 2):
        label = f'{a}-{b}'
        pairs[label] = (getters[a], getters[b], a == b)

    return pairs


def get_concentrations(frames):
    """Return average mole fractions {element: c} across all frames."""
    counts = {}
    for frame in frames:
        for el, n in zip(*np.unique(frame['elements'], return_counts=True)):
            counts[el] = counts.get(el, 0) + n
    total = sum(counts.values())
    return {el: n / total for el, n in counts.items()}



# =============================================================================
# Correlation-function conventions — see CONVENTIONS in the module docstring
# =============================================================================

# Weight definitions, for the printed table.  f = 2 - delta_AB throughout.
NORMALIZATION_EQUATIONS = {
    'unity':    'w = f c_A c_B',
    'FZ':       'w = f c_A c_B b_A b_B / <b>²',
    'absolute': 'w = f c_A c_B b_A b_B / 100',
    'formula':  'w = n f c_A c_B b_A b_B / 100   (n = atoms per formula unit)',
}

NORMALIZATION_UNITS = {'unity': '', 'FZ': '', 'absolute': 'barn/sr/atom',
                       'formula': 'barn/sr/formula-unit'}

FUNCTION_EQUATIONS = {
    'g': 'Σ w g_AB',
    'h': 'Σ w [g_AB − 1]',
    'D': '4πrρ Σ w [g_AB − 1]',
    'T': '4πrρ Σ w g_AB',
}


def parse_keys(value, valid, name):
    """
    Parse a semicolon-separated convention list, failing loudly on unknown keys.

    The separator is ';' and not ',' because submit_pipeline.sh passes settings
    through `sbatch --export`, which is itself comma-delimited and silently
    truncates a value at the first embedded comma — a comma-separated list would
    arrive on the cluster as its first entry alone, with no error.  Same reason
    bad_freud.py takes ELEMENTS as "Si;O;H".
    """
    if ',' in value:
        raise ValueError(
            f"{name}={value!r} uses ',' but the separator is ';' — a comma would be "
            f"truncated by `sbatch --export` in submit_pipeline.sh. Write it as "
            f"{value.replace(',', ';')!r}."
        )
    keys = [k.strip() for k in value.split(';') if k.strip()]
    if not keys:
        raise ValueError(f"{name} is empty; choose from {list(valid)}.")
    unknown = [k for k in keys if k not in valid]
    if unknown:
        raise ValueError(f"Unknown {name} {unknown}; choose from {list(valid)}.")
    return keys


def pair_weights(normalization, elements, concentrations, b_override=None):
    """
    Return {pair_label: w_AB} for one normalization, or None if it is unusable.

    b_override supplies effective scattering lengths for isotope mixtures — e.g.
    {'H': 0.64 * (-3.7406) + 0.36 * 6.671} for null water — and falls back to
    NEUTRON_SCATTERING_LENGTHS for every element it does not name.
    """
    pairs = list(itertools.combinations_with_replacement(sorted(elements), 2))

    if normalization == 'unity':
        # b = 1 for every element, so no scattering lengths are consulted at all.
        return {f'{a}-{bl}': (1 if a == bl else 2) * concentrations[a] * concentrations[bl]
                for a, bl in pairs}

    b = dict(b_override or {})
    missing = [el for el in elements if el not in b and el not in NEUTRON_SCATTERING_LENGTHS]
    if missing:
        print(f"Warning: no scattering length for {missing}; '{normalization}' columns skipped.")
        return None
    for el in elements:
        b.setdefault(el, NEUTRON_SCATTERING_LENGTHS[el])

    b_mean = sum(concentrations[el] * b[el] for el in elements)

    if normalization == 'formula':
        # Per formula unit rather than per atom. Wright's eq. (10) writes the
        # baseline as T0 = 4*pi*r*rho0*(sum_j b_j)^2 with rho0 in formula units
        # per A^3 and sum_j b_j the TOTAL scattering length of one unit. Since
        # sum_j b_j = n<b> and rho0 = rho_atom/n, that is
        #     (rho_a/n)(n<b>)^2 = n * rho_a * <b>^2
        # i.e. exactly n times the per-atom result, for any composition. Folding
        # the n into the weights keeps rho per atom everywhere else in this
        # script, so only this one factor changes.
        if RDF_ATOMS_PER_FORMULA_UNIT <= 0:
            raise SystemExit(
                "rdf_freud.py: RDF_NORMALIZATION=formula needs RDF_ATOMS_PER_FORMULA_UNIT,\n"
                "  the number of atoms in the formula unit the paper quotes (SiO2 -> 3).\n"
                f"  This system's concentrations are "
                + ", ".join(f"{el}={concentrations[el]:.4f}" for el in sorted(elements))
                + ",\n  so the smallest integer formula implies n = "
                + f"{1.0/min(concentrations[el] for el in elements):.2f} — set it explicitly."
            )
        denom = 100.0 / RDF_ATOMS_PER_FORMULA_UNIT
    elif normalization == 'FZ':
        # FZ divides out <b>², which is a nearly-cancelling sum for H-rich or
        # null samples; the resulting amplification is numerical, not physical.
        if abs(b_mean) < 1e-9:
            print("Warning: <b> = 0 (null mixture); 'FZ' columns skipped — "
                  "use 'absolute', which never divides by <b>.")
            return None
        if abs(b_mean) < 1.0:
            print(f"Warning: <b> = {b_mean:.4f} fm is small, so 'FZ' divides by "
                  f"<b>² = {b_mean ** 2:.4f} fm² and amplifies a nearly-cancelling "
                  f"sum. Prefer the 'absolute' columns for this composition.")
        denom = b_mean ** 2
    else:                                    # 'absolute': fm² -> barn, no <b>
        denom = 100.0

    return {f'{a}-{bl}': (1 if a == bl else 2) * concentrations[a] * concentrations[bl]
                         * b[a] * b[bl] / denom
            for a, bl in pairs}


def apply_function(func, r, weights, partial_results, rho_mean):
    """Build one convention curve from pair weights and the partial g_AB(r)."""
    y = np.zeros(len(r))
    for label, w in weights.items():
        g = partial_results[label][1]
        y += w * (g - 1.0) if func in ('h', 'D') else w * g
    if func in ('D', 'T'):
        y = y * 4 * np.pi * r * rho_mean
    return y


def _column_units(func, norm):
    base = NORMALIZATION_UNITS[norm]
    if func in ('D', 'T'):
        return f'{base} Å⁻²'.strip() if base else 'Å⁻²'
    return base or 'dimensionless'


def build_conventions(partial_results, elements, concentrations, rho_mean,
                      normalizations, functions, b_override=None):
    """
    Build every requested (function, normalization) column from the partials.

    Returns (results, meta) where results is {label: (r, y)} and meta is
    {label: {...}} carrying the defining equation, Σw, asymptotic limits and
    units.  The limits are printed, not drawn: the plotting pipeline draws only
    curves, and a printed limit checked against one is the stronger test.
    """
    r = next(iter(partial_results.values()))[0]
    results, meta = {}, {}

    for norm in normalizations:
        weights = pair_weights(norm, elements, concentrations, b_override)
        if weights is None:
            continue
        sum_w = sum(weights.values())

        for func in functions:
            label = f'{func}_{norm}'
            results[label] = (r, apply_function(func, r, weights, partial_results, rho_mean))
            meta[label] = {
                'equation':  f'{FUNCTION_EQUATIONS[func]},  {NORMALIZATION_EQUATIONS[norm]}',
                'sum_w':     sum_w,
                'units':     _column_units(func, norm),
                'limit_0':   {'g': '0', 'h': f'{-sum_w:.4f}', 'D': '0', 'T': '0'}[func],
                'limit_inf': {'g': f'{sum_w:.4f}', 'h': '0', 'D': '0 (oscillates)',
                              'T': f'4πrρ·{sum_w:.4f}'}[func],
            }

    if not results:
        raise ValueError(
            "No convention columns could be built. Every requested normalization "
            "was skipped — see the warnings above."
        )
    return results, meta


def lorch_peak_function(r, delta_r, q_max):
    """
    Wright's real-space peak function P(r), the cosine transform of the Q-space
    modification function M(Q) = sin(dr*Q)/(dr*Q) truncated at q_max.

    The integral has a closed form in sine integrals:

        int_0^Qmax  sin(dr Q) cos(r Q) / (dr Q) dQ
            = [ Si((dr + r) Qmax) + Si((dr - r) Qmax) ] / (2 dr)

    so no numerical quadrature is needed. Returned unnormalized; the caller
    normalizes it to unit area, which is what makes it a resolution kernel
    rather than a weighted one (Wright carries b_j b_k in P_jk; here those
    weights already live in the convention's w_AB).
    """
    from scipy.special import sici
    return (sici((delta_r + r) * q_max)[0] + sici((delta_r - r) * q_max)[0]) / (2 * delta_r)


def apply_lorch(r, y, delta_r, q_max, chunk=512):
    """
    Convolve with the Lorch peak function, including the reflected term:

        y'(r) = int_0^inf y(r') [ P(r - r') - P(r + r') ] dr'

    The P(r + r') term enforces the odd symmetry of the correlation function
    about the origin. It is negligible beyond a few kernel widths but is what
    keeps the excluded-volume region correct at small r, which is exactly the
    region used to check the normalization.
    """
    dr_grid = r[1] - r[0]
    norm = np.trapezoid(lorch_peak_function(np.linspace(-2.0, 2.0, 20001), delta_r, q_max),
                        np.linspace(-2.0, 2.0, 20001))
    out = np.empty_like(y)
    for start in range(0, len(r), chunk):
        stop = min(start + chunk, len(r))
        block = r[start:stop, None]                       # (chunk, 1)
        kernel = (lorch_peak_function(block - r[None, :], delta_r, q_max)
                  - lorch_peak_function(block + r[None, :], delta_r, q_max))
        out[start:stop] = kernel @ y * dr_grid / norm
    return out


def modified_lorch_kernel(r_block, r_grid, delta):
    """
    The radial kernel K(r, r') for Soper's modified Lorch function, eq. (60):

        L'(r, Δ) = 3/(4πΔ³)   |r| <= Δ
                 = 0          |r| >  Δ

    a uniform sphere of radius Δ and unit volume integral. Unlike the standard
    Lorch function this is applied in r space, not as a window in Q space, so it
    is a genuine 3D convolution rather than a 1D one.

    IN PRACTICE THIS IS APPLIED TO h(r), the baseline-subtracted correlation
    function, and the derivation below is written with that in mind. h is not an
    arbitrary choice of column: it is the one this kernel can be evaluated on
    honestly against a truncated grid. See apply_modified_lorch() for the
    measurement that makes the point.

    For spherically symmetric h and L' the 3D convolution collapses to a single
    radial integral,

        (h * L')(r) = (2π/r) ∫ dr' r' h(r') ∫_{|r-r'|}^{r+r'} du u L'(u)

    and with L' constant out to Δ the inner integral is elementary:

        r (h * L')(r) = 3/(4Δ³) ∫ dr' [r' h(r')] K(r, r')
        K(r, r')      = min(r + r', Δ)² − min(|r − r'|, Δ)²

    The two min()s are what carry the geometry: when |r − r'| >= Δ the shells do
    not overlap and K is identically zero, so no separate cutoff is needed. This
    was checked against direct 3D quadrature (agreement to ~1e-10) and against
    the requirement that a constant convolve to itself.
    """
    return (np.minimum(r_block + r_grid, delta) ** 2
            - np.minimum(np.abs(r_block - r_grid), delta) ** 2)


def apply_modified_lorch(r, y, delta, r_weighted, chunk=512):
    """
    Convolve one column with the uniform sphere of Soper eq. (60).

    MEANT FOR h(r). h → 0 at large r, and that is what makes this kernel usable
    on a finite grid: the integral wants r' out to r + Δ, and past R_MAX there
    is nothing there. For h the missing shell contributes nothing, because h is
    already zero out there. For a column carrying a baseline it contributes the
    baseline, and losing it shows.

    Measured on an R_MAX = 8 Å grid with Δ = 0.2, structure decaying to its
    far-field value:

        h convolved, last Δ of the grid : 0.00000        (truth 0)
        g convolved, last Δ of the grid : 0.50 .. 1.00   (truth 1)

    Analytically the two carry the same information — the kernel has unit volume
    integral, so it maps the constant 1 to itself and (g * L') = (h * L') + 1
    exactly. On a truncated grid that identity holds in the bulk (agreement to
    1.6e-4 here) and fails at the edges (0.50 discrepancy in the last 2Δ),
    because the shell beyond R_MAX is precisely the part that would have
    supplied the baseline. Subtracting the baseline first is what removes the
    artifact rather than hiding it.

    Nothing here refuses another column, and the r_weighted flag below is what a
    different one would need. It says whether the column already carries a
    factor of r, because the convolution acts on the underlying 3D radial
    function and not on the plotted curve: h is that function, so it is
    multiplied by r going in and divided by r coming out, while D and T are
    4πrρ × it and pass through as they stand. Getting that backwards changes the
    answer without changing its shape enough to notice, which is why it is an
    explicit argument rather than a guess — but note that D and T also carry the
    r-weighting that makes the edge loss above worse, not better.
    """
    if delta <= 0:
        raise ValueError("modified Lorch needs a positive Δ.")
    dr_grid = r[1] - r[0]
    g = y if r_weighted else r * y
    out = np.empty_like(y)
    for start in range(0, len(r), chunk):
        stop = min(start + chunk, len(r))
        kernel = modified_lorch_kernel(r[start:stop, None], r[None, :], delta)
        out[start:stop] = kernel @ g * dr_grid
    out *= 3.0 / (4.0 * delta ** 3)
    return out if r_weighted else out / r


def modified_lorch_fwhm(delta):
    """
    Real-space resolution of eq. (60), for checking against a quoted number.

    Far from the origin the sphere smears a shell by its own projection onto one
    axis, p(x) = 3(Δ² − x²)/(4Δ³), a parabola on [−Δ, Δ]. Half maximum sits at
    x = Δ/√2, so

        FWHM = √2 Δ ≈ 1.4142 Δ

    Confirmed numerically against a delta shell put through the full radial
    convolution. Δ is therefore NOT the resolution — it is about 1.41x smaller,
    the same trap as the standard Lorch function's Δr.
    """
    return np.sqrt(2.0) * delta


def lorch_fwhm(delta_r, q_max):
    """
    FWHM of the resulting peak function — the number papers actually quote as
    their real-space resolution, so it is what you check against the text.
    Also reports whether the kernel came out double-humped, which happens when
    delta_r*q_max > pi and is the signature of having fed in the quoted
    resolution instead of the M(Q) parameter.
    """
    x = np.linspace(0.0, 20.0 / q_max, 20001)
    p = lorch_peak_function(x, delta_r, q_max)
    p = p / p.max()
    half = x[p >= 0.5]
    return 2 * half.max(), bool(x[np.argmax(p)] > 1e-3)


def broaden(results, meta, sigma, dr, mode='gaussian', delta_r=0.0, q_max=0.0,
            delta=0.0):
    """
    Add a *_broadened twin of every convention column, so modeled peaks are not
    sharper than measured ones purely for instrumental reasons.  Mutates both
    dicts.

    'gaussian'        the generic stand-in (Soper used ~0.1 Å).
    'lorch'           the modification function neutron glass diffraction
                      actually uses: a window in Q space, cosine-transformed to
                      a 1D peak function. Matches a Gaussian in width but has
                      negative side lobes a Gaussian cannot produce.
    'modified_lorch'  Soper eq. (60): a uniform sphere of radius Δ convolved in
                      r space instead of a window applied in Q space, and meant
                      for h — the baseline-subtracted form is the one whose
                      convolution survives the grid ending at R_MAX. Being a 3D
                      top hat it is strictly non-negative — no side lobes at all
                      — so it cannot push a correlation function negative
                      between peaks the way 'lorch' can. The cost is that it has
                      no Q-space counterpart, so it does not correspond to any
                      particular measurement's truncation.
    """
    if mode == 'gaussian' and sigma <= 0:
        return
    for label in list(results):
        r, y = results[label]
        if mode == 'modified_lorch':
            # h is the intended column and is not r-weighted; D and T carry an
            # explicit factor of r. See apply_modified_lorch().
            r_weighted = label.split('_')[0] in ('D', 'T')
            wide = apply_modified_lorch(r, y, delta, r_weighted)
            note = (f"⊗ modified Lorch (Soper eq. 60) Δ={delta:.4f} Å "
                    f"-> resolution FWHM {modified_lorch_fwhm(delta):.4f} Å")
        elif mode == 'lorch':
            wide = apply_lorch(r, y, delta_r, q_max)
            fwhm, double = lorch_fwhm(delta_r, q_max)
            note = (f"⊗ Lorch Δr={delta_r:.4f} Å, Q_max={q_max} Å⁻¹ "
                    f"-> resolution FWHM {fwhm:.4f} Å"
                    + ("  [WARNING: kernel is double-humped, Δr*Q_max > pi — Δr is "
                       "probably the paper's quoted resolution, not its M(Q) parameter]"
                       if double else ""))
        else:
            wide = gaussian_filter1d(y, sigma / dr, mode='nearest')
            note = f"⊗ Gaussian σ={sigma} Å"
        results[f'{label}_broadened'] = (r, wide)
        meta[f'{label}_broadened'] = dict(
            meta[label], equation=f"{meta[label]['equation']}, {note}"
        )


def print_convention_table(meta):
    """
    Print each column's defining equation, Σw, and limits.  These are the QC
    numbers: a printed limit checked against the curve beats a label, and for
    h_absolute the Σw must reproduce Soper's Table 1 column sum.
    """
    width = max(len(name) for name in meta)
    print("\nConvention columns (check these limits against the curves):")
    print(f"  {'column'.ljust(width)}  {'Σw':>10}  {'r→0':>10}  {'r→∞':>16}  units")
    for name, m in meta.items():
        print(f"  {name.ljust(width)}  {m['sum_w']:>10.4f}  {m['limit_0']:>10}  "
              f"{m['limit_inf']:>16}  {m['units']}")
        print(f"  {' '.ljust(width)}    {m['equation']}")


def warn_hydrogen_is_deuterium(elements):
    """Say out loud that H is being weighted as D — the sign of b_H rides on it."""
    if 'H' in elements:
        print("Note: 'H' is weighted as DEUTERIUM (b = 6.671 fm). Protium is "
              "b = -3.7406 fm, so a light sample would flip the sign of every "
              "H term. All neutron columns below are for a deuterated sample.")


def build_coordination(nr_raw, concentrations):
    """
    Name the coordination numbers by direction, and add the reverse of each.

    freud's n_r counts SYSTEM points around each QUERY point, and compute_rdf
    passes A as the system and B as the query for label 'A-B' — so the raw curve
    is A around each B.  Every A-B pair within r is counted once from each side,
    N_B·n(A around B) = N_A·n(B around A), so the reverse direction follows from
    the concentration ratio exactly (given the constant atom count assumed here)
    and needs no second freud pass.
    """
    out = {}
    for label, (r, nr) in nr_raw.items():
        a, b = label.split('-')
        out[f'{a}_around_{b}'] = (r, nr)
        if a != b:
            out[f'{b}_around_{a}'] = (r, nr * concentrations[b] / concentrations[a])
    return out


def build_wright(partial_results, elements, concentrations, rho_mean, q_max, n_formula):
    """
    T(r) on a published neutron-diffraction footing, ready to overlay.

    Three choices, all of which have to agree with the paper or the curves will
    not lie on top of each other:

      function       T(r) = 4*pi*r*rho*sum w g_AB, which is Wright's
                     T(r) = D(r) + T0(r) — it oscillates about the straight
                     baseline T0 rather than about zero.
      normalization  per FORMULA UNIT, so a factor n above the per-atom result
                     (see the 'formula' entry in CONVENTIONS).
      broadening     the Lorch modification function with dr = pi/q_max, applied
                     as the full P(r-r') - P(r+r') convolution.

    Returns (columns, info) where columns is {label: array} and info carries the
    numbers worth checking against the paper: the T0 slope, Sum(w), and the
    resolution FWHM.
    """
    r = next(iter(partial_results.values()))[0]
    weights = pair_weights('formula', elements, concentrations)
    if weights is None:
        raise SystemExit("rdf_freud.py: RDF_WRIGHT=yes could not build weights — see above.")

    sum_w = sum(weights.values())
    delta_r = np.pi / q_max
    fwhm, double_humped = lorch_fwhm(delta_r, q_max)

    t_raw = apply_function('T', r, weights, partial_results, rho_mean)
    t_broad = apply_lorch(r, t_raw, delta_r, q_max)
    baseline = 4 * np.pi * r * rho_mean * sum_w        # Wright's T0(r)

    info = {
        'sum_w': sum_w, 'delta_r': delta_r, 'fwhm': fwhm, 'q_max': q_max,
        'n_formula': n_formula, 'rho': rho_mean,
        't0_slope': 4 * np.pi * rho_mean * sum_w, 'double_humped': double_humped,
    }
    return {'T_wright': t_broad, 'T_wright_unbroadened': t_raw, 'T0_baseline': baseline}, info


def report_wright(info):
    """Print the numbers that decide whether the comparison is set up right."""
    print("\nWright-comparison output (overlay T_wright on the published T(r)):")
    print(f"  normalization  per formula unit, n = {info['n_formula']:g} atoms   "
          f"Σw = {info['sum_w']:.4f} barn/sr/formula-unit")
    print(f"  broadening     Lorch Δr = π/Q_max = {info['delta_r']:.4f} Å at "
          f"Q_max = {info['q_max']:g} Å⁻¹  ->  resolution FWHM {info['fwhm']:.4f} Å")
    print(f"                 ^ compare that FWHM against the resolution the paper quotes")
    if info['double_humped']:
        print("  WARNING: the kernel is double-humped — check Q_max")
    print(f"  density        ρ = {info['rho']:.6f} atoms/Å³ (per atom; the n is in the weights)")
    print(f"  T0 slope       4πρΣw = {info['t0_slope']:.4f} barn/sr/formula-unit/Å³")
    print(f"                 ^ THE check: measure the slope of the paper's average-density")
    print(f"                   line and it must equal this. It is fixed by composition and")
    print(f"                   density alone, so it needs no fitting.")


def save_csv(results, filename):
    """Save results dict {label: (r, y)} to CSV with one column per label."""
    r = next(iter(results.values()))[0]
    header = 'r_Angstrom,' + ','.join(results.keys())
    data = np.column_stack([r] + [y for _, y in results.values()])
    np.savetxt(filename, data, delimiter=',', header=header, comments='', fmt='%.6f')
    print(f"Data table saved to {filename}")


if __name__ == '__main__':
    # Validate the convention keys before reading the trajectory, so a typo costs
    # a second rather than a full parse of a multi-GB dump.
    normalizations = parse_keys(RDF_NORMALIZATION, NORMALIZATION_EQUATIONS, 'RDF_NORMALIZATION')
    functions      = parse_keys(RDF_FUNCTIONS, FUNCTION_EQUATIONS, 'RDF_FUNCTIONS')
    print(f"RDF_NORMALIZATION: {normalizations}")
    print(f"RDF_FUNCTIONS:     {functions}")

    print(f"Reading trajectory: {DUMP_FILE}")
    _t_parse_start = time.time()
    frames = read_lammps_dump(DUMP_FILE)
    t_parse = time.time() - _t_parse_start

    atom_counts = [len(f['positions']) for f in frames]
    print(f"Frames read:  {len(frames)}")
    print(f"Atoms range:  {min(atom_counts)} – {max(atom_counts)}")
    print(f"Atoms mean:   {np.mean(atom_counts):.1f}")

    elements = sorted(set(frames[0]['elements'].tolist()))
    print(f"Elements:     {elements}")
    warn_hydrogen_is_deuterium(elements)

    concentrations = get_concentrations(frames)
    print("Concentrations: " + ", ".join(f"{el}={c:.3f}" for el, c in concentrations.items()))

    rho = np.array([f['number_density'] for f in frames])   # atoms/Å³, shape (n_frames,)
    print(f"Number density: {rho.mean():.6e} atoms/Å³")

    min_box_dim = min(min(f['box'].Lx, f['box'].Ly, f['box'].Lz) for f in frames)
    print(f"Shortest box edge (min over frames): {min_box_dim:.3f} Å")
    if R_MAX >= min_box_dim / 2:
        raise ValueError(
            f"R_MAX={R_MAX} is too large: freud requires r_max < half the shortest "
            f"box dimension ({min_box_dim / 2:.3f} Å). Lower R_MAX in the config section."
        )

    pairs = build_pair_getters(elements)
    gr_results = {}
    nr_raw     = {}
    sampling   = []
    _t_compute_start = time.time()
    for name, (get_a, get_b, is_self) in pairs.items():
        print(f"Computing RDF: {name}...")
        r, g, nr, counts = compute_rdf(frames, get_a, get_b, self_pair=is_self)
        gr_results[name] = (r, g)
        nr_raw[name]     = (r, nr)

        # M at the first peak: the bin people actually read numbers off. The
        # 0.5 A floor skips the empty excluded-volume bins, whose g(r) = 0
        # would otherwise win the argmax on a noisy curve.
        physical = r > 0.5
        if physical.any() and counts.sum() > 0:
            peak = np.flatnonzero(physical)[np.argmax(g[physical])]
            sampling.append((name, f"first peak {r[peak]:.2f} Å", float(counts[peak])))
        else:
            sampling.append((name, "no pairs found", 0.0))
    t_compute = time.time() - _t_compute_start

    nr_results = build_coordination(nr_raw, concentrations)

    conventions, meta = build_conventions(
        gr_results, elements, concentrations, rho.mean(), normalizations, functions
    )
    r_grid = next(iter(gr_results.values()))[0]
    broaden(conventions, meta, RDF_RESOLUTION_SIGMA, r_grid[1] - r_grid[0],
            mode=RDF_RESOLUTION_MODE, delta_r=RDF_LORCH_DR, q_max=RDF_LORCH_QMAX,
            delta=RDF_MODIFIED_LORCH_DELTA)
    print_convention_table(meta)
    gr_results.update(conventions)

    report_sampling(sampling)
    report_cost(
        t_parse, t_compute,
        threads=freud.parallel.get_num_threads() or os.cpu_count(),
        m_limiting=min((m for _, _, m in sampling), default=0.0),
        scaling=(f"cost ~ frames x N x R_MAX^3 (measured 8.2x per R_MAX doubling at 8->16 Å); "
                 f"R_MAX={R_MAX} here"),
    )

    if OUTPUT_CSV is not None:
        save_csv(gr_results, OUTPUT_CSV)
    if OUTPUT_NR_CSV is not None:
        save_csv(nr_results, OUTPUT_NR_CSV)


    if RDF_WRIGHT == 'yes':
        wright_cols, wright_info = build_wright(
            {k: v for k, v in gr_results.items() if '-' in k},   # partials only
            elements, concentrations, rho.mean(),
            RDF_WRIGHT_QMAX, RDF_ATOMS_PER_FORMULA_UNIT,
        )
        report_wright(wright_info)
        if OUTPUT_WRIGHT_CSV is not None:
            save_csv({k: (r_grid, v) for k, v in wright_cols.items()}, OUTPUT_WRIGHT_CSV)
