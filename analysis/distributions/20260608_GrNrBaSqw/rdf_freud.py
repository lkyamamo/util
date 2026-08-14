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
- rdfs.png  — subplot grid: one panel per element pair, plus one per convention column
- rdfs.csv  — r (Å), partial g_AB(r), one column per requested convention
                                                              (set OUTPUT_CSV=None to skip)
- nrs.png   — cumulative coordination number n(r) plots       (set OUTPUT_NR_PLOT=None to skip)
- nrs.csv   — r (Å), n(r) both directions per pair            (set OUTPUT_NR_CSV=None to skip)

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
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# =============================================================================
# CONFIGURATION — edit these variables between runs
# =============================================================================

# Input trajectory file
DUMP_FILE = os.environ.get("TRAJ", "../int_dump.lammpstrj")

# Output plot file
OUTPUT_PLOT = "rdfs.png"

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

# Plot layout: how many columns in the subplot grid
PLOT_NCOLS = 2

# DPI for saved plot
PLOT_DPI = 150

# n(r) output files; set to None to skip
OUTPUT_NR_PLOT = "nrs.png"
OUTPUT_NR_CSV  = "nrs.csv"

# =============================================================================
# END CONFIGURATION
# =============================================================================

# Prepend today's date (YYYYMMDD_) to every output filename.
def _dated(filename):
    return None if filename is None else f"{date.today():%Y%m%d}_{filename}"

OUTPUT_PLOT    = _dated(OUTPUT_PLOT)
OUTPUT_CSV     = _dated(OUTPUT_CSV)
OUTPUT_NR_PLOT = _dated(OUTPUT_NR_PLOT)
OUTPUT_NR_CSV  = _dated(OUTPUT_NR_CSV)

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
}

NORMALIZATION_UNITS = {'unity': '', 'FZ': '', 'absolute': 'barn/sr/atom'}

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

    if normalization == 'FZ':
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
    {label: {...}} carrying the defining equation, Σw, asymptotic limits, units,
    and the y value of the plot's reference line (None for no line).
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
                'reference': {'g': sum_w, 'h': 0.0, 'D': 0.0, 'T': None}[func],
            }

    if not results:
        raise ValueError(
            "No convention columns could be built. Every requested normalization "
            "was skipped — see the warnings above."
        )
    return results, meta


def broaden(results, meta, sigma, dr):
    """
    Add a *_broadened twin of every convention column: Gaussian resolution
    matching, so modeled peaks are not sharper than measured ones purely for
    instrumental reasons (Soper used ~0.1 Å).  Mutates both dicts.
    """
    if sigma <= 0:
        return
    for label in list(results):
        r, y = results[label]
        results[f'{label}_broadened'] = (r, gaussian_filter1d(y, sigma / dr, mode='nearest'))
        meta[f'{label}_broadened'] = dict(
            meta[label], equation=f"{meta[label]['equation']}, ⊗ Gaussian σ={sigma} Å"
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


def save_csv(results, filename):
    """Save results dict {label: (r, y)} to CSV with one column per label."""
    r = next(iter(results.values()))[0]
    header = 'r_Angstrom,' + ','.join(results.keys())
    data = np.column_stack([r] + [y for _, y in results.values()])
    np.savetxt(filename, data, delimiter=',', header=header, comments='', fmt='%.6f')
    print(f"Data table saved to {filename}")


def plot_rdfs(results, meta):
    """
    Plot partials and convention columns.  meta carries per-column units and the
    reference level, since these panels no longer share one y axis: a partial is
    a dimensionless g(r) about 1, while h_absolute is barn/sr/atom about 0.
    """
    n = len(results)
    ncols = PLOT_NCOLS
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4 * nrows), squeeze=False)
    axes = axes.flatten()

    for ax, (name, (r, g)) in zip(axes, results.items()):
        ax.plot(r, g)
        m         = meta.get(name)
        reference = 1.0 if m is None else m['reference']       # partials: g(r) -> 1
        if reference is not None:
            ax.axhline(reference, color='gray', linestyle='--', linewidth=0.8)
        ax.set_xlabel('r (Å)')
        ax.set_ylabel('g(r)' if m is None else f"{name} ({m['units']})")
        ax.set_title(name)

    # hide any unused subplots
    for ax in axes[n:]:
        ax.set_visible(False)

    fig.tight_layout()
    if OUTPUT_PLOT is not None:
        fig.savefig(OUTPUT_PLOT, dpi=PLOT_DPI)
        print(f"Plot saved to {OUTPUT_PLOT}")
    plt.close(fig)


def plot_nrs(nr_results):
    n = len(nr_results)
    ncols = PLOT_NCOLS
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4 * nrows), squeeze=False)
    axes = axes.flatten()

    for ax, (name, (r, nr)) in zip(axes, nr_results.items()):
        ax.plot(r, nr)
        ax.set_xlabel('r (Å)')
        ax.set_ylabel('n(r)')
        ax.set_title(name)

    for ax in axes[n:]:
        ax.set_visible(False)

    fig.tight_layout()
    if OUTPUT_NR_PLOT is not None:
        fig.savefig(OUTPUT_NR_PLOT, dpi=PLOT_DPI)
        print(f"n(r) plot saved to {OUTPUT_NR_PLOT}")
    plt.close(fig)


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
    broaden(conventions, meta, RDF_RESOLUTION_SIGMA, r_grid[1] - r_grid[0])
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

    plot_rdfs(gr_results, meta)
    if OUTPUT_NR_PLOT is not None:
        plot_nrs(nr_results)