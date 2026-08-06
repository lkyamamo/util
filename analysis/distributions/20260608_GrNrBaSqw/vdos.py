"""
vdos.py — Vibrational Density of States from LAMMPS dump trajectories.

Aligned with the existing analysis/dynamics/src/msd.cpp (qmd_msd.exe) DOS
calculation: same input parameter names/units, and — by default — the same
windowing and normalization conventions, so the two tools are directly
comparable. See METHOD below for what "aligned" means and what the
alternative (this script's own FFT-periodogram approach) trades off.

QUICK START
-----------
1. Set DUMP_FILE to your LAMMPS custom dump trajectory (element, vx, vy, vz
   columns required). Defaults to dynamics.lammpstrj (env var DYNAMICS_TRAJ)
   — a separate, higher-frequency trajectory from dump.lammpstrj (used by
   rdf_freud.py/bad_freud.py), dumped specifically for dsf.py/vdos.py since
   resolving vibrational frequencies needs much finer time sampling than
   structural analysis does; see OH-therm.input/b-SiO-therm.input's second
   dump block.
2. Set the required environment variables listed under CONFIGURATION below.
3. Run:  python vdos.py

CONFIGURATION
-------------
Each setting comes from exactly one environment variable — no fallback names,
and no silent defaults for anything that changes the result.

REQUIRED (unset or empty is a hard error; the script will not guess):
  DYNAMICS_DT            fs between dumped frames — the ONLY variable vdos.py
                         and msd.py share, since it describes
                         dynamics.lammpstrj rather than either analysis
  VDOS_CORR_LENGTH       fs; VACF max lag / Welch-segment length. Sets the
                         frequency resolution: dnu ~ 1/CORR_LENGTH
  VDOS_CORR_INTERVAL     fs; spacing between VACF reference frames. Sets the
                         noise level only
  VDOS_MAX_FREQUENCY_EV  eV; upper limit of the output DOS grid

Optional — sampling defaults are the identity choice (read the whole
trajectory), the rest are documented algorithm conventions, not physics:
  VDOS_N_FRAMES      max frames to read; 0 = all   (default 0)
  VDOS_STRIDE        read every Nth frame          (default 1)
  VDOS_NUM_GRIDS     frequency grid points         (default 5000)
  VDOS_METHOD        see METHOD below              (default vacf_cosine_transform)
  VDOS_WINDOW        see METHOD below              (default cosine_lag / hann)
  VDOS_NORMALIZATION see NORMALIZATION below       (default phonon)
  VDOS_THREADS       FFT workers, 0 = all cores    (default 0)
  VDOS_PLOT_XUNIT    meV | THz | cm-1 | eV         (default meV / THz)

The VDOS_ prefix is the point: vdos.py and msd.py read the same
dynamics.lammpstrj but want different settings — VDOS needs a short lag for
frequency resolution, MSD a long one to reach the diffusive regime — so a name
that reached both would make one of them silently wrong. DYNAMICS_DT is the
sole shared name: dt is a fact about the trajectory, so one value must reach
both scripts.

CORR_LENGTH previously defaulted to 75% of the trajectory when unset, and
MAX_FREQUENCY_EV to 0.1 eV (806 cm^-1, which silently truncates most spectra).
Neither defaults any more.

METHOD
------
Two interchangeable ways to get VDOS(ν), selected via METHOD:

'vacf_cosine_transform' (default) — replicates qmd_msd.exe's algorithm:
  1. Build the velocity autocorrelation function (VACF) per element using
     multiple time origins: a new reference frame is taken every
     CORR_INTERVAL fs, correlated against every frame up to CORR_LENGTH fs
     later, normalized so C(0)=1 (same as msd.cpp's vac[type][t] /= vac0[type]).
  2. Taper it with a quarter-cosine lag window (WINDOW='cosine_lag'):
     Zt(t) = C(t)*cos(pi*t/(2*CORR_LENGTH)) — same taper as msd.cpp's Zt.
  3. Get DOS via a direct (non-FFT) cosine transform onto an explicit
     frequency grid up to MAX_FREQUENCY_EV with NUM_GRIDS points — same as
     msd.cpp's Gw. Vectorized here as one matrix multiply per element rather
     than msd.cpp's nested loop, but numerically the same calculation.
  4. Combine elements into a total using msd.cpp's phonon-DOS convention
     (NORMALIZATION='phonon'): Total = sum_el (6/pi)*(N_el/N_total)*Gw_el.

'fft_periodogram' — this script's own approach (see prior revisions): skips
  the explicit VACF and computes VDOS directly as the batched FFT periodogram
  |FFT(v(t))|^2 (Wiener-Khinchin), which is O(N log N) instead of msd.cpp's
  O(N^2) direct VACF + direct cosine transform. CORR_LENGTH/CORR_INTERVAL are
  reused here as the Welch-segment length/spacing (in fs) rather than a VACF
  window. Its natural window is a Hann taper on the raw velocity trace
  (WINDOW='hann'), and MAX_FREQUENCY_EV/NUM_GRIDS are ignored — the frequency
  grid is fixed by CORR_LENGTH and TIME_UNIT (FFT bin spacing/Nyquist), not
  freely chosen.

NORMALIZATION applies regardless of METHOD:
  'phonon' (default) — matches msd.cpp: partial curves are left as computed
     (C(0)=1-normalized VACF cosine-transform, or — for fft_periodogram,
     which has no intrinsic physical scale — individually unit-area-rescaled
     first as a stand-in), and the total is the mole-fraction-weighted
     6/pi-scaled sum, targeting integral(g(nu))dnu ~ 3 per atom.
  'unit_area' — this script's original convention: every curve (each partial
     and the total) is independently rescaled so integral(curve)dnu = 1.

OUTPUT
------
- vdos.csv — freq_meV, freq_THz, freq_cm-1, freq_eV, then DoS(<element>) per
             element and DoS(Total) — column naming matches msd.cpp's dos.dat
             (Freq(meV), DoS(<el>), DoS(Total)), with THz/cm-1/eV as extras.
- vdos.png — all curves overlaid on one axes  (set OUTPUT_PLOT=None to skip)

Not replicated from msd.cpp: real-space MSD (a different quantity; use a
dedicated script if needed) and the charge-weighted current-current spectrum
JDoS (requires a per-element CHARGE map specific to msd.cpp's ionic systems,
not part of this generic pipeline).

DEPENDENCIES
------------
  pip install numpy matplotlib
  pip install scipy   # optional: multi-threaded FFT via VDOS_THREADS (fft_periodogram only)

COLUMN LAYOUT
-------------
Expects LAMMPS custom dump format with at minimum columns:
id element vx vy vz (positions are not read/needed). Column positions are
read automatically from the ITEM: ATOMS header line. Atom order must be
consistent frame-to-frame (dump_modify ... sort id in the LAMMPS input,
already used throughout this pipeline).

vx/vy/vz are read as-is and expected to be in Angstrom/ps — LAMMPS "units
metal" velocity units, which is what OH-therm.input/b-SiO-therm.input use.
This doesn't need to match TIME_UNIT's fs: peak positions and each curve's
normalized shape depend only on TIME_UNIT (the real time between frames),
never on the velocity values' own units, since a global rescaling of v
cancels out of the normalization either method uses.
"""

import os
from datetime import date

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# =============================================================================
# CONFIGURATION — edit these variables between runs
# =============================================================================

# Input trajectory file
DUMP_FILE = os.environ.get("DYNAMICS_TRAJ", "dynamics.lammpstrj")

# Every setting below is read from exactly one environment variable, all
# VDOS_-prefixed except DYNAMICS_DT — no fallback names. vdos.py and msd.py
# read the same dynamics.lammpstrj but want different settings (VDOS needs a
# short lag for frequency resolution, MSD a long one to reach the diffusive
# regime), so a name that reached both would make one of them silently wrong.
# Required settings have no default: an unset or empty value is a hard error,
# never a silently substituted number. Missing ones are collected so a single
# run reports every one of them at once instead of failing one at a time.
_MISSING = []

def _require(name, description):
    """Value of `name`, or None after recording it as missing."""
    value = os.environ.get(name, "")
    if value == "":
        _MISSING.append(f"  {name:<22} {description}")
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
    where the default cannot silently distort the result — reading every frame,
    or an algorithm choice documented in this file's header."""
    value = os.environ.get(name, "")
    return value if value != "" else default


# Trajectory sampling (frame-reading only — msd.cpp has no equivalent, it
# always reads every frame of its input file).
N_FRAMES = int(_env("VDOS_N_FRAMES", "0"))   # max frames to read; 0 = all
STRIDE   = int(_env("VDOS_STRIDE", "1"))    # read every Nth frame

# --- Required. DYNAMICS_DT is the one variable vdos.py and msd.py share: it is
# a property of dynamics.lammpstrj itself, not of either analysis, so both
# scripts must read the same value. The rest are VDOS's alone.
_DYNAMICS_DT_ENV = _require("DYNAMICS_DT", "fs between dumped frames")

# Which algorithm to run — see METHOD in the module docstring.
METHOD = _env("VDOS_METHOD", "vacf_cosine_transform")
if METHOD not in ("vacf_cosine_transform", "fft_periodogram"):
    raise ValueError(f"Unknown VDOS_METHOD={METHOD!r}; use 'vacf_cosine_transform' or 'fft_periodogram'.")

# Maximum VACF time lag / Welch-segment length, in fs — same name/units as
# msd.cpp's corr_length_fs argument (its "correlation length"). Sets the
# frequency resolution of the output: dnu ~ 1/CORR_LENGTH.
_CORR_LENGTH_ENV = _require("VDOS_CORR_LENGTH", "fs; VACF max lag / Welch-segment length")

# Spacing between VACF reference frames / Welch-segment starts, in fs — same
# name/units as msd.cpp's corr_interval_fs argument. Sets how many origins are
# averaged, i.e. the noise level, and nothing else.
_CORR_INTERVAL_ENV = _require("VDOS_CORR_INTERVAL", "fs; spacing between VACF reference frames")

# Upper frequency limit of the output DOS grid, in eV — same name/units as
# msd.cpp's max_frequency_ev argument. Only used by METHOD='vacf_cosine_transform';
# fft_periodogram's frequency range is fixed by CORR_LENGTH/TIME_UNIT instead.
_MAX_FREQUENCY_EV_ENV = _require("VDOS_MAX_FREQUENCY_EV", "eV; upper limit of the output DOS grid")
_check_required("vdos.py")

TIME_UNIT = float(_DYNAMICS_DT_ENV)
MAX_FREQUENCY_EV = float(_MAX_FREQUENCY_EV_ENV)

# Number of frequency grid points. Same value as msd.cpp's hardcoded
# num_grids=5000 (not a msd.cpp CLI argument, but its only value in practice).
# Only used by METHOD='vacf_cosine_transform'.
NUM_GRIDS = int(_env("VDOS_NUM_GRIDS", "5000"))

# Window applied before the frequency transform. Valid choices depend on
# METHOD (the two methods window physically different things — a VACF vs. a
# raw velocity trace) — see module docstring. Defaults to whichever matches
# msd.cpp when METHOD='vacf_cosine_transform', or this script's own Hann
# choice under 'fft_periodogram'.
_DEFAULT_WINDOW = {'vacf_cosine_transform': 'cosine_lag', 'fft_periodogram': 'hann'}[METHOD]
WINDOW = _env("VDOS_WINDOW", _DEFAULT_WINDOW)
_VALID_WINDOWS = {'vacf_cosine_transform': ('cosine_lag', 'none'), 'fft_periodogram': ('hann', 'none')}
if WINDOW not in _VALID_WINDOWS[METHOD]:
    raise ValueError(
        f"WINDOW={WINDOW!r} is not valid for METHOD={METHOD!r}; use one of {_VALID_WINDOWS[METHOD]}."
    )

# How partial/total curves are combined and scaled — see module docstring.
NORMALIZATION = _env("VDOS_NORMALIZATION", "phonon")
if NORMALIZATION not in ("phonon", "unit_area"):
    raise ValueError(f"Unknown VDOS_NORMALIZATION={NORMALIZATION!r}; use 'phonon' or 'unit_area'.")

# FFT threading — only used by METHOD='fft_periodogram' and only if scipy is
# installed; 0 = all available cores.
N_THREADS = int(_env("VDOS_THREADS", "0"))

# x-axis unit for vdos.png — 'meV', 'THz', 'cm-1', or 'eV'. The CSV always
# contains all four regardless of this setting. Defaults to meV (msd.cpp's
# own dos.dat convention) under vacf_cosine_transform, THz otherwise.
PLOT_XUNIT = _env("VDOS_PLOT_XUNIT", "meV" if METHOD == "vacf_cosine_transform" else "THz")

# Output files (set to None to skip writing)
OUTPUT_CSV  = "vdos.csv"
OUTPUT_PLOT = "vdos.png"

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

try:
    import scipy.fft as _fft
    _HAVE_SCIPY = True
except ImportError:
    import numpy.fft as _fft
    _HAVE_SCIPY = False

THZ_PER_CM1     = 0.0299792458   # 1 cm^-1 = 0.0299792458 THz
EV_PER_THZ      = 0.0041356677   # 1 THz = h * 1e12 Hz = 4.1356677e-3 eV (E = h*nu)
FS_PER_EV_CYCLE = 1.0 / 4.1356677  # eV -> cycles/fs; matches msd.cpp's freq_unit_ev2fs

FREQ_UNIT_LABELS = {'meV': 'E (meV)', 'THz': 'ν (THz)', 'cm-1': 'ν (cm⁻¹)', 'eV': 'E (eV)'}


def read_lammps_dump(filename, n_frames=0, stride=1):
    """
    Read only what VDOS needs: per-frame element labels and velocities.
    Positions/box are not parsed. Returns (elements, velocities) where
    elements is a (n_atoms,) array (from the first frame) and velocities is
    a (n_kept_frames, n_atoms, 3) array.
    """
    elements = None
    vel_frames = []
    frame_idx = 0

    with open(filename) as f:
        while True:
            line = f.readline()
            if not line:
                break  # EOF

            # TIMESTEP
            f.readline()

            # NUMBER OF ATOMS
            f.readline()
            n_atoms = int(f.readline().strip())

            # BOX BOUNDS (header + 3 dims) — not needed for VDOS
            for _ in range(4):
                f.readline()

            # ATOMS header — parse column positions dynamically
            header = f.readline().split()  # ['ITEM:', 'ATOMS', 'id', 'type', 'element', ...]
            cols = header[2:]
            missing = [c for c in ('element', 'vx', 'vy', 'vz') if c not in cols]
            if missing:
                raise ValueError(
                    f"Dump is missing required column(s) {missing} — VDOS needs "
                    f"velocities (dump ... vx vy vz) in addition to element."
                )
            col_element = cols.index('element')
            col_vx      = cols.index('vx')
            col_vy      = cols.index('vy')
            col_vz      = cols.index('vz')

            keep = (frame_idx % stride == 0)
            if keep:
                frame_elements = []
                velocities = np.empty((n_atoms, 3))
                for i in range(n_atoms):
                    parts = f.readline().split()
                    frame_elements.append(parts[col_element])
                    velocities[i, 0] = float(parts[col_vx])
                    velocities[i, 1] = float(parts[col_vy])
                    velocities[i, 2] = float(parts[col_vz])

                frame_elements = np.array(frame_elements)
                if elements is None:
                    elements = frame_elements
                elif not np.array_equal(elements, frame_elements):
                    raise ValueError(
                        "Atom element order changed between frames — VDOS requires "
                        "consistent atom identity/order across the trajectory "
                        "(check 'dump_modify ... sort id' is set in the LAMMPS input)."
                    )
                vel_frames.append(velocities)
            else:
                for _ in range(n_atoms):
                    f.readline()

            frame_idx += 1
            if n_frames and len(vel_frames) >= n_frames:
                break

    if not vel_frames:
        raise ValueError(f"No frames read from {filename} (empty file or STRIDE too large).")

    return elements, np.stack(vel_frames, axis=0)


# ---------------------------------------------------------------------------
# METHOD='vacf_cosine_transform' — replicates msd.cpp
# ---------------------------------------------------------------------------

def compute_vacf_multi_origin(velocities, elements, corr_length_frames, corr_interval_frames):
    """
    Per-element VACF via multiple time origins, same accumulation as
    msd.cpp's get_msd(): a reference frame every corr_interval_frames,
    correlated against every frame up to corr_length_frames later; C(0)=1
    normalization uses the sum of |v_ref|^2 over ALL references combined
    (not re-averaged per reference), matching msd.cpp exactly.

    Returns (vac, n_refs) where vac is {element: (corr_length_frames,) array}.
    """
    n_frames = velocities.shape[0]
    unique_els = sorted(set(elements.tolist()))
    masks = {el: (elements == el) for el in unique_els}

    vac_sum  = {el: np.zeros(corr_length_frames) for el in unique_els}
    vac0_sum = {el: 0.0 for el in unique_els}

    refs = range(0, n_frames - corr_length_frames + 1, corr_interval_frames)
    n_refs = 0
    for r in refs:
        v_ref   = velocities[r]                                  # (n_atoms, 3)
        segment = velocities[r:r + corr_length_frames]            # (corr_length_frames, n_atoms, 3)
        dot     = np.einsum('tik,ik->ti', segment, v_ref)         # (corr_length_frames, n_atoms)
        for el in unique_els:
            mask = masks[el]
            vac_sum[el]  += dot[:, mask].sum(axis=1)
            vac0_sum[el] += (v_ref[mask] ** 2).sum()
        n_refs += 1

    if n_refs == 0:
        raise ValueError(
            "No valid VACF reference frames: CORR_LENGTH is longer than the trajectory "
            "— shorten CORR_LENGTH or provide more frames."
        )

    vac = {el: (vac_sum[el] / vac0_sum[el] if vac0_sum[el] > 0 else vac_sum[el]) for el in unique_els}
    return vac, n_refs


def apply_cosine_lag_window(vac, corr_length_frames, window):
    """Quarter-cosine taper on the VACF: Zt(t) = C(t)*cos(pi*t/(2*corr_length)),
    same as msd.cpp's Zt. window='none' leaves the VACF untapered."""
    if window == 'none':
        return {el: c.copy() for el, c in vac.items()}
    t = np.arange(corr_length_frames)
    taper = np.cos(0.5 * np.pi * t / corr_length_frames)
    return {el: c * taper for el, c in vac.items()}


def cosine_transform_dos(Zt, time_unit_fs, corr_length_frames, max_frequency_ev, num_grids):
    """
    Direct (non-FFT) cosine transform of the tapered VACF onto an explicit
    frequency grid, vectorized as one matrix multiply per element — same
    calculation as msd.cpp's nested-loop Gw, just not looped in Python.

    Returns (freq_eV, Gw) where Gw is {element: (num_grids,) array}.
    """
    freq_eV      = max_frequency_ev * np.arange(num_grids) / num_grids       # (num_grids,)
    freq_cyc_fs  = freq_eV * FS_PER_EV_CYCLE
    t_fs         = np.arange(corr_length_frames) * time_unit_fs
    cos_matrix   = np.cos(2 * np.pi * np.outer(freq_cyc_fs, t_fs))          # (num_grids, corr_length_frames)

    Gw = {el: (cos_matrix @ z) * time_unit_fs for el, z in Zt.items()}
    return freq_eV, Gw


# ---------------------------------------------------------------------------
# METHOD='fft_periodogram' — this script's own approach
# ---------------------------------------------------------------------------

def _make_window(length, window):
    if window == 'hann':
        w = np.hanning(length)
        return w / np.sqrt(np.mean(w ** 2))   # keep mean power unbiased by windowing
    if window == 'none':
        return np.ones(length)
    raise ValueError(f"Unknown WINDOW={window!r}; use 'hann' or 'none'.")


def _segment_starts(n_frames, seg_len, step_frames):
    """Start indices of non-truncated segments of length seg_len, stepping by
    step_frames. Always returns at least one segment."""
    step = max(1, int(step_frames))
    starts = list(range(0, n_frames - seg_len + 1, step))
    return starts or [0]


def compute_power_spectrum(velocities, dt_fs, window='hann', seg_len_frames=None, step_frames=None):
    """
    Periodogram VDOS: batched real FFT over all atoms and x/y/z components at
    once, no explicit VACF needed (Wiener-Khinchin theorem).

    The trajectory is split into windows of seg_len_frames frames (clipped to
    n_frames if larger; None = whole trajectory, one window) starting every
    step_frames, and the periodograms from each window are averaged (Welch's
    method) to reduce per-frequency noise at the cost of resolution.

    velocities : (n_frames, n_atoms, 3)
    Returns (freq_THz, power, n_segments) where power is (n_atoms, n_freq) —
    squared FFT magnitude summed over x/y/z and averaged over segments, per
    atom, unnormalized.
    """
    n_frames = velocities.shape[0]
    seg_len = n_frames if seg_len_frames is None else min(seg_len_frames, n_frames)
    step = 1 if step_frames is None else max(1, step_frames)
    starts = _segment_starts(n_frames, seg_len, step)

    w = _make_window(seg_len, window)

    fft_kwargs = {}
    if _HAVE_SCIPY:
        fft_kwargs['workers'] = N_THREADS if N_THREADS > 0 else -1

    n_atoms = velocities.shape[1]
    n_freq = seg_len // 2 + 1
    power = np.zeros((n_atoms, n_freq))

    for start in starts:
        segment = velocities[start:start + seg_len] * w[:, None, None]
        spectrum = np.asarray(_fft.rfft(segment, axis=0, **fft_kwargs))   # (n_freq, n_atoms, 3), complex
        power += (np.abs(spectrum) ** 2).sum(axis=2).T   # (n_atoms, n_freq)

    power /= len(starts)

    freq_cycles_per_fs = np.fft.rfftfreq(seg_len, d=dt_fs)
    freq_THz = freq_cycles_per_fs * 1000.0   # 1 cycle/fs = 1000 THz

    return freq_THz, power, len(starts)


# ---------------------------------------------------------------------------
# Shared: normalization/combination, frequency units, output
# ---------------------------------------------------------------------------

def combine_results(raw_curves, elements, freq, normalization, prenormalize_partials):
    """
    Combine per-element raw curves into a results dict (partials + 'total').

    normalization='phonon' (msd.cpp's convention): partials are used as-is
      (or, if prenormalize_partials, first individually rescaled to unit
      area — needed for fft_periodogram, whose raw partials have no
      intrinsic physical scale the way a C(0)=1 VACF does), and the total is
      the mole-fraction-weighted sum sum_el (6/pi)*(N_el/N_total)*curve_el.
    normalization='unit_area': every partial and the total are independently
      rescaled so integral(curve) dnu = 1.

    Returns dict {element_or_'total': (n_freq,) array}.
    """
    counts = {el: int((elements == el).sum()) for el in raw_curves}
    n_total = sum(counts.values())
    d_nu = freq[1] - freq[0]

    def unit_area(curve):
        area = curve.sum() * d_nu
        return curve / area if area > 0 else curve

    if normalization == 'phonon':
        partials = ({el: unit_area(c) for el, c in raw_curves.items()}
                    if prenormalize_partials else dict(raw_curves))
        total = np.zeros_like(freq)
        for el, c in partials.items():
            total += (6.0 / np.pi) * (counts[el] / n_total) * c
    else:  # 'unit_area'
        partials = {el: unit_area(c) for el, c in raw_curves.items()}
        total = unit_area(sum(raw_curves.values()))

    results = dict(partials)
    results['total'] = total
    return results


def _freq_all_units(freq_eV):
    """Expand a native eV frequency grid into all four supported units."""
    freq_THz = freq_eV / EV_PER_THZ
    return {
        'eV':    freq_eV,
        'meV':   freq_eV * 1000.0,
        'THz':   freq_THz,
        'cm-1':  freq_THz / THZ_PER_CM1,
    }


def save_csv(results, freq_by_unit, filename):
    """Save results dict {element_or_'total': array} to CSV, with columns
    named to match msd.cpp's dos.dat (Freq(meV), DoS(<el>), DoS(Total))."""
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


def plot_vdos(results, freq_by_unit, filename, xunit=PLOT_XUNIT):
    if xunit not in FREQ_UNIT_LABELS:
        raise ValueError(f"Unknown PLOT_XUNIT={xunit!r}; use one of {list(FREQ_UNIT_LABELS)}.")
    x = freq_by_unit[xunit]

    fig, ax = plt.subplots(figsize=(8, 5))
    for label, curve in results.items():
        ax.plot(x, curve, label=label, linewidth=1.5 if label == 'total' else 1.0)
    ax.set_xlabel(FREQ_UNIT_LABELS[xunit])
    ax.set_ylabel('DOS (phonon-normalized)' if NORMALIZATION == 'phonon' else 'VDOS (unit-area normalized)')
    ax.legend()
    fig.tight_layout()
    fig.savefig(filename, dpi=PLOT_DPI)
    plt.close(fig)
    print(f"Plot saved to {filename}")


if __name__ == '__main__':
    print(f"Reading trajectory: {DUMP_FILE}")
    print(f"  N_FRAMES={N_FRAMES or 'all'}, STRIDE={STRIDE}, TIME_UNIT={TIME_UNIT} fs")
    print(f"  METHOD={METHOD}, WINDOW={WINDOW}, NORMALIZATION={NORMALIZATION}")

    import time
    t0 = time.time()
    elements, velocities = read_lammps_dump(DUMP_FILE, n_frames=N_FRAMES, stride=STRIDE)
    t1 = time.time()
    n_frames_read = velocities.shape[0]
    print(f"  Frames read: {n_frames_read}, atoms: {velocities.shape[1]} ({t1 - t0:.2f}s)")

    unique_els = sorted(set(elements.tolist()))
    print(f"  Elements: {unique_els}")

    CORR_LENGTH = float(_CORR_LENGTH_ENV)
    CORR_INTERVAL = float(_CORR_INTERVAL_ENV)
    print(f"  CORR_LENGTH={CORR_LENGTH:.1f} fs, CORR_INTERVAL={CORR_INTERVAL:.1f} fs")

    corr_length_frames = min(round(CORR_LENGTH / TIME_UNIT), n_frames_read)
    corr_interval_frames = max(1, round(CORR_INTERVAL / TIME_UNIT))
    if corr_length_frames < round(CORR_LENGTH / TIME_UNIT):
        print(f"  Warning: CORR_LENGTH clipped to the trajectory length "
              f"({corr_length_frames} frames, ~{corr_length_frames * TIME_UNIT:.1f} fs).")

    print(f"  FFT backend: {'scipy.fft (workers=' + str(N_THREADS or -1) + ')' if _HAVE_SCIPY else 'numpy.fft (single-threaded)'}"
          if METHOD == 'fft_periodogram' else "  (fft_periodogram FFT backend not used by this METHOD)")

    print(f"Computing VDOS ({METHOD})...")
    t2 = time.time()

    if METHOD == 'vacf_cosine_transform':
        vac, n_refs = compute_vacf_multi_origin(velocities, elements, corr_length_frames, corr_interval_frames)
        Zt = apply_cosine_lag_window(vac, corr_length_frames, WINDOW)
        freq_eV, Gw = cosine_transform_dos(Zt, TIME_UNIT, corr_length_frames, MAX_FREQUENCY_EV, NUM_GRIDS)
        results = combine_results(Gw, elements, freq_eV, NORMALIZATION, prenormalize_partials=False)
        freq_by_unit = _freq_all_units(freq_eV)
        print(f"  VACF reference frames used: {n_refs}")
    else:  # 'fft_periodogram'
        freq_THz, power, n_segments = compute_power_spectrum(
            velocities, TIME_UNIT, window=WINDOW,
            seg_len_frames=corr_length_frames, step_frames=corr_interval_frames,
        )
        raw_curves = {el: power[elements == el].sum(axis=0) for el in unique_els}
        freq_eV = freq_THz * EV_PER_THZ
        results = combine_results(raw_curves, elements, freq_eV, NORMALIZATION, prenormalize_partials=True)
        freq_by_unit = _freq_all_units(freq_eV)
        print(f"  Segments averaged: {n_segments}")

    t3 = time.time()
    print(f"  Compute: {t3 - t2:.2f}s")
    print(f"  Total: {t3 - t0:.2f}s")

    if OUTPUT_CSV is not None:
        save_csv(results, freq_by_unit, OUTPUT_CSV)
    if OUTPUT_PLOT is not None:
        plot_vdos(results, freq_by_unit, OUTPUT_PLOT)
