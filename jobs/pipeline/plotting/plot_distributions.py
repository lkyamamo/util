"""
plot_distributions.py — turn an analysis directory's CSVs into figures.

This is the whole plotting stage. The analysis scripts in
analysis/distributions/20260608_GrNrBaSqw/ write CSVs and nothing else; every
figure in the pipeline is produced here, from those CSVs. Nothing is recomputed
— if a curve is not in a CSV it cannot be drawn, which is deliberate: a figure
that disagrees with the table it came from is the failure mode this structure
removes.

QUICK START
-----------
Do not run this directly; run plot_pipeline.sh, which resolves the config and
sets the environment this reads:

    cd <analysis directory>
    /path/to/jobs/pipeline/plotting/plot_pipeline.sh

OUTPUT
------
One PNG per quantity, in <date>_<calc>/ next to the CSV it came from:

    <date>_rdf/           each partial, each convention column and its
                          _broadened twin, each n(r), the wright curves
    <date>_bad/           one per A-B-C triplet
    <date>_dsf/           S(q) partials/total/neutron, the S(q,w) heatmaps
    <date>_vdos/          one per species and per weighted total
    <date>_msd/           one per species
    <date>_vdos_dynmat/   DOS curves, stretch/bend/rock, participation ratio,
                          reduced DOS

THE DATE COMES FROM THE CSV, NOT FROM TODAY. Plotting 20260813_bads.csv writes
20260813_bad/, so figures stay married to the data that produced them and
re-plotting an old run does not stamp it with the current date.

NO REFERENCE LINES are drawn. An individual plot carries its curve and nothing
else — no g(r) = 1 line, no baseline under T(r). A composite is written only
where the combination is itself the result (the Wright overlay, the mode
character decomposition, the species overlays), and is named so it reads as one.
"""

import glob
import os
import re
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# =============================================================================
# Configuration — read from the environment, which plot_pipeline.sh fills in
# from plot_pipeline.conf. Defaults here match the shipped conf, so this module
# is still runnable if something calls it directly.
# =============================================================================

def _env(name, default):
    v = os.environ.get(name)
    return default if v is None or v == "" else v


def _flag(name, default="1"):
    v = _env(name, default).strip().lower()
    if v not in ("0", "1", "yes", "no", "true", "false"):
        raise SystemExit(f"plot_distributions.py: {name}={v!r} must be 0 or 1.")
    return v in ("1", "yes", "true")


TARGET_DIR   = _env("PLOT_TARGET_DIR", ".")
PLOT_STYLE   = _env("PLOT_STYLE", "analysis")
PLOT_DPI     = int(_env("PLOT_DPI", "150"))
PLOT_FORMAT  = _env("PLOT_FORMAT", "png")

# Which run to plot when a directory holds several dates: 'latest' (default),
# 'all', or an explicit YYYYMMDD.
PLOT_DATE    = _env("PLOT_DATE", "latest")

# Per-calculation switches.
PLOT_RDF         = _flag("PLOT_RDF")
PLOT_BAD         = _flag("PLOT_BAD")
PLOT_DSF         = _flag("PLOT_DSF")
PLOT_VDOS        = _flag("PLOT_VDOS")
PLOT_MSD         = _flag("PLOT_MSD")
PLOT_VDOS_DYNMAT = _flag("PLOT_VDOS_DYNMAT")

# Write the composite figures alongside the individual curves.
PLOT_COMPOSITES  = _flag("PLOT_COMPOSITES")

# x axis for the frequency plots: meV | THz | cm-1 | eV. The CSVs carry all
# four, so this only chooses which column becomes the axis.
PLOT_FREQ_UNIT   = _env("PLOT_FREQ_UNIT", "meV")

# Colormap for the S(q,w) heatmaps.
PLOT_HEATMAP_CMAP = _env("PLOT_HEATMAP_CMAP", "inferno")


PLOT_STYLES = {
    # 'analysis'    — titled and fully labelled, for reading a run.
    # 'publication' — heavy lines and spines, large bold labels, no y ticks.
    #                 This is what the old rdf_freud_plot.py / bad_freud_plot.py
    #                 second pass produced; it is a config choice now, and it
    #                 applies to all six calculations rather than two.
    'analysis':    dict(figsize=(7.0, 4.5), linewidth=1.5, color='C0', spine_lw=0.8,
                        weight='normal', label_fs=12, tick_fs=10,
                        tick_len=4, tick_w=1.0, yticks=True, titles=True),
    'publication': dict(figsize=(4.0, 3.0), linewidth=3.0, color='steelblue', spine_lw=2.0,
                        weight='bold', label_fs=20, tick_fs=14,
                        tick_len=6, tick_w=2.0, yticks=False, titles=False),
}

FREQ_COLUMNS = {'meV': 'freq_meV', 'THz': 'freq_THz',
                'cm-1': 'freq_cm-1', 'eV': 'freq_eV'}
FREQ_LABELS  = {'meV': 'Energy (meV)', 'THz': 'Frequency (THz)',
                'cm-1': 'Wavenumber (cm⁻¹)', 'eV': 'Energy (eV)'}

# rdfs.csv mixes units: the partials are dimensionless g_AB(r) about 1, while
# the combined columns are named <function>_<normalization> and carry their own
# units. See CONVENTIONS in rdf_freud.py for the defining equations.
NORMALIZATION_UNITS = {'unity': '', 'FZ': '', 'absolute': 'barn/sr/atom',
                       'formula': 'barn/sr/formula-unit'}
FUNCTION_KEYS = ('g', 'h', 'D', 'T')

# The three orthogonal directions a bridging atom can move in, as vdos_dynmat.py
# names them. Their DOS columns are prefixed on output so they group together and
# cannot be mistaken for an element's partial DOS sitting in the same directory.
CHARACTER_NAMES = ('stretch', 'bend', 'rock')


# =============================================================================
# Plot primitives
# =============================================================================

def style():
    if PLOT_STYLE not in PLOT_STYLES:
        raise SystemExit(
            f"plot_distributions.py: PLOT_STYLE={PLOT_STYLE!r}; "
            f"use one of {list(PLOT_STYLES)}.")
    return PLOT_STYLES[PLOT_STYLE]


def new_plot():
    """A single-axes figure in the configured style."""
    st = style()
    return plt.subplots(figsize=st['figsize']) + (st,)


def save_plot(fig, ax, out_dir, name, xlabel, ylabel, title=None, legend=False):
    """
    Finish one figure and write it as out_dir/<name>.<format>.

    `name` becomes the filename, so anything a path cannot carry is substituted
    rather than left to mangle the path silently.
    """
    st = style()
    ax.set_xlabel(xlabel, fontsize=st['label_fs'], fontweight=st['weight'])
    ax.set_ylabel(ylabel, fontsize=st['label_fs'], fontweight=st['weight'])
    if title and st['titles']:
        ax.set_title(title)
    ax.tick_params(axis='x', labelsize=st['tick_fs'],
                   length=st['tick_len'], width=st['tick_w'])
    if st['yticks']:
        ax.tick_params(axis='y', labelsize=st['tick_fs'],
                       length=st['tick_len'], width=st['tick_w'])
    else:
        ax.yaxis.set_ticks([])
    if legend:
        ax.legend(fontsize=8 if st['titles'] else 10)
    for spine in ax.spines.values():
        spine.set_linewidth(st['spine_lw'])
    fig.tight_layout()
    safe = re.sub(r'[^A-Za-z0-9._+-]', '_', name)
    path = os.path.join(out_dir, f"{safe}.{PLOT_FORMAT}")
    fig.savefig(path, dpi=PLOT_DPI)
    plt.close(fig)
    return path


def curve(out_dir, x, y, name, xlabel, ylabel, title=None, xlim=None):
    """One quantity, one file, nothing else on the axes."""
    fig, ax, st = new_plot()
    ax.plot(x, y, color=st['color'], linewidth=st['linewidth'])
    if xlim is not None:
        ax.set_xlim(*xlim)
    return save_plot(fig, ax, out_dir, name, xlabel, ylabel, title=title)


def overlay(out_dir, x, series, name, xlabel, ylabel, title=None, emphasize=()):
    """
    A composite: several curves on shared axes, with a legend.

    Written only where the comparison between the curves is itself the result —
    a partial DOS read against the total it sums into, one species outrunning
    another. `emphasize` names the curves drawn heavier.
    """
    fig, ax, st = new_plot()
    for label, y in series.items():
        ax.plot(x, y, label=label,
                linewidth=1.5 if label in emphasize else 1.0)
    return save_plot(fig, ax, out_dir, name, xlabel, ylabel,
                     title=title, legend=True)


# =============================================================================
# CSV handling
# =============================================================================

def read_csv(path):
    """
    Read one of the analysis CSVs.

    Returns (labels, columns) with labels[0]/columns[0] the x axis. np.loadtxt
    rather than a csv reader because these are dense numeric tables and some of
    them (dsf.csv especially) are large.
    """
    with open(path) as f:
        labels = f.readline().strip().split(',')
    data = np.loadtxt(path, delimiter=',', skiprows=1, ndmin=2)
    return labels, [data[:, i] for i in range(data.shape[1])]


def rdf_column_label(name):
    """
    (ylabel, title) for one rdfs.csv column.

    A bare pair like 'O-Si' is a partial; anything matching
    <function>_<normalization>[_broadened] is a combined column that carries its
    own units, which is exactly why these cannot share an axis.
    """
    parts = name.split('_')
    if len(parts) < 2 or parts[0] not in FUNCTION_KEYS or parts[1] not in NORMALIZATION_UNITS:
        return 'g(r)', name                       # a partial, e.g. 'O-Si'
    units = NORMALIZATION_UNITS[parts[1]]
    if parts[0] in ('D', 'T'):
        units = f'{units} Å⁻²'.strip() if units else 'Å⁻²'
    return (f'{name} ({units})' if units else name), name


def discover(target_dir):
    """
    Group the CSVs in target_dir by their YYYYMMDD_ prefix.

    Returns {date: {kind: path}}. The prefix is what every analysis script
    stamps on its output via _dated(), so it is the only thing that reliably
    ties one run's files together.
    """
    kinds = {
        'rdfs': 'rdfs', 'nrs': 'nrs', 'wright': 'wright', 'bads': 'bads',
        'sq': 'sq', 'dsf': 'dsf', 'vdos': 'vdos', 'msd': 'msd',
    }
    runs = {}
    for path in sorted(glob.glob(os.path.join(target_dir, '[0-9]' * 8 + '_*.csv'))):
        base = os.path.basename(path)
        date, _, rest = base.partition('_')
        stem = rest[:-4]                                   # drop '.csv'
        # vdos_dynmat's outputs are named from VDOS_DYNMAT_OUTPUT, and its
        # modes table shares that stem, so match those two before the fixed set
        # (otherwise 'vdos_dynmat' would never be distinguished from 'vdos').
        if stem.endswith('_modes'):
            kind = 'dynmat_modes'
        elif stem in kinds:
            kind = kinds[stem]
        elif stem == 'diffusion':
            # msd.py's table of D values: a summary for the temperature-sweep
            # aggregator, not a curve. Without this it falls through to 'dynmat'.
            continue
        elif stem not in ('vdos', 'msd'):
            kind = 'dynmat'                                # vdos_dynmat.csv
        else:
            continue
        runs.setdefault(date, {})[kind] = path
    return runs


def select_dates(runs):
    """Which run(s) to plot, per PLOT_DATE."""
    dates = sorted(runs)
    if not dates:
        return []
    if PLOT_DATE == 'latest':
        return [dates[-1]]
    if PLOT_DATE == 'all':
        return dates
    if PLOT_DATE not in runs:
        raise SystemExit(
            f"plot_distributions.py: PLOT_DATE={PLOT_DATE} matches no CSVs in "
            f"{TARGET_DIR}. Present: {', '.join(dates) or '(none)'}")
    return [PLOT_DATE]


def out_dir_for(target_dir, date, calc):
    """<date>_<calc>/ — created here, since this script owns its own output."""
    path = os.path.join(target_dir, f'{date}_{calc}')
    os.makedirs(path, exist_ok=True)
    return path


# =============================================================================
# Per-calculation plotting
# =============================================================================

def plot_rdf(target_dir, date, files):
    out = out_dir_for(target_dir, date, 'rdf')
    n = 0

    if 'rdfs' in files:
        labels, cols = read_csv(files['rdfs'])
        r = cols[0]
        for name, y in zip(labels[1:], cols[1:]):
            ylabel, title = rdf_column_label(name)
            curve(out, r, y, name, 'r (Å)', ylabel, title=title)
            n += 1

    if 'nrs' in files:
        labels, cols = read_csv(files['nrs'])
        r = cols[0]
        for name, y in zip(labels[1:], cols[1:]):
            # prefixed so an n(r) cannot be confused with the g(r) of the same
            # pair sitting in the same directory
            curve(out, r, y, f'nr_{name}', 'r (Å)', 'n(r)', title=f'n(r)  {name}')
            n += 1

    if 'wright' in files:
        labels, cols = read_csv(files['wright'])
        r = cols[0]
        by = dict(zip(labels[1:], cols[1:]))
        ylabel = r'T(r)  (barn sr$^{-1}$ formula-unit$^{-1}$ Å$^{-2}$)'
        names = {'T_wright': 'wright',
                 'T_wright_unbroadened': 'wright_unbroadened',
                 'T0_baseline': 'wright_T0_baseline'}
        for col, name in names.items():
            if col in by:
                curve(out, r, by[col], name, 'r (Å)', ylabel, title=name)
                n += 1
        if PLOT_COMPOSITES and 'T_wright' in by:
            # The figure a published T(r) is actually overlaid on: seeing the
            # broadened curve, the raw curve and the baseline together is the
            # whole point of the comparison.
            fig, ax, st = new_plot()
            ax.plot(r, by['T_wright'], color='C0', linewidth=1.6,
                    label='T(r), Lorch-broadened')
            if 'T_wright_unbroadened' in by:
                ax.plot(r, by['T_wright_unbroadened'], color='C0',
                        linewidth=0.8, alpha=0.45, label='T(r), unbroadened')
            if 'T0_baseline' in by:
                ax.plot(r, by['T0_baseline'], color='0.4', linestyle='--',
                        linewidth=1.0,
                        label=r'$T^0(r)=4\pi r\rho^0\langle b\rangle^2$')
            save_plot(fig, ax, out, 'wright_composite', 'r (Å)', ylabel,
                      title='Wright comparison', legend=True)
            n += 1
    return out, n


def plot_bad(target_dir, date, files):
    out = out_dir_for(target_dir, date, 'bad')
    labels, cols = read_csv(files['bads'])
    angles = cols[0]
    for name, y in zip(labels[1:], cols[1:]):
        # 0-180 is kept on every one: a bond angle distribution read against an
        # auto-scaled axis invites comparing two triplets that do not share a
        # range.
        curve(out, angles, y, name, 'Angle (degrees)', 'P(θ)',
              title=name, xlim=(0, 180))
    return out, len(labels) - 1


def plot_dsf(target_dir, date, files):
    out = out_dir_for(target_dir, date, 'dsf')
    n = 0

    if 'sq' in files:
        labels, cols = read_csv(files['sq'])
        q = cols[0]
        series = {}
        for name, y in zip(labels[1:], cols[1:]):
            # dynasor names the partials Sq_A_B; the pair reads as A-B, and that
            # is what the rdf side calls the same pair.
            short = name[3:] if name.startswith('Sq_') else name
            if name.startswith('Sq_') and short.count('_') == 1:
                short = short.replace('_', '-')
            curve(out, q, y, f'sq_{short}', 'q (Å⁻¹)', 'S(q)', title=f'S(q)  {short}')
            series[short] = y
            n += 1
        if PLOT_COMPOSITES and series:
            overlay(out, q, series, 'sq_all_curves', 'q (Å⁻¹)', 'S(q)',
                    title='S(q), all curves', emphasize=('total', 'neutron'))
            n += 1

    if 'dsf' in files:
        # dsf.csv is long-format: one row per (q, omega) pair, q varying
        # slowest. Reshaping back to the (n_q, n_omega) grid is what makes a
        # heatmap possible without recomputing anything.
        labels, cols = read_csv(files['dsf'])
        q_flat, w_flat = cols[0], cols[1]
        q = np.unique(q_flat)
        w = np.unique(w_flat)
        if len(q) * len(w) != len(q_flat):
            raise SystemExit(
                f"plot_distributions.py: {os.path.basename(files['dsf'])} is not a "
                f"complete (q, omega) grid: {len(q)} x {len(w)} != {len(q_flat)} rows. "
                f"It should be one row per pair, q varying slowest — see "
                f"save_csv_dsf() in dsf.py.")
        if len(w) < 2:
            # One frequency bin is not a spectrum. dynasor truncates the ACF when
            # the window outruns the trajectory, and this is what that looks like
            # downstream, so say so rather than drawing a one-pixel-wide heatmap.
            print(f"  note: {os.path.basename(files['dsf'])} has a single omega value "
                  f"— no S(q,omega) heatmap to draw. Check DSF_WINDOW_SIZE against "
                  f"the trajectory length.")
            return out, n
        st = style()
        for name, y in zip(labels[2:], cols[2:]):
            grid = y.reshape(len(q), len(w))
            # A heatmap carries a colorbar rather than a y scale, so it does not
            # go through save_plot(): stripping its ticks for the publication
            # style would make it unreadable.
            fig, ax = plt.subplots(figsize=(7, 5))
            im = ax.imshow(grid, origin='lower', aspect='auto',
                           extent=[w[0], w[-1], q[0], q[-1]],
                           vmin=0, vmax=np.nanpercentile(grid, 99),
                           cmap=PLOT_HEATMAP_CMAP)
            ax.set_xlabel('ω (THz)', fontsize=st['label_fs'], fontweight=st['weight'])
            ax.set_ylabel('q (Å⁻¹)', fontsize=st['label_fs'], fontweight=st['weight'])
            if st['titles']:
                ax.set_title(name)
            fig.colorbar(im, ax=ax, label=name)
            fig.tight_layout()
            short = name[4:] if name.startswith('Sqw_') else name
            if name.startswith('Sqw_') and short.count('_') == 1:
                short = short.replace('_', '-')          # Sqw_O_Si -> O-Si
            safe = re.sub(r'[^A-Za-z0-9._+-]', '_', f'sqw_{short}')
            fig.savefig(os.path.join(out, f'{safe}.{PLOT_FORMAT}'), dpi=PLOT_DPI)
            plt.close(fig)
            n += 1
    return out, n


def _dos_columns(labels, cols):
    """Split a vdos-style CSV into (x, xlabel, {name: curve}) for the chosen unit."""
    if PLOT_FREQ_UNIT not in FREQ_COLUMNS:
        raise SystemExit(f"plot_distributions.py: PLOT_FREQ_UNIT={PLOT_FREQ_UNIT!r}; "
                         f"use one of {list(FREQ_COLUMNS)}.")
    xcol = FREQ_COLUMNS[PLOT_FREQ_UNIT]
    if xcol not in labels:
        raise SystemExit(f"plot_distributions.py: no {xcol} column; found {labels[:4]}")
    x = cols[labels.index(xcol)]
    series = {name: y for name, y in zip(labels, cols)
              if not name.startswith('freq_')}
    return x, FREQ_LABELS[PLOT_FREQ_UNIT], series


def plot_vdos(target_dir, date, files, calc='vdos'):
    out = out_dir_for(target_dir, date, calc)
    labels, cols = read_csv(files['vdos' if calc == 'vdos' else 'dynmat'])
    x, xlabel, series = _dos_columns(labels, cols)

    dos, extra = {}, {}
    for name, y in series.items():
        # DoS(...) columns are the density of states proper; g/nu^2 is the
        # reduced DOS, a different quantity that must not join their overlay.
        (extra if not name.startswith('DoS(') else dos)[name] = y

    n = 0
    for name, y in dos.items():
        short = name[4:-1] if name.endswith(')') else name
        label = f'character_{short}' if short in CHARACTER_NAMES else short
        curve(out, x, y, label, xlabel, 'DOS', title=f'VDOS  {short}')
        n += 1
    if PLOT_COMPOSITES and dos:
        plain = {k[4:-1]: v for k, v in dos.items()}
        overlay(out, x, plain, 'all_curves', xlabel, 'DOS',
                title='VDOS, all curves',
                emphasize=[k for k in plain if k.lower().startswith('total')])
        n += 1

    for name, y in extra.items():
        if name == 'g/nu^2':
            curve(out, x, y, 'reduced_dos', xlabel, f'g / {PLOT_FREQ_UNIT}²',
                  title='Reduced DOS g(ν)/ν² — a peak here is the boson peak')
        else:
            curve(out, x, y, name, xlabel, name, title=name)
        n += 1
    return out, n


def plot_dynmat_modes(target_dir, date, files, out):
    """
    The per-mode table: one point per mode, not a curve.

    The scatter is the right form because the spread of participation ratio at
    a given frequency is itself the information — a tight low band means every
    mode there is equally extended, a wide one means they are not.
    """
    labels, cols = read_csv(files['dynmat_modes'])
    x, xlabel, _ = _dos_columns(labels, cols)
    by = dict(zip(labels, cols))
    if 'participation_ratio' not in by:
        return 0
    st = style()
    fig, ax, _ = new_plot()
    ax.plot(x, by['participation_ratio'], '.', markersize=2, alpha=0.4)
    ax.set_ylim(0, 1)
    save_plot(fig, ax, out, 'participation_ratio', xlabel, 'participation ratio',
              title='Localization (1 = every atom moves, 1/N = one atom moves)')
    return 1


def plot_dynmat(target_dir, date, files):
    out, n = plot_vdos(target_dir, date, files, calc='vdos_dynmat')

    # The character decomposition lives in the same CSV as DoS(stretch/bend/
    # rock); it has already been drawn per curve above, so only the composite
    # is left. A band assignment is made by reading the three against each
    # other and against the total, which is why it earns a composite.
    if PLOT_COMPOSITES:
        labels, cols = read_csv(files['dynmat'])
        x, xlabel, series = _dos_columns(labels, cols)
        wanted = ['DoS(Total_unity)', 'DoS(stretch)', 'DoS(bend)', 'DoS(rock)']
        present = {k[4:-1]: series[k] for k in wanted if k in series}
        if len(present) > 1:
            overlay(out, x, present, 'character_composite', xlabel, 'DOS',
                    title='Character-resolved DOS (bridging-atom motion)',
                    emphasize=('Total_unity',))
            n += 1

    if 'dynmat_modes' in files:
        n += plot_dynmat_modes(target_dir, date, files, out)
    return out, n


def plot_msd(target_dir, date, files):
    out = out_dir_for(target_dir, date, 'msd')
    labels, cols = read_csv(files['msd'])
    t = cols[0]
    series = {}
    for name, y in zip(labels[1:], cols[1:]):
        short = name[4:] if name.startswith('MSD_') else name
        curve(out, t, y, short, 't (fs)', 'MSD (Å²)', title=f'MSD  {short}')
        series[short] = y
    n = len(series)
    if PLOT_COMPOSITES and series:
        # The comparison between species IS the result here — a light element
        # outrunning a heavy one is the thing being read.
        overlay(out, t, series, 'all_species', 't (fs)', 'MSD (Å²)',
                title='MSD, all species', emphasize=('total',))
        n += 1
    return out, n


# =============================================================================
# Main
# =============================================================================

CALCULATIONS = [
    # (name, enabled, the CSV kind that must be present, plotting function)
    ('rdf',         PLOT_RDF,         ('rdfs', 'nrs', 'wright'), plot_rdf),
    ('bad',         PLOT_BAD,         ('bads',),                 plot_bad),
    ('dsf',         PLOT_DSF,         ('sq', 'dsf'),             plot_dsf),
    ('vdos',        PLOT_VDOS,        ('vdos',),                 plot_vdos),
    ('msd',         PLOT_MSD,         ('msd',),                  plot_msd),
    ('vdos_dynmat', PLOT_VDOS_DYNMAT, ('dynmat',),               plot_dynmat),
]


if __name__ == '__main__':
    target = os.path.abspath(TARGET_DIR)
    if not os.path.isdir(target):
        raise SystemExit(f"plot_distributions.py: no such directory: {target}")

    runs = discover(target)
    dates = select_dates(runs)
    if not dates:
        raise SystemExit(
            f"plot_distributions.py: no analysis CSVs in {target}.\n"
            f"  Expected files named <YYYYMMDD>_<name>.csv, written by the "
            f"scripts in analysis/distributions/20260608_GrNrBaSqw/.\n"
            f"  Run the analysis first, or point at the directory that holds it.")

    print(f"Plotting from: {target}")
    print(f"  style={PLOT_STYLE}  freq-unit={PLOT_FREQ_UNIT}  composites="
          f"{'on' if PLOT_COMPOSITES else 'off'}  dpi={PLOT_DPI}")
    print(f"  runs found: {', '.join(sorted(runs))}"
          f"   plotting: {', '.join(dates)}")

    total, skipped = 0, []
    for date in dates:
        files = runs[date]
        for name, enabled, needs, fn in CALCULATIONS:
            present = [k for k in needs if k in files]
            if not present:
                continue
            if not enabled:
                skipped.append(f"{name} ({date}): disabled in config")
                continue
            out, n = fn(target, date, files)
            print(f"  {os.path.basename(out)}/  {n} plot(s)"
                  f"   <- {', '.join(os.path.basename(files[k]) for k in present)}")
            total += n

    for line in skipped:
        print(f"  skipped: {line}")
    if total == 0:
        print("  nothing plotted — every calculation present was disabled in the config")
    else:
        print(f"Done. {total} plot(s) written.")
