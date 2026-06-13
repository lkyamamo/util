"""
rdf_freud.py — Radial Distribution Function and Coordination Number calculator using freud

QUICK START
-----------
1. Set DUMP_FILE to your LAMMPS custom dump trajectory.
2. Verify COL_ELEMENT / COL_X / COL_Y / COL_Z match your dump's ITEM: ATOMS column order.
3. Set R_MAX to less than half the shortest box dimension.
4. Run:  python rdf_freud.py

OUTPUT
------
- rdfs.png  — g(r) subplot grid, one panel per element pair + total + neutron-weighted
- rdfs.csv  — r (Å), g(r) per pair, total, neutron-weighted  (set OUTPUT_CSV=None to skip)
- nrs.png   — cumulative coordination number n(r) plots       (set OUTPUT_NR_PLOT=None to skip)
- nrs.csv   — r (Å), n(r) per pair                            (set OUTPUT_NR_CSV=None to skip)

COLUMN LAYOUT
-------------
Adjust COL_ELEMENT / COL_X / COL_Y / COL_Z if your dump's ITEM: ATOMS columns differ
from the default layout:  id  element  x  y  z  vx  vy  vz
                          0   1        2  3  4  5   6   7
"""

import itertools
import os

import numpy as np
import freud
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
R_MAX = 20.0        # maximum r in Angstroms; must be < half shortest box dimension
BINS  = 2000         # number of bins

# Column indices in the ITEM: ATOMS line (0-indexed)
# Default matches: id element x y z vx vy vz
COL_ELEMENT = 1
COL_X       = 2
COL_Y       = 3
COL_Z       = 4

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

# Coherent neutron scattering lengths (fm).  Add elements as needed.
# Values from NIST: https://www.ncnr.nist.gov/resources/n-lengths/
NEUTRON_SCATTERING_LENGTHS = {
    'H':   6.671,    # deuterium (D); protium b = -3.7406
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

            # BOX BOUNDS
            f.readline()
            xlo, xhi = map(float, f.readline().split())
            ylo, yhi = map(float, f.readline().split())
            zlo, zhi = map(float, f.readline().split())

            # ATOMS header line
            f.readline()

            elements, positions = [], []
            for _ in range(n_atoms):
                parts = f.readline().split()
                elements.append(parts[COL_ELEMENT])
                positions.append([
                    float(parts[COL_X]),
                    float(parts[COL_Y]),
                    float(parts[COL_Z]),
                ])

            box = freud.box.Box(
                Lx=xhi - xlo,
                Ly=yhi - ylo,
                Lz=zhi - zlo,
            )
            positions = np.array(positions)

            # wrap positions into [-L/2, L/2] as freud expects
            center = np.array([(xlo + xhi) / 2, (ylo + yhi) / 2, (zlo + zhi) / 2])
            positions -= center

            frames.append({
                'timestep':       timestep,
                'box':            box,
                'positions':      positions,
                'elements':       np.array(elements),
                'number_density': n_atoms / box.volume,   # atoms / Å³
            })

    return frames


def compute_rdf(frames, get_a, get_b, self_pair=False):
    """
    Compute the frame-averaged RDF and n(r) using freud's built-in accumulation.
    freud accumulates across compute() calls when reset=False, so rdf.rdf and
    rdf.n_r at the end are already the properly normalized frame averages.

    Returns: r, mean_g, mean_nr
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
        return rdf.bin_centers, empty, empty

    return rdf.bin_centers, rdf.rdf, rdf.n_r


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



def _neutron_weights(elements, concentrations):
    """Return (b, b_mean, w_neutron dict) or None if any element is missing."""
    missing = [el for el in elements if el not in NEUTRON_SCATTERING_LENGTHS]
    if missing:
        print(f"Warning: no scattering length for {missing}; neutron g(r) skipped.")
        return None
    b = {el: NEUTRON_SCATTERING_LENGTHS[el] for el in elements}
    b_mean = sum(concentrations[el] * b[el] for el in elements)
    weights = {}
    for a, bl in itertools.combinations_with_replacement(sorted(elements), 2):
        factor = 1 if a == bl else 2
        w_total = factor * concentrations[a] * concentrations[bl]
        weights[f'{a}-{bl}'] = w_total * b[a] * b[bl] / b_mean ** 2
    return weights


def combine_partials(partial_results, elements, concentrations, rho):
    """
    Combine partial g(r)s into total and neutron-weighted g(r) and t(r).

    Total g(r):   w_αβ = c_α c_β  (×2 for cross-pairs)
    Neutron g(r): w_αβ = c_α c_β b_α b_β / <b>²  (×2 for cross-pairs)
    t(r):         g_neutron(r) · 4π r · mean(ρ)

    rho : per-frame number density array (atoms / Å³), shape (n_frames,)

    Returns a dict {label: (r, y)} with keys 'total', and optionally 'neutron' and 't'.
    """
    r = next(iter(partial_results.values()))[0]

    g_total   = np.zeros(len(r))
    g_neutron = np.zeros(len(r))

    nw = _neutron_weights(elements, concentrations)

    for a, bl in itertools.combinations_with_replacement(sorted(elements), 2):
        label  = f'{a}-{bl}'
        factor = 1 if a == bl else 2
        _, g, *_ = partial_results[label]

        g_total += factor * concentrations[a] * concentrations[bl] * g
        if nw is not None:
            g_neutron += nw[label] * g

    out = {'total': (r, g_total)}
    if nw is not None:
        t = g_neutron * 4 * np.pi * r * rho.mean()
        out['neutron'] = (r, g_neutron)
        out['t']       = (r, t)
    return out


def save_csv(results, filename):
    """Save results dict {label: (r, y)} to CSV with one column per label."""
    r = next(iter(results.values()))[0]
    header = 'r_Angstrom,' + ','.join(results.keys())
    data = np.column_stack([r] + [y for _, y in results.values()])
    np.savetxt(filename, data, delimiter=',', header=header, comments='', fmt='%.6f')
    print(f"Data table saved to {filename}")


def plot_rdfs(results):
    n = len(results)
    ncols = PLOT_NCOLS
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(6 * ncols, 4 * nrows), squeeze=False)
    axes = axes.flatten()

    for ax, (name, (r, g)) in zip(axes, results.items()):
        ax.plot(r, g)
        ax.axhline(1.0, color='gray', linestyle='--', linewidth=0.8)
        ax.set_xlabel('r (Å)')
        ax.set_ylabel('g(r)')
        ax.set_title(name)

    # hide any unused subplots
    for ax in axes[n:]:
        ax.set_visible(False)

    fig.tight_layout()
    if OUTPUT_PLOT is not None:
        fig.savefig(OUTPUT_PLOT, dpi=PLOT_DPI)
        print(f"Plot saved to {OUTPUT_PLOT}")
    plt.show()


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
    plt.show()


if __name__ == '__main__':
    print(f"Reading trajectory: {DUMP_FILE}")
    frames = read_lammps_dump(DUMP_FILE)

    atom_counts = [len(f['positions']) for f in frames]
    print(f"Frames read:  {len(frames)}")
    print(f"Atoms range:  {min(atom_counts)} – {max(atom_counts)}")
    print(f"Atoms mean:   {np.mean(atom_counts):.1f}")

    elements = sorted(set(frames[0]['elements'].tolist()))
    print(f"Elements:     {elements}")

    concentrations = get_concentrations(frames)
    print("Concentrations: " + ", ".join(f"{el}={c:.3f}" for el, c in concentrations.items()))

    rho = np.array([f['number_density'] for f in frames])   # atoms/Å³, shape (n_frames,)
    print(f"Number density: {rho.mean():.6e} atoms/Å³")

    pairs = build_pair_getters(elements)
    gr_results = {}
    nr_results = {}
    for name, (get_a, get_b, is_self) in pairs.items():
        print(f"Computing RDF: {name}...")
        r, g, nr = compute_rdf(frames, get_a, get_b, self_pair=is_self)
        gr_results[name] = (r, g)
        nr_results[name] = (r, nr)

    gr_results.update(combine_partials(gr_results, elements, concentrations, rho))

    if OUTPUT_CSV is not None:
        save_csv(gr_results, OUTPUT_CSV)
    if OUTPUT_NR_CSV is not None:
        save_csv(nr_results, OUTPUT_NR_CSV)

    plot_rdfs(gr_results)
    if OUTPUT_NR_PLOT is not None:
        plot_nrs(nr_results)