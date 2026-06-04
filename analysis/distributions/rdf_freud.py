import numpy as np
import freud
import matplotlib.pyplot as plt

# =============================================================================
# CONFIGURATION — edit these variables between runs
# =============================================================================

# Input trajectory file
DUMP_FILE = "dump.lammpstrj"

# Output plot file
OUTPUT_PLOT = "rdfs.png"

# Output data table (CSV); set to None to skip
OUTPUT_CSV = "rdfs.csv"

# RDF parameters
R_MAX = 10.0        # maximum r in Angstroms; must be < half shortest box dimension
BINS  = 200         # number of bins

# Oxygen classification toggle
# If True: compute additional RDFs for bridging O, silanol O, and water O
# If False: only compute bulk element-pair RDFs (faster)
CLASSIFY_OXYGENS = True

# Cutoffs used for oxygen classification (Angstroms)
# Tune these to the first minimum of the Si-O and O-H RDFs respectively
SI_O_CUTOFF = 2.0
O_H_CUTOFF  = 1.2

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

# =============================================================================
# END CONFIGURATION
# =============================================================================


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
                'timestep': timestep,
                'box':      box,
                'positions': positions,
                'elements': np.array(elements),
            })

    return frames


def classify_oxygens(box, si_pos, o_pos, h_pos):
    n_o = len(o_pos)
    si_count = np.zeros(n_o, dtype=int)
    h_count  = np.zeros(n_o, dtype=int)

    if len(si_pos) > 0:
        si_query = freud.locality.AABBQuery(box, si_pos)
        si_nlist = si_query.query(
            o_pos, {'r_max': SI_O_CUTOFF, 'exclude_ii': False}
        ).toNeighborList()
        si_count = np.bincount(si_nlist[:, 0], minlength=n_o)

    if len(h_pos) > 0:
        h_query = freud.locality.AABBQuery(box, h_pos)
        h_nlist = h_query.query(
            o_pos, {'r_max': O_H_CUTOFF, 'exclude_ii': False}
        ).toNeighborList()
        h_count = np.bincount(h_nlist[:, 0], minlength=n_o)

    bridging = (si_count == 2) & (h_count == 0)
    silanol  = (si_count == 1) & (h_count == 1)
    water    = (si_count == 0) & (h_count == 2)
    other    = ~(bridging | silanol | water)

    return bridging, silanol, water, other


_o_classification_cache = {}

def get_classified_o(frame, kind):
    """Return oxygen positions of a given kind for a single frame."""
    ts = frame['timestep']
    if ts not in _o_classification_cache:
        pos, el = frame['positions'], frame['elements']
        b, s, w, _ = classify_oxygens(
            frame['box'],
            pos[el == 'Si'],
            pos[el == 'O'],
            pos[el == 'H'],
        )
        o_pos = pos[el == 'O']
        _o_classification_cache[ts] = {
            'bridging': o_pos[b],
            'silanol':  o_pos[s],
            'water':    o_pos[w],
        }
    return _o_classification_cache[ts][kind]


def compute_rdf(frames, get_a, get_b, self_pair=False):
    """
    Compute per-frame RDF between position sets returned by get_a and get_b,
    then return the simple mean across frames.
    For self-pairs (O-O, H-H, Si-Si), set self_pair=True so freud excludes
    the i=j distance-zero contribution automatically.
    """
    rdf = freud.density.RDF(bins=BINS, r_max=R_MAX)
    all_rdf = []

    for frame in frames:
        pos_a = get_a(frame)
        pos_b = get_b(frame)

        if len(pos_a) == 0 or len(pos_b) == 0:
            continue

        if self_pair:
            rdf.compute((frame['box'], pos_a), reset=True)
        else:
            rdf.compute((frame['box'], pos_a), query_points=pos_b, reset=True)

        all_rdf.append(rdf.rdf.copy())

    if not all_rdf:
        return rdf.bin_centers, np.zeros(BINS)

    return rdf.bin_centers, np.mean(all_rdf, axis=0)


def build_pair_getters(classify):
    """
    Return a dict of {label: (get_a, get_b, self_pair)} for all pairs to compute.
    classify controls whether oxygen-type pairs are included.
    self_pair=True omits query_points in freud so the i=j contribution is excluded.
    """
    def get_si(f): return f['positions'][f['elements'] == 'Si']
    def get_o(f):  return f['positions'][f['elements'] == 'O']
    def get_h(f):  return f['positions'][f['elements'] == 'H']

    # (get_a, get_b, self_pair)
    pairs = {
        'Si-Si': (get_si, get_si, True),
        'Si-O':  (get_si, get_o,  False),
        'Si-H':  (get_si, get_h,  False),
        'O-O':   (get_o,  get_o,  True),
        'O-H':   (get_o,  get_h,  False),
        'H-H':   (get_h,  get_h,  True),
    }

    if classify:
        def get_o_bridging(f): return get_classified_o(f, 'bridging')
        def get_o_silanol(f):  return get_classified_o(f, 'silanol')
        def get_o_water(f):    return get_classified_o(f, 'water')

        pairs.update({
            'Si-O_bridging': (get_si,      get_o_bridging, False),
            'Si-O_silanol':  (get_si,      get_o_silanol,  False),
            'O_water-H':     (get_o_water, get_h,          False),
        })

    return pairs


def save_csv(results, filename):
    """Save all RDF results to a single CSV with r and one g(r) column per pair."""
    r = next(iter(results.values()))[0]
    header = 'r_Angstrom,' + ','.join(results.keys())
    data = np.column_stack([r] + [g for _, g in results.values()])
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

    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=PLOT_DPI)
    plt.show()
    print(f"Plot saved to {OUTPUT_PLOT}")


if __name__ == '__main__':
    print(f"Reading trajectory: {DUMP_FILE}")
    frames = read_lammps_dump(DUMP_FILE)

    # sanity check on atom counts
    atom_counts = [len(f['positions']) for f in frames]
    print(f"Frames read:  {len(frames)}")
    print(f"Atoms range:  {min(atom_counts)} – {max(atom_counts)}")
    print(f"Atoms mean:   {np.mean(atom_counts):.1f}")

    # oxygen classification sanity check (first frame only)
    if CLASSIFY_OXYGENS:
        f0  = frames[0]
        pos = f0['positions']
        el  = f0['elements']
        b, s, w, other = classify_oxygens(
            f0['box'], pos[el == 'Si'], pos[el == 'O'], pos[el == 'H']
        )
        print(f"\nOxygen classification (frame 0):")
        print(f"  Bridging:  {b.sum()}")
        print(f"  Silanol:   {s.sum()}")
        print(f"  Water:     {w.sum()}")
        print(f"  Other:     {other.sum()}  ← should be ~0; if not, tune cutoffs")

    # compute RDFs
    pairs = build_pair_getters(CLASSIFY_OXYGENS)
    results = {}
    for name, (get_a, get_b, is_self) in pairs.items():
        print(f"Computing RDF: {name}...")
        r, g = compute_rdf(frames, get_a, get_b, self_pair=is_self)
        results[name] = (r, g)

    if OUTPUT_CSV is not None:
        save_csv(results, OUTPUT_CSV)

    plot_rdfs(results)