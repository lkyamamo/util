"""
bad_freud.py — Bond Angle Distribution calculator using freud

QUICK START
-----------
1. Set DUMP_FILE to your LAMMPS custom dump trajectory.
2. Set ELEMENTS to the list of element symbols present in the simulation.
3. Set R_CUTOFF and R_MINCUT pairwise cutoff radii for every element pair.
4. Populate TRIPLET_CUTOFFS with the BADs to compute and their per-arm cutoffs.
5. Run:  python bad_freud.py

OUTPUT
------
- A PNG plot of all BAD curves  (OUTPUT_PLOT)
- A CSV table of angle vs P(θ)  (OUTPUT_CSV; set to None to skip)

TRIPLET_CUTOFFS
---------------
Each entry in TRIPLET_CUTOFFS is a dict specifying one BAD:
    'triplet' : (el_a, el_b, el_c)  — B is the central atom; wings sorted alphabetically
    'label'   : str                 — column name in the CSV / plot title (unique per triplet)
    'r_max_ab', 'r_min_ab'          — radial cutoffs for the A-B bond
    'r_max_cb', 'r_min_cb'          — radial cutoffs for the C-B bond
                                      (any omitted value falls back to R_CUTOFF/R_MINCUT)

The same triplet may appear more than once under different labels to compute
the BAD with different cutoff windows; each entry is an independent calculation.
Two entries sharing the same (triplet, label) raise a ValueError at startup.

COLUMN LAYOUT
-------------
Expects LAMMPS custom dump format with at minimum columns: id type element x y z
Column positions are read automatically from the ITEM: ATOMS header line.
"""

import os

import numpy as np
import freud
import matplotlib.pyplot as plt

# =============================================================================
# CONFIGURATION — edit these variables between runs
# =============================================================================

# Input trajectory file
DUMP_FILE = os.environ.get("TRAJ", "../OH.lammpstrj")

# Output plot file
OUTPUT_PLOT = "bads.png"

# Output data table (CSV); set to None to skip
OUTPUT_CSV = "bads.csv"

# Elements present in the simulation.
# The code trusts this list and never auto-detects species from the trajectory.
# Order does not matter — code sorts internally.
ELEMENTS = ['Si', 'O', 'H']

# Pairwise upper neighbor cutoff radii in Angstroms.
# Keys must be "El1-El2" with elements sorted alphabetically.
R_CUTOFF = {
    'H-H':   2.0,
    'H-O':   1.4,
    'H-Si':  2.0,
    'O-O':   2.8,
    'O-Si':  2.2,
    'Si-Si': 3.2,
}

# Pairwise lower neighbor cutoff radii in Angstroms.
# Pairs closer than this are excluded (removes self-interactions and unphysical contacts).
R_MINCUT = {
    'H-H':   0.5,
    'H-O':   0.5,
    'H-Si':  0.5,
    'O-O':   0.5,
    'O-Si':  0.5,
    'Si-Si': 0.5,
}

# List of BADs to compute. Each entry is a dict with:
#   'triplet' : (el_a, el_b, el_c)  — B is the central atom; wings sorted alphabetically
#   'label'   : str                 — output label (must be unique per triplet)
#   'r_max_ab', 'r_min_ab'          — cutoffs for the A-B bond (fall back to R_CUTOFF/R_MINCUT)
#   'r_max_cb', 'r_min_cb'          — cutoffs for the C-B bond (fall back to R_CUTOFF/R_MINCUT)
#
# The same triplet may appear more than once with different labels and cutoffs;
# each entry is treated as an independent BAD.
# Two entries with identical (triplet, label) raise a ValueError at startup.
TRIPLET_CUTOFFS = [
    {'triplet': ('O', 'Si', 'O'), 'label': 'O-Si-O',
     'r_max_ab': 2.2, 'r_min_ab': 0.5, 'r_max_cb': 2.2, 'r_min_cb': 0.5},
]

# Number of bins spanning 0–180°
BINS = 180

# Plot layout
PLOT_NCOLS = 3
PLOT_DPI   = 150

# =============================================================================
# END CONFIGURATION
# =============================================================================


def _pair_key(el1, el2):
    """Return canonical pairwise key with elements sorted alphabetically."""
    return '-'.join(sorted([el1, el2]))


def _validate_triplet_cutoffs(entries):
    """Raise ValueError if any (triplet, label) pair appears more than once."""
    seen = set()
    for entry in entries:
        key = (tuple(entry['triplet']), entry['label'])
        if key in seen:
            raise ValueError(
                f"Duplicate (triplet, label) in TRIPLET_CUTOFFS: "
                f"triplet={entry['triplet']}, label={entry['label']!r}"
            )
        seen.add(key)


def _query_neighbors(aq, query_pos, r_max, r_min):
    """
    Query wing-atom positions against a pre-built AABBQuery over B atoms.

    Parameters
    ----------
    aq        : freud.locality.AABBQuery built over pos_b
    query_pos : positions of A or C atoms (shape N_wing x 3)
    r_max     : upper cutoff for this pair (from R_CUTOFF)
    r_min     : lower cutoff for this pair (from R_MINCUT)

    Returns
    -------
    neighbors : list of length len(pos_b); neighbors[b_idx] is a list of
                indices into query_pos that are within [r_min, r_max] of
                that B atom.
    """
    n_b = len(aq.points)
    if len(query_pos) == 0:
        return [[] for _ in range(n_b)]

    result = aq.query(query_pos, {'r_max': r_max, 'r_min': r_min, 'exclude_ii': False})
    neighbors = [[] for _ in range(n_b)]
    for bond in result:
        # bond = (query_point_index, point_index, distance)
        # bond[0] → index into query_pos (A or C atoms)
        # bond[1] → index into pos_b (B reference atoms)
        neighbors[bond[1]].append(bond[0])
    return neighbors


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
                'timestep':  timestep,
                'box':       box,
                'positions': positions,
                'elements':  np.array(elements),
            })

    return frames


def compute_bad(frames, el_a, el_b, el_c, r_max_ab, r_min_ab, r_max_cb, r_min_cb):
    """
    Compute the bond angle distribution for A-B-C triplets (B is the central atom).

    Only A-B and C-B distances within their respective pairwise [r_min, r_max]
    are included. Each frame's histogram is normalized to unit-area probability
    density, then the mean is taken across frames.

    Parameters
    ----------
    r_max_ab, r_min_ab : upper/lower cutoff for the A-B bond
    r_max_cb, r_min_cb : upper/lower cutoff for the C-B bond

    Returns
    -------
    bin_centers : ndarray, shape (BINS,) — angles in degrees
    histogram   : ndarray, shape (BINS,) — mean P(θ), integrates to 1
    """
    # Use the symmetric same-wing branch only when both element and cutoffs match.
    # If the cutoffs differ (e.g. O-H--O), fall through to the different-wing
    # branch so each arm is queried independently.
    same_wing = (
        el_a == el_c
        and r_max_ab == r_max_cb
        and r_min_ab == r_min_cb
    )

    bin_edges     = np.linspace(0.0, 180.0, BINS + 1)
    bin_width     = bin_edges[1] - bin_edges[0]
    hist_accum    = np.zeros(BINS)
    n_frames_used = 0

    for frame in frames:
        els = frame['elements']
        pos = frame['positions']
        box = frame['box']

        pos_b = pos[els == el_b]
        pos_a = pos[els == el_a]
        pos_c = pos[els == el_c]

        if len(pos_b) == 0 or len(pos_a) == 0 or len(pos_c) == 0:
            continue

        aq = freud.locality.AABBQuery(box, pos_b)

        # ---- same-wing (A == C, same cutoffs, e.g. O-Si-O) ------------------
        if same_wing:
            a_neighbors = _query_neighbors(aq, pos_a, r_max_ab, r_min_ab)

            frame_angles = []
            for b_idx in range(len(pos_b)):
                a_nbrs = a_neighbors[b_idx]
                if len(a_nbrs) < 2:
                    continue

                v_a = box.wrap(pos_a[a_nbrs] - pos_b[b_idx])  # (n_a, 3) MIC vectors

                i, j = np.triu_indices(len(a_nbrs), k=1)      # unique pairs, no self-pairs
                v1, v2 = v_a[i], v_a[j]

                cos_theta = (
                    np.einsum('ij,ij->i', v1, v2)
                    / (np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1))
                )
                frame_angles.append(
                    np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
                )

        # ---- different-wing (A != C, or same element with different cutoffs) -
        else:
            a_neighbors = _query_neighbors(aq, pos_a, r_max_ab, r_min_ab)
            c_neighbors = _query_neighbors(aq, pos_c, r_max_cb, r_min_cb)

            frame_angles = []
            for b_idx in range(len(pos_b)):
                a_nbrs = np.array(a_neighbors[b_idx])
                c_nbrs = np.array(c_neighbors[b_idx])
                if len(a_nbrs) == 0 or len(c_nbrs) == 0:
                    continue

                v_a = box.wrap(pos_a[a_nbrs] - pos_b[b_idx])  # (n_a, 3) MIC vectors
                v_c = box.wrap(pos_c[c_nbrs] - pos_b[b_idx])  # (n_c, 3) MIC vectors

                ii = np.repeat(np.arange(len(a_nbrs)), len(c_nbrs))
                jj = np.tile(np.arange(len(c_nbrs)), len(a_nbrs))

                # When wings are the same element, exclude self-pairs (same atom
                # appearing on both arms gives a spurious 0° angle).
                if el_a == el_c:
                    mask = a_nbrs[ii] != c_nbrs[jj]
                    ii, jj = ii[mask], jj[mask]

                if len(ii) == 0:
                    continue

                v1, v2 = v_a[ii], v_c[jj]                     # A x C Cartesian product

                cos_theta = (
                    np.einsum('ij,ij->i', v1, v2)
                    / (np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1))
                )
                frame_angles.append(
                    np.degrees(np.arccos(np.clip(cos_theta, -1.0, 1.0)))
                )

        if not frame_angles:
            continue

        all_angles = np.concatenate(frame_angles)
        h, _ = np.histogram(all_angles, bins=bin_edges)
        total = h.sum()
        if total > 0:
            # normalize so that Σ P(θᵢ) · Δθ = 1
            hist_accum += h / (total * bin_width)
            n_frames_used += 1

    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    if n_frames_used == 0:
        return bin_centers, np.zeros(BINS)
    return bin_centers, hist_accum / n_frames_used


def save_csv(results, filename):
    """Save results dict {label: (angles, hist)} to CSV with one column per label."""
    angles = next(iter(results.values()))[0]
    header = 'angle_deg,' + ','.join(results.keys())
    data   = np.column_stack([angles] + [h for _, h in results.values()])
    np.savetxt(filename, data, delimiter=',', header=header, comments='', fmt='%.6f')
    print(f"Data table saved to {filename}")


def plot_bads(results):
    n     = len(results)
    ncols = PLOT_NCOLS
    nrows = (n + ncols - 1) // ncols

    fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 3.5 * nrows), squeeze=False)
    axes = axes.flatten()

    for ax, (name, (angles, hist)) in zip(axes, results.items()):
        ax.plot(angles, hist)
        ax.set_xlabel('Angle (degrees)')
        ax.set_ylabel('P(θ)')
        ax.set_title(name)
        ax.set_xlim(0, 180)

    for ax in axes[n:]:
        ax.set_visible(False)

    plt.tight_layout()
    plt.savefig(OUTPUT_PLOT, dpi=PLOT_DPI)
    plt.show()
    print(f"Plot saved to {OUTPUT_PLOT}")


if __name__ == '__main__':
    _validate_triplet_cutoffs(TRIPLET_CUTOFFS)

    print(f"Reading trajectory: {DUMP_FILE}")
    frames = read_lammps_dump(DUMP_FILE)

    atom_counts = [len(f['positions']) for f in frames]
    print(f"Frames read: {len(frames)}")
    print(f"Atoms range: {min(atom_counts)} – {max(atom_counts)}")
    print(f"Atoms mean:  {np.mean(atom_counts):.1f}")

    print(f"BADs to compute: {len(TRIPLET_CUTOFFS)}")

    results = {}
    for entry in TRIPLET_CUTOFFS:
        el_a, el_b, el_c = entry['triplet']
        label    = entry['label']
        key_ab   = _pair_key(el_a, el_b)
        key_cb   = _pair_key(el_c, el_b)
        r_max_ab = entry.get('r_max_ab', R_CUTOFF[key_ab])
        r_min_ab = entry.get('r_min_ab', R_MINCUT[key_ab])
        r_max_cb = entry.get('r_max_cb', R_CUTOFF[key_cb])
        r_min_cb = entry.get('r_min_cb', R_MINCUT[key_cb])
        print(f"Computing BAD: {label}...")
        angles, hist = compute_bad(frames, el_a, el_b, el_c,
                                   r_max_ab, r_min_ab, r_max_cb, r_min_cb)
        results[label] = (angles, hist)

    if OUTPUT_CSV is not None:
        save_csv(results, OUTPUT_CSV)

    plot_bads(results)
