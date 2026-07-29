"""
bad_freud.py — Bond Angle Distribution calculator using freud

QUICK START
-----------
1. Set DUMP_FILE to your LAMMPS custom dump trajectory.
2. Set ELEMENTS to the list of element symbols present in the simulation.
3. Set R_CUTOFF and R_MINCUT pairwise cutoff radii for every element pair.
4. Optionally add entries to TRIPLET_CUTOFFS for BADs that need custom
   per-arm cutoffs; every other symmetric-unique A-B-C triplet from ELEMENTS
   is computed automatically with the default R_CUTOFF/R_MINCUT.
5. Run:  python bad_freud.py

OUTPUT
------
- A PNG plot of all BAD curves  (OUTPUT_PLOT)
- A CSV table of angle vs P(θ)  (OUTPUT_CSV; set to None to skip)

DEFAULT TRIPLET SWEEP
----------------------
For every central element B in ELEMENTS and every symmetric-unique pair of
wing elements (A, C) — A-B-C and C-B-A are the same angle, so wings are only
enumerated once — a BAD is computed using the default R_CUTOFF/R_MINCUT for
the A-B and C-B pairs, labeled plainly as 'A-B-C'.

TRIPLET_CUTOFFS
---------------
Each entry in TRIPLET_CUTOFFS is a dict specifying one additional, independently
cutoff BAD:
    'triplet' : (el_a, el_b, el_c)  — B is the central atom; wings sorted alphabetically
    'label'   : str                 — column name in the CSV / plot title (unique per triplet)
    'r_max_ab', 'r_min_ab'          — radial cutoffs for the A-B bond
    'r_max_cb', 'r_min_cb'          — radial cutoffs for the C-B bond
                                      (any omitted value falls back to R_CUTOFF/R_MINCUT)

The same triplet may appear more than once under different labels to compute
the BAD with different cutoff windows; each entry is an independent calculation.
Two entries sharing the same (triplet, label) raise a ValueError at startup.
If a TRIPLET_CUTOFFS entry's triplet collides with one from the default sweep,
its label gets a '-specific' suffix so the two calculations stay distinct.

COLUMN LAYOUT
-------------
Expects LAMMPS custom dump format with at minimum columns: id type element x y z
Column positions are read automatically from the ITEM: ATOMS header line.
"""

import itertools
import os
from datetime import date
from typing import NamedTuple

import numpy as np
import freud
import matplotlib
matplotlib.use('Agg')
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

# R_CUTOFF/R_MINCUT/TRIPLET_CUTOFFS are read from plain delimited strings
# (not JSON) because sbatch --export=... is comma-delimited and silently
# truncates any value containing a literal comma.

def _parse_pair_dict(env_name, default):
    """Parse 'H-H:2.0;H-O:1.4' into {'H-H': 2.0, 'H-O': 1.4}, or return default if unset.
    Keys are normalized to alphabetical element order (e.g. 'Si-O' and 'O-Si' both
    become 'O-Si'), matching _pair_key()'s convention used for triplet cutoff lookups.
    """
    raw = os.environ.get(env_name)
    if raw is None:
        return default
    result = {}
    for entry in raw.split(";"):
        key, value = entry.split(":")
        key = "-".join(sorted(key.split("-")))
        result[key] = float(value)
    return result


def _parse_triplet_cutoffs(env_name, default):
    """Parse 'label:elA-elB-elC:r_max_ab:r_min_ab:r_max_cb:r_min_cb|...' into the
    TRIPLET_CUTOFFS list-of-dicts structure, or return default if unset.
    A trailing cutoff field left empty (e.g. 'label:elA-elB-elC:::: ') omits
    that key so it falls back to R_CUTOFF/R_MINCUT, same as a dict literal
    simply not including it.
    """
    raw = os.environ.get(env_name)
    if raw is None:
        return default
    entries = []
    for chunk in raw.split("|"):
        label, triplet_str, r_max_ab, r_min_ab, r_max_cb, r_min_cb = chunk.split(":")
        entry = {'triplet': tuple(triplet_str.split("-")), 'label': label}
        if r_max_ab:
            entry['r_max_ab'] = float(r_max_ab)
        if r_min_ab:
            entry['r_min_ab'] = float(r_min_ab)
        if r_max_cb:
            entry['r_max_cb'] = float(r_max_cb)
        if r_min_cb:
            entry['r_min_cb'] = float(r_min_cb)
        entries.append(entry)
    return entries


# Elements present in the simulation.
# The code trusts this list and never auto-detects species from the trajectory.
# Order does not matter — code sorts internally.
# Semicolon-separated (not comma) since sbatch --export=... is comma-delimited.
ELEMENTS = [e.strip() for e in os.environ.get("ELEMENTS", "Si;O;H").split(";")]

# Pairwise upper neighbor cutoff radii in Angstroms.
# Keys must be "El1-El2" with elements sorted alphabetically.
R_CUTOFF = _parse_pair_dict("R_CUTOFF", {
    'H-H':   2.0,
    'H-O':   1.4,
    'H-Si':  2.0,
    'O-O':   2.8,
    'O-Si':  2.2,
    'Si-Si': 3.2,
})

# Pairwise lower neighbor cutoff radii in Angstroms.
# Pairs closer than this are excluded (removes self-interactions and unphysical contacts).
R_MINCUT = _parse_pair_dict("R_MINCUT", {
    'H-H':   0.5,
    'H-O':   0.5,
    'H-Si':  0.5,
    'O-O':   0.5,
    'O-Si':  0.5,
    'Si-Si': 0.5,
})

# List of BADs to compute. Each entry is a dict with:
#   'triplet' : (el_a, el_b, el_c)  — B is the central atom; wings sorted alphabetically
#   'label'   : str                 — output label (must be unique per triplet)
#   'r_max_ab', 'r_min_ab'          — cutoffs for the A-B bond (fall back to R_CUTOFF/R_MINCUT)
#   'r_max_cb', 'r_min_cb'          — cutoffs for the C-B bond (fall back to R_CUTOFF/R_MINCUT)
#
# The same triplet may appear more than once with different labels and cutoffs;
# each entry is treated as an independent BAD.
# Two entries with identical (triplet, label) raise a ValueError at startup.
TRIPLET_CUTOFFS = _parse_triplet_cutoffs("TRIPLET_CUTOFFS", [
    {'triplet': ('O', 'Si', 'O'), 'label': 'O-Si-O',
     'r_max_ab': 2.2, 'r_min_ab': 0.5, 'r_max_cb': 2.2, 'r_min_cb': 0.5},
])

# Number of bins spanning 0–180°
BINS = int(os.environ.get("BAD_BINS", "180"))

# Plot layout
PLOT_NCOLS = 3
PLOT_DPI   = 150

# =============================================================================
# END CONFIGURATION
# =============================================================================

# Prepend today's date (YYYYMMDD_) to every output filename.
def _dated(filename):
    return None if filename is None else f"{date.today():%Y%m%d}_{filename}"

OUTPUT_PLOT = _dated(OUTPUT_PLOT)
OUTPUT_CSV  = _dated(OUTPUT_CSV)


def _pair_key(el1, el2):
    """Return canonical pairwise key with elements sorted alphabetically."""
    return '-'.join(sorted([el1, el2]))


def _canonical_triplet(el_a, el_b, el_c):
    """Sort wing elements so A-B-C and C-B-A map to the same key (same angle)."""
    wing1, wing2 = sorted([el_a, el_c])
    return (wing1, el_b, wing2)


def _default_triplet_cutoffs(elements):
    """
    Generate every symmetric-unique A-B-C triplet (B central) from `elements`,
    using the default R_CUTOFF/R_MINCUT pairwise radii (no explicit r_max_*/r_min_*
    keys — the main loop falls back to R_CUTOFF/R_MINCUT for those).
    """
    entries = []
    for el_b in elements:
        for el_a, el_c in itertools.combinations_with_replacement(sorted(elements), 2):
            entries.append({'triplet': (el_a, el_b, el_c), 'label': f'{el_a}-{el_b}-{el_c}'})
    return entries


def _build_triplet_cutoffs():
    """
    Combine the auto-generated default-cutoff sweep with the user-specified
    TRIPLET_CUTOFFS. If a TRIPLET_CUTOFFS entry's triplet collides with one
    from the default sweep, its label gets a '-specific' suffix so both stay
    distinct in the results.
    """
    defaults = _default_triplet_cutoffs(ELEMENTS)
    default_keys = {tuple(e['triplet']) for e in defaults}

    specifics = []
    for entry in TRIPLET_CUTOFFS:
        entry = dict(entry)
        if _canonical_triplet(*entry['triplet']) in default_keys:
            entry['label'] = f"{entry['label']}-specific"
        specifics.append(entry)

    return defaults + specifics


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


class _Q(NamedTuple):
    """
    One frame's worth of bonds for a single unique neighbor query.

    Bonds are stored flat and grouped by central atom: bonds
    ``start[i] : start[i] + cnt[i]`` all belong to central atom i.

    u     : (n_bonds, 3) float32 — unit vectors pointing central → wing
    w     : (n_bonds,)   intp    — wing atom index, in the wing element's own index space
    start : (n_b,)       intp    — first bond offset for each central atom
    cnt   : (n_b,)       intp    — number of bonds for each central atom
    """
    u:     np.ndarray
    w:     np.ndarray
    start: np.ndarray
    cnt:   np.ndarray


def _run_query(box, pos_central, pos_wing, r_max, r_min):
    """
    Evaluate one neighbor query and return it as a central-atom-grouped _Q.

    The AABBQuery is built over the *wing* atoms and probed with the *central*
    atoms — the reverse of the intuitive direction. Distance is symmetric so the
    resulting bond set is identical either way, but this orientation makes freud
    return the bonds already segmented by central atom (segments/neighbor_counts),
    which is exactly the grouping the angle kernel needs. Querying the other way
    would require re-bucketing every bond in Python.

    freud's vectors are already minimum-image corrected, so no box.wrap is needed.
    """
    n_b = len(pos_central)

    # No wing atoms → no bonds, but still return per-central arrays of the right
    # length so callers need no special case.
    if len(pos_wing) == 0 or n_b == 0:
        return _Q(u=np.empty((0, 3), dtype=np.float32),
                  w=np.empty(0, dtype=np.intp),
                  start=np.zeros(n_b, dtype=np.intp),
                  cnt=np.zeros(n_b, dtype=np.intp))

    nl = freud.locality.AABBQuery(box, pos_wing).query(
        pos_central, {'r_max': r_max, 'r_min': r_min, 'exclude_ii': False}
    ).toNeighborList()

    d = nl.distances
    # r_min > 0 excludes coincident atoms, so this should be unreachable; assert
    # rather than let a zero distance produce nan angles that silently vanish
    # from the histogram.
    if len(d) and d.min() <= 0.0:
        raise ValueError(f"zero-length bond in query r_max={r_max}, r_min={r_min}")

    # freud hands back uint32; cast now because every downstream use feeds integer
    # arithmetic that goes negative mid-expression and would wrap silently.
    return _Q(u=nl.vectors / d[:, None],
              w=nl.point_indices.astype(np.intp),
              start=nl.segments.astype(np.intp),
              cnt=nl.neighbor_counts.astype(np.intp))


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


def _resolve_cutoffs(entry):
    """
    Resolve an entry's four radial cutoffs, falling back to the pairwise
    R_CUTOFF/R_MINCUT defaults for any the entry leaves unspecified.
    """
    el_a, el_b, el_c = entry['triplet']
    key_ab = _pair_key(el_a, el_b)
    key_cb = _pair_key(el_c, el_b)
    # NOTE: dict.get(key, default_expr) always evaluates default_expr eagerly,
    # even when key is already present — so R_CUTOFF[key_ab] would raise
    # KeyError even for entries that fully specify their own cutoffs. Use an
    # explicit conditional instead so R_CUTOFF/R_MINCUT are only consulted
    # when actually needed.
    return (
        entry['r_max_ab'] if 'r_max_ab' in entry else R_CUTOFF[key_ab],
        entry['r_min_ab'] if 'r_min_ab' in entry else R_MINCUT[key_ab],
        entry['r_max_cb'] if 'r_max_cb' in entry else R_CUTOFF[key_cb],
        entry['r_min_cb'] if 'r_min_cb' in entry else R_MINCUT[key_cb],
    )


def _plan_queries(entries):
    """
    Resolve every entry's cutoffs up front and reduce its two arms to neighbor
    query keys, collecting the unique set of queries the whole run needs.

    A query is identified by (central element, wing element, r_max, r_min), so
    arms agreeing on all four are evaluated once and shared by every triplet
    that asked for them: O-H-O and O-H-H both need the O neighbors of every H
    under the same cutoffs, and the default sweep builds O(E^3) triplets out of
    only O(E^2) distinct arms.

    Returns
    -------
    plans   : list of dicts — 'triplet', 'label', 'key_ab', 'key_cb', 'same_wing'
    queries : list of unique query keys, in first-use order
    """
    plans   = []
    queries = {}
    for entry in entries:
        el_a, el_b, el_c = entry['triplet']
        r_max_ab, r_min_ab, r_max_cb, r_min_cb = _resolve_cutoffs(entry)

        key_ab = (el_b, el_a, r_max_ab, r_min_ab)
        key_cb = (el_b, el_c, r_max_cb, r_min_cb)
        queries[key_ab] = None
        queries[key_cb] = None

        plans.append({
            'triplet': tuple(entry['triplet']),
            'label':   entry['label'],
            'key_ab':  key_ab,
            'key_cb':  key_cb,
            # The symmetric branch is valid exactly when both arms reduce to the
            # same query — same wing element AND same cutoffs. If the cutoffs
            # differ (e.g. O-H--O), the keys differ and each arm is queried
            # independently.
            'same_wing': key_ab == key_cb,
        })
    return plans, list(queries)


# Max pair slots expanded at once, to bound peak memory. A dense structure with
# generous cutoffs can otherwise produce tens of millions of pairs in one frame.
PAIR_CHUNK = 2_000_000


def _triplet_angles(qa, qc, same_wing, same_wing_element):
    """
    Every A-B-C angle in one frame, computed without looping over central atoms.

    For each central atom the candidate pairs form an (n_a x n_c) grid: every
    arm-A neighbor against every arm-C neighbor. Rather than build those grids
    one atom at a time, reserve one slot for every pair in the whole frame and
    recover, by arithmetic on a slot's position alone, which central atom it
    belongs to and which cell of that atom's grid it is — division gives the row,
    remainder gives the column. Grids stay ragged; nothing is iterated.

    Parameters
    ----------
    qa, qc            : _Q records for the two arms (same central element, so
                        their per-central arrays are the same length)
    same_wing         : both arms are literally the same query
    same_wing_element : el_a == el_c (wing indices are comparable between arms)

    Returns
    -------
    list of angle arrays in degrees (possibly empty)
    """
    n_a, n_c = qa.cnt, qc.cnt
    npair_all = n_a * n_c                       # pair slots per central atom
    n_b = len(n_a)
    if n_b == 0 or npair_all.sum() == 0:
        return []

    # Chunk boundaries over central atoms, sized so no chunk expands past
    # PAIR_CHUNK slots. An atom whose own npair exceeds it still lands in one
    # chunk, which is fine — per-atom counts are bounded by coordination number.
    cum = np.cumsum(npair_all)
    edges = np.searchsorted(cum, np.arange(PAIR_CHUNK, cum[-1] + PAIR_CHUNK, PAIR_CHUNK))
    bounds = np.unique(np.concatenate(([0], np.clip(edges, 0, n_b), [n_b])))

    out = []
    for b0, b1 in zip(bounds[:-1], bounds[1:]):
        npair = npair_all[b0:b1]
        tot = int(npair.sum())
        if tot == 0:
            continue

        b    = np.repeat(np.arange(b0, b1, dtype=np.intp), npair)
        off  = np.cumsum(npair) - npair                                # exclusive prefix
        ordn = np.arange(tot, dtype=np.intp) - np.repeat(off, npair)   # restarts per central

        # Unravel the flat ordinal into that central atom's (n_a x n_c) grid.
        # Dividing by arm C's count makes l_a the row and l_c the column.
        # No division by zero: a central with n_c == 0 has npair == 0, so
        # np.repeat drops it before it can appear in b.
        n_c_b = n_c[b]
        l_a   = ordn // n_c_b
        l_c   = ordn - l_a * n_c_b

        i_a = qa.start[b] + l_a
        i_c = qc.start[b] + l_c

        if same_wing:
            # Both arms are the same neighbor list, so the grid holds each pair
            # twice plus the diagonal; keeping the upper triangle is exactly
            # np.triu_indices(n, k=1).
            keep = l_a < l_c
            i_a, i_c = i_a[keep], i_c[keep]
        elif same_wing_element:
            # Same element on both arms but different cutoffs: drop the cases
            # where both arms picked the same physical atom (spurious 0° angle).
            keep = qa.w[i_a] != qc.w[i_c]
            i_a, i_c = i_a[keep], i_c[keep]
        # Different elements: no atom can appear on both arms, so nothing to drop.

        if len(i_a) == 0:
            continue

        # Unit vectors, so the dot product is the cosine directly — no norms.
        dot = np.einsum('ij,ij->i', qa.u[i_a], qc.u[i_c])
        out.append(np.degrees(np.arccos(np.clip(dot, -1.0, 1.0))))

    return out


def compute_bads(frames, plans, queries):
    """
    Compute every planned bond angle distribution in a single pass over `frames`.

    Only A-B and C-B distances within their respective [r_min, r_max] are
    included. Each frame's histogram is normalized to unit-area probability
    density, then the mean is taken across the frames that contributed to it.

    Frames are the outer loop so per-frame work is shared across triplets: each
    element's positions are masked out once, and each unique query in `queries`
    is evaluated once — including its unit bond vectors — then reused by every
    arm pointing at it. Angles for a triplet are computed for all central atoms
    at once (see _triplet_angles).

    Returns
    -------
    dict {label: (bin_centers, histogram)} — histogram integrates to 1
    """
    bin_edges   = np.linspace(0.0, 180.0, BINS + 1)
    bin_width   = bin_edges[1] - bin_edges[0]
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    # Accumulate by plan index rather than label: two plans may share a label,
    # and keying by label would silently sum them into one curve.
    hist_accum    = [np.zeros(BINS) for _ in plans]
    n_frames_used = [0] * len(plans)

    elements = {key[0] for key in queries} | {key[1] for key in queries}

    for frame in frames:
        els = frame['elements']
        pos = frame['positions']
        box = frame['box']

        pos_by_el = {el: pos[els == el] for el in elements}

        # Evaluate each unique query once for this frame. Bond vectors depend
        # only on the query, not on which triplet consumes it, so unit vectors
        # are built here too rather than once per triplet.
        neighbors = {
            key: _run_query(box, pos_by_el[key[0]], pos_by_el[key[1]], key[2], key[3])
            for key in queries
        }

        for idx, plan in enumerate(plans):
            el_a, _, el_c = plan['triplet']

            frame_angles = _triplet_angles(
                neighbors[plan['key_ab']],
                neighbors[plan['key_cb']],
                same_wing=plan['same_wing'],
                same_wing_element=(el_a == el_c),
            )

            if not frame_angles:
                continue

            h, _ = np.histogram(np.concatenate(frame_angles), bins=bin_edges)
            total = h.sum()
            if total > 0:
                # normalize so that Σ P(θᵢ) · Δθ = 1
                hist_accum[idx] += h / (total * bin_width)
                n_frames_used[idx] += 1

    results = {}
    for idx, plan in enumerate(plans):
        n = n_frames_used[idx]
        results[plan['label']] = (
            bin_centers,
            hist_accum[idx] / n if n else np.zeros(BINS),
        )
    return results


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
    plt.close(fig)
    print(f"Plot saved to {OUTPUT_PLOT}")


if __name__ == '__main__':
    if not ELEMENTS:
        raise ValueError("ELEMENTS is empty — set it to the element symbols present in the simulation.")

    all_triplet_cutoffs = _build_triplet_cutoffs()
    _validate_triplet_cutoffs(all_triplet_cutoffs)

    if not all_triplet_cutoffs:
        raise ValueError("No BADs configured: ELEMENTS and TRIPLET_CUTOFFS are both empty.")

    print(f"Reading trajectory: {DUMP_FILE}")
    frames = read_lammps_dump(DUMP_FILE)

    atom_counts = [len(f['positions']) for f in frames]
    print(f"Frames read: {len(frames)}")
    print(f"Atoms range: {min(atom_counts)} – {max(atom_counts)}")
    print(f"Atoms mean:  {np.mean(atom_counts):.1f}")

    n_default = len(all_triplet_cutoffs) - len(TRIPLET_CUTOFFS)
    print(f"BADs to compute: {len(all_triplet_cutoffs)} "
          f"({n_default} default + {len(TRIPLET_CUTOFFS)} specific)")

    plans, queries = _plan_queries(all_triplet_cutoffs)
    print(f"Unique neighbor queries per frame: {len(queries)} "
          f"(from {2 * len(plans)} triplet arms)")

    print("Computing BADs...")
    results = compute_bads(frames, plans, queries)

    if OUTPUT_CSV is not None:
        save_csv(results, OUTPUT_CSV)

    plot_bads(results)
