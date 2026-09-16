#!/usr/bin/env python3
"""
ring_geometry.py

Geometry of a single Si-O ring, for studies that walk something through it.

The ring files this reads are the per-ring data files written by
silica_rings_to_lammps.py: 2n atoms (n Si, n O) in ring order, with image flags
chosen so the ring is contiguous once unwrapped. Everything here works on the
unwrapped ring, so no box or minimum-image convention is involved - a ring that
straddles a periodic boundary is already whole.

What the aperture means
-----------------------
The atoms are treated as points. The aperture of a line through the ring is the
distance from that line to the nearest atom *centre*,

    aperture(line) = min_i dist(p_i, line)

and the ring's aperture is that maximised over the lines considered. No radii
enter it, so it is a pure statement about where the nuclei are and nothing has
to be assumed about how big an atom is.

Contact radii have not gone away, they have just stopped deciding the answer.
Every ring is additionally checked against CONTACT_RADII and reports how many of
its atoms lie closer to the axis than their own radius - the cases a hard-sphere
picture would have called an overlap. That is a flag to look at, not a number
the aperture is measured with, and changing the radii changes the flags without
moving a single aperture.

The axis
--------
By default the path axis is the simplest defensible one: the line through the
ring's centroid along its best-fit plane normal. No optimisation, nothing to
tune, and identical for two people who implement it independently.

`--optimise-axis` turns on a search for a better line - first sliding it across
the plane (aperture_plane), then tilting it within a cone (aperture_cyl). The
cone is load-bearing: an isolated ring has nothing on its flanks, so an
unconstrained "largest cylinder" just goes around the outside. Both optimised
apertures are >= the centroid one by construction. Useful for judging how much
the simple axis gives up on a puckered ring; not the default.

"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize

# LAMMPS atom types in the ring files, and the contact radii used to flag close
# approaches. These are deliberately smaller than vdW radii. They do not enter
# any aperture - apertures are measured between the axis and the atom centres -
# so changing them moves the flags and nothing else.
SI_TYPE, O_TYPE = 1, 2
CONTACT_RADII = {SI_TYPE: 1.5, O_TYPE: 0.9}

# How far the optimised axis may tilt away from the best-fit plane normal, in
# degrees.
#
# This constraint is load-bearing physics, not a numerical convenience. An
# isolated ring has nothing on its flanks, so the largest cylinder that avoids
# every ring atom is unbounded - tip the axis far enough and it simply goes
# around the outside of the ring instead of through it. "Through" is what the
# cone defines, and without it the optimum is meaningless.
#
# 15 degrees comfortably covers the pucker of a realistic glass ring while
# staying well short of sideways escape. It is a modelling choice: rings that
# optimise to the wall are counted and reported, so the choice can be audited
# rather than assumed.
MAX_TILT_DEG = 15.0


@dataclass
class Ring:
    """One ring, unwrapped. Arrays are in ring order: Si, O, Si, O, ..."""
    ids: np.ndarray          # LAMMPS atom ids
    types: np.ndarray        # LAMMPS atom types
    positions: np.ndarray    # (2n, 3) unwrapped coordinates
    box: dict                # the original box, carried through for writers
    masses: dict             # type -> mass, carried through for writers
    source: str              # path this was read from

    @property
    def n(self) -> int:
        """Ring size: the number of silicons."""
        return int((self.types == SI_TYPE).sum())

    def radii(self, radii: dict[int, float] | None = None) -> np.ndarray:
        """Per-atom contact radii, for flagging close approaches only."""
        table = CONTACT_RADII if radii is None else radii
        return np.array([table[int(t)] for t in self.types], dtype=float)

    @property
    def silicon_ids(self) -> np.ndarray:
        return self.ids[self.types == SI_TYPE]

    @property
    def oxygen_ids(self) -> np.ndarray:
        return self.ids[self.types == O_TYPE]


# -- reading ------------------------------------------------------------------

SECTION_NAMES = {"Atoms", "Velocities", "Masses"}


def read_ring(path: str) -> Ring:
    """
    Read one per-ring data file and unwrap it with its image flags.

    This is a deliberately small reader rather than lammps_data.load: these files
    are written by silica_rings_to_lammps.py in one fixed shape (atom style
    atomic, always image flags, a handful of atoms), and there are thousands of
    them to scan. The cluster builder, which writes files and has to get every
    column right, uses lammps_data instead.
    """
    box: dict[str, float] = {}
    masses: dict[int, float] = {}
    rows: list[list[str]] = []
    section = None

    for raw in open(path):
        line = raw.split("#")[0].strip()
        if not line:
            continue
        if line in SECTION_NAMES:
            section = line
            continue
        words = line.split()
        if section == "Atoms":
            rows.append(words)
        elif section == "Masses":
            masses[int(words[0])] = float(words[1])
        elif words[-2:] in (["xlo", "xhi"], ["ylo", "yhi"], ["zlo", "zhi"]):
            box[words[-2]], box[words[-1]] = float(words[0]), float(words[1])

    if not rows:
        raise ValueError(f"{path}: no Atoms section found")
    if len(rows[0]) < 8:
        raise ValueError(
            f"{path}: expected atom style atomic with image flags "
            f"(id type x y z ix iy iz), found {len(rows[0])} columns"
        )

    ids = np.array([int(r[0]) for r in rows])
    types = np.array([int(r[1]) for r in rows])
    wrapped = np.array([[float(r[2]), float(r[3]), float(r[4])] for r in rows])
    images = np.array([[int(r[5]), int(r[6]), int(r[7])] for r in rows])

    lengths = np.array([box["xhi"] - box["xlo"],
                        box["yhi"] - box["ylo"],
                        box["zhi"] - box["zlo"]])
    return Ring(ids=ids, types=types, positions=wrapped + images * lengths,
                box=box, masses=masses, source=path)


@dataclass
class Structure:
    """A whole source structure, wrapped into its box. Orthogonal cells only."""
    ids: np.ndarray
    types: np.ndarray
    positions: np.ndarray
    lengths: np.ndarray
    origin: np.ndarray

    def index_of(self) -> dict[int, int]:
        return {int(a): i for i, a in enumerate(self.ids)}


def read_structure(path: str) -> Structure:
    """
    Read a whole LAMMPS data file as bare arrays, wrapped into the box.

    Same reasoning as read_ring: this is for measuring, not editing, and the
    structures are large enough that a light reader is worth having. Anything
    that *writes* a data file should go through lammps_data instead.
    """
    box: dict[str, float] = {}
    rows, section = [], None
    for raw in open(path):
        line = raw.split("#")[0].strip()
        if not line:
            continue
        if line in SECTION_NAMES:
            section = line
            continue
        words = line.split()
        if section == "Atoms":
            rows.append(words)
        elif section is None and words[-2:] in (["xlo", "xhi"], ["ylo", "yhi"],
                                                ["zlo", "zhi"]):
            box[words[-2]], box[words[-1]] = float(words[0]), float(words[1])

    if not rows:
        raise ValueError(f"{path}: no Atoms section found")
    origin = np.array([box["xlo"], box["ylo"], box["zlo"]])
    lengths = np.array([box["xhi"] - box["xlo"], box["yhi"] - box["ylo"],
                        box["zhi"] - box["zlo"]])
    positions = np.array([[float(r[2]), float(r[3]), float(r[4])] for r in rows])
    return Structure(
        ids=np.array([int(r[0]) for r in rows]),
        types=np.array([int(r[1]) for r in rows]),
        positions=origin + np.mod(positions - origin, lengths),
        lengths=lengths, origin=origin)


def si_o_neighbours(structure: Structure, cutoff: float = 1.9) -> dict[int, list[int]]:
    """
    Si id -> the ids of the oxygens bonded to it, under PBC.

    The cutoff must match the one the ring finder used, or the two will disagree
    about what is bonded and a ring's own oxygens may not come back as
    neighbours of its own silicons.
    """
    from scipy.spatial import cKDTree

    silicon = np.where(structure.types == SI_TYPE)[0]
    oxygen = np.where(structure.types == O_TYPE)[0]
    folded = structure.positions - structure.origin
    tree = cKDTree(folded[oxygen], boxsize=structure.lengths)
    return {int(structure.ids[s]): [int(structure.ids[oxygen[h]]) for h in hits]
            for s, hits in zip(silicon, tree.query_ball_point(folded[silicon], cutoff))}


def minimum_image(delta: np.ndarray, lengths: np.ndarray) -> np.ndarray:
    """Wrap a displacement into [-L/2, L/2) per axis. Orthogonal cells only."""
    return delta - lengths * np.round(delta / lengths)


def nearest_distance_fn(positions: np.ndarray, lengths: np.ndarray):
    """
    Build a fast "how close is this point to anything" function for one
    structure: given (M, 3) query points it returns (M,) minimum-image distances
    to the nearest of `positions`.

    lammps_data.nearest_existing answers the same question and also says *which*
    atom, which is what a per-frame report needs. This is for the places that
    ask hundreds of thousands of times and only want the number - solving a path
    length, or scanning a path for its tightest point - where going through
    pandas for each query dominates the runtime.

    Orthogonal cells only, like everything else here.
    """
    positions = np.asarray(positions, dtype=float)
    lengths = np.asarray(lengths, dtype=float)

    def nearest(points: np.ndarray) -> np.ndarray:
        points = np.atleast_2d(np.asarray(points, dtype=float))
        delta = positions[None, :, :] - points[:, None, :]
        delta -= lengths * np.round(delta / lengths)
        return np.linalg.norm(delta, axis=2).min(axis=1)

    return nearest


def clearance_along(points: np.ndarray, origin: np.ndarray,
                    axis: np.ndarray) -> float:
    """
    Public form of the aperture measure: the distance from the given line to
    the nearest of the given atom centres. Used to re-measure a ring's opening
    once the charge-balancing oxygens have been added to it.
    """
    return _clearance(points, origin, np.asarray(axis, dtype=float))


# -- plane and axis -----------------------------------------------------------

def best_fit_plane(points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Least-squares plane through the points.

    Returns (centroid, normal, e1, e2), where e1 and e2 span the plane and are
    ordered by how much the points spread along them, so e1 is the ring's long
    axis. The normal is sign-fixed by its largest component, so two runs on the
    same ring agree on which way is "through".
    """
    centroid = points.mean(axis=0)
    _, _, vt = np.linalg.svd(points - centroid)
    e1, e2, normal = vt[0], vt[1], vt[2]
    if normal[np.argmax(np.abs(normal))] < 0:
        normal = -normal
    return centroid, normal, e1, e2


def planarity_rmsd(points: np.ndarray) -> float:
    """RMS distance of the atoms from their best-fit plane: how puckered the ring is."""
    centroid, normal, _, _ = best_fit_plane(points)
    return float(np.sqrt((((points - centroid) @ normal) ** 2).mean()))


def eccentricity(points: np.ndarray) -> float:
    """
    How far from circular the ring is in projection: the ratio of the two
    in-plane principal spreads. 1.0 is a circle; a larger value is a slot.
    """
    centroid, _, e1, e2 = best_fit_plane(points)
    local = points - centroid
    in_plane = np.column_stack([local @ e1, local @ e2])
    spreads = np.sort(np.linalg.eigvalsh(np.cov(in_plane.T)))[::-1]
    return float(np.sqrt(spreads[0] / max(spreads[1], 1e-12)))


def distances_to_line(points: np.ndarray, origin: np.ndarray,
                      axis: np.ndarray) -> np.ndarray:
    """Perpendicular distance from each atom centre to the line (origin, axis)."""
    local = points - np.asarray(origin, dtype=float)
    axis = np.asarray(axis, dtype=float)
    return np.linalg.norm(local - np.outer(local @ axis, axis), axis=1)


def _clearance(points: np.ndarray, origin: np.ndarray, axis: np.ndarray) -> float:
    """
    The aperture of a line: the distance from it to the nearest atom centre.

    Atoms are points here. Nothing is assumed about their size, so this is a
    statement about the positions alone; see contact_violations for the
    hard-sphere reading of the same geometry.
    """
    return float(distances_to_line(points, origin, axis).min())


def contact_violations(points: np.ndarray, types: np.ndarray,
                       origin: np.ndarray, axis: np.ndarray,
                       radii: dict[int, float] | None = None) -> dict:
    """
    Which atoms sit closer to the line than their own contact radius - the
    approaches a hard-sphere picture would have called an overlap.

    Reported rather than acted on. The aperture is measured between the axis and
    the atom centres and does not know these radii exist; this is the separate
    note that says where a finite atom size would have started to matter.
    """
    table = CONTACT_RADII if radii is None else radii
    distances = distances_to_line(points, origin, axis)
    limits = np.array([table[int(t)] for t in types], dtype=float)
    overlap = limits - distances          # positive where the atom is inside its radius
    inside = overlap > 0
    return {
        "n_inside_contact_radius": int(inside.sum()),
        "n_si_inside": int((inside & (types == SI_TYPE)).sum()),
        "n_o_inside": int((inside & (types == O_TYPE)).sum()),
        "worst_contact_overlap": float(overlap.max()) if inside.any() else 0.0,
        "nearest_atom_type": int(types[int(distances.argmin())]),
    }


def _optimise(objective, starts: np.ndarray, bounds: list[tuple[float, float]]
              ) -> tuple[float, np.ndarray]:
    """
    Maximise a min-of-distances objective over a box.

    The objective is continuous but not smooth - it is a minimum over atoms, so
    it creases wherever the closest atom changes - which rules out gradient
    methods and leaves Nelder-Mead from several starts.

    The bounds are not optional. Unbounded, the simplex walks off to infinity in
    the flat region outside the feasible set, the coordinates overflow, and the
    objective evaluates to inf, which then wins: the optimiser returns a
    "perfect" aperture for a line nowhere near the ring. Bounding the box keeps
    every evaluation finite and meaningful.
    """
    best_value, best_params = -np.inf, starts[0]
    for guess in starts:
        guess = np.clip(guess, [b[0] for b in bounds], [b[1] for b in bounds])
        result = minimize(lambda p: -objective(p), guess, method="Nelder-Mead",
                          bounds=bounds,
                          options={"xatol": 1e-4, "fatol": 1e-6, "maxiter": 2000})
        value = float(-result.fun)
        if np.isfinite(value) and value > best_value:
            best_value, best_params = value, result.x
    if not np.isfinite(best_value):
        raise ValueError("aperture optimisation produced no finite result")
    return best_value, best_params


def centroid_aperture(ring: Ring) -> tuple[float, np.ndarray, np.ndarray]:
    """
    The largest probe that passes along the line through the ring's centroid in
    the direction of its best-fit plane normal. No optimisation.

    Returns (aperture, center, axis). This is the default path axis: it is
    reproducible, has nothing to tune, and is what every other aperture here is
    measured against.
    """
    points = ring.positions
    centroid, normal, _, _ = best_fit_plane(points)
    return _clearance(points, centroid, normal), centroid, normal


def plane_aperture(ring: Ring) -> tuple[float, np.ndarray, np.ndarray]:
    """
    The largest probe that passes along the best-fit plane normal, optimising
    only where the line crosses the plane.

    Returns (aperture, center, axis). The axis is the plane normal.
    """
    points = ring.positions
    centroid, normal, e1, e2 = best_fit_plane(points)
    spread = float(np.linalg.norm(points - centroid, axis=1).max())

    def objective(p):
        return _clearance(points, centroid + p[0] * e1 + p[1] * e2, normal)

    starts = np.array([[0.0, 0.0], [0.3, 0.0], [-0.3, 0.0],
                       [0.0, 0.3], [0.0, -0.3]]) * spread
    bounds = [(-spread, spread)] * 2
    value, params = _optimise(objective, starts, bounds)
    return value, centroid + params[0] * e1 + params[1] * e2, normal


def cylinder_aperture(ring: Ring, max_tilt_deg: float = MAX_TILT_DEG
                      ) -> tuple[float, np.ndarray, np.ndarray, bool]:
    """
    The largest probe that passes through the ring along any straight line whose
    direction lies within `max_tilt_deg` of the plane normal: the maximal
    inscribed cylinder, cone-constrained (see MAX_TILT_DEG for why the cone has
    to be there).

    Returns (aperture, center, axis, at_wall), where center is the point on the
    axis closest to the ring's centroid - the natural origin for a path
    coordinate - and at_wall says the optimum sits on the cone boundary, so the
    ring wanted to tilt further than it was allowed.

    Seeded from the plane-normal optimum and never returns less than it, so
    aperture_cyl >= aperture_plane holds by construction rather than by luck.
    """
    points = ring.positions
    centroid, normal, e1, e2 = best_fit_plane(points)
    spread = float(np.linalg.norm(points - centroid, axis=1).max())
    max_tilt = float(np.tan(np.radians(max_tilt_deg)))

    def unpack(p):
        # Project the tilt onto the disc so the reachable set is a true cone.
        # Clipping the components independently would make it a square, with a
        # diagonal reaching sqrt(2) times as far as the stated half-angle.
        tilt = np.asarray(p[:2], dtype=float)
        length = np.linalg.norm(tilt)
        if length > max_tilt:
            tilt = tilt * (max_tilt / length)
        axis = normal + tilt[0] * e1 + tilt[1] * e2
        return centroid + p[2] * e1 + p[3] * e2, axis / np.linalg.norm(axis), length

    def objective(p):
        origin, axis, _ = unpack(p)
        return _clearance(points, origin, axis)

    flat_value, flat_center, _ = plane_aperture(ring)
    offset = np.array([np.dot(flat_center - centroid, e1),
                       np.dot(flat_center - centroid, e2)])

    starts = np.zeros((7, 4))
    starts[:, 2:] = offset                      # every start begins at the plane optimum
    starts[1:3, 0] = [0.5 * max_tilt, -0.5 * max_tilt]
    starts[3:5, 1] = [0.5 * max_tilt, -0.5 * max_tilt]
    starts[5, 2] += 0.25 * spread
    starts[6, 3] += 0.25 * spread
    bounds = [(-max_tilt, max_tilt)] * 2 + [(-spread, spread)] * 2

    value, params = _optimise(objective, starts, bounds)
    if value < flat_value:
        # The plane normal is inside the cone, so it is always available; if the
        # search somehow did worse, keep the answer we already had.
        return flat_value, flat_center, normal, False

    origin, axis, tilt_length = unpack(params)
    # Slide the origin to the point on the axis nearest the centroid, so the
    # path coordinate is zero at the ring rather than at an arbitrary offset.
    center = origin + np.dot(centroid - origin, axis) * axis
    return value, center, axis, bool(tilt_length >= 0.999 * max_tilt)


# -- the whole picture --------------------------------------------------------

AXIS_CHOICES = ("centroid", "slid", "tilted")


def ring_metrics(ring: Ring, radii: dict[int, float] | None = None,
                 axis_choice: str = "centroid") -> dict:
    """
    Every geometric number this module knows how to compute, as a flat dict.

    All three apertures are always reported, because the gap between them is
    itself the interesting number - it says how far the ring's atom-average sits
    from its actual opening. `axis_choice` only decides which line is written
    out as `center`/`axis` for a path generator to use.

    `radii` reaches only the contact flags, which are measured on whichever axis
    `axis_choice` selected. No aperture depends on it.
    """
    if axis_choice not in AXIS_CHOICES:
        raise ValueError(f"axis_choice must be one of {AXIS_CHOICES}, not {axis_choice!r}")

    points = ring.positions
    at_centroid, centroid, normal = centroid_aperture(ring)
    slid, slid_center, _ = plane_aperture(ring)
    tilted, tilted_center, tilted_axis, at_wall = cylinder_aperture(ring)

    # Each search contains the previous one's answer among its candidates, so
    # this ordering cannot fail unless an optimiser is broken. One has been, once.
    if slid < at_centroid - 1e-6 or tilted < slid - 1e-6:
        raise ValueError(
            f"{ring.source}: apertures are not ordered (centroid {at_centroid:.4f}, "
            f"slid {slid:.4f}, tilted {tilted:.4f})")

    center, axis = {
        "centroid": (centroid, normal),
        "slid": (slid_center, normal),
        "tilted": (tilted_center, tilted_axis),
    }[axis_choice]
    aperture = {"centroid": at_centroid, "slid": slid, "tilted": tilted}[axis_choice]

    return {
        "n": ring.n,
        "n_atoms": len(ring.ids),
        "axis_choice": axis_choice,
        "aperture": aperture,
        "aperture_centroid": at_centroid,
        "aperture_slid": slid,
        "aperture_tilted": tilted,
        "hole_offset": float(np.linalg.norm(slid_center - centroid)),
        "tilt_deg": float(np.degrees(np.arccos(np.clip(abs(np.dot(tilted_axis, normal)), 0, 1)))),
        "axis_at_cone_wall": int(at_wall),
        "planarity_rmsd": planarity_rmsd(points),
        "eccentricity": eccentricity(points),
        # Two radii, because "the radius of the ring" is ambiguous and the
        # difference is not random. mean_radius is the plain 3-D distance from
        # the centroid; mean_radius_inplane projects onto the best-fit plane
        # first, which is what "radius" means for a ring and what keeps pucker
        # out of it. The 3-D version runs larger by 0.005 A at n=3 and 0.085 A
        # at n=8 - a bias that grows with ring size, so prefer the in-plane one
        # for anything plotted or correlated against size.
        #
        # Both average silicon and oxygen together. They are not at the same
        # radius: silicon sits 0.26-0.38 A further out, so either number is
        # halfway between the two shells rather than on one of them.
        "mean_radius": float(np.linalg.norm(points - centroid, axis=1).mean()),
        "mean_radius_inplane": float(np.linalg.norm(
            (points - centroid) - np.outer((points - centroid) @ normal, normal),
            axis=1).mean()),
        "mean_radius_si": float(np.linalg.norm(
            points[ring.types == SI_TYPE] - centroid, axis=1).mean()),
        "mean_radius_o": float(np.linalg.norm(
            points[ring.types == O_TYPE] - centroid, axis=1).mean()),
        **contact_violations(points, ring.types, center, axis, radii),
        "centroid_x": centroid[0], "centroid_y": centroid[1], "centroid_z": centroid[2],
        "normal_x": normal[0], "normal_y": normal[1], "normal_z": normal[2],
        "center_x": center[0], "center_y": center[1], "center_z": center[2],
        "axis_x": axis[0], "axis_y": axis[1], "axis_z": axis[2],
    }
