#!/usr/bin/env python3
"""
water_geometry.py

The geometry of one rigid water molecule: building its two O->H offsets, aiming
them along a path, and rolling them about that path.

Nothing here knows what a file is. These are vectors in and vectors out, which
is why they are shared rather than duplicated - generate_water_path.py places
water into LAMMPS data files and VASP POSCARs through structure_io, while
ring_path places it into LAMMPS clusters through lammps_data, and both need the
same molecule built the same way.

The format-dependent half - reading a template molecule out of a file, resolving
atom types, writing a frame, reporting what the water nearly hit - lives with
whichever writer is doing it, not here.
"""

from __future__ import annotations

import numpy as np

# The three atoms a frame adds, in the order they are written and reported.
WATER_LABELS = ("O", "H1", "H2")

def normalize(vector: np.ndarray) -> np.ndarray:
    vector = np.asarray(vector, dtype=float)
    length = np.linalg.norm(vector)
    if length < 1e-9:
        raise ValueError("cannot normalize a zero-length vector")
    return vector / length


def perpendicular_to(direction: np.ndarray, hint: np.ndarray | None = None) -> np.ndarray:
    """
    A unit vector perpendicular to `direction`, deterministically chosen so that
    two runs with the same inputs produce the same water orientation. `hint`
    picks which perpendicular, when the caller cares; a hint parallel to
    `direction` is ignored.
    """
    direction = normalize(direction)
    if hint is not None:
        residual = np.asarray(hint, dtype=float)
        residual = residual - np.dot(residual, direction) * direction
        if np.linalg.norm(residual) > 1e-6:
            return normalize(residual)
    axis = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(direction, axis)) > 0.9:
        axis = np.array([0.0, 1.0, 0.0])
    return normalize(np.cross(direction, axis))


def canonical_offsets(oh_length: float, hoh_angle_deg: float) -> tuple[np.ndarray, np.ndarray]:
    """Two O->H offsets in a local frame, bisector along +z and both H in the xz plane."""
    half = np.radians(hoh_angle_deg) / 2.0
    bisector = np.array([0.0, 0.0, 1.0])
    spread = np.array([1.0, 0.0, 0.0])
    h1 = oh_length * (np.cos(half) * bisector + np.sin(half) * spread)
    h2 = oh_length * (np.cos(half) * bisector - np.sin(half) * spread)
    return h1, h2


def orient_offsets(h1: np.ndarray, h2: np.ndarray, bisector: np.ndarray,
                   in_plane_hint: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """
    Rigidly rotate a pair of O->H offsets so their bisector points along
    `bisector`, preserving both bond lengths and the H-O-H angle exactly. The
    molecular plane is fixed by `in_plane_hint`, whose component perpendicular
    to the bisector becomes the direction the two hydrogens spread along.
    """
    h1 = np.asarray(h1, dtype=float)
    h2 = np.asarray(h2, dtype=float)

    source_bisector = normalize(h1) + normalize(h2)
    if np.linalg.norm(source_bisector) < 1e-6:
        raise ValueError("the two O-H offsets are anti-parallel; there is no bisector")
    source_bisector = normalize(source_bisector)
    source_spread = normalize(h1 - np.dot(h1, source_bisector) * source_bisector)

    target_bisector = normalize(bisector)
    target_spread = perpendicular_to(target_bisector, in_plane_hint)

    return tuple(
        np.dot(offset, source_bisector) * target_bisector
        + np.dot(offset, source_spread) * target_spread
        for offset in (h1, h2)
    )


def rotate_about_bisector(h1: np.ndarray, h2: np.ndarray,
                          degrees: float) -> tuple[np.ndarray, np.ndarray]:
    """
    Spin the two O->H offsets about their own H-O-H bisector by `degrees`.

    The bisector is the rotation axis, so it does not move: the oxygen keeps
    facing wherever the orientation aimed it, and what turns is the plane the
    two hydrogens lie in. Under --orientation away/toward the bisector is the
    line of approach, so this is the free parameter that nothing else fixes -
    which way the water is rolled around it.

    The rotation is counterclockwise about the bisector, seen looking down the
    bisector toward the oxygen (equivalently, the right-hand rule with the thumb
    along the bisector). It is always counterclockwise: the angle is reduced
    into [0, 360) first, so a negative value becomes the counterclockwise turn
    that lands in the same place rather than a clockwise one.

    Bond lengths and the H-O-H angle are untouched - this is a rigid rotation.
    """
    angle = np.radians(float(degrees) % 360.0)
    axis = normalize(normalize(h1) + normalize(h2))

    def rotated(vector):
        # Rodrigues' rotation formula about the unit `axis`.
        return (vector * np.cos(angle)
                + np.cross(axis, vector) * np.sin(angle)
                + axis * np.dot(axis, vector) * (1.0 - np.cos(angle)))

    return rotated(np.asarray(h1, dtype=float)), rotated(np.asarray(h2, dtype=float))


def describe_offsets(h1: np.ndarray, h2: np.ndarray) -> str:
    """The bond lengths and angle a pair of offsets actually came out at."""
    angle = np.degrees(np.arccos(np.dot(normalize(h1), normalize(h2))))
    return (f"O-H {np.linalg.norm(h1):.4f} / {np.linalg.norm(h2):.4f} A, "
            f"H-O-H {angle:.2f} deg")


def resolve_offsets(approach: np.ndarray, *,
                    orientation: str = "away",
                    h_rotation: float = 0.0,
                    oh_length: float = 0.9572,
                    hoh_angle: float = 104.52,
                    offsets: tuple[np.ndarray, np.ndarray] | None = None
                    ) -> tuple[np.ndarray, np.ndarray, str]:
    """
    The two rigid O->H offsets for a whole path. `approach` is the unit vector
    pointing the way the oxygen travels.

    `offsets` supplies a molecule measured from somewhere else - a template file,
    say - in place of building one from `oh_length` and `hoh_angle`. Reading that
    file is the caller's job, because the file format is not this module's
    business.

    Returns the offsets and a sentence describing how they were arrived at, for
    the caller to print.
    """
    if offsets is not None:
        h1, h2 = np.asarray(offsets[0], dtype=float), np.asarray(offsets[1], dtype=float)
        source = "from the given template"
    else:
        h1, h2 = canonical_offsets(oh_length, hoh_angle)
        source = f"built from --oh-length {oh_length} and --hoh-angle {hoh_angle}"

    if orientation == "away":
        # Bisector anti-parallel to the approach: the hydrogens trail behind the
        # oxygen and its lone pairs face forward, which is the geometry a
        # nucleophilic attack starts from.
        h1, h2 = orient_offsets(h1, h2, -approach)
        source += ", hydrogens pointing away from the target"
    elif orientation == "toward":
        # Bisector along the approach: the hydrogens lead and arrive before the
        # oxygen does. Use this when the proton, not the oxygen, should arrive
        # first.
        h1, h2 = orient_offsets(h1, h2, approach)
        source += ", hydrogens pointing toward the target"
    else:
        # Molecular plane perpendicular to the approach: both hydrogens stay off
        # the line of travel, which keeps them clear of a tight channel.
        bisector = perpendicular_to(approach)
        h1, h2 = orient_offsets(h1, h2, bisector,
                                in_plane_hint=np.cross(approach, bisector))
        source += ", molecular plane perpendicular to the path"

    # Applied last, and to every source of offsets, so a roll always means the
    # same thing: a turn about the bisector the orientation just set. The
    # bisector is the axis, so this cannot undo the orientation above.
    rotation = float(h_rotation) % 360.0
    if rotation:
        h1, h2 = rotate_about_bisector(h1, h2, rotation)
        source += f", rolled {rotation:g} deg counterclockwise about the bisector"

    return h1, h2, source
