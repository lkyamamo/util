#!/usr/bin/env python3
"""
bond_events.py — bond breaking/forming events in bulk water, with a
flicker-immune criterion, a mechanistic classification of every event, and an
OVITO-ingestible annotated trajectory.

QUICK START
-----------
1. Run the analysis:
       python bond_events.py --traj dynamics.lammpstrj --dt 2.5
2. Open output/bond_events.lammpstrj in OVITO and add an Expression Selection:
       EventType != 0

The default cutoffs are stated below with their justification. Deriving them
automatically from the run's own g(r) is parked as future work — see
calibrate_cutoffs.py, which is a sketch and is NOT wired into this tool.

WHY NOT A PLAIN CUTOFF
----------------------
The O-H stretch is ~3700 cm^-1, period ~9 fs. dynamics.lammpstrj is dumped
every 10 steps of a 0.25 fs timestep, i.e. dt = 2.5 fs — about 3.6 samples per
stretch period. A single hard cutoff has no hysteresis, so any pair whose
vibration straddles it "breaks" and "re-forms" every couple of frames. That
does not matter for an equilibrium average (g(r), bond angles, silanol counts),
where the flicker averages out; it completely dominates an *event count*, which
would then scale with dump frequency rather than with chemistry.

So bonds are tracked with a Schmitt trigger plus a persistence requirement:

    BONDED   -> UNBONDED  when r > r_break for >= tau_persist consecutive frames
    UNBONDED -> BONDED    when r < r_form  for >= tau_persist consecutive frames

with r_form < r_break. Excursions shorter than tau_persist are counted
separately as `recrossings` — that number is exactly how much a naive
single-cutoff analysis would have over-reported, and it is reported as a
first-class result rather than silently discarded.

An event is *timestamped* at the frame the threshold was first crossed but only
*emitted* once persistence is confirmed, so the tool runs tau_persist frames
behind. Crossings still unconfirmed at end-of-trajectory are reported as
`pending_at_eof` rather than dropped silently.

MECHANISM ("why did it break?")
-------------------------------
Every confirmed break is classified by what the H does within --transfer-window:

  transfer      H bonds to a *different* O — Grotthuss proton hop. Expected to
                dominate overwhelmingly in equilibrium bulk water.
  recross       H re-bonds to the *same* O after a confirmed unbonded period.
  dissociation  H bonds to no O for the whole window. Rare; a large fraction
                means the cutoffs are wrong, not that the water is reacting.
  truncated     the window ran past the end of the trajectory.

and carries descriptors that say *why* it happened:

  r_at_break            O-H distance at the crossing frame
  v_radial              relative O/H velocity along the bond (A/ps); large and
                        positive means the bond was kinetically driven apart
  ke_h, ke_h_over_kT    kinetic energy of the H and its ratio to (3/2)kT — are
                        breaks coming from the hot tail of the distribution?
  coord_o_donor_before/after, coord_o_acceptor_before/after
                        donor 2->1 with acceptor 2->3 is the OH- / H3O+ signature
  n_o_neighbors_donor   O...O neighbours of the donor within --hbond-cutoff:
                        was the donor over-solvated?
  oo_distance_at_break, ohO_angle_at_break, oo_distance_at_form
                        the classic Grotthuss reaction coordinate. A hop at
                        O...O ~ 2.4 A with a near-linear O-H...O angle is
                        textbook; one at 3.2 A means the cutoffs are wrong.

OUTPUT (all under output/, or --output-dir)
-------------------------------------------
  bond_events.lammpstrj    input trajectory + per-atom annotation columns
                           (Coordination EventType EventAge PartnerID Species)
  events.jsonl             one JSON object per event, every descriptor
  events_per_frame.csv     per-frame break/form counts and species populations
  summary.txt              aggregate report
  bond_events.png          diagnostic panels

COLUMN LAYOUT
-------------
Expects a LAMMPS custom dump with at least `id x y z` and either `element` or
`type` (with --type-map). Column positions are read from the ITEM: ATOMS header,
so column order does not matter. `vx vy vz` are optional; without them the
velocity-derived descriptors are omitted. Atom identity is taken from the `id`
column rather than from row order.

DEPENDENCIES
------------
  pip install numpy freud matplotlib
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from collections import deque
from dataclasses import dataclass
from itertools import islice
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

import numpy as np
import freud

# --------------------------------------------------------------------------
# Physical constants (LAMMPS metal units: distance A, time ps, velocity A/ps)
# --------------------------------------------------------------------------

# 1 (g/mol) * (A/ps)^2 = 1.0364269e-4 eV
AMU_ANG2_PS2_TO_EV = 1.0364269e-4
BOLTZMANN_EV_PER_K = 8.617333262e-5

DEFAULT_MASSES = {"H": 1.00784, "O": 15.9994, "Si": 28.0855}

# Per-atom EventType codes written into the annotated dump.
EV_NONE = 0
EV_O_BREAK = 1
EV_H_BREAK = 2
EV_O_FORM = 3
EV_H_FORM = 4

# Per-atom Species codes written into the annotated dump.
SP_WATER = 0
SP_HYDRONIUM = 1
SP_HYDROXIDE = 2
SP_FREE_H = 3
SP_FREE_O = 4
SP_OTHER = 5

# ==========================================================================
# Trajectory reading
# ==========================================================================


@dataclass
class Frame:
    """One trajectory frame. `index` counts kept frames, not file frames."""

    index: int
    timestep: int
    pos: np.ndarray                  # (n_atoms, 3)
    vel: Optional[np.ndarray]        # (n_atoms, 3) or None
    box_L: np.ndarray                # (3,) edge lengths
    box_lo: np.ndarray               # (3,) lower bounds


@dataclass
class Topology:
    """Atom identity, fixed for the whole trajectory."""

    ids: np.ndarray                  # (n_atoms,) int64
    elements: np.ndarray             # (n_atoms,) unicode
    o_rows: np.ndarray               # (n_o,) indices of oxygens
    h_rows: np.ndarray               # (n_h,) indices of hydrogens
    has_velocities: bool

    @property
    def n_atoms(self) -> int:
        return len(self.ids)


def _parse_type_map(raw: Optional[str]) -> Dict[int, str]:
    """Parse '1:O;2:H' into {1: 'O', 2: 'H'}.

    Semicolon-separated rather than comma-separated so the same string survives
    `sbatch --export=...`, which is comma-delimited — the convention already
    used by bad_freud.py's _parse_pair_dict.
    """
    if not raw:
        return {}
    out: Dict[int, str] = {}
    for entry in raw.split(";"):
        entry = entry.strip()
        if not entry:
            continue
        key, value = entry.split(":")
        out[int(key)] = value.strip()
    return out


def _resolve_columns(header: Sequence[str], type_map: Dict[int, str]) -> Dict[str, int]:
    """Map required column names to their positions in the ITEM: ATOMS header."""
    cols = list(header[2:])  # drop 'ITEM:' 'ATOMS'
    index = {name: i for i, name in enumerate(cols)}

    missing = [c for c in ("id", "x", "y", "z") if c not in index]
    if missing:
        raise ValueError(
            f"Dump is missing required column(s) {missing}. bond_events.py needs "
            f"wrapped positions and atom ids (dump ... id element x y z [vx vy vz]). "
            f"Found columns: {cols}"
        )
    if "element" not in index and "type" not in index:
        raise ValueError(
            f"Dump has neither an 'element' nor a 'type' column; one is needed to "
            f"tell O from H. Found columns: {cols}"
        )
    if "element" not in index and not type_map:
        raise ValueError(
            "Dump has no 'element' column, only 'type'. Pass --type-map to say "
            "which numeric type is which element, e.g. --type-map '1:O;2:H'."
        )

    resolved = {
        "id": index["id"],
        "x": index["x"],
        "y": index["y"],
        "z": index["z"],
        "n_cols": len(cols),
    }
    resolved["species"] = index.get("element", index.get("type"))
    resolved["species_is_type"] = "element" not in index
    if all(c in index for c in ("vx", "vy", "vz")):
        resolved["vx"] = index["vx"]
        resolved["vy"] = index["vy"]
        resolved["vz"] = index["vz"]
    return resolved


def read_frames(
    path: Path,
    *,
    start: int = 0,
    stop: Optional[int] = None,
    stride: int = 1,
    type_map: Optional[Dict[int, str]] = None,
) -> Iterator[Tuple[Frame, Topology]]:
    """
    Stream frames from a LAMMPS custom dump.

    Yields (Frame, Topology). Topology is the same object every time; it is
    built from the first frame and validated against every later one, since
    keying bonds by atom id only works if the id set is stable.

    Atom blocks are tokenized in one shot and converted column-wise with numpy
    rather than looped over in Python — at 26k atoms x 5000 frames the per-atom
    loop used by msd.py would dominate the runtime.
    """
    type_map = type_map or {}
    topology: Optional[Topology] = None
    columns: Optional[Dict[str, int]] = None
    file_frame = 0
    kept = 0

    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                raise ValueError(
                    f"Expected 'ITEM: TIMESTEP' at frame {file_frame} of {path}, got {line!r}"
                )
            timestep = int(f.readline())

            f.readline()                       # ITEM: NUMBER OF ATOMS
            n_atoms = int(f.readline())

            f.readline()                       # ITEM: BOX BOUNDS ...
            lo = np.empty(3)
            hi = np.empty(3)
            for d in range(3):
                bounds = f.readline().split()
                lo[d], hi[d] = float(bounds[0]), float(bounds[1])

            header = f.readline().split()      # ITEM: ATOMS ...
            if columns is None:
                columns = _resolve_columns(header, type_map)

            in_range = file_frame >= start and (stop is None or file_frame < stop)
            keep = in_range and ((file_frame - start) % stride == 0)

            if not keep:
                # Consume the atom block without parsing it.
                for _ in islice(f, n_atoms):
                    pass
                file_frame += 1
                if stop is not None and file_frame >= stop:
                    break
                continue

            block = "".join(islice(f, n_atoms))
            tokens = np.array(block.split()).reshape(n_atoms, columns["n_cols"])

            ids = tokens[:, columns["id"]].astype(np.int64)
            pos = tokens[:, [columns["x"], columns["y"], columns["z"]]].astype(np.float64)
            if "vx" in columns:
                vel = tokens[:, [columns["vx"], columns["vy"], columns["vz"]]].astype(np.float64)
            else:
                vel = None

            if topology is None:
                raw_species = tokens[:, columns["species"]]
                if columns["species_is_type"]:
                    types = raw_species.astype(np.int64)
                    unknown = sorted(set(types.tolist()) - set(type_map))
                    if unknown:
                        raise ValueError(
                            f"Dump contains atom type(s) {unknown} not covered by "
                            f"--type-map {type_map}."
                        )
                    elements = np.array([type_map[t] for t in types])
                else:
                    elements = raw_species.astype("<U4")
                topology = Topology(
                    ids=ids,
                    elements=elements,
                    o_rows=np.flatnonzero(elements == "O"),
                    h_rows=np.flatnonzero(elements == "H"),
                    has_velocities=vel is not None,
                )
                if len(topology.o_rows) == 0 or len(topology.h_rows) == 0:
                    raise ValueError(
                        f"Found {len(topology.o_rows)} O and {len(topology.h_rows)} H "
                        f"in {path}. bond_events.py analyses O-H bonds in water."
                    )
            elif not np.array_equal(topology.ids, ids):
                raise ValueError(
                    "Atom id order changed between frames. Bond tracking keys pairs by "
                    "atom id and requires a stable ordering — check 'dump_modify ... "
                    "sort id' is set in the LAMMPS input."
                )

            yield Frame(
                index=kept,
                timestep=timestep,
                pos=pos,
                vel=vel,
                box_L=hi - lo,
                box_lo=lo,
            ), topology

            kept += 1
            file_frame += 1
            if stop is not None and file_frame >= stop:
                break

    if topology is None:
        raise ValueError(f"No frames read from {path} (empty file, or --frames/--stride excluded everything).")


# ==========================================================================
# Geometry helpers
# ==========================================================================


def minimum_image(delta: np.ndarray, box_L: np.ndarray) -> np.ndarray:
    """Minimum-image a displacement vector (or array of them).

    Same one-liner as msd.py:234; np.round handles any number of wraps at once.
    """
    return delta - box_L * np.round(delta / box_L)


def _freud_box_and_points(frame: Frame) -> Tuple[freud.box.Box, np.ndarray]:
    """
    Build a freud box and recentre positions onto it.

    freud boxes are centred on the origin, so positions have to be shifted by
    the box centre before querying — the same recentring bad_freud.py does.
    Using freud (rather than a bare cKDTree) also means periodicity is applied
    by construction; voxel_analysis.py:294 builds a cKDTree with no `boxsize=`
    and is silently non-periodic as a result.
    """
    box = freud.box.Box(Lx=frame.box_L[0], Ly=frame.box_L[1], Lz=frame.box_L[2])
    centred = frame.pos - (frame.box_lo + frame.box_L / 2.0)
    return box, box.wrap(centred.astype(np.float32))


def _oh_neighbors(
    box: freud.box.Box,
    points: np.ndarray,
    topology: Topology,
    r_query: float,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    All O-H pairs within r_query. Returns (o_rows, h_rows, distances) as
    global atom-row indices.

    The AABBQuery is built over the H atoms and probed with the O atoms — the
    same orientation _run_query uses in bad_freud.py. exclude_ii must be False
    because the two point sets are different arrays, so matching indices are
    unrelated atoms.
    """
    o_pts = points[topology.o_rows]
    h_pts = points[topology.h_rows]
    nlist = (
        freud.locality.AABBQuery(box, h_pts)
        .query(o_pts, {"r_max": r_query, "exclude_ii": False})
        .toNeighborList()
    )
    return (
        topology.o_rows[nlist.query_point_indices],
        topology.h_rows[nlist.point_indices],
        np.asarray(nlist.distances, dtype=np.float64),
    )


def _oo_neighbor_counts(
    box: freud.box.Box,
    points: np.ndarray,
    topology: Topology,
    hbond_cutoff: float,
) -> np.ndarray:
    """Number of O...O neighbours within hbond_cutoff, indexed by O position in
    topology.o_rows."""
    o_pts = points[topology.o_rows]
    nlist = (
        freud.locality.AABBQuery(box, o_pts)
        .query(o_pts, {"r_max": hbond_cutoff, "exclude_ii": True})
        .toNeighborList()
    )
    return np.bincount(nlist.query_point_indices, minlength=len(topology.o_rows))


# ==========================================================================
# Pair state machine
# ==========================================================================


@dataclass
class PairState:
    """Hysteretic state of one O-H pair."""

    bonded: bool
    pending_start: int = -1        # frame the current candidate crossing began
    pending_r: float = 0.0         # O-H distance at that frame
    last_seen: int = -1            # last frame this pair was inside r_query


@dataclass
class BreakEvent:
    """A confirmed bond break, awaiting mechanistic classification."""

    frame: int
    timestep: int
    o_row: int
    h_row: int
    r_at_break: float
    descriptors: dict
    classification: str = "pending"
    acceptor_row: int = -1
    acceptor_coord_before: int = -1
    acceptor_coord_after: int = -1
    form_frame: int = -1


class BondTracker:
    """
    Schmitt-trigger + persistence bond state machine over a streamed trajectory.

    Only pairs that have come within r_query are tracked, so the state dict
    stays proportional to the number of O-H contacts (~2.2 per O in water),
    not to N^2.
    """

    def __init__(self, args: argparse.Namespace, topology: Topology):
        self.args = args
        self.topology = topology
        self.r_form = args.r_form
        self.r_break = args.r_break
        # The neighbour search only has to decide which side of r_break a pair
        # is on: anything outside r_query is by definition beyond r_break and is
        # handled by the fallback below. So the radius is kept just above
        # r_break rather than out at the hydrogen-bond shell — at ~1.6 A each O
        # has ~2 neighbours instead of ~6, and the per-pair state machine is the
        # runtime bottleneck. Distances used in descriptors are recomputed from
        # positions, so this costs no accuracy.
        self.r_query = args.r_break * 1.2
        self.tau = args.tau_persist
        self.transfer_window = args.transfer_window

        self.pairs: Dict[Tuple[int, int], PairState] = {}
        # Bond set as of frame 0, used by pass 2 to replay coordination. Burn-in
        # transitions mutate it instead of emitting events, so it always reflects
        # the state the first reportable frame actually starts from.
        self.initial_bonds: set = set()
        self.events: List[dict] = []
        self.pending_breaks: List[BreakEvent] = []

        # Atom row -> position within topology.o_rows. Fixed for the whole
        # trajectory, so it is built once rather than per frame.
        self.o_index: Dict[int, int] = {int(row): i for i, row in enumerate(topology.o_rows)}

        # Confirmed coordination number per atom, kept in step with the state
        # machine so every event record can carry the donor/acceptor before and
        # after counts without cross-referencing the annotated trajectory.
        self.coord = np.zeros(topology.n_atoms, dtype=np.int32)

        self.n_recrossings = 0        # sub-tau excursions above r_break while bonded
        self.n_failed_forms = 0       # sub-tau excursions below r_form while unbonded
        self.n_suppressed_burnin = 0
        self.per_frame: List[dict] = []

        # Ring buffer of recent frames, deep enough to reach back from a
        # transfer confirmation to the originating break frame.
        self.buffer: deque = deque(maxlen=self.tau + self.transfer_window + 2)
        self._initialised = False

    # -- helpers ---------------------------------------------------------

    def _frame_at(self, index: int) -> Optional[Frame]:
        for frame in reversed(self.buffer):
            if frame.index == index:
                return frame
            if frame.index < index:
                break
        return None

    def _kt(self) -> float:
        return BOLTZMANN_EV_PER_K * self.args.temperature

    def _break_descriptors(
        self,
        frame: Frame,
        o_row: int,
        h_row: int,
        oo_counts: np.ndarray,
    ) -> dict:
        """Local environment of an O-H pair at the frame the bond broke."""
        d = minimum_image(frame.pos[h_row] - frame.pos[o_row], frame.box_L)
        r = float(np.linalg.norm(d))
        desc = {
            "r_at_break": round(r, 4),
            "n_o_neighbors_donor": int(oo_counts[self.o_index[o_row]]),
        }
        if frame.vel is not None:
            unit = d / r
            dv = frame.vel[h_row] - frame.vel[o_row]
            ke_h = 0.5 * DEFAULT_MASSES["H"] * float(np.dot(frame.vel[h_row], frame.vel[h_row]))
            ke_h *= AMU_ANG2_PS2_TO_EV
            desc["v_radial"] = round(float(np.dot(dv, unit)), 4)
            desc["ke_h"] = round(ke_h, 6)
            desc["ke_h_over_kT"] = round(ke_h / (1.5 * self._kt()), 4)
        return desc

    # -- main step -------------------------------------------------------

    def step(self, frame: Frame) -> None:
        self.buffer.append(frame)
        box, points = _freud_box_and_points(frame)
        o_rows, h_rows, dists = _oh_neighbors(box, points, self.topology, self.r_query)
        oo_counts = _oo_neighbor_counts(box, points, self.topology, self.args.hbond_cutoff)

        t = frame.index
        breaks_this_frame: List[Tuple[int, int, int]] = []
        forms_this_frame: List[Tuple[int, int, int]] = []

        if not self._initialised:
            # Seed the state machine: a pair starts bonded only if it is
            # unambiguously inside r_form. Pairs sitting in the hysteresis band
            # start unbonded and will form if they persist. Events during the
            # first tau frames are suppressed so this seeding cannot masquerade
            # as chemistry.
            for o_row, h_row, r in zip(o_rows, h_rows, dists):
                key = (int(o_row), int(h_row))
                bonded = bool(r < self.r_form)
                self.pairs[key] = PairState(bonded=bonded, last_seen=t)
                if bonded:
                    self.initial_bonds.add(key)
                    self.coord[o_row] += 1
                    self.coord[h_row] += 1
            self._initialised = True
            self._record_frame(frame, breaks_this_frame, forms_this_frame)
            self._advance_pending(t)
            return

        seen = set()
        for o_row, h_row, r in zip(o_rows, h_rows, dists):
            key = (int(o_row), int(h_row))
            seen.add(key)
            state = self.pairs.get(key)
            if state is None:
                # New contact: appears for the first time already inside r_query.
                state = PairState(bonded=False, last_seen=t)
                self.pairs[key] = state
            state.last_seen = t
            self._update(key, state, float(r), t, frame, breaks_this_frame,
                         forms_this_frame, oo_counts)

        # Pairs that were bonded but have left r_query entirely: their distance
        # is by definition > r_break, so they must keep advancing towards a break.
        for key, state in self.pairs.items():
            if key in seen or not state.bonded:
                continue
            self._update(key, state, float("inf"), t, frame, breaks_this_frame,
                         forms_this_frame, oo_counts)

        self._record_frame(frame, breaks_this_frame, forms_this_frame)
        self._advance_pending(t)

        if t % 200 == 0:
            self._prune(t)

    def _update(
        self,
        key: Tuple[int, int],
        state: PairState,
        r: float,
        t: int,
        frame: Frame,
        breaks_this_frame: List[Tuple[int, int, int]],
        forms_this_frame: List[Tuple[int, int, int]],
        oo_counts: np.ndarray,
    ) -> None:
        """Advance one pair's Schmitt trigger by one frame."""
        o_row, h_row = key
        if state.bonded:
            if r > self.r_break:
                if state.pending_start < 0:
                    state.pending_start = t
                    state.pending_r = r
                elif t - state.pending_start + 1 >= self.tau:
                    crossing = state.pending_start
                    state.bonded = False
                    state.pending_start = -1
                    if crossing < self.tau:
                        self.n_suppressed_burnin += 1
                        self.initial_bonds.discard(key)
                        self.coord[o_row] -= 1
                        self.coord[h_row] -= 1
                    else:
                        breaks_this_frame.append((crossing, o_row, h_row))
                        self._emit_break(crossing, o_row, h_row, oo_counts)
            elif state.pending_start >= 0:
                # Went back below r_break before persisting: vibrational
                # excursion, not a chemical event. This is the count a plain
                # cutoff analysis would have reported as a break.
                self.n_recrossings += 1
                state.pending_start = -1
        else:
            if r < self.r_form:
                if state.pending_start < 0:
                    state.pending_start = t
                    state.pending_r = r
                elif t - state.pending_start + 1 >= self.tau:
                    crossing = state.pending_start
                    state.bonded = True
                    state.pending_start = -1
                    if crossing < self.tau:
                        self.n_suppressed_burnin += 1
                        self.initial_bonds.add(key)
                        self.coord[o_row] += 1
                        self.coord[h_row] += 1
                    else:
                        forms_this_frame.append((crossing, o_row, h_row))
                        self._emit_form(crossing, o_row, h_row)
            elif state.pending_start >= 0:
                self.n_failed_forms += 1
                state.pending_start = -1

    def _emit_break(
        self,
        crossing: int,
        o_row: int,
        h_row: int,
        oo_counts: np.ndarray,
    ) -> None:
        frame = self._frame_at(crossing)
        if frame is None:                      # should not happen; buffer is deep enough
            return
        desc = self._break_descriptors(frame, o_row, h_row, oo_counts)
        desc["coord_o_donor_before"] = int(self.coord[o_row])
        desc["coord_o_donor_after"] = int(self.coord[o_row]) - 1
        self.coord[o_row] -= 1
        self.coord[h_row] -= 1
        self.pending_breaks.append(
            BreakEvent(
                frame=crossing,
                timestep=frame.timestep,
                o_row=o_row,
                h_row=h_row,
                r_at_break=desc["r_at_break"],
                descriptors=desc,
            )
        )

    def _emit_form(self, crossing: int, o_row: int, h_row: int) -> None:
        frame = self._frame_at(crossing)
        timestep = frame.timestep if frame is not None else -1
        acceptor_coord_before = int(self.coord[o_row])
        self.coord[o_row] += 1
        self.coord[h_row] += 1
        self.events.append(
            {
                "kind": "form",
                "frame": crossing,
                "timestep": timestep,
                "o_id": int(self.topology.ids[o_row]),
                "h_id": int(self.topology.ids[h_row]),
                "o_row": o_row,
                "h_row": h_row,
                "coord_o_before": acceptor_coord_before,
                "coord_o_after": acceptor_coord_before + 1,
            }
        )
        # A form may complete a pending break for the same H.
        for pending in self.pending_breaks:
            if pending.classification != "pending" or pending.h_row != h_row:
                continue
            if crossing < pending.frame or crossing - pending.frame > self.transfer_window:
                continue
            pending.classification = "recross" if o_row == pending.o_row else "transfer"
            pending.acceptor_row = o_row
            pending.acceptor_coord_before = acceptor_coord_before
            pending.acceptor_coord_after = acceptor_coord_before + 1
            pending.form_frame = crossing
            break

    def _advance_pending(self, t: int) -> None:
        """Finalise breaks whose classification window has fully elapsed."""
        deadline = self.transfer_window + self.tau
        still_pending: List[BreakEvent] = []
        for pending in self.pending_breaks:
            resolved = pending.classification != "pending"
            expired = t - pending.frame >= deadline
            if resolved or expired:
                if not resolved:
                    pending.classification = "dissociation"
                self.events.append(self._finalise(pending))
            else:
                still_pending.append(pending)
        self.pending_breaks = still_pending

    def _finalise(self, pending: BreakEvent, truncated: bool = False) -> dict:
        ids = self.topology.ids
        record = {
            "kind": "break",
            "frame": pending.frame,
            "timestep": pending.timestep,
            "o_id": int(ids[pending.o_row]),
            "h_id": int(ids[pending.h_row]),
            "o_row": pending.o_row,
            "h_row": pending.h_row,
            "classification": "truncated" if truncated else pending.classification,
        }
        record.update(pending.descriptors)
        if pending.acceptor_row >= 0:
            record["acceptor_o_id"] = int(ids[pending.acceptor_row])
            record["form_frame"] = pending.form_frame
            record["coord_o_acceptor_before"] = pending.acceptor_coord_before
            record["coord_o_acceptor_after"] = pending.acceptor_coord_after
            record.update(self._transfer_geometry(pending))
        return record

    def _transfer_geometry(self, pending: BreakEvent) -> dict:
        """O...O distance and O-H...O angle — the Grotthuss reaction coordinate."""
        out: dict = {}
        at_break = self._frame_at(pending.frame)
        if at_break is not None:
            d_oo = minimum_image(
                at_break.pos[pending.acceptor_row] - at_break.pos[pending.o_row],
                at_break.box_L,
            )
            d_oh = minimum_image(
                at_break.pos[pending.h_row] - at_break.pos[pending.o_row],
                at_break.box_L,
            )
            d_ho = minimum_image(
                at_break.pos[pending.acceptor_row] - at_break.pos[pending.h_row],
                at_break.box_L,
            )
            out["oo_distance_at_break"] = round(float(np.linalg.norm(d_oo)), 4)
            # Angle at the H: donor-O -- H ... acceptor-O. 180 deg is a linear,
            # textbook proton transfer.
            v1, v2 = -d_oh, d_ho
            cos = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
            out["ohO_angle_at_break"] = round(float(np.degrees(np.arccos(np.clip(cos, -1.0, 1.0)))), 3)
        at_form = self._frame_at(pending.form_frame)
        if at_form is not None:
            d_oo = minimum_image(
                at_form.pos[pending.acceptor_row] - at_form.pos[pending.o_row],
                at_form.box_L,
            )
            out["oo_distance_at_form"] = round(float(np.linalg.norm(d_oo)), 4)
        return out

    def _record_frame(
        self,
        frame: Frame,
        breaks_this_frame: List[Tuple[int, int, int]],
        forms_this_frame: List[Tuple[int, int, int]],
    ) -> None:
        self.per_frame.append(
            {
                "frame": frame.index,
                "timestep": frame.timestep,
                "n_breaks_confirmed": len(breaks_this_frame),
                "n_forms_confirmed": len(forms_this_frame),
            }
        )

    def _prune(self, t: int) -> None:
        """Drop unbonded pairs that have drifted well outside the query shell."""
        stale = t - 500
        for key in [k for k, s in self.pairs.items()
                    if not s.bonded and s.pending_start < 0 and s.last_seen < stale]:
            del self.pairs[key]

    def finish(self) -> Tuple[int, int]:
        """Flush unresolved state at end of trajectory.

        Returns (n_pending_at_eof, n_unconfirmed_crossings) — crossings that ran
        out of trajectory before they could be confirmed or classified. Reported
        rather than dropped, since silently discarding them would bias the rate.
        """
        n_pending = len(self.pending_breaks)
        for pending in self.pending_breaks:
            self.events.append(self._finalise(pending, truncated=True))
        self.pending_breaks = []
        n_unconfirmed = sum(1 for s in self.pairs.values() if s.pending_start >= 0)
        self.events.sort(key=lambda e: (e["frame"], e["kind"]))
        return n_pending, n_unconfirmed


# ==========================================================================
# Pass 2 — annotated trajectory for OVITO
# ==========================================================================


def write_annotated_trajectory(
    args: argparse.Namespace,
    topology: Topology,
    initial_bonds: Sequence[Tuple[int, int]],
    events: Sequence[dict],
    out_path: Path,
) -> List[dict]:
    """
    Re-emit the trajectory with per-atom annotation columns OVITO reads natively.

    The bond set is *replayed* from the initial bonds plus the event list rather
    than recomputed, so pass 2 needs no neighbour search and is guaranteed
    consistent with the events in events.jsonl.

    Columns appended: Coordination EventType EventAge PartnerID Species.
    These become user particle properties in OVITO, so a single frame-independent
    expression such as `EventType != 0` selects events across the whole
    trajectory — unlike silanol.py:507-519, which emits a per-frame OR-chain of
    (ParticleIdentifier == N) that has to be re-pasted at every frame.

    Returns per-frame species populations for the summary.
    """
    n = topology.n_atoms
    is_o = topology.elements == "O"

    coord = np.zeros(n, dtype=np.int32)
    h_partner = np.full(n, -1, dtype=np.int64)
    for o_row, h_row in initial_bonds:
        coord[o_row] += 1
        coord[h_row] += 1
        h_partner[h_row] = o_row

    by_frame: Dict[int, List[dict]] = {}
    for event in events:
        by_frame.setdefault(event["frame"], []).append(event)

    event_type = np.zeros(n, dtype=np.int32)
    event_age = np.full(n, -1, dtype=np.int64)
    partner_id = np.full(n, -1, dtype=np.int64)

    id_str = np.char.mod("%d", topology.ids)
    populations: List[dict] = []

    with open(out_path, "w") as out:
        for frame, _ in read_frames(
            args.traj,
            start=args.frames[0] if args.frames else 0,
            stop=args.frames[1] if args.frames else None,
            stride=args.stride,
            type_map=_parse_type_map(args.type_map),
        ):
            event_type[:] = EV_NONE
            live = event_age >= 0
            event_age[live] += 1

            for event in by_frame.get(frame.index, []):
                o_row, h_row = event["o_row"], event["h_row"]
                if event["kind"] == "break":
                    coord[o_row] -= 1
                    coord[h_row] -= 1
                    if h_partner[h_row] == o_row:
                        h_partner[h_row] = -1
                    event_type[o_row], event_type[h_row] = EV_O_BREAK, EV_H_BREAK
                else:
                    coord[o_row] += 1
                    coord[h_row] += 1
                    h_partner[h_row] = o_row
                    event_type[o_row], event_type[h_row] = EV_O_FORM, EV_H_FORM
                event_age[o_row] = event_age[h_row] = 0
                partner_id[o_row] = topology.ids[h_row]
                partner_id[h_row] = topology.ids[o_row]

            species = _species_from_coordination(coord, h_partner, is_o)
            populations.append(_population_counts(frame, species))
            _write_frame(out, frame, topology, id_str, coord, event_type,
                         event_age, partner_id, species)

    return populations


def _species_from_coordination(
    coord: np.ndarray, h_partner: np.ndarray, is_o: np.ndarray
) -> np.ndarray:
    """Classify every atom from its confirmed coordination number.

    Coordination comes from the hysteretic bond state, not from an instantaneous
    cutoff, so these labels do not flicker with the O-H stretch.
    """
    species = np.full(len(coord), SP_OTHER, dtype=np.int32)

    o_coord = np.where(is_o, coord, -1)
    species[o_coord == 2] = SP_WATER
    species[o_coord == 3] = SP_HYDRONIUM
    species[o_coord == 1] = SP_HYDROXIDE
    species[o_coord == 0] = SP_FREE_O

    is_h = ~is_o
    # An H inherits the label of the O it belongs to, so a whole H3O+ lights up
    # in OVITO under `Species == 1`, not just its oxygen.
    bonded_h = is_h & (coord == 1) & (h_partner >= 0)
    species[bonded_h] = species[h_partner[bonded_h]]
    species[is_h & (coord == 0)] = SP_FREE_H
    return species


def _population_counts(frame: Frame, species: np.ndarray) -> dict:
    counts = np.bincount(species, minlength=SP_OTHER + 1)
    return {
        "frame": frame.index,
        "timestep": frame.timestep,
        "n_hydronium": int(counts[SP_HYDRONIUM]) ,
        "n_hydroxide": int(counts[SP_HYDROXIDE]),
        "n_free_h": int(counts[SP_FREE_H]),
        "n_free_o": int(counts[SP_FREE_O]),
        "n_other": int(counts[SP_OTHER]),
    }


def _write_frame(
    out,
    frame: Frame,
    topology: Topology,
    id_str: np.ndarray,
    coord: np.ndarray,
    event_type: np.ndarray,
    event_age: np.ndarray,
    partner_id: np.ndarray,
    species: np.ndarray,
) -> None:
    """Write one annotated frame. Columns are formatted with vectorised numpy
    string ops; a per-atom Python loop would dominate the runtime at 26k atoms
    x thousands of frames."""
    has_v = frame.vel is not None
    names = ["id", "element", "x", "y", "z"]
    parts = [
        id_str,
        topology.elements,
        np.char.mod("%.5f", frame.pos[:, 0]),
        np.char.mod("%.5f", frame.pos[:, 1]),
        np.char.mod("%.5f", frame.pos[:, 2]),
    ]
    if has_v:
        names += ["vx", "vy", "vz"]
        parts += [np.char.mod("%.5f", frame.vel[:, d]) for d in range(3)]
    names += ["Coordination", "EventType", "EventAge", "PartnerID", "Species"]
    parts += [
        np.char.mod("%d", coord),
        np.char.mod("%d", event_type),
        np.char.mod("%d", event_age),
        np.char.mod("%d", partner_id),
        np.char.mod("%d", species),
    ]

    joined = parts[0]
    for part in parts[1:]:
        joined = np.char.add(np.char.add(joined, " "), part)

    lo, hi = frame.box_lo, frame.box_lo + frame.box_L
    out.write("ITEM: TIMESTEP\n%d\n" % frame.timestep)
    out.write("ITEM: NUMBER OF ATOMS\n%d\n" % topology.n_atoms)
    out.write("ITEM: BOX BOUNDS pp pp pp\n")
    for d in range(3):
        out.write("%.10e %.10e\n" % (lo[d], hi[d]))
    out.write("ITEM: ATOMS " + " ".join(names) + "\n")
    out.write("\n".join(joined.tolist()))
    out.write("\n")


# ==========================================================================
# Reporting
# ==========================================================================


def write_outputs(
    args: argparse.Namespace,
    topology: Topology,
    tracker: BondTracker,
    populations: List[dict],
    n_frames: int,
    eof: Tuple[int, int],
    output_dir: Path,
) -> None:
    events = tracker.events
    breaks = [e for e in events if e["kind"] == "break"]
    forms = [e for e in events if e["kind"] == "form"]

    with open(output_dir / "events.jsonl", "w") as f:
        for event in events:
            f.write(json.dumps(event) + "\n")

    _write_per_frame_csv(output_dir / "events_per_frame.csv", tracker.per_frame, populations)

    classes: Dict[str, int] = {}
    for event in breaks:
        classes[event["classification"]] = classes.get(event["classification"], 0) + 1
    (output_dir / "summary.json").write_text(
        json.dumps(
            {
                "n_frames": n_frames,
                "n_breaks": len(breaks),
                "n_forms": len(forms),
                "classifications": classes,
                "n_recrossings": tracker.n_recrossings,
                "n_failed_forms": tracker.n_failed_forms,
                "n_suppressed_burnin": tracker.n_suppressed_burnin,
                "n_pending_at_eof": eof[0],
                "n_unconfirmed_at_eof": eof[1],
            },
            indent=2,
        )
        + "\n"
    )

    summary = _summary_text(args, topology, tracker, breaks, forms, populations, n_frames, eof)
    (output_dir / "summary.txt").write_text(summary)
    print(summary)
    _plot(args, breaks, tracker.per_frame, populations, output_dir / "bond_events.png")


def _write_per_frame_csv(path: Path, per_frame: List[dict], populations: List[dict]) -> None:
    pops = {p["frame"]: p for p in populations}
    header = ("frame,timestep,n_breaks_confirmed,n_forms_confirmed,"
              "n_hydronium,n_hydroxide,n_free_h,n_free_o,n_other")
    lines = [header]
    for row in per_frame:
        pop = pops.get(row["frame"], {})
        lines.append(
            f"{row['frame']},{row['timestep']},{row['n_breaks_confirmed']},"
            f"{row['n_forms_confirmed']},{pop.get('n_hydronium', '')},"
            f"{pop.get('n_hydroxide', '')},{pop.get('n_free_h', '')},"
            f"{pop.get('n_free_o', '')},{pop.get('n_other', '')}"
        )
    path.write_text("\n".join(lines) + "\n")


def _summary_text(
    args: argparse.Namespace,
    topology: Topology,
    tracker: BondTracker,
    breaks: List[dict],
    forms: List[dict],
    populations: List[dict],
    n_frames: int,
    eof: Tuple[int, int],
) -> str:
    n_pending, n_unconfirmed = eof
    n_o = len(topology.o_rows)
    duration_ps = n_frames * args.dt / 1000.0

    classes: Dict[str, int] = {}
    for event in breaks:
        classes[event["classification"]] = classes.get(event["classification"], 0) + 1

    lines = [
        "bond_events.py summary",
        "=" * 60,
        f"trajectory        {args.traj}",
        f"frames            {n_frames}  (dt = {args.dt} fs, {duration_ps:.3f} ps total)",
        f"atoms             {topology.n_atoms}  ({n_o} O, {len(topology.h_rows)} H)",
        f"velocities        {'yes' if topology.has_velocities else 'no (kinetic descriptors omitted)'}",
        "",
        "criterion",
        "-" * 60,
        f"r_form            {args.r_form} A",
        f"r_break           {args.r_break} A",
        f"tau_persist       {args.tau_persist} frames = {args.tau_persist * args.dt:.1f} fs",
        f"transfer_window   {args.transfer_window} frames = {args.transfer_window * args.dt:.1f} fs",
        f"hbond_cutoff      {args.hbond_cutoff} A",
        "",
        "events",
        "-" * 60,
        f"breaks confirmed  {len(breaks)}",
        f"forms  confirmed  {len(forms)}",
    ]
    for name in ("transfer", "recross", "dissociation", "truncated"):
        count = classes.get(name, 0)
        share = 100.0 * count / len(breaks) if breaks else 0.0
        lines.append(f"  {name:<14} {count:>8}  ({share:5.1f}%)")

    if duration_ps > 0 and n_o:
        rate = len(breaks) / duration_ps / n_o
        lines.append(f"break rate        {rate:.4g} per ps per O")

    lines += [
        "",
        "flicker suppressed by the persistence criterion",
        "-" * 60,
        f"recrossings       {tracker.n_recrossings}",
        "  (excursions above r_break that did not persist tau_persist frames —",
        "   a plain single-cutoff analysis would have reported every one of these",
        "   as a bond break)",
        f"failed forms      {tracker.n_failed_forms}",
    ]
    if tracker.n_recrossings and breaks:
        lines.append(
            f"  ratio           {tracker.n_recrossings / len(breaks):.1f}x more flicker than real breaks"
        )

    lines += [
        "",
        "boundary effects (reported, not discarded)",
        "-" * 60,
        f"suppressed burn-in  {tracker.n_suppressed_burnin}  (first {args.tau_persist} frames, seeding the state machine)",
        f"pending at EOF      {n_pending}  (breaks classified 'truncated')",
        f"unconfirmed at EOF  {n_unconfirmed}  (crossings that ran out of trajectory)",
    ]

    if breaks:
        r_vals = np.array([e["r_at_break"] for e in breaks])
        lines += ["", "descriptors at the break frame", "-" * 60,
                  f"r_at_break        mean {r_vals.mean():.3f} A, max {r_vals.max():.3f} A"]
        v_vals = np.array([e["v_radial"] for e in breaks if "v_radial" in e])
        if len(v_vals):
            frac = 100.0 * float((v_vals > 0).mean())
            lines.append(f"v_radial          mean {v_vals.mean():+.3f} A/ps, {frac:.0f}% separating")
        ke_vals = np.array([e["ke_h_over_kT"] for e in breaks if "ke_h_over_kT" in e])
        if len(ke_vals):
            lines.append(f"ke_h / (3/2)kT    mean {ke_vals.mean():.2f} (1.0 = thermal average)")
        oo_vals = np.array([e["oo_distance_at_break"] for e in breaks if "oo_distance_at_break" in e])
        if len(oo_vals):
            lines.append(f"O...O at transfer mean {oo_vals.mean():.3f} A")
        ang_vals = np.array([e["ohO_angle_at_break"] for e in breaks if "ohO_angle_at_break" in e])
        if len(ang_vals):
            lines.append(f"O-H...O angle     mean {ang_vals.mean():.1f} deg (180 = linear hop)")

    if populations:
        h3o = np.array([p["n_hydronium"] for p in populations])
        oh = np.array([p["n_hydroxide"] for p in populations])
        lines += ["", "species populations", "-" * 60,
                  f"H3O+              mean {h3o.mean():.2f}, max {h3o.max()}",
                  f"OH-               mean {oh.mean():.2f}, max {oh.max()}",
                  f"charge imbalance  mean {float((h3o - oh).mean()):+.3f}"]
        lines.append("  (in equilibrium bulk water each H3O+ should be paired with an OH-;")
        lines.append("   a drifting imbalance means the cutoffs are wrong)")

    return "\n".join(lines) + "\n"


def _plot(
    args: argparse.Namespace,
    breaks: List[dict],
    per_frame: List[dict],
    populations: List[dict],
    path: Path,
) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))

    ax = axes[0, 0]
    if populations:
        t = np.array([p["frame"] for p in populations]) * args.dt / 1000.0
        ax.plot(t, [p["n_hydronium"] for p in populations], label="H3O+", lw=1.0)
        ax.plot(t, [p["n_hydroxide"] for p in populations], label="OH-", lw=1.0)
        ax.plot(t, [p["n_free_h"] for p in populations], label="free H", lw=1.0)
        ax.legend(fontsize=8)
    ax.set_xlabel("t (ps)")
    ax.set_ylabel("count")
    ax.set_title("Species populations")

    ax = axes[0, 1]
    if breaks:
        ax.hist([e["r_at_break"] for e in breaks], bins=40, color="C3")
        ax.axvline(args.r_break, color="k", ls="--", lw=1.0, label="r_break")
        ax.legend(fontsize=8)
    ax.set_xlabel("O-H distance at break (A)")
    ax.set_ylabel("events")
    ax.set_title("Where bonds break")

    ax = axes[1, 0]
    v_vals = [e["v_radial"] for e in breaks if "v_radial" in e]
    if v_vals:
        ax.hist(v_vals, bins=40, color="C0")
        ax.axvline(0.0, color="k", ls="--", lw=1.0)
        ax.set_title("Radial velocity at break (>0 = separating)")
    else:
        ax.set_title("Radial velocity — no velocities in dump")
    ax.set_xlabel("v_radial (A/ps)")
    ax.set_ylabel("events")

    ax = axes[1, 1]
    oo = [e["oo_distance_at_break"] for e in breaks if "oo_distance_at_break" in e]
    ang = [e["ohO_angle_at_break"] for e in breaks if "ohO_angle_at_break" in e]
    if len(oo) >= 5:
        ax.hist2d(oo, ang, bins=(30, 30), cmap="viridis")
    elif oo:
        ax.scatter(oo, ang, s=12, color="C2")
    ax.set_xlabel("O...O distance (A)")
    ax.set_ylabel("O-H...O angle (deg)")
    ax.set_title("Proton-transfer reaction coordinate")

    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"Plot saved to {path}")


# ==========================================================================
# CLI
# ==========================================================================


def _parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Bond breaking/forming events in bulk water, with an OVITO-ready output.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--traj", type=Path, required=True,
                   help="LAMMPS custom dump trajectory (id, x, y, z, element or type)")
    p.add_argument("--dt", type=float, default=2.5,
                   help="time between consecutive dumped frames, in fs "
                        "(OH-therm.input: 0.25 fs timestep x dynamics_dump_frequency 10)")
    p.add_argument("--r-form", type=float, default=1.15,
                   help="O-H distance below which a bond forms, in A")
    p.add_argument("--r-break", type=float, default=1.35,
                   help="O-H distance above which a bond breaks, in A (must exceed --r-form)")
    p.add_argument("--tau-persist", type=int, default=4,
                   help="consecutive frames a state change must persist to count "
                        "(default 4 x 2.5 fs = 10 fs, about one O-H stretch period)")
    p.add_argument("--transfer-window", type=int, default=20,
                   help="frames to watch a freed H before classifying the break")
    p.add_argument("--hbond-cutoff", type=float, default=3.5,
                   help="O...O distance defining the first solvation shell, in A")
    p.add_argument("--temperature", type=float, default=300.0,
                   help="run temperature in K, for the ke_h/kT descriptor")
    p.add_argument("--frames", nargs=2, type=int, default=None, metavar=("START", "STOP"),
                   help="half-open range of file frames to read")
    p.add_argument("--stride", type=int, default=1,
                   help="read every Nth frame (changes the effective dt — see warning)")
    p.add_argument("--type-map", default=None,
                   help="numeric type to element map for dumps without an 'element' "
                        "column, e.g. '1:O;2:H' (semicolon-separated)")
    p.add_argument("--output-dir", type=Path, default=None,
                   help="where to write results (default: output/ next to this script)")
    p.add_argument("--no-ovito-dump", action="store_true",
                   help="skip the annotated trajectory (pass 2), which is the slow part")
    return p.parse_args(argv)


def _validate(args: argparse.Namespace) -> None:
    if args.r_break <= args.r_form:
        raise ValueError(
            f"--r-break ({args.r_break}) must exceed --r-form ({args.r_form}); the gap "
            f"between them is the hysteresis that suppresses vibrational flicker. Equal "
            f"cutoffs reduce this to a plain threshold crossing."
        )
    if args.tau_persist < 1:
        raise ValueError("--tau-persist must be at least 1 frame.")
    if args.transfer_window < 1:
        raise ValueError("--transfer-window must be at least 1 frame.")
    if args.stride < 1:
        raise ValueError("--stride must be at least 1.")
    if not args.traj.exists():
        raise FileNotFoundError(f"Trajectory not found: {args.traj}")


def main(argv: Optional[List[str]] = None) -> int:
    args = _parse_args(argv)
    _validate(args)

    if args.stride > 1:
        print(
            f"WARNING: --stride {args.stride} multiplies the effective frame spacing to "
            f"{args.dt * args.stride:.1f} fs. --tau-persist is counted in *kept* frames, so "
            f"the persistence requirement is now {args.tau_persist * args.dt * args.stride:.1f} fs. "
            f"Striding a trajectory that already samples the O-H stretch at only ~3.6 points "
            f"per period will lose real events.",
            file=sys.stderr,
        )

    output_dir = args.output_dir or (Path(__file__).resolve().parent / "output")
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Reading trajectory: {args.traj}")
    t0 = time.time()
    tracker: Optional[BondTracker] = None
    topology: Optional[Topology] = None
    n_frames = 0

    for frame, topo in read_frames(
        args.traj,
        start=args.frames[0] if args.frames else 0,
        stop=args.frames[1] if args.frames else None,
        stride=args.stride,
        type_map=_parse_type_map(args.type_map),
    ):
        if tracker is None:
            topology = topo
            tracker = BondTracker(args, topo)
            print(f"  {topo.n_atoms} atoms ({len(topo.o_rows)} O, {len(topo.h_rows)} H)")
        tracker.step(frame)
        n_frames += 1

    if tracker is None or topology is None:
        raise ValueError(f"No frames read from {args.traj}.")

    eof = tracker.finish()
    t1 = time.time()
    print(f"  {n_frames} frames, pass 1 in {t1 - t0:.2f}s")

    populations: List[dict] = []
    if not args.no_ovito_dump:
        out_path = output_dir / "bond_events.lammpstrj"
        populations = write_annotated_trajectory(
            args, topology, tracker.initial_bonds, tracker.events, out_path
        )
        print(f"  pass 2 in {time.time() - t1:.2f}s -> {out_path}")

    write_outputs(args, topology, tracker, populations, n_frames, eof, output_dir)
    print(f"output_dir={output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
