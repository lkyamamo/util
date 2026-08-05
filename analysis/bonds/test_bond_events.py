#!/usr/bin/env python3
"""
test_bond_events.py — end-to-end validation of bond_events.py on synthetic
trajectories with known ground truth.

Run:  python test_bond_events.py

Same strategy as analysis/shock/etc/test_pipeline.py: build a LAMMPS dump whose
answer is known by construction, subprocess the real script, and compare.
Nothing here is checked against a reference output file, so the tests stay
meaningful if the implementation is rewritten.

The repo has no multi-frame water trajectory to use as a fixture —
analysis/structure/SiOH.0.lammpstrj is a single frame with no velocities — so
the trajectories are synthesised here.

CASES
-----
1. flicker       An O-H pair oscillating across r_break every few frames, never
                 staying out for tau_persist. Ground truth: ZERO breaks, one
                 recrossing per cycle. A plain single-cutoff implementation
                 would report a break on every cycle. This test is the whole
                 argument for the hysteresis criterion, made executable.
2. transfer      An H migrating from one O to another. Ground truth: one break
                 classified 'transfer', one form, donor 2->1 (OH-), acceptor
                 2->3 (H3O+), O...O = 2.6 A.
3. dissociation  An H that leaves and never comes back. Ground truth: one break
                 classified 'dissociation'.
4. pbc           Case 2 translated so the cluster straddles a periodic face.
                 Ground truth: results IDENTICAL to case 2. Catches the class of
                 bug in voxel_analysis.py:294, which builds a cKDTree with no
                 boxsize= and is silently non-periodic.
5. eof           Case 3 truncated shortly after the break confirms, so the
                 classification window runs off the end. Ground truth: the break
                 is reported as 'truncated' and counted in n_pending_at_eof
                 rather than silently dropped.
6. no_velocities Case 2 written without vx vy vz columns. Ground truth:
                 detection identical to case 2, since the criterion is purely
                 geometric; only v_radial/ke_h/ke_h_over_kT are absent, and they
                 must be absent rather than zero.
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

SCRIPT = Path(__file__).resolve().parent / "bond_events.py"

BOX_LO = np.array([0.0, 0.0, 0.0])
BOX_HI = np.array([20.0, 20.0, 20.0])
BOX_L = BOX_HI - BOX_LO

# Criterion used by every case. r_break - r_form is the hysteresis band.
R_FORM = 1.15
R_BREAK = 1.35
TAU = 4
TRANSFER_WINDOW = 20

_failures: List[str] = []


def chk(label: str, got, expected, tol: float = 1e-3) -> None:
    if isinstance(expected, float) or isinstance(got, float):
        ok = got is not None and abs(float(got) - float(expected)) <= tol
    else:
        ok = got == expected
    status = "PASS" if ok else "FAIL"
    print(f"  [{status}] {label}: got {got!r}, expected {expected!r}")
    if not ok:
        _failures.append(label)


# --------------------------------------------------------------------------
# Synthetic trajectory construction
# --------------------------------------------------------------------------


def write_dump(
    path: Path,
    elements: Sequence[str],
    positions: np.ndarray,     # (n_frames, n_atoms, 3), unwrapped
    velocities: Optional[np.ndarray] = None,
    with_velocity_columns: bool = True,
) -> None:
    """Write a LAMMPS custom dump matching OH-therm.input's column layout.

    with_velocity_columns=False omits vx vy vz entirely, as a dump written
    without them would — not zeros, but no columns at all.
    """
    n_frames, n_atoms, _ = positions.shape
    ids = np.arange(1, n_atoms + 1)
    with open(path, "w") as f:
        for t in range(n_frames):
            wrapped = (positions[t] - BOX_LO) % BOX_L + BOX_LO
            f.write(f"ITEM: TIMESTEP\n{t * 10}\n")
            f.write(f"ITEM: NUMBER OF ATOMS\n{n_atoms}\n")
            f.write("ITEM: BOX BOUNDS pp pp pp\n")
            for d in range(3):
                f.write(f"{BOX_LO[d]:.10e} {BOX_HI[d]:.10e}\n")
            if with_velocity_columns:
                f.write("ITEM: ATOMS id element x y z vx vy vz\n")
            else:
                f.write("ITEM: ATOMS id element x y z\n")
            for i in range(n_atoms):
                row = (
                    f"{ids[i]} {elements[i]} "
                    f"{wrapped[i, 0]:.6f} {wrapped[i, 1]:.6f} {wrapped[i, 2]:.6f}"
                )
                if with_velocity_columns:
                    v = velocities[t, i] if velocities is not None else (0.0, 0.0, 0.0)
                    row += f" {v[0]:.6f} {v[1]:.6f} {v[2]:.6f}"
                f.write(row + "\n")


def flicker_trajectory(n_cycles: int = 10) -> Tuple[List[str], np.ndarray]:
    """One O, two H. H2's O-H distance cycles across r_break but never stays
    beyond it for tau_persist frames.

    The cycle is [1.00, 1.10, 1.40, 1.45, 1.40, 1.10] A: three consecutive
    frames above r_break = 1.35, one short of the four required.
    """
    cycle = [1.00, 1.10, 1.40, 1.45, 1.40, 1.10]
    radii = cycle * n_cycles
    n_frames = len(radii)
    elements = ["O", "H", "H"]
    pos = np.zeros((n_frames, 3, 3))
    o = np.array([10.0, 10.0, 10.0])
    for t, r in enumerate(radii):
        pos[t, 0] = o
        pos[t, 1] = o + np.array([0.0, 1.00, 0.0])   # spectator, always bonded
        pos[t, 2] = o + np.array([r, 0.0, 0.0])      # the flickering one
    return elements, pos


def transfer_trajectory(shift: np.ndarray = np.zeros(3)) -> Tuple[List[str], np.ndarray]:
    """
    Two O 2.6 A apart; one H migrates from O1 to O2 between frames 10 and 19.

    Ground truth by construction, with r_form=1.15, r_break=1.35, tau=4:
      x(f) = 11.0 + 0.6*(f-10)/9 for f in [10, 19]
      r(H,O1) = x-10 first exceeds 1.35 at f=16  -> break timestamped frame 16,
                confirmed at frame 19
      r(H,O2) = 12.6-x first drops below 1.15 at f=17 -> form timestamped 17
      break and form are 1 frame apart, well inside the transfer window
      -> exactly one break, classified 'transfer'
    """
    n_frames = 50
    elements = ["O", "O", "H", "H", "H", "H"]
    o1 = np.array([10.0, 10.0, 10.0])
    o2 = np.array([12.6, 10.0, 10.0])
    pos = np.zeros((n_frames, 6, 3))
    for t in range(n_frames):
        if t < 10:
            x = 11.0
        elif t < 19:
            x = 11.0 + 0.6 * (t - 10) / 9.0
        else:
            x = 11.6
        pos[t, 0] = o1
        pos[t, 1] = o2
        pos[t, 2] = o1 + np.array([0.0, 1.0, 0.0])    # O1 spectator
        pos[t, 3] = o2 + np.array([0.0, 1.0, 0.0])    # O2 spectator
        pos[t, 4] = o2 + np.array([0.0, -1.0, 0.0])   # O2 spectator
        pos[t, 5] = np.array([x, 10.0, 10.0])         # the migrating H
    return elements, pos + shift


def dissociation_trajectory(n_frames: int = 60) -> Tuple[List[str], np.ndarray]:
    """One O, two H; H2 walks out to 5 A between frames 10 and 24 and stays.

    r first exceeds r_break = 1.35 at the frame where 1.0 + 4.0*(f-10)/14 > 1.35,
    i.e. (f-10) > 1.225 -> f = 12. Break timestamped frame 12.
    """
    elements = ["O", "H", "H"]
    o = np.array([10.0, 10.0, 10.0])
    pos = np.zeros((n_frames, 3, 3))
    for t in range(n_frames):
        if t < 10:
            r = 1.0
        elif t < 24:
            r = 1.0 + 4.0 * (t - 10) / 14.0
        else:
            r = 5.0
        pos[t, 0] = o
        pos[t, 1] = o + np.array([0.0, 1.0, 0.0])
        pos[t, 2] = o + np.array([r, 0.0, 0.0])
    return elements, pos


# --------------------------------------------------------------------------
# Runner
# --------------------------------------------------------------------------


def run_case(
    workdir: Path,
    name: str,
    elements: Sequence[str],
    positions: np.ndarray,
    extra_args: Sequence[str] = (),
    with_velocity_columns: bool = True,
) -> Tuple[dict, List[dict], Path]:
    """Write the trajectory, run bond_events.py on it, return (summary, events, outdir)."""
    traj = workdir / f"{name}.lammpstrj"
    outdir = workdir / f"{name}_output"
    write_dump(traj, elements, positions, with_velocity_columns=with_velocity_columns)

    cmd = [
        sys.executable, str(SCRIPT),
        "--traj", str(traj),
        "--output-dir", str(outdir),
        "--dt", "2.5",
        "--r-form", str(R_FORM),
        "--r-break", str(R_BREAK),
        "--tau-persist", str(TAU),
        "--transfer-window", str(TRANSFER_WINDOW),
        *extra_args,
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        print(proc.stdout)
        print(proc.stderr, file=sys.stderr)
        raise SystemExit(f"bond_events.py failed on case {name!r}")

    summary = json.loads((outdir / "summary.json").read_text())
    events = [json.loads(line) for line in (outdir / "events.jsonl").read_text().splitlines() if line]
    return summary, events, outdir


def read_annotated(path: Path) -> List[Dict[str, np.ndarray]]:
    """Parse the annotated dump back into per-frame column dicts."""
    frames: List[Dict[str, np.ndarray]] = []
    with open(path) as f:
        while True:
            line = f.readline()
            if not line:
                break
            f.readline()                      # timestep
            f.readline()
            n_atoms = int(f.readline())
            f.readline()
            for _ in range(3):
                f.readline()
            names = f.readline().split()[2:]
            rows = [f.readline().split() for _ in range(n_atoms)]
            table = np.array(rows)
            frame = {}
            for i, name in enumerate(names):
                col = table[:, i]
                frame[name] = col if name == "element" else col.astype(np.float64)
            frames.append(frame)
    return frames


# --------------------------------------------------------------------------
# Cases
# --------------------------------------------------------------------------


def case_flicker(workdir: Path) -> None:
    print("\n1. flicker — vibration across the cutoff must NOT count as breaking")
    n_cycles = 10
    elements, pos = flicker_trajectory(n_cycles)
    summary, events, _ = run_case(workdir, "flicker", elements, pos)

    chk("breaks", summary["n_breaks"], 0)
    chk("forms", summary["n_forms"], 0)
    chk("events written", len(events), 0)
    # Exactly one excursion above r_break per cycle. Recrossings are a
    # diagnostic counter rather than a confirmed event, so unlike a real break
    # they are counted during the burn-in frames too.
    chk("recrossings", summary["n_recrossings"], n_cycles)
    print(f"       a plain cutoff would have reported ~{n_cycles} breaks here")


def case_transfer(workdir: Path) -> dict:
    print("\n2. transfer — H migrating between two O")
    elements, pos = transfer_trajectory()
    summary, events, outdir = run_case(workdir, "transfer", elements, pos)

    chk("breaks", summary["n_breaks"], 1)
    chk("forms", summary["n_forms"], 1)
    chk("classified transfer", summary["classifications"].get("transfer", 0), 1)

    break_event = next(e for e in events if e["kind"] == "break")
    chk("break frame", break_event["frame"], 16)
    chk("donor O id", break_event["o_id"], 1)
    chk("migrating H id", break_event["h_id"], 6)
    chk("acceptor O id", break_event.get("acceptor_o_id"), 2)
    chk("donor coord before/after",
        (break_event.get("coord_o_donor_before"), break_event.get("coord_o_donor_after")), (2, 1))
    chk("acceptor coord before/after",
        (break_event.get("coord_o_acceptor_before"), break_event.get("coord_o_acceptor_after")), (2, 3))
    chk("O...O at break", break_event.get("oo_distance_at_break"), 2.6, tol=1e-3)
    chk("O-H...O angle", break_event.get("ohO_angle_at_break"), 180.0, tol=0.5)

    form_event = next(e for e in events if e["kind"] == "form")
    chk("form frame", form_event["frame"], 17)
    chk("form acceptor", form_event["o_id"], 2)

    # Coordination in the annotated trajectory: donor 2->1 (OH-), acceptor 2->3 (H3O+).
    frames = read_annotated(outdir / "bond_events.lammpstrj")
    chk("annotated frame count", len(frames), 50)
    chk("donor coord before", frames[15]["Coordination"][0], 2.0)
    chk("donor coord after", frames[20]["Coordination"][0], 1.0)
    chk("acceptor coord before", frames[15]["Coordination"][1], 2.0)
    chk("acceptor coord after", frames[20]["Coordination"][1], 3.0)
    chk("donor species after", frames[20]["Species"][0], 2.0)     # OH-
    chk("acceptor species after", frames[20]["Species"][1], 1.0)  # H3O+
    chk("H inherits H3O+ label", frames[20]["Species"][5], 1.0)

    # EventType/EventAge drive the OVITO expressions, so check them directly.
    chk("EventType at break frame (O)", frames[16]["EventType"][0], 1.0)
    chk("EventType at break frame (H)", frames[16]["EventType"][5], 2.0)
    chk("EventType away from events", float(frames[5]["EventType"].max()), 0.0)
    chk("EventAge before any event", frames[5]["EventAge"][0], -1.0)
    chk("EventAge at break frame", frames[16]["EventAge"][0], 0.0)
    chk("EventAge 3 frames later", frames[19]["EventAge"][0], 3.0)
    chk("PartnerID at break", frames[16]["PartnerID"][0], 6.0)
    return summary


def case_dissociation(workdir: Path) -> None:
    print("\n3. dissociation — H leaves and never returns")
    elements, pos = dissociation_trajectory()
    summary, events, _ = run_case(workdir, "dissociation", elements, pos)

    chk("breaks", summary["n_breaks"], 1)
    chk("forms", summary["n_forms"], 0)
    chk("classified dissociation", summary["classifications"].get("dissociation", 0), 1)
    break_event = next(e for e in events if e["kind"] == "break")
    chk("break frame", break_event["frame"], 12)
    chk("no acceptor recorded", "acceptor_o_id" in break_event, False)


def case_pbc(workdir: Path, reference: dict) -> None:
    print("\n4. pbc — the same event across a periodic face must be identical")
    # Translate so the O...O pair straddles x = 0 after wrapping.
    elements, pos = transfer_trajectory(shift=np.array([-11.3, 0.0, 0.0]))
    summary, events, _ = run_case(workdir, "pbc", elements, pos)

    chk("breaks", summary["n_breaks"], reference["n_breaks"])
    chk("forms", summary["n_forms"], reference["n_forms"])
    chk("classifications", summary["classifications"], reference["classifications"])
    break_event = next(e for e in events if e["kind"] == "break")
    chk("break frame", break_event["frame"], 16)
    chk("O...O across boundary", break_event.get("oo_distance_at_break"), 2.6, tol=1e-3)
    chk("O-H...O angle across boundary", break_event.get("ohO_angle_at_break"), 180.0, tol=0.5)


def case_eof(workdir: Path) -> None:
    print("\n5. eof — a break whose classification window runs off the end")
    # The break is timestamped frame 12 and confirms at frame 15; ending at
    # frame 20 leaves the 20-frame classification window incomplete.
    elements, pos = dissociation_trajectory(n_frames=20)
    summary, events, _ = run_case(workdir, "eof", elements, pos)

    chk("breaks", summary["n_breaks"], 1)
    chk("pending at EOF", summary["n_pending_at_eof"], 1)
    chk("classified truncated", summary["classifications"].get("truncated", 0), 1)
    break_event = next(e for e in events if e["kind"] == "break")
    chk("break frame", break_event["frame"], 12)


def case_no_velocities(workdir: Path, reference: dict) -> None:
    print("\n6. no velocities — a dump without vx vy vz must still work")
    # dump.lammpstrj carries velocities, but plenty of dumps do not (e.g.
    # analysis/structure/SiOH.0.lammpstrj, or dielectric-therm.input's
    # 'id type x y z'). Detection is purely geometric, so everything except the
    # three kinetic descriptors must be unchanged.
    elements, pos = transfer_trajectory()
    summary, events, outdir = run_case(
        workdir, "no_velocities", elements, pos, with_velocity_columns=False
    )

    chk("breaks", summary["n_breaks"], reference["n_breaks"])
    chk("forms", summary["n_forms"], reference["n_forms"])
    chk("classifications", summary["classifications"], reference["classifications"])

    break_event = next(e for e in events if e["kind"] == "break")
    chk("break frame unchanged", break_event["frame"], 16)
    chk("acceptor unchanged", break_event.get("acceptor_o_id"), 2)
    chk("O...O unchanged", break_event.get("oo_distance_at_break"), 2.6, tol=1e-3)
    chk("angle unchanged", break_event.get("ohO_angle_at_break"), 180.0, tol=0.5)
    chk("donor coord unchanged",
        (break_event.get("coord_o_donor_before"), break_event.get("coord_o_donor_after")), (2, 1))

    # The velocity-derived descriptors must be absent, not zero — a zero
    # v_radial would silently pollute the summary statistics.
    for key in ("v_radial", "ke_h", "ke_h_over_kT"):
        chk(f"{key} omitted", key in break_event, False)

    # The annotated trajectory must drop the velocity columns rather than
    # emit placeholders, so OVITO does not show a bogus zero Velocity property.
    with open(outdir / "bond_events.lammpstrj") as f:
        header = [f.readline() for _ in range(9)][8].split()[2:]
    chk("annotated columns", header,
        ["id", "element", "x", "y", "z",
         "Coordination", "EventType", "EventAge", "PartnerID", "Species"])


def main() -> int:
    if not SCRIPT.exists():
        raise SystemExit(f"bond_events.py not found next to this test: {SCRIPT}")

    print("bond_events.py — synthetic ground-truth tests")
    print("=" * 62)
    with tempfile.TemporaryDirectory(prefix="bond_events_test_") as tmp:
        workdir = Path(tmp)
        case_flicker(workdir)
        transfer_summary = case_transfer(workdir)
        case_dissociation(workdir)
        case_pbc(workdir, transfer_summary)
        case_eof(workdir)
        case_no_velocities(workdir, transfer_summary)

    print("\n" + "=" * 62)
    if _failures:
        print(f"{len(_failures)} FAILED: " + ", ".join(_failures))
        return 1
    print("all checks PASSED")
    return 0


if __name__ == "__main__":
    sys.exit(main())
