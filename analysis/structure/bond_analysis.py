"""
bond_analysis.py

Identifies and shifts surface-dangling (bondless) atoms at the top and bottom
of a non-periodic slab, exports the free-atom subsets as LAMMPS dump files,
and writes a shifted LAMMPS data file with the surface atoms moved inward.

Outputs (written to the directory of the input file):
  <stem>_shifted.data       – modified structure
  free_atoms_bottom.dump    – bondless atoms at the bottom surface (pre-shift)
  free_atoms_bottom_ids.txt – OVITO expression-selection filter string
  free_atoms_top.dump       – bondless atoms at the top surface (post-bottom-shift)
  free_atoms_top_ids.txt    – OVITO expression-selection filter string
"""

from pathlib import Path

import numpy as np
from ovito.io import import_file, export_file
from ovito.modifiers import (
    CreateBondsModifier,
    DeleteSelectedModifier,
    InvertSelectionModifier,
)

# ── Configuration ─────────────────────────────────────────────────────────────
FILE_PATH = "/Users/loganyamamoto/Desktop/Research/grants/geo_sciences/finalized-bubble-collapse/structures/water/0013_19x1x1.data"
SURFACE_DIST = 1.2   # Å from top/bottom surfaces
BOND_CUTOFF  = 1.1   # Å pairwise O–H cutoff
REFIT_BOX    = True  # Shift atoms so min == 0 and resize box to atom extent

# Non-periodic direction: 'x', 'y', or 'z' (slab normal / confined axis)
DIRECTION = "x"

# If True, atoms that are bondless even under full 3-D PBC are detected on the
# original structure and excluded from the surface-shifting candidates.
EXCLUDE_PREEXISTING_ORPHANS = True
# ─────────────────────────────────────────────────────────────────────────────


def main():
    file_path   = Path(FILE_PATH)
    shifted_path = file_path.parent / f"{file_path.stem}_shifted{file_path.suffix}"
    axis = {"x": 0, "y": 1, "z": 2}[DIRECTION]

    # ── Helper modifiers ──────────────────────────────────────────────────────

    def set_pbc(frame, data):
        """Non-periodic along the slab normal; periodic in the other two."""
        pbc = [True, True, True]
        pbc[axis] = False
        data.cell_.pbc = tuple(pbc)

    def compute_coordination_from_bonds(frame, data):
        """Count bonds per particle and store as 'Coordination' property."""
        bonds = data.particles.bonds
        if bonds is None:
            raise RuntimeError("No bonds — ensure CreateBondsModifier runs first.")
        topology = bonds.topology.array
        n = data.particles.count
        coordination = np.zeros(n, dtype=int)
        np.add.at(coordination, topology[:, 0], 1)
        np.add.at(coordination, topology[:, 1], 1)
        data.particles_.create_property("Coordination", data=coordination)

    def make_bond_modifier(cutoff=BOND_CUTOFF):
        m = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)
        m.set_pairwise_cutoff(1, 2, cutoff)
        return m

    # ── Export a subset of particles ─────────────────────────────────────────

    def export_subset(pipeline, mask, filename):
        """Temporarily add select/invert/delete modifiers, export, then remove."""
        mask_copy = mask.copy()

        def select_subset(frame, data):
            data.particles_.create_property("Selection", data=mask_copy.astype(int))

        pipeline.modifiers.append(select_subset)
        pipeline.modifiers.append(InvertSelectionModifier())
        pipeline.modifiers.append(DeleteSelectedModifier())
        try:
            export_file(
                pipeline, filename, "lammps/dump",
                columns=["Particle Identifier", "Particle Type",
                         "Position.X", "Position.Y", "Position.Z"],
            )
        finally:
            for _ in range(3):
                pipeline.modifiers.pop()

    # ── Pre-existing orphan detection (full 3-D PBC) ──────────────────────────
    if EXCLUDE_PREEXISTING_ORPHANS:
        pbc_pip = import_file(str(file_path))
        pbc_pip.modifiers.append(lambda f, d: setattr(d.cell_, "pbc", (True, True, True)))
        pbc_pip.modifiers.append(make_bond_modifier())
        pbc_pip.modifiers.append(compute_coordination_from_bonds)
        d_pbc = pbc_pip.compute()
        coord_pbc = d_pbc.particles["Coordination"].array
        preexisting_orphan_ids = set(
            d_pbc.particles["Particle Identifier"].array[coord_pbc == 0].tolist()
        )
        print(f"Pre-existing orphans (full PBC, excluded from shifting): "
              f"{len(preexisting_orphan_ids)}")
    else:
        preexisting_orphan_ids = set()
        print("Pre-existing orphan exclusion disabled.")

    # ── Build the main pipeline ───────────────────────────────────────────────
    pipeline = import_file(str(file_path))
    pipeline.modifiers.append(set_pbc)
    pipeline.modifiers.append(make_bond_modifier())
    pipeline.modifiers.append(compute_coordination_from_bonds)

    # ── Pass 1: find bondless atoms near the bottom surface ──────────────────
    data = pipeline.compute()

    positions    = data.particles.positions
    coord        = data.particles["Coordination"].array
    particle_ids = data.particles["Particle Identifier"].array
    axis_coords  = positions[:, axis]

    axis_min      = data.cell[axis, 3]
    axis_max      = data.cell[axis, 3] + data.cell[axis, axis]
    axis_mid      = 0.5 * (axis_min + axis_max)
    axis_box_size = axis_max - axis_min
    print(f"{DIRECTION}_min={axis_min:.4f}  {DIRECTION}_max={axis_max:.4f}")

    preexisting_mask = np.isin(particle_ids, list(preexisting_orphan_ids))
    free_mask = (coord == 0) & ~preexisting_mask

    bottom_mask = (
        free_mask
        & (axis_coords <= axis_mid)
        & (axis_coords <= axis_min + SURFACE_DIST)
    )

    out_dir = file_path.parent
    export_subset(pipeline, bottom_mask, str(out_dir / "free_atoms_bottom.dump"))
    bottom_ids = particle_ids[bottom_mask]
    (out_dir / "free_atoms_bottom_ids.txt").write_text(
        " || ".join(f"ParticleIdentifier == {pid}" for pid in bottom_ids)
    )
    print(f"Pass 1: {bottom_mask.sum()} bottom free atoms saved.")

    # ── Shift bottom free atoms out of the box (+axis_box_size) ──────────────
    bottom_mask_copy = bottom_mask.copy()

    def shift_bottom_along_axis(frame, data):
        pos = np.copy(data.particles.positions)
        pos[bottom_mask_copy, axis] += axis_box_size
        data.particles_.positions_ = pos

    pipeline.modifiers.append(shift_bottom_along_axis)

    # Re-run bond + coordination after the shift
    pipeline.modifiers.append(make_bond_modifier())
    pipeline.modifiers.append(compute_coordination_from_bonds)

    # ── Pass 2: find bondless atoms near the top surface ─────────────────────
    data2 = pipeline.compute()

    coord2       = data2.particles["Coordination"].array
    axis_coords2 = data2.particles.positions[:, axis]
    free_mask2   = (coord2 == 0)

    top_mask = (
        free_mask2
        & (axis_coords2 > axis_mid)
        & (axis_coords2 >= axis_max - SURFACE_DIST)
        & ~bottom_mask_copy
    )

    export_subset(pipeline, top_mask, str(out_dir / "free_atoms_top.dump"))
    top_ids = data2.particles["Particle Identifier"].array[top_mask]
    (out_dir / "free_atoms_top_ids.txt").write_text(
        " || ".join(f"ParticleIdentifier == {pid}" for pid in top_ids)
    )
    print(f"Pass 2: {top_mask.sum()} top free atoms saved.")

    # ── Shift top atoms inward and export final structure ─────────────────────
    top_mask_copy = top_mask.copy()

    def shift_top_along_axis(frame, data):
        pos = np.copy(data.particles.positions)
        pos[top_mask_copy, axis] -= axis_box_size
        data.particles_.positions_ = pos

    def shift_and_resize_box(frame, data):
        pos = np.copy(data.particles.positions)
        axis_min_val = pos[:, axis].min()
        pos[:, axis] -= axis_min_val
        axis_range = pos[:, axis].max()
        data.particles_.positions_ = pos
        data.cell_[axis, 3] = 0.0
        data.cell_[axis, axis] = axis_range

    pipeline.modifiers.append(shift_top_along_axis)
    n_extra = 1
    if REFIT_BOX:
        pipeline.modifiers.append(shift_and_resize_box)
        n_extra += 1

    export_file(pipeline, str(shifted_path), "lammps/data")
    for _ in range(n_extra):
        pipeline.modifiers.pop()
    print(f"Saved modified structure → {shifted_path}")

    # ── Assertions ────────────────────────────────────────────────────────────
    assert bottom_mask.sum() == len(bottom_ids), "Bottom mask / saved IDs mismatch"
    assert top_mask.sum()    == len(top_ids),    "Top mask / saved IDs mismatch"

    # ── Validation: no new orphans introduced by shifting ────────────────────
    val_pip = import_file(str(shifted_path))
    val_pip.modifiers.append(set_pbc)
    val_pip.modifiers.append(make_bond_modifier())
    val_pip.modifiers.append(compute_coordination_from_bonds)
    data_final = val_pip.compute()

    coord_final      = data_final.particles["Coordination"].array
    final_orphan_ids = set(
        data_final.particles["Particle Identifier"].array[coord_final == 0].tolist()
    )
    new_orphan_ids = final_orphan_ids - preexisting_orphan_ids
    n_pre = len(final_orphan_ids & preexisting_orphan_ids)

    print(f"\nValidation — orphan atoms after shifts (coord==0): {len(final_orphan_ids)}")
    print(f"  Pre-existing (full PBC, expected) : {n_pre}")
    print(f"  Newly introduced by shifting      : {len(new_orphan_ids)}  [must be 0]")

    if new_orphan_ids:
        print("Newly-introduced orphan IDs (first 20):", sorted(new_orphan_ids)[:20])

    assert len(new_orphan_ids) == 0, (
        f"Shifting introduced {len(new_orphan_ids)} new orphan atoms."
    )

    # ── Geometry check on shifted file ───────────────────────────────────────
    axis_all   = data_final.particles.positions[:, axis]
    box_origin = data_final.cell[axis, 3]
    box_size   = data_final.cell[axis, axis]
    print(f"\nGeometry check on {shifted_path.name}:")
    print(f"  {DIRECTION} range (particles): {axis_all.min():.6f} → {axis_all.max():.6f}")
    print(f"  Box origin={box_origin:.6f}  size={box_size:.6f}")
    print(f"  {DIRECTION} axis length: {axis_all.max() - axis_all.min():.6f} Å")

    # ── Double-shift guard ────────────────────────────────────────────────────
    double_shifted = bottom_mask & top_mask
    n_double = int(double_shifted.sum())
    print(f"\nDouble-shift check: {n_double} atoms in both masks (expected 0)")
    assert n_double == 0, f"{n_double} atoms were shifted twice."

    print("\nAll checks passed.")


if __name__ == "__main__":
    main()
