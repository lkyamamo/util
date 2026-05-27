"""Build OVITO pipelines with per-type sphere radii and pairwise bond cutoffs."""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path


def merge_symmetric_pairs(
    pair_cutoffs: Mapping[tuple[str, str], float],
    *,
    symmetric: bool = True,
) -> dict[tuple[str, str], float]:
    """Mirror (a,b) → (b,a) so OVITO Pairwise tables get both directions."""
    merged: dict[tuple[str, str], float] = {}
    for (a, b), dist in pair_cutoffs.items():
        ta, tb = str(a).strip(), str(b).strip()
        merged[(ta, tb)] = float(dist)
        if symmetric and ta != tb:
            merged[(tb, ta)] = float(dist)
    return merged


def apply_lammps_type_id_to_symbol(
    pipeline,
    id_to_symbol: Mapping[int, str] | None,
) -> None:
    """Rename ParticleType.name from LAMMPS numeric IDs to symbols (must run first)."""

    if not id_to_symbol:
        return

    mapping = {int(k): str(v).strip() for k, v in id_to_symbol.items()}

    def modify(frame: int, data):  # noqa: ARG001 — OVITO signature
        pts = data.particles_
        if pts is None:
            return
        ptypes = pts.particle_types_
        if ptypes is None:
            return
        for element_type in ptypes.types_:
            tid_i = int(element_type.id)
            if tid_i in mapping:
                element_type.name = mapping[tid_i]

    pipeline.modifiers.append(modify)


def _make_atom_radius_assigner(atom_display_radii: dict[str, float]):
    normalized = {k.strip(): float(v) for k, v in atom_display_radii.items()}

    def modify(frame: int, data):  # noqa: ARG001 — OVITO signature
        pts = data.particles_
        if pts is None:
            return
        ptypes = pts.particle_types_
        if ptypes is None:
            return
        # Use types_ (mutable) — .types is read-only and won't persist radius changes.
        for element_type in ptypes.types_:
            key = str(element_type.name).strip()
            if key in normalized:
                element_type.radius = normalized[key]

    return modify


def apply_atom_display_radii(
    pipeline,
    atom_display_radii: Mapping[str, float] | None,
) -> None:
    """Append a modifier that assigns ParticleType.radius for rendered sphere sizes."""
    if not atom_display_radii:
        return
    pipeline.modifiers.append(_make_atom_radius_assigner(dict(atom_display_radii)))


def apply_pairwise_cutoff_bonds(
    pipeline,
    pair_cutoffs: Mapping[tuple[str, str], float] | None,
    *,
    bond_visual_radius: float | None = None,
    lower_cutoff: float = 0.0,
    symmetric_pairs: bool = True,
) -> None:
    """
    Draw bonds when separation is below the per-type-pair cutoff (Å, usual case).

    OVITO CreateBonds Pairwise mode requires at least one positive cutoff that matches
    present particle types. Type names are resolved from frame 0; unknown names are skipped.
    If nothing remains, the modifier is omitted (no bonds).
    """
    if not pair_cutoffs:
        return

    from ovito.modifiers import CreateBondsModifier

    merged = merge_symmetric_pairs(pair_cutoffs, symmetric=symmetric_pairs)

    data = pipeline.compute(0)
    if data.particles is None:
        return
    type_prop = data.particles.particle_types
    available = {str(t.name).strip() for t in type_prop.types}

    mod = CreateBondsModifier()
    mod.mode = CreateBondsModifier.Mode.Pairwise
    mod.lower_cutoff = float(lower_cutoff)

    n_set = 0
    for (a, b), dist in merged.items():
        ta, tb = str(a).strip(), str(b).strip()
        if dist <= 0.0:
            continue
        if ta not in available or tb not in available:
            continue
        mod.set_pairwise_cutoff(ta, tb, float(dist))
        n_set += 1

    if n_set == 0:
        print(
            "ovito_visual_pipeline: skipping CreateBonds — no positive pair cutoffs apply "
            f"to loaded types {sorted(available)} (check BOND_PAIR_CUTOFFS names vs dump)."
        )
        return

    if bond_visual_radius is not None and bond_visual_radius > 0:
        mod.vis.radius = float(bond_visual_radius)

    pipeline.modifiers.append(mod)


def load_styled_pipeline(
    dump_path: str | Path,
    *,
    atom_display_radii: Mapping[str, float] | None,
    pairwise_bond_cutoffs: Mapping[tuple[str, str], float] | None,
    bond_visual_radius: float | None = None,
    symmetric_bond_pairs: bool = True,
    lammps_type_id_to_symbol: Mapping[int, str] | None = None,
):
    """import_file + optional LAMMPS type renaming + display radii + pairwise CreateBonds."""
    from ovito.io import import_file

    pipeline = import_file(str(dump_path))
    apply_lammps_type_id_to_symbol(pipeline, lammps_type_id_to_symbol)
    apply_atom_display_radii(pipeline, atom_display_radii)
    apply_pairwise_cutoff_bonds(
        pipeline,
        pairwise_bond_cutoffs,
        bond_visual_radius=bond_visual_radius,
        symmetric_pairs=symmetric_bond_pairs,
    )
    return pipeline
