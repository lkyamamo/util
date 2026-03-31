#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path
from typing import List, Optional, Tuple

import numpy as np
import ovito
from ovito.io import export_file, import_file
from ovito.modifiers import (
    CreateBondsModifier,
    RadialDistributionFunctionModifier,
    TimeAveragingModifier,
)


def _ensure_bonds(
    pipeline,
    *,
    si_o_cutoff: float = 1.9,
    o_h_cutoff: float = 1.2,
) -> None:
    """Ensure the pipeline has a CreateBondsModifier with desired cutoffs."""
    cb = next((m for m in pipeline.modifiers if isinstance(m, CreateBondsModifier)), None)
    if cb is None:
        cb = CreateBondsModifier(mode=CreateBondsModifier.Mode.Pairwise)
        pipeline.modifiers.append(cb)

    # OVITO API: set_pairwise_cutoff(typeA, typeB, cutoff)
    cb.set_pairwise_cutoff("Si", "O", si_o_cutoff)
    cb.set_pairwise_cutoff("O", "H", o_h_cutoff)


def _type_id(data, name: str) -> int:
    """Return OVITO particle-type id for an element name like 'Si', 'O', 'H'."""
    for t in data.particles.particle_types.types:
        if t.name == name:
            return t.id
    available = [(t.id, t.name) for t in data.particles.particle_types.types]
    raise KeyError(f"Particle type '{name}' not found. Available types: {available}")


def average_count_silanols(
    pipeline,
    *,
    frames: Optional[Tuple[int, int]] = None,
    si_o_cutoff: float = 1.9,
    o_h_cutoff: float = 1.2,
    allow_other_neighbors: bool = True,
    return_identifiers_if_available: bool = True,
) -> Tuple[float, List[List[int]]]:
    """Return (avg_count, ids_per_frame) for O bonded to exactly 1 Si (and any # of H).

    ids are OVITO particle identifiers if present, otherwise particle indices (0-based).
    """
    _ensure_bonds(pipeline, si_o_cutoff=si_o_cutoff, o_h_cutoff=o_h_cutoff)

    total_silanols = 0
    ids: List[List[int]] = []

    start_frame = frames[0] if frames else 0
    end_frame = frames[1] if frames else len(pipeline.frames)

    for frame in range(start_frame, end_frame):
        data = pipeline.compute(frame)
        ptype = data.particles["Particle Type"].array
        bonds = data.particles.bonds

        t_Si = _type_id(data, "Si")
        t_O = _type_id(data, "O")
        t_H = _type_id(data, "H")

        neighbors: List[List[int]] = [[] for _ in range(data.particles.count)]
        for i, j in bonds.topology:
            neighbors[i].append(j)
            neighbors[j].append(i)

        matched_indices: List[int] = []
        for i in range(data.particles.count):
            if ptype[i] != t_O:
                continue

            nb = neighbors[i]
            n_si = sum(ptype[j] == t_Si for j in nb)
            if n_si != 1:
                continue

            if not allow_other_neighbors:
                if any((ptype[j] != t_Si and ptype[j] != t_H) for j in nb):
                    continue

            matched_indices.append(i)

        if return_identifiers_if_available and "Particle Identifier" in data.particles:
            pid = data.particles["Particle Identifier"].array
            matched_ids = [int(pid[i]) for i in matched_indices]
        else:
            matched_ids = matched_indices

        total_silanols += len(matched_ids)
        ids.append(matched_ids)

    return total_silanols / (end_frame - start_frame), ids


def _parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute silanol count and RDF (OVITO).")
    p.add_argument(
        "--base-dir",
        type=Path,
        default=None,
        help="Directory containing '2.0_64.data' (defaults to this script's directory).",
    )
    return p.parse_args(argv)


def main(argv: Optional[List[str]] = None) -> None:
    args = _parse_args(argv)
    base_dir = args.base_dir.resolve() if args.base_dir is not None else Path(__file__).resolve().parent
    output_dir = base_dir / "output"
    output_dir.mkdir(parents=True, exist_ok=True)

    vasp_xdatcar_path = base_dir / "2.0_64.data"

    pipeline = import_file(str(vasp_xdatcar_path))
    pipeline.modifiers.append(
        RadialDistributionFunctionModifier(
            number_of_bins=300,
            cutoff=10,
        )
    )

    count, ids = average_count_silanols(pipeline, frames=(150, 198))
    print((count, ids[:20]))

    pipeline.modifiers.append(TimeAveragingModifier(operate_on="table:coordination-rdf"))

    total_rdf = pipeline.compute().tables["coordination-rdf[average]"].xy()
    np.savetxt(str(output_dir / "rdf.txt"), total_rdf)
    export_file(pipeline, str(output_dir / "rdf.txt"), "txt/table", key="coordination-rdf[average]")


if __name__ == "__main__":
    main()

