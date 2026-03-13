"""
compare_structures_pymatgen.py

Compares an original structure and a CG-relaxed structure, both in LAMMPS
data file format, using pymatgen. Computes:
  1. ΔV/V  — volumetric change
  2. Lattice parameter comparison (a, b, c, α, β, γ)
  3. MSD in real coordinates (Å), MIC via original lattice
  4. COM shift — reported as diagnostic only, not subtracted
  5. StructureMatcher — symmetry-aware structural similarity
  6. Per-species MSD breakdown
  7. Top 5 most displaced atoms

Usage:
  python compare_structures_pymatgen.py                         # default output file
  python compare_structures_pymatgen.py orig.lammps cg.lammps
  python compare_structures_pymatgen.py --output my_results.txt
  python compare_structures_pymatgen.py --no-export             # print only
"""

import argparse
import numpy as np
from pymatgen.io.lammps.data import LammpsData
from pymatgen.analysis.structure_matcher import StructureMatcher


# =============================================================================
# ARGUMENT PARSING
# =============================================================================

def parse_args():
    parser = argparse.ArgumentParser(
        description="Compare original and CG-relaxed LAMMPS structures."
    )
    parser.add_argument(
        "orig_file",
        nargs="?",
        default="original.lammps",
        help="Path to original LAMMPS data file (default: original.lammps)"
    )
    parser.add_argument(
        "cg_file",
        nargs="?",
        default="cg_relaxed.lammps",
        help="Path to CG-relaxed LAMMPS data file (default: cg_relaxed.lammps)"
    )
    parser.add_argument(
        "--output", "-o",
        default="structure_comparison.txt",
        help="Output file path (default: structure_comparison.txt)"
    )
    parser.add_argument(
        "--no-export",
        action="store_true",
        default=False,
        help="Disable file export, print to terminal only"
    )
    return parser.parse_args()


# =============================================================================
# OUTPUT HANDLER
# =============================================================================

class OutputHandler:
    def __init__(self, export: bool, filepath: str):
        self.export   = export
        self.filepath = filepath
        self.file     = None

    def __enter__(self):
        if self.export:
            self.file = open(self.filepath, "w")
            print(f"Exporting results to: {self.filepath}")
        return self

    def __exit__(self, *args):
        if self.file:
            self.file.close()

    def print(self, *args, **kwargs):
        print(*args, **kwargs)
        if self.file:
            print(*args, **kwargs, file=self.file)


# =============================================================================
# DISPLACEMENT HELPER
# =============================================================================
# Displacements are computed in Cartesian space (Å).
#
# The MIC is applied by converting the raw Cartesian delta into fractional
# coordinates of the ORIGINAL lattice, applying np.round() to map each
# component into (-0.5, 0.5], then converting back to Cartesian Å.
#
# Using the original lattice as the reference frame is correct because:
#   - It defines the periodic boundaries the atoms started in
#   - The CG lattice has a different basis after box relaxation, so converting
#     through it would give displacements in an inconsistent frame
#
# The COM shift is computed and returned for diagnostic reporting but is NOT
# subtracted from the displacements. For a two-frame CG comparison, any
# apparent COM drift is an artifact of box remapping or floating point
# accumulation during minimization — not a physical rigid-body translation.
# Subtracting it would be correcting for something that was never physically
# real. It is reported so the user can verify it is small and flag anything
# unexpected.

def compute_displacements(cart_orig, cart_cg, lattice_orig, masses):
    """
    Compute per-atom displacements in Cartesian Å with MIC applied via the
    original lattice. COM shift is computed as a diagnostic but not removed.

    Parameters
    ----------
    cart_orig    : (N, 3) Cartesian coordinates of original structure (Å)
    cart_cg      : (N, 3) Cartesian coordinates of CG structure (Å)
    lattice_orig : (3, 3) lattice vectors of original structure (rows, Å)
    masses       : (N,)  atomic masses, used only for COM diagnostic

    Returns
    -------
    per_atom_disp : (N,)  scalar displacement magnitudes (Å)
    com_shift     : (3,)  mass-weighted mean displacement vector (Å),
                          returned for reporting only
    """
    delta = cart_cg - cart_orig

    # MIC: convert to fractional of original lattice, wrap, convert back
    inv_lat    = np.linalg.inv(lattice_orig)
    delta_frac = delta @ inv_lat
    delta_frac -= np.round(delta_frac)
    delta      = delta_frac @ lattice_orig

    # COM shift — diagnostic only, not subtracted
    com_shift = np.average(delta, axis=0, weights=masses)

    per_atom_disp = np.sqrt(np.sum(delta**2, axis=1))
    return per_atom_disp, com_shift


# =============================================================================
# MAIN ANALYSIS
# =============================================================================

def run_analysis(orig_file: str, cg_file: str, out: OutputHandler):

    # -------------------------------------------------------------------------
    # Load structures
    # -------------------------------------------------------------------------
    out.print("=" * 55)
    out.print("  Structure Comparison: Original vs CG Relaxed")
    out.print("=" * 55)

    orig_lmp = LammpsData.from_file(orig_file, atom_style="atomic")
    cg_lmp   = LammpsData.from_file(cg_file,   atom_style="atomic")

    orig = orig_lmp.structure
    cg   = cg_lmp.structure

    N = len(orig)
    assert len(orig) == len(cg), \
        f"Atom count mismatch: {len(orig)} vs {len(cg)}"

    out.print(f"\nOriginal file   : {orig_file}")
    out.print(f"CG file         : {cg_file}")
    out.print(f"Number of atoms : {N}")
    out.print(f"Species present : {set(str(s.specie) for s in orig)}")

    masses         = np.array([site.specie.atomic_mass for site in orig])
    species        = [str(site.specie) for site in orig]
    unique_species = sorted(set(species))
    lattice_orig   = orig.lattice.matrix    # (3,3), rows = lattice vectors, Å

    # -------------------------------------------------------------------------
    # Volume
    # -------------------------------------------------------------------------
    V_orig    = orig.volume
    V_cg      = cg.volume
    dV_over_V = (V_cg - V_orig) / V_orig

    out.print(f"\n--- Volume ---")
    out.print(f"  Original : {V_orig:.4f} Å³")
    out.print(f"  CG       : {V_cg:.4f} Å³")
    out.print(f"  ΔV/V     : {dV_over_V:+.4f}  ({dV_over_V * 100:+.2f}%)")

    # -------------------------------------------------------------------------
    # Lattice parameters
    # -------------------------------------------------------------------------
    lo = orig.lattice
    lc = cg.lattice

    out.print(f"\n--- Lattice Parameters ---")
    out.print(f"  {'Param':<8} {'Original':>12} {'CG':>12} {'Δ':>12}")
    out.print(f"  {'-'*46}")
    for name, vo, vc in [
        ("a (Å)",   lo.a,     lc.a),
        ("b (Å)",   lo.b,     lc.b),
        ("c (Å)",   lo.c,     lc.c),
        ("α (°)",   lo.alpha, lc.alpha),
        ("β (°)",   lo.beta,  lc.beta),
        ("γ (°)",   lo.gamma, lc.gamma),
    ]:
        out.print(f"  {name:<8} {vo:>12.5f} {vc:>12.5f} {vc - vo:>+12.5f}")

    # -------------------------------------------------------------------------
    # Displacements
    # -------------------------------------------------------------------------
    cart_orig = np.array([site.coords for site in orig])
    cart_cg   = np.array([site.coords for site in cg])

    per_atom_disp, com_shift = compute_displacements(
        cart_orig, cart_cg, lattice_orig, masses
    )

    msd       = np.mean(per_atom_disp**2)
    mean_disp = np.mean(per_atom_disp)

    # COM shift — diagnostic only
    out.print(f"\n--- Center of Mass Shift (diagnostic only) ---")
    out.print(f"  shift     : {com_shift} Å")
    out.print(f"  magnitude : {np.linalg.norm(com_shift):.6f} Å")
    out.print(f"  note      : not subtracted — expected to be numerical artifact")

    out.print(f"\n--- Cartesian Displacements (MIC via original lattice) ---")
    out.print(f"  MSD (Å²)              : {msd:.6f}")
    out.print(f"  sqrt(MSD) (Å)         : {np.sqrt(msd):.6f}")
    out.print(f"  Mean displacement (Å) : {mean_disp:.6f}")

    # -------------------------------------------------------------------------
    # Per-species MSD
    # -------------------------------------------------------------------------
    if len(unique_species) > 1:
        out.print(f"\n--- Per-Species MSD ---")
        out.print(f"  {'Species':<8} {'MSD (Å²)':>14} {'sqrt(MSD) (Å)':>16} {'n':>6}")
        out.print(f"  {'-'*48}")
        for sp in unique_species:
            idx    = [i for i, s in enumerate(species) if s == sp]
            sp_msd = np.mean(per_atom_disp[idx]**2)
            out.print(
                f"  {sp:<8} {sp_msd:>14.6f} {np.sqrt(sp_msd):>16.6f} {len(idx):>6}"
            )

    # -------------------------------------------------------------------------
    # StructureMatcher
    # -------------------------------------------------------------------------
    out.print(f"\n--- StructureMatcher (symmetry-aware) ---")
    matcher = StructureMatcher()

    fit = matcher.fit(orig, cg)
    out.print(f"  Structures match : {fit}")

    rms_result = matcher.get_rms_dist(orig, cg)
    if rms_result is not None:
        rms, max_dist = rms_result
        out.print(f"  RMS distance (normalized)  : {rms:.6f}")
        out.print(f"  Max distance (normalized)  : {max_dist:.6f}")
    else:
        out.print("  RMS distance : could not be computed (structures too dissimilar)")

    # -------------------------------------------------------------------------
    # Top 5 most displaced atoms
    # -------------------------------------------------------------------------
    top5_idx = np.argsort(per_atom_disp)[::-1][:5]

    out.print(f"\n--- Top 5 Most Displaced Atoms ---")
    out.print(f"  {'Index':>6}  {'Species':>8}  {'|Δr| (Å)':>10}")
    for idx in top5_idx:
        out.print(f"  {idx:>6}  {species[idx]:>8}  {per_atom_disp[idx]:>10.6f}")

    out.print("\nDone.")


# =============================================================================
# ENTRY POINT
# =============================================================================

if __name__ == "__main__":
    args = parse_args()

    with OutputHandler(
        export   = not args.no_export,
        filepath = args.output
    ) as out:
        run_analysis(args.orig_file, args.cg_file, out)