# Shock voxel analysis

Per-voxel properties of a bubble-collapse trajectory: density, temperature,
pressure, speeds, species and the silica surface, then a merged trajectory with
the jet tip and crater geometry.

| File | What it is |
|---|---|
| `voxel_analysis.py` | One dump frame → one h5 of per-voxel properties (NaN-aware 7-voxel smoothing). `voxel_analysis.slurm` runs it per frame. |
| `merge_h5.py` | Per-frame h5 files → `trajectory.h5`, with the crater reference, depth map and sphere fit. `merge_h5.slurm` runs it. |
| `plot_trajectory.py`, `shock_vis_pyvista.py` | Plots and 3D views of `trajectory.h5`. |
| `submit_voxel_pipeline.sh` | Submits the per-frame jobs and the merge. |
| `etc/` | The pipeline test plus design notes and a hand-driven notebook. |

## Verification

`etc/test_pipeline.py` tests `voxel_analysis.py` and `merge_h5.py` together: it
writes two synthetic frames on a 5×6×6 voxel grid (10 Å voxels), runs both
scripts, and compares every property of the five populated voxels to values
computed independently in the test.

| Test | Method | Checked | Tolerance |
|---|---|---|---|
| `test_smoothed_voxel_properties` (both frames) | known answer, independent implementation | Raw density, temperature (3N−3 dof), kinetic + virial pressure, speeds and number density computed by hand per voxel, then the same NaN-aware `uniform_filter` smoothing re-implemented in the test; voxel types exact | 1e-3 in each property's units |
| `test_frame_scalars` (both frames) | known answer | No false jet tip on a sparse low-density grid (`jet_tip_x` NaN); Si surface at x = 5 then 8 Å; one hydronium; exactly one surface voxel | exact; surface 1e-3 Å |
| `test_crater_datasets` | known answer | Crater dataset shapes; reference x ≈ 5 and 8 Å; no sphere fit when no bin is deeper than the cutoff | 0.1 Å (voxel-binned) |
| `test_jet_tip_detection` | known answer | `detect_jet_tip` on a hand-built grid: excludes outside the rod and below the density threshold, picks the smallest qualifying x (25 Å) | 1e-3 Å |

- **Invariants:** none tested yet. Translation of the frame by a whole number of
  voxels should shift every grid; worth adding when this code next changes.
- **Regression:** none. The tracked `test_output/` folder here is a stale
  artifact of the old test (debt), not a baseline.
- `plot_trajectory.py` and `shock_vis_pyvista.py` draw figures and have no tests.
- `submit_voxel_pipeline.sh` only submits the SLURM jobs; it has no test beyond the gate's `bash -n`.

### Last verified

<!-- Rewritten by tools/repo_check.py on every run that verifies this folder; never appended. -->
- Fast: 2026-11-13 on top of `00e7f05`: `python -m pytest -q -p no:anyio` in `analysis/`: 20 passed (9.0 s)
- Slow: none defined
