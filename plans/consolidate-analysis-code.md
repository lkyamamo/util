# Move all analysis Python into analysis/, with verification

**State:** draft (approve after plan 1 lands; it depends on plan 1's testkit, pytest setup and gate)

## Goal

Every Python script that turns simulation output into a number you report lives in `analysis/`, so the analysis
verification rules from plan 1 cover all of it. Each moved script leaves with tests (all four methods) and a
README `## Verification` section. The two diverged copies of the dielectric code become one. Every run records
which commit of the analysis code it used.

## Scope

**Moves (`git mv`, then the edits below):**

| From | To |
|---|---|
| `simulation/lammps/20260617_dielectric_multi_traj/1.calc_mpi.py` | `analysis/dielectric/calc_mpi.py` |
| `simulation/lammps/20260617_dielectric_multi_traj/2.dipole_std.py` | `analysis/dielectric/dipole_std.py` |
| `simulation/lammps/20260929_dielectric_chunked/2.calc_mpi.py`, `3.dipole_std.py` | deleted (merged into the above) |
| `jobs/pipeline/lammps_to_temperature_sweep/aggregate_dielectric_vs_temperature.py` | `analysis/dielectric/` |
| `jobs/pipeline/lammps_to_temperature_sweep/aggregate_diffusion_vs_temperature.py` | `analysis/dynamics/` |
| `plotting/lammps/thermodynamics.py`, `plotting/lammps/csv_combine.py` | `analysis/thermo/` |
| `plotting/vasp/parse_oszicar.py` | `analysis/energetics/` |
| `simulation/lammps/structure_comparison/struct_compare.py` | `analysis/structure/` |
| `simulation/lammps/msd_code/lammps_custom_to_msd_input_xyz.py` | merged into `analysis/dynamics/lammps_custom_to_input_xyz.py` (6 lines differ) |
| `simulation/lammps/msd_code/lammps_to_msd_input_xyz.py` | deleted if `analysis/dynamics/1.robust_lammps_to_input_xyz.py` covers it; otherwise to `analysis/dynamics/` |

`plotting/`, `simulation/lammps/msd_code/` and `simulation/lammps/structure_comparison/` are then empty and
removed. The `msd_code` README content goes into `analysis/dynamics/README.md`.

**Callers updated:**

- `jobs/pipeline/lammps_to_temperature_sweep/submit_temperature_sweep.sh` and its `.conf.example`: the default
  dielectric template dir becomes `analysis/dielectric/`, and the `AGGREGATE_*_SCRIPT` paths change. Copied
  scripts lose their `1.`/`2.` prefixes.
- `jobs/pipeline/lammps_to_temperature_sweep/dielectric-production.input`: the comments naming `1.calc_mpi.py`.
- `jobs/slurm/lammps_submit.slurm`: `THERMO_SCRIPT`.
- `jobs/slurm/lammps_dielectric_submit.slurm`
- `simulation/lammps/20260617_dielectric_multi_traj/*.slurm`
- `simulation/lammps/20260929_dielectric_chunked/*.slurm` and its README: copy the scripts from
  `analysis/dielectric/` at submit, as below.

**Provenance:** wherever a pipeline or run-folder SLURM script copies analysis code into a run directory, it also
writes `analysis_commit.txt` (`git rev-parse HEAD` plus a dirty flag). The ledger can attach it later; that is
out of scope here.

**New tests, testkit and docs:**

- `analysis/dielectric/test_calc_mpi.py`
- `analysis/dielectric/test_dipole_std.py`
- `analysis/dielectric/test_aggregate_dielectric.py`
- `analysis/dynamics/test_aggregate_diffusion.py`
- `analysis/dynamics/test_lammps_custom_to_input_xyz.py`
- `analysis/thermo/test_thermodynamics.py`
- `analysis/thermo/test_csv_combine.py`
- `analysis/energetics/test_parse_oszicar.py`
- `analysis/structure/test_struct_compare.py`
- `analysis/testkit/` additions: a water-box generator with known molecular dipoles, and log.lammps, OSZICAR and
  sweep-summary writers.
- `analysis/testkit/data/` (new): small real fixtures from endeavour, each under 1 MB:
  - one dielectric dump subset, plus LAMMPS's own `compute dipole` output for it
  - one `log.lammps`
  - one OSZICAR
  - one sweep's summaries
- READMEs with `## Verification`: `analysis/dielectric/README.md`, `analysis/thermo/README.md`,
  `analysis/energetics/README.md`, `analysis/structure/README.md`, and an updated `analysis/dynamics/README.md`.
- `analysis/pyproject.toml`: add mpi4py, polars, log-lammps-reader and pymatgen.
- `CLAUDE.md`: the repo map, with `plotting/` gone and the new `analysis/` folders listed.

## Approach

1. **Fetch the fixtures.** Pull the regression and independent-implementation fixtures from endeavour over SSH,
   and trim them to under 1 MB. You pick the source runs (see Verification).
2. **Write the tests against the old paths first,** and run them on the unmoved code. They must pass before
   anything moves, so the move itself is then verified by the same tests.
3. **Merge the dielectric copies.**
   - The two `calc_mpi.py` copies differ only in usage and error text that names their own SLURM file. The merged
     version names neither.
   - The two `dipole_std.py` copies are identical.
   - The tests run on both old copies and the merged one, and must give identical outputs.
4. **Move and update callers in one commit per destination folder.** Each commit carries its tests, README and
   `Plan:` trailer. A pure `git mv` counts as a rename under the plan 1 gate; the content merges carry tests.
5. **Coordinate the cluster.** `lammps_submit.slurm` runs `THERMO_SCRIPT` at the end of each job, so a job queued
   with the old path breaks if the cluster checkout is pulled first. Before pulling on endeavour:
   - check `squeue -u lkyamamo` for queued jobs that use these scripts
   - pull only when none are queued, or let them finish first

   Existing run directories hold their own copies and are unaffected.

## Verification

### Known answer

- **calc_mpi:** water boxes built with set positions and charges, so the total dipole M = Σ qᵢrᵢ is known exactly:
  - molecules split across the periodic boundary (tests the unwrapping and molecule assignment)
  - a box of molecules all aligned (|M| = N·μ)
  - random orientations
- **dipole_std:** dipole series drawn with a set ⟨M²⟩ − ⟨M⟩². The expected ε comes from the fluctuation formula
  ε = 1 + (⟨M²⟩ − ⟨M⟩²)/(3Vε₀k_BT), with conducting boundaries, worked out by hand in the test.
- **Aggregators:** fake summary files at set temperatures. The output CSV rows must match exactly.
- **thermodynamics:** a synthetic `log.lammps` with thermo columns of known mean over the second half.
- **csv_combine:** numbered CSVs with known rows.
- **parse_oszicar:** a synthetic OSZICAR in the documented line format, with set T/E/F/E0/EK/SP/SK.
- **lammps_custom_to_input_xyz:** a custom dump with known positions and velocities, converted and compared field
  by field.
- **struct_compare:** a structure and a copy with a known uniform strain and known per-atom displacements. ΔV/V,
  the lattice parameters and the MSD must match by hand.

### Invariants

- **calc_mpi:**
  - M does not change under a whole-box translation (the system is neutral) or a permutation of atom order
  - the output is identical with 1 and 2 MPI ranks
  - resuming from checkpoints gives the same result as a straight run
- **dipole_std:**
  - ε does not change when a constant is added to M
  - ε scales as 1/T at a fixed ⟨δM²⟩
- **Aggregators:** the input order of the entries does not change the output.
- **thermodynamics:** the averages do not depend on how the log is split across `log.lammps.resumeN`.
- **struct_compare:** comparing a structure with itself gives zero for everything. Swapping the two files flips
  the sign of ΔV/V.

### Independent implementation

- **calc_mpi** against LAMMPS's own `compute dipole` on the fixture dump.
  `dielectric-production.input` already writes it "for direct comparison".
- **thermodynamics** against simledger's LAMMPS parser (`jobs/ledger/simledger/parsers/lammps.py`), which
  independently takes second-half means of T, P and density from the same log.
- **parse_oszicar** against simledger's VASP parser (`jobs/ledger/simledger/parsers/vasp.py`) for the measured MD
  temperature.
- **lammps_custom_to_input_xyz** against ASE/OVITO reading the same dump.
- **struct_compare** against ASE (volume, cell parameters) on the same pair.
- **dipole_std:** n/a. No second fluctuation-formula implementation exists, so the known-answer case is its
  check.
- **Aggregators and csv_combine:** n/a. They are pure parsing with no second implementation, so known-answer plus
  regression cover them.

### Regression

- Each script's output on its real fixture is pinned in `analysis/<folder>/baselines/<script>.json`, for example
  ε_total and its components from one real sweep temperature.
- A change fails the test until the baseline is updated in the same commit, with the reason in the commit message.
  A baseline change on its own means a result changed, not that the code is wrong.

### Tolerances and sources

- **calc_mpi M:** 1e-10 relative (exact arithmetic in float64).
- **Against LAMMPS `compute dipole`:** 1e-6 relative, since the dump is written with limited precision. The
  achieved error is reported in the README's `### Last verified` block (latest run only).
- **dipole_std ε:** analytic expectation. The tolerance is 5 standard errors of the variance estimate,
  √(2/N)·⟨δM²⟩, with N = 10⁵ samples and a fixed seed.
- **thermodynamics against simledger:** 1e-9 relative. Both read the same numbers; only the summation order
  differs.
- **struct_compare:** 1e-8 Å.
- **Regression baselines:** 1e-12 relative. Deterministic code on the same input should reproduce exactly; the
  tolerance absorbs only floating-point reordering.

### Slow checks

- **calc_mpi rank independence at 4 ranks** on a 5184-atom box (`@pytest.mark.slow`; needs MPI).
- **dipole_std convergence:** ε against trajectory length on the fixture, to show the estimate is converged at
  production length. It runs nightly and records the curve.

### Done when

- All the new tests pass on the old copies before the move, and on the moved copies after.
- The merged and old dielectric code give identical outputs on every fixture.
- `python3 tools/repo_check.py --all` shows no analysis code outside `analysis/`, and no run-by-test or README
  errors for any moved file.
- After the cluster pull, one sweep stage and one `lammps_submit.slurm` job run on endeavour using the new paths
  and write `analysis_commit.txt`.

## Out of scope

- Behavior changes to any moved script beyond the merge and path updates.
- Tests for `analysis/distributions/20260608_GrNrBaSqw/`: plan 3.
- Setup-side scripts that stay put: `simulation/vasp/*` (INCAR/POSCAR helpers), `simulation/lammps/input_description.py`,
  `simulation/lammps/plot_template.py`, and `jobs/pipeline/plotting/` (it draws figures, it doesn't produce numbers).
- Teaching simledger to read `analysis_commit.txt`.

## Open choices for you

- **Fixture runs:** which sweep temperature, `log.lammps` and OSZICAR on endeavour to use (default: the most recent
  completed runs in `20260825_locked`).
- **The small converter:** if `lammps_to_msd_input_xyz.py` turns out not to be covered by the robust converter,
  keep it as a separate script or fold it into the robust one.
