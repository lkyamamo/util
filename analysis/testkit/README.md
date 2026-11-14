# testkit

Shared helpers for the analysis tests. Not analysis code: nothing here
produces a result, so the analysis rules (run-by-test, README per script)
apply to the tests that use it, not to it.

| Module | What it gives a test |
|---|---|
| `systems.py` | Synthetic systems with analytic answers: `ideal_gas` (g(r) = 1), `random_walk` (MSD = 6Dt), `harmonic_oscillators` (VDOS peak at ω), `dipole_series` (set ⟨M²⟩ − ⟨M⟩²), `lattice_from_data` (a crystal from `starting-structures/`) |
| `lammps_io.py` | `write_dump`, `read_dump`, `write_data` |
| `scripts.py` | `run_script` (subprocess with this interpreter, failing with its output) and `load_script` (import from a file in a non-package folder) |
| `checks.py` | `Checker`, behind the `checker` fixture in `analysis/conftest.py`: record every check, fail once with all of them |
| `data/` | Small committed reference inputs (under 1 MB each), for regression and cross-checks |

## Verification

`test_testkit.py` checks each generator against the answer it claims, since
every test built on a wrong generator inherits the error.

| Test | Method | Checked | Tolerance |
|---|---|---|---|
| `test_random_walk_msd_is_6Dt` | known answer | MSD(t) = 6Dt at every lag | 5 standard errors, √(2/(3N)) relative |
| `test_ideal_gas_is_uniform` | invariant | All positions in [0, L); flat histogram | each bin within 5σ (Poisson) |
| `test_harmonic_oscillator_spectrum_peaks_at_omega` | known answer, invariant | Velocity spectrum peak at ω/2π; v = dx/dt | one frequency bin; the central-difference truncation error Aω(ωΔt)²/6 |
| `test_dipole_series_fluctuation` | known answer | ⟨\|M\|²⟩ − \|⟨M⟩\|² and ⟨M⟩ as set | 5 standard errors |
| `test_dump_round_trip` | known answer | `write_dump` → `read_dump` keeps timesteps, box, elements and wrapped positions | 1e-6 (written precision) |
| `test_data_file_read_by_ase` | independent implementation | ASE reads `write_data` output back to the same positions, cell and types | float equality |
| `test_alpha_cristobalite_lattice` | known answer | 4 Si + 8 O; nearest Si–O 1.60–1.61 Å (Downs & Palmer, Am. Mineral. 79, 9 (1994)) | 1.59–1.63 Å |

### Last verified

<!-- Rewritten by tools/repo_check.py on every run that verifies this folder; never appended. -->
- Fast: 2026-11-13 on top of `875f945`: `python -m pytest -q -p no:anyio` in `analysis/`: 20 passed (9.1 s)
- Slow: none defined
