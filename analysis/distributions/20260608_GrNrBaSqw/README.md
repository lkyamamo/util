# Structural Analysis Pipeline

Scripts for computing structural and dynamical properties from a LAMMPS MD trajectory.
Each script produces independent output and can be run in any order or simultaneously.
`rdf_freud.py`/`bad_freud.py` read a structural trajectory (`dump.lammpstrj`);
`dsf.py`/`vdos.py`/`msd.py` read a separate, higher-frequency trajectory
(`dynamics.lammpstrj`) — resolving vibrational frequencies and diffusion
needs much finer time sampling than structural analysis does. See "All
scripts — trajectory input" below.

---

## Scripts

| Script | What it computes | Output files |
|---|---|---|
| `rdf_freud.py` | Radial distribution function g(r), coordination number n(r), and selectable neutron correlation functions for all element pairs | `rdfs.csv`, `rdfs.png`, `nrs.csv`, `nrs.png` |
| `bad_freud.py` | Bond angle distribution P(θ) for all A-B-C triplets | `bads.csv`, `bads.png` |
| `dsf.py` | Static structure factor S(q) and dynamic structure factor S(q,ω) | `sq.csv`, `sq.png`, `dsf.csv`, `dsf.png` |
| `vdos.py` | Vibrational density of states (aligned with `analysis/dynamics/src/msd.cpp` by default) | `vdos.csv`, `vdos.png` |
| `msd.py` | Mean square displacement and self-diffusion coefficient (10⁻⁵ cm²/s) per element | `msd.csv`, `msd.png` |

---

## How to Run

### On the cluster (SLURM)

```bash
sbatch distribution_submit.slurm
```

The slurm script runs all five Python scripts in sequence under the same job.
Toggle which scripts run at the top of `distribution_submit.slurm`:

```bash
RUN_DSF=1   # 1 = run, 0 = skip
RUN_RDF=1
RUN_BAD=1
RUN_VDOS=1
RUN_MSD=1
```

### Locally (all cores)

```bash
./distribution_run.sh          # all cores
./distribution_run.sh 8        # 8 threads
```

`distribution_run.sh` runs all five scripts, mirroring `distribution_submit.slurm`'s
`RUN_DSF`/`RUN_RDF`/`RUN_BAD`/`RUN_VDOS`/`RUN_MSD` flags (edit them at the top of
the script to skip any). To run a single script instead, call it directly:

```bash
python rdf_freud.py
python bad_freud.py
python vdos.py
python msd.py
```

---

## What to Modify Before Running

### All scripts — trajectory input

There are two standardized trajectory files, each read by env var (`TRAJ` /
`DYNAMICS_TRAJ`) into each script's own `DUMP_FILE` variable:

| Script | Trajectory | Env var | Default |
|---|---|---|---|
| `rdf_freud.py` | Structural | `TRAJ` | `"dump.lammpstrj"` |
| `bad_freud.py` | Structural | `TRAJ` | `"dump.lammpstrj"` |
| `dsf.py` | Dynamics (higher-frequency) | `DYNAMICS_TRAJ` | `"dynamics.lammpstrj"` |
| `vdos.py` | Dynamics (higher-frequency) | `DYNAMICS_TRAJ` | `"dynamics.lammpstrj"` |
| `msd.py` | Dynamics (higher-frequency) | `DYNAMICS_TRAJ` | `"dynamics.lammpstrj"` |

`dsf.py`/`vdos.py`/`msd.py` need `dynamics.lammpstrj` sampled much more often
than `dump.lammpstrj` — resolving vibrational frequencies and diffusion
requires a time step between frames well below the period of the fastest
mode you care about (Nyquist limit), whereas structural analysis (g(r), P(θ))
only needs occasional snapshots. `OH-therm.input`/`b-SiO-therm.input` write
both dumps from the same production run, at independent frequencies
(`dump_frequency` vs. `dynamics_dump_frequency`).

---

### `rdf_freud.py` — required changes

| Variable | What it controls | Notes |
|---|---|---|
| `DUMP_FILE` | Trajectory path | `$TRAJ` |
| `R_MAX` | Max r in Å for g(r) | `$R_MAX`; must be < half the shortest box dimension |

Column positions are read from the `ITEM: ATOMS` header line, so there is nothing to set by hand.
Wrapped coordinates (`x y z`) and a constant atom count are assumed; a triclinic dump aborts with an
explanation rather than being silently mis-wrapped.

Optional:

| Variable | What it controls | Default |
|---|---|---|
| `RDF_BINS` | Number of r-bins | `2000` |
| `RDF_NORMALIZATION` | Pair-weight convention(s): `unity`, `FZ`, `absolute` | `"FZ;absolute;unity"` |
| `RDF_FUNCTIONS` | Correlation function(s): `g`, `h`, `D`, `T` | `"g;h;D"` |
| `RDF_RESOLUTION_SIGMA` | Gaussian resolution broadening in Å, written as `*_broadened` twins; `0` disables | `0.1` |
| `OUTPUT_CSV` | g(r) CSV path; `None` to skip | `"rdfs.csv"` |
| `OUTPUT_NR_CSV` | n(r) CSV path; `None` to skip | `"nrs.csv"` |

The two list-valued keys are **semicolon-separated**, not comma-separated: `submit_pipeline.sh` passes
settings through `sbatch --export`, which is itself comma-delimited and silently truncates a value at the
first embedded comma — so a comma-separated list would reach the cluster as its first entry alone, with no
error. `bad_freud.py` takes `ELEMENTS="Si;O;H"` for the same reason. A comma in either key is rejected with
the corrected string in the message.

Every key is `RDF_`-prefixed, env key and Python constant alike, so the analysis it configures is explicit:
`NORMALIZATION` on its own does not say which distribution it belongs to, and `vdos.py` has its own
unrelated `VDOS_NORMALIZATION` (`phonon` | `unit_area`) in the same shared environment.

**Conventions.** `RDF_NORMALIZATION` × `RDF_FUNCTIONS` is emitted as a cross product, one column per combination
named `<function>_<normalization>`. The script prints every column's defining equation, weight sum Σw, and
asymptotic limits at startup; check a printed limit against the curve rather than trusting the symbol
(Keen, *J. Appl. Cryst.* **34**, 172 (2001) tabulates why this matters). Columns written before this
became configurable map as `total → g_unity`, `neutron → g_FZ`, `t → T_FZ`.

*Normalization* — how much each partial contributes, via the pair weight `w_AB` (`f = 2 − δ_AB`):

| key | `w_AB` | Σw | units | what it is |
|---|---|---|---|---|
| `unity` | `f c_A c_B` | 1 | — | Every element scatters identically (b = 1) — the name refers to the scattering lengths, not to Σw, since FZ also sums to 1. Not measurable; it's the composition-averaged structure and the b-free baseline the neutron curves depart from. Formerly `total`. |
| `FZ` | `f c_A c_B b_A b_B / ⟨b⟩²` | 1 | — | Faber–Ziman. Tends to 1 at large r like a partial does, and being dimensionless it superimposes across compositions — at the cost of dividing by a nearly-cancelling sum. |
| `absolute` | `f c_A c_B b_A b_B / 100` | ⟨b⟩²/100 | barn/sr/atom | No division: the weighted sum in the units a measured differential cross-section carries. The scale is physical, so the excluded-volume plateau lands at −Σw. Stays well conditioned as ⟨b⟩ → 0, so use it for light water (⟨b⟩² = 0.0031 barn) or a null mixture (⟨b⟩ = 0), where FZ is useless. |

*Function* — what gets built from those weights and the partials `g_AB`:

| key | definition | what it is |
|---|---|---|
| `g` | `Σ w g_AB` | Weighted pair distribution; baseline at Σw. |
| `h` | `Σ w [g_AB − 1]` | Total correlation function: baseline subtracted first, so peaks sit on zero and the excluded-volume region reads −Σw. **`h_absolute` is Soper's eq. (20)**, the neutron G_n(r). |
| `D` | `4πrρ Σ w [g_AB − 1]` | Differential correlation function. The 4πr factor offsets the decay of peak amplitude with distance, keeping far-field oscillations legible; oscillates about 0. **`D_FZ` is the PDF community's G(r)**. |
| `T` | `4πrρ Σ w g_AB` | Total radial distribution function — `D` keeping the bulk baseline, so it climbs as 4πrρΣw. Area under a peak is a coordination number. Not in the default `RDF_FUNCTIONS`. |

**Two things worth knowing.** `'H'` in `NEUTRON_SCATTERING_LENGTHS` is *deuterium* (b = 6.671 fm), since
LAMMPS labels both isotopes `H`; every neutron column is therefore for a deuterated sample, and the script
says so at startup. And `n(r)` columns are named by direction — `O_around_Si` ≈ 4 for silica,
`Si_around_O` its reciprocal — because freud counts system points around query points, which makes a bare
`Si-O` label ambiguous.

---

### `bad_freud.py` — required changes

| Variable | What it controls | Notes |
|---|---|---|
| `DUMP_FILE` | Trajectory path | |
| `ELEMENTS` | List of element symbols in the simulation | e.g. `['Si', 'O', 'H']` — code never auto-detects |
| `R_CUTOFF` | Upper neighbor cutoff per pair (Å) | Keys must be `"El1-El2"` sorted alphabetically; must include every pair in `ELEMENTS` |
| `R_MINCUT` | Lower neighbor cutoff per pair (Å) | Same key format; excludes unphysical close contacts |
| `COL_ELEMENT` | Column index of element symbol | Default `1` for layout `id element x y z` |
| `COL_X`, `COL_Y`, `COL_Z` | Column indices of coordinates | Default `2, 3, 4` |

Optional:

| Variable | What it controls | Default |
|---|---|---|
| `TRIPLETS` | Restrict which A-B-C triplets to compute; `None` = all | `None` |
| `BINS` | Number of angle bins (0–180°) | `180` |

---

### `dsf.py` — required changes

| Variable | What it controls | Notes |
|---|---|---|
| `DUMP_FILE` | Dynamics trajectory path (env var `DYNAMICS_TRAJ`, default `dynamics.lammpstrj`) | Requires `element` column in dump (not numeric type) |
| `DT` | Time between consecutive dumped frames in **femtoseconds** | e.g. if LAMMPS dumps every 100 steps at 0.5 fs/step → `DT = 50.0` |
| `N_FRAMES` | Max frames to read | Controls how much of the trajectory is used |
| `WINDOW_SIZE` | Number of time lags for F(q,t) | Sets frequency resolution: Δω ∝ 1/(WINDOW_SIZE × DT); set equal to `N_FRAMES` to use full trajectory |

Optional:

| Variable | What it controls | Default |
|---|---|---|
| `STRIDE` | Read every Nth frame | `1` |
| `Q_MAX` | Max q in Å⁻¹ | `20.0` |
| `N_Q_BINS` | Radial q-bins after spherical averaging | `200` |
| `COMPUTE_STATIC` | Compute S(q) | `True` |
| `COMPUTE_DYNAMIC` | Compute S(q,ω) | `True` |
| `COMPUTE_SELF` | Compute incoherent/self part (slow) | `False` |
| `N_THREADS` | numba thread count; `0` = all cores | `0` |

---

### `vdos.py` — required changes

| Variable | What it controls | Notes |
|---|---|---|
| `DUMP_FILE` | Dynamics trajectory path (env var `DYNAMICS_TRAJ`, default `dynamics.lammpstrj`) | Requires `element vx vy vz` columns (positions not read) |
| `TIME_UNIT` | Time between consecutive dumped frames in **femtoseconds** | Falls back to `DT` if unset — same quantity dsf.py calls `DT`. **The only env var vdos.py and msd.py share**; every other vdos.py setting is read from a `VDOS_`-prefixed name |

Optional (see the module docstring for the full METHOD/WINDOW/NORMALIZATION
explanation):

| Variable (env var) | What it controls | Default |
|---|---|---|
| `METHOD` (`VDOS_METHOD`) | `'vacf_cosine_transform'` (matches `analysis/dynamics/src/msd.cpp`) or `'fft_periodogram'` (faster, this script's own approach) | `'vacf_cosine_transform'` |
| `CORR_LENGTH` (`VDOS_CORR_LENGTH`) | VACF max lag / Welch-segment length, in fs | 75% of total trajectory duration |
| `CORR_INTERVAL` (`VDOS_CORR_INTERVAL`) | Spacing between VACF reference frames / segment starts, in fs | 10% of `CORR_LENGTH` |
| `MAX_FREQUENCY_EV` (`VDOS_MAX_FREQUENCY_EV`) | Upper frequency limit of the output grid, in eV (`vacf_cosine_transform` only) | `0.1` |
| `NUM_GRIDS` (`VDOS_NUM_GRIDS`) | Frequency grid points (`vacf_cosine_transform` only) | `5000` |
| `WINDOW` (`VDOS_WINDOW`) | `'cosine_lag'`/`'none'` under `vacf_cosine_transform`; `'hann'`/`'none'` under `fft_periodogram` | Matches `METHOD` |
| `VDOS_NORMALIZATION` | Sum rule: `'phonon'` (∫ ≈ 3 per atom, matches msd.cpp) or `'unit_area'` (∫ = 1) | `'phonon'` |
| `VDOS_WEIGHTING` | Species weighting, semicolon-separated: `unity`, `coherent`, `incoherent`, `total` | `"unity"` |
| `N_FRAMES`, `STRIDE` (`VDOS_N_FRAMES`, `VDOS_STRIDE`) | Max frames to read / read every Nth frame | `0` (all), `1` |
| `VDOS_THREADS` | scipy FFT thread count (`fft_periodogram` only); `0` = all cores | `0` |

**Weighting.** `VDOS_WEIGHTING` decides how much each element contributes to a total; one
`DoS(Total_<weighting>)` column is emitted per entry. It is a *different axis* from
`VDOS_NORMALIZATION`, which sets the sum rule — weighting is relative species contribution, normalization
is overall scale.

| key | `w_el` | what it is |
|---|---|---|
| `unity` | `c_el` | Mole fractions — no scattering physics. Exactly what this script produced before weighting existed. |
| `coherent` | `c_el · σ_coh,el / m_el` | |
| `incoherent` | `c_el · σ_inc,el / m_el` | |
| `total` | `c_el · (σ_coh,el + σ_inc,el) / m_el` | What a chopper spectrometer collects; the usual generalized-DOS weighting. |

Weights are normalized to Σw = 1, so the `phonon` 3-per-atom sum rule survives and every total stays
comparable to the unity-weighted one. The resolved per-element weights are printed at startup — check
them against the curves rather than trusting the column name. For SiO₂: `unity` gives Si 0.333 / O 0.667,
while `coherent` and `total` give Si 0.127 / O 0.873.

**Why σ/m and not `b`.** Inelastic scattering measures the generalized DOS, in which the one-phonon
incoherent cross-section carries a factor σ/m per species. That is a different quantity from the coherent
scattering length that weights diffraction, so `rdf_freud.py`'s and `dsf.py`'s `b`-weighting does not
transfer here — σ_inc cannot be derived from b_coh at all.

**Hydrogen is refused.** With `H` or `D` present, any weighting other than `unity` exits with an
explanation instead of a number: protium and deuterium differ by ~21× in σ/m, a LAMMPS dump labels both
`H`, and under `incoherent` weighting H would carry >99.9% of the weight — so the result would be set
almost entirely by the species whose treatment is undecided. `unity` still works on those systems.

**Debye–Waller is not applied**, and the startup output says so. At fixed Q it would be a pure per-species
re-weighting, but Q and energy transfer are kinematically coupled in a real spectrometer while this DOS is
not Q-resolved, and the harmonic fixed-site assumption behind ⟨u²⟩ fails for diffusing species.

The semicolon separator is required for the same reason as `RDF_NORMALIZATION` — `sbatch --export` is
comma-delimited and would truncate the value. A comma is rejected with the corrected string.

> **Behavior change:** under `VDOS_NORMALIZATION=unit_area` the total was previously an unweighted sum of
> the raw per-element curves, which gave every element equal weight regardless of atom count — inconsistent
> with the `phonon` branch's mole fractions. It is now the same weighted sum as `phonon`. On a 1:2 Si:O
> test this moved the O:Si peak ratio from 1.00 to 1.96 (up to 32% change in the curve). `phonon` output is
> unchanged.

---

### `msd.py` — required changes

| Variable | What it controls | Notes |
|---|---|---|
| `DUMP_FILE` | Dynamics trajectory path (env var `DYNAMICS_TRAJ`, default `dynamics.lammpstrj`) | Requires `element x y z` columns (wrapped, not `xu yu zu`; velocities not read) |
| `TIME_UNIT` | Time between consecutive dumped frames in **femtoseconds** | Falls back to `DT` if unset — same quantity dsf.py/vdos.py use. **The only env var msd.py and vdos.py share**; every other msd.py setting is read from an `MSD_`-prefixed name |

Optional (see the module docstring for the full METHOD explanation — no
unwrapped coordinates needed, minimum-image correction on each
reference-to-current displacement instead):

| Variable (env var) | What it controls | Default |
|---|---|---|
| `CORR_LENGTH` (`MSD_CORR_LENGTH`) | Max time lag / reference-frame spacing window, in fs — same meaning as, but set independently of, vdos.py's `CORR_LENGTH` | 75% of total trajectory duration |
| `CORR_INTERVAL` (`MSD_CORR_INTERVAL`) | Spacing between reference frames, in fs — same meaning as, but set independently of, vdos.py's `CORR_INTERVAL` | 10% of `CORR_LENGTH` |
| `MSD_FIT_FRACTION` | Fraction of the tail of the correlation window used for the diffusion-coefficient linear fit (excludes the early ballistic regime) | `0.5` |
| `N_FRAMES`, `STRIDE` (`MSD_N_FRAMES`, `MSD_STRIDE`) | Max frames to read / read every Nth frame | `0` (all), `1` |

Prints a self-diffusion coefficient per element (and total) to the console,
in units of 10⁻⁵ cm²/s (the standard unit for reporting liquid/solid
self-diffusion coefficients), via the 3D Einstein relation MSD(t) = 6Dt.

---

### `distribution_submit.slurm` — required changes

| Item | Location | Notes |
|---|---|---|
| `RUN_DSF` / `RUN_RDF` / `RUN_BAD` / `RUN_VDOS` / `RUN_MSD` | Lines 38–42 (run flags) | Set `1` to run, `0` to skip |
| `--time` | SBATCH header | Increase for large trajectories |
| `--cpus-per-task` | SBATCH header | Sets thread count; passed to `OMP_NUM_THREADS` and `VDOS_THREADS` |
| venv path | `source /home1/lkyamamo/venv/struc_analysis/bin/activate` | Update if environment moves |

---

### `distribution_run.sh` — required changes

No changes required for HPC use — activates the same venv as `distribution_submit.slurm`:
`/home1/lkyamamo/venv/struc_analysis/bin/activate`

---

## Environment

### HPC (cluster)

The shared venv on the HPC contains all dependencies:

```bash
source /home1/lkyamamo/venv/struc_analysis/bin/activate
```

This is what `distribution_submit.slurm` activates automatically.

### Local install

```bash
pip install dynasor matplotlib freud
pip install icc_rt    # optional: 5–10× speedup for dynasor's numba backend
pip install scipy     # optional: multi-threaded FFT for vdos.py's fft_periodogram method
```

---

## Dump Format Requirements

| Script | Element column | Coordinate columns | Notes |
|---|---|---|---|
| `rdf_freud.py` | `element` (symbol) — **required** | `x y z` (wrapped, real, Å) | Column layout auto-detected from `ITEM: ATOMS` header; orthogonal box and constant atom count required |
| `bad_freud.py` | `element` (symbol) or any string | `x y z` (real, Å) | Column indices set manually via `COL_*` |
| `dsf.py` | `element` (symbol) — **required** | `x y z` or `xs ys zs` or `xu yu zu` | Reads `dynamics.lammpstrj`; column layout auto-detected from `ITEM: ATOMS` header |
| `vdos.py` | `element` (symbol) — **required** | `vx vy vz` (velocities) — **required**; positions not read | Reads `dynamics.lammpstrj`; column layout auto-detected from `ITEM: ATOMS` header |
| `msd.py` | `element` (symbol) — **required** | `x y z` (wrapped, real, Å) — **required**; velocities not read | Reads `dynamics.lammpstrj`; column layout auto-detected from `ITEM: ATOMS` header |
