# Structural Analysis Pipeline

Scripts for computing structural and dynamical properties from a LAMMPS MD trajectory.
Each script produces independent output and can be run in any order or simultaneously.
`rdf_freud.py`/`bad_freud.py` read a structural trajectory (`dump.lammpstrj`);
`dsf.py`/`vdos.py`/`msd.py` read a separate, higher-frequency trajectory
(`dynamics.lammpstrj`) — resolving vibrational frequencies and diffusion
needs much finer time sampling than structural analysis does. See "All
scripts — trajectory input" below. `vdos_dynmat.py` is the exception: it reads
no trajectory at all, only the dynamical matrix LAMMPS wrote (plus the first
frame of `dump.lammpstrj`, for element labels).

---

## Scripts

| Script | What it computes | Output files |
|---|---|---|
| `rdf_freud.py` | Radial distribution function g(r), coordination number n(r), and selectable neutron correlation functions for all element pairs | `rdfs.csv`, `rdfs.png`, `nrs.csv`, `nrs.png` |
| `bad_freud.py` | Bond angle distribution P(θ) for all A-B-C triplets | `bads.csv`, `bads.png` |
| `dsf.py` | Static structure factor S(q) and dynamic structure factor S(q,ω) | `sq.csv`, `sq.png`, `dsf.csv`, `dsf.png` |
| `vdos.py` | Vibrational density of states (aligned with `analysis/dynamics/src/msd.cpp` by default) | `vdos.csv`, `vdos.png` |
| `msd.py` | Mean square displacement and self-diffusion coefficient (10⁻⁵ cm²/s) per element | `msd.csv`, `msd.png` |
| `vdos_dynmat.py` | Vibrational density of states from the LAMMPS dynamical matrix — harmonic, 0 K, no trajectory; optionally the stretch/bend/rock band assignment, participation ratio and boson peak | `vdos_dynmat.csv`, `vdos_dynmat.png`, `vdos_dynmat_modes.csv`, `vdos_dynmat_character.png` |

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
RUN_VDOS_DYNMAT=0   # defaults off: needs a dynmat.dat from the LAMMPS stage
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
python vdos_dynmat.py
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
| `DYNAMICS_DT` | Time between consecutive dumped frames in **femtoseconds** | The *same* key `vdos.py` and `msd.py` read — all three analyse `dynamics.lammpstrj`, so its dt is one number. Optional here (defaults to 2.0) and used only on the dynamic path; static S(q) has no time axis |
| `DSF_N_FRAMES` | `frame_stop`, an **index** into the dump — not a count | Frames used = `DSF_N_FRAMES / DSF_STRIDE` |
| `DSF_WINDOW_SIZE` | Number of time lags for F(q,t) | Sets frequency resolution: Δν = 1/(2 × WINDOW_SIZE × DT × STRIDE); must cover several periods of the slowest mode |

Optional:

| Variable | What it controls | Default |
|---|---|---|
| `DSF_STRIDE` | Read every Nth frame | `600` |
| `DSF_Q_MAX` | Max q in Å⁻¹ (static) | `20.0` |
| `DSF_N_Q_BINS` | Radial q-bins after spherical averaging | `130` |
| `DSF_WINDOW_STEP`, `DSF_Q_MAX_DYN`, `DSF_N_Q_BINS_DYN`, `DSF_MAX_Q_POINTS_DYN` | Dynamic-only q/window settings | `1`, `4.0`, `25`, `25000` |
| `COMPUTE_STATIC` | Compute S(q) | `True` |
| `COMPUTE_DYNAMIC` | Compute S(q,ω) | `True` |
| `COMPUTE_SELF` | Compute incoherent/self part (slow) | `False` |
| `DSF_NEUTRON_WEIGHTING` | Emit the neutron-weighted S(q) columns: `yes` or `no` | `yes` |
| `N_THREADS` | numba thread count; `0` = all cores | `0` |

Every `dsf.py` setting is `DSF_`-prefixed, like `vdos.py`'s `VDOS_*` and `msd.py`'s `MSD_*`. It previously
read bare `DT`, `N_FRAMES`, `STRIDE`, `Q_MAX`, `N_Q_BINS`, `WINDOW_SIZE`, `WINDOW_STEP`, `Q_MAX_DYN`,
`N_Q_BINS_DYN` and `MAX_Q_POINTS_DYN` — generic enough to collide with anything else in a shared pipeline
environment, and requiring `submit_pipeline.sh` to translate them on the way in. Setting one of the old
names now **aborts with a rename notice** rather than being silently ignored.

The one exception is the frame spacing, which comes from `DYNAMICS_DT` rather than a `DSF_` key: `dsf.py`,
`vdos.py` and `msd.py` all read `dynamics.lammpstrj`, so its dt is a property of the trajectory and giving
each script its own copy would only let them disagree about one physical number — which would shift the
frequency axis of S(q,ω) relative to a VDOS built from the very same frames. The short-lived `DSF_DT` aborts
with the same notice.

**Which normalization `Sq_neutron` is.** dynasor weights partials as `S_AB → f_A f_B S_AB` with `f = b_coh`
and sums, with **no division by ⟨b⟩²**. So `Sq_neutron` is the unnormalized weighted sum in fm² — the
reciprocal-space counterpart of `rdf_freud.py`'s `absolute` convention, *not* of its `FZ`. Compare it
against `g_absolute`/`h_absolute`, not `g_FZ`. There is no convention selector here yet, because `dsf.py`
delegates the weighting arithmetic to dynasor rather than owning it.

**Hydrogen is refused**, as in `vdos.py` and `vdos_dynmat.py`, but for a sharper reason: dynasor weights by
**natural abundance**, so its `H` is protium at `b_coh = −3.7406 fm`, while `rdf_freud.py` treats `H` as
deuterium at `+6.671 fm`. Those have opposite signs, so S(q) and g(r) — Fourier transform pairs that should
describe the same sample — would silently disagree on every H-containing term. Run H-bearing systems with
`DSF_NEUTRON_WEIGHTING=no` to get the unweighted partials and total.

**Every other element agrees**, and `dsf.py` now checks this at runtime against `REFERENCE_B_COH` (a copy of
`rdf_freud.py`'s table), warning on any drift above 2% — so a future dynasor update cannot change the
weights underneath you unnoticed. Verified across the shared table: Al, C, Na, P, S, O, Si, N, Mg and Cl
agree to ≤0.1%; Ni (10.300 vs 10.332) and Zr (7.160 vs 7.119) differ by <1% because NIST tabulates one
value per element while dynasor sums over isotopes at natural abundance, which is why the tolerance is 2%
rather than exact.

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
| `vdos_dynmat.py` | `element` (symbol) — **required** | not read | Reads only frame 0 of `dump.lammpstrj`, for the element of each matrix row; the physics comes from `dynmat.dat` |

---

## `vdos_dynmat.py` — VDOS from the dynamical matrix

A second, independent route to the same quantity `vdos.py` produces. `vdos.py`
measures what the atoms actually did at the simulation temperature; this measures
the curvature of the potential energy surface at a single minimum.

| | `vdos.py` | `vdos_dynmat.py` |
|---|---|---|
| Source | velocity autocorrelation of `dynamics.lammpstrj` | eigenvalues of the force-constant matrix |
| Temperature | the MD temperature | 0 K (harmonic) |
| Linewidth | thermal + `1/CORR_LENGTH` resolution limit | discrete modes; width is whatever you set with `SMEARING` |
| Anharmonicity | included | excluded — peaks sit slightly higher |
| Cost | reading a long trajectory | 6N force evaluations + a (3N)² diagonalization |

Both write the same CSV column layout, so the two can be overlaid directly. Peak
positions should agree. Absolute heights will not quite, under `phonon`
normalization — the two reach "≈3 per atom" by different routes — so use
`unit_area` on both if you want the heights to line up too.

### Two stages

Enabled with a single pipeline flag, `--run-vdos-dynmat 1`, which drives both:

1. **LAMMPS.** Sets `RUN_DYNMAT=1`, which activates a block at the end of
   `OH-therm.input` / `b-SiO-therm.input`: minimize the configuration the MD run
   just finished with, then `dynamical_matrix all regular ${DYNMAT_DISPLACEMENT}
   file ${DYNMAT_FILE}`. Requires a LAMMPS build with the **PHONON** package —
   which is why `jobs/slurm/lammps_submit.slurm` now defaults to
   `lmp_mpi_phonon_2019` (override with `--lmp-bin`).
2. **Analysis.** `vdos_dynmat.py` reshapes the matrix to (3N, 3N), symmetrizes,
   diagonalizes, converts eigenvalues to frequencies, and histograms them.

### Why no reference dump is needed

`dynamical_matrix` iterates atoms by global ID, so matrix row `3i+α` is atom
`i+1`. `dump.lammpstrj` is written with `dump_modify sort id`, so frame 0 lists
atoms in that same order — and because the matrix is built on the final
configuration of the run that wrote that dump (same atoms, same IDs, no
`replicate` in between), the dump's `element` column indexes the matrix rows
directly. The script checks `9N²` against the dump's atom count and refuses to
run if they disagree, and checks that the IDs really are `1..N` in order.

### Units

The matrix is mass-weighted, so with `units metal`:

| `dynamical_matrix` style | eigenvalue λ | ν [THz] |
|---|---|---|
| `regular` (default) | eV/(Å²·amu) | `sqrt(λ) * 15.6333042` |
| `eskm` | 1/ps² | `sqrt(λ) / 2π` |

These agree: LAMMPS's own eskm factor is 9648.5, and `sqrt(9648.53)/2π = 15.6333`.
`VDOS_DYNMAT_MATRIX_STYLE` must match whichever style the `.input` file used —
it sets this conversion, so a mismatch rescales the whole spectrum.

### Variables

`DYNMAT_*` drive the LAMMPS stage (they reach `in.input` as `-var`);
`VDOS_DYNMAT_*` drive the Python. All are optional except one.

| Variable | Default | Meaning |
|---|---|---|
| `DYNMAT_MIN_STYLE` | `cg` | `min_style` for the pre-dynmat minimization |
| `DYNMAT_MIN_ETOL` | `1.0e-12` | `minimize` etol |
| `DYNMAT_MIN_FTOL` | `1.0e-12` | `minimize` ftol — the one that matters; residual forces become spurious imaginary modes |
| `DYNMAT_MIN_MAXITER` | `100000` | `minimize` maxiter |
| `DYNMAT_MIN_MAXEVAL` | `1000000` | `minimize` maxeval |
| `DYNMAT_DISPLACEMENT` | `0.0001` | finite-difference step, Å |
| `DYNMAT_FILE` | `dynmat.dat` | matrix filename; one value reaches both stages |
| `DYNMAT_BINARY` | `no` | `yes` = raw float64 (half the size, no text parse) |
| `VDOS_DYNMAT_MAX_FREQUENCY` | **required** | DOS grid upper limit, in `XUNIT` |
| `VDOS_DYNMAT_XUNIT` | `meV` | `meV` \| `THz` \| `cm-1` \| `eV`; sets the plot axis and the unit the two above are read in |
| `VDOS_DYNMAT_BINS` | `500` | frequency grid points |
| `VDOS_DYNMAT_SMEARING` | `0` | Gaussian FWHM in `XUNIT`; 0 = plain histogram |
| `VDOS_DYNMAT_MATRIX_STYLE` | `regular` | must match the `.input` file |
| `VDOS_DYNMAT_NORMALIZATION` | `phonon` | Sum rule: `phonon` (∫ = 3 per atom) or `unit_area` |
| `VDOS_DYNMAT_WEIGHTING` | `unity` | Species weighting, semicolon-separated: `unity`, `coherent`, `incoherent`, `total` |
| `VDOS_DYNMAT_PARTIAL` | `yes` | `no` skips per-element curves, uses `eigvalsh`, halves memory and runtime |
| `VDOS_DYNMAT_ASR` | `none` | `simple` imposes the acoustic sum rule |
| `VDOS_DYNMAT_THREADS` | unset | BLAS threads; the pipeline's `OMP_NUM_THREADS` already covers this |
| `VDOS_DYNMAT_CHARACTER` | `no` | `yes` adds the mode-character analysis below |
| `DYNMAT_REF_TRAJ` | `dynmat_ref.lammpstrj` | minimized coordinates; one value reaches both stages |
| `VDOS_DYNMAT_BRIDGE_ELEMENT` | `O` | the bridging atom |
| `VDOS_DYNMAT_NEIGHBOR_ELEMENT` | `Si` | its two neighbours |
| `VDOS_DYNMAT_BOND_CUTOFF` | `2.2` | bridge–neighbour max distance, Å |
| `VDOS_DYNMAT_OUTPUT` | `vdos_dynmat` | output basename |

**Weighting.** `VDOS_DYNMAT_WEIGHTING` works exactly like `vdos.py`'s `VDOS_WEIGHTING` and emits the same
`DoS(Total_<weighting>)` columns, with the same σ/m physics, the same Σ-normalization, the same refusal of
H/D-bearing systems, and the same omission of the Debye–Waller factor. Everything but `unity` needs
`VDOS_DYNMAT_PARTIAL=yes`, since without eigenvectors there are no per-element participations to weight.
The mass in σ/m comes from `ELEMENT_MASSES` — the masses the dynamical matrix was mass-weighted with — so
adding an element for weighting means adding it to both `ELEMENT_MASSES` and `NEUTRON_CROSS_SECTIONS`.

The two scripts print the same per-element *shares* for a given composition (SiO₂: `unity` Si 0.333 /
O 0.667, `coherent` and `total` Si 0.127 / O 0.873), which is the cross-check that they agree — but they
reach it differently, and the difference is easy to get wrong when reading the code. `vdos.py`'s partials
are per-element *shapes*, so its weight carries the concentration (`w_el = c_el·σ/m`). `vdos_dynmat.py`'s
partials already carry it — each integrates to 3·c_el — so its factor is `k_el = σ/m` normalized so
`Σ k_el·c_el = 1`. Multiplying by another `c_el` here would count concentration twice.

> **Output rename:** the `DoS(Total)` column is now `DoS(Total_unity)`, matching `vdos.py`. Values are
> unchanged; the reduced DOS `g/ν²` and the mode-character plot still use the unweighted total.

`MAX_FREQUENCY` is required, like `vdos.py`'s `VDOS_MAX_FREQUENCY_EV`, because a
default that is too low truncates the spectrum silently rather than failing.
`BINS` only changes smoothness, so it has one.

### Cost

The finite-difference loop is 6N force evaluations and the matrix is (3N)². The
cell is whatever the MD run used, so at `replicate 6 6 6` (N=5184) that is 31104
force evaluations and a 15552×15552 matrix: ~1.9 GB in memory, ~2.9 GB as text,
peaking near 3× that during the diagonalization. Budget the trajectory job's
`--time` for it and consider `DYNMAT_BINARY="yes"`.

Parallelization: the LAMMPS half is MPI parallel like any other force
computation, so it scales across the trajectory job's ranks (though LAMMPS
gathers 9N doubles on *every* rank, so per-rank memory does not fall). The Python
half is single-node and thread parallel — `np.linalg.eigh` is a threaded LAPACK
call driven by `OMP_NUM_THREADS`. There is no distributed diagonalization.

### Mode character — what kind of motion each band is

`VDOS_DYNMAT_CHARACTER=yes` adds the standard amorphous-silica band assignment,
following Bell & Dean and Taraskin & Elliott. At every **bridging** oxygen (one
with exactly two Si neighbours) the two bond directions r̂₁, r̂₂ define three
mutually orthogonal directions — two in the Si–O–Si plane, one normal to it:

| direction | definition | motion | band |
|---|---|---|---|
| **stretch** | `norm(r̂₁ − r̂₂)`, along Si···Si | one Si–O lengthens as the other shortens | ~1050–1200 cm⁻¹ (130–150 meV) |
| **bend** | `norm(r̂₁ + r̂₂)`, along the bisector | the Si–O–Si angle opens and closes | ~800 cm⁻¹ (~100 meV) |
| **rock** | `norm(r̂₁ × r̂₂)`, ⊥ to the plane | O moves out of the Si–O–Si plane | ~400–500 cm⁻¹ (50–60 meV) |

The two in-plane directions are perpendicular for free: `(r̂₁−r̂₂)·(r̂₁+r̂₂) =
|r̂₁|² − |r̂₂|² = 0` because both are unit vectors — the diagonals of a rhombus.
So the three form a *complete* basis and each oxygen's displacement splits
exactly, `|u|² = (u·ŝ)² + (u·b̂)² + (u·r̂)²`.

**Nothing is classified.** Every mode gets three fractions summing to 1, e.g.
`(0.62, 0.21, 0.17)` — no thresholds, no labels. The bands appear when the DOS is
weighted by those fractions, which is why `DoS(stretch)+DoS(bend)+DoS(rock)`
equals `DoS(Total)` exactly.

Three things to know about what the fractions mean:

- Only **bridging-oxygen** motion is in the denominator. Si motion and
  non-bridging-O motion contribute nothing, so a fraction reads "of the bridging-O
  motion in this mode, how much is stretch" — not "of the whole mode". The
  element-partial DOS covers the rest.
- A 0.5/0.5 mode is genuinely ambiguous between "half the oxygens rocking, half
  stretching" and "every oxygen at 45°". The per-mode table narrows this; only
  looking at the eigenvector settles it.
- Displacements are `u = e/√m`, **not** the eigenvectors. The eigenvectors of a
  mass-weighted matrix are not displacements, and using them directly would
  misweight oxygen against silicon throughout.

A **linear** bridge (180°) has no plane: `r̂₁+r̂₂` and `r̂₁×r̂₂` both vanish and
bend/rock become physically degenerate. Ideal β-cristobalite is exactly this, so
the script substitutes an arbitrary perpendicular pair and reports how many
bridges are within 5° of linear. Their *sum* (transverse motion) is still
meaningful; the split between bend and rock is not. Stretch is unaffected.

Alongside the decomposition you also get:

- **Participation ratio**, `PR(m) = 1/(N Σᵢ|eᵢ|⁴)`, from `1/N` (all motion on one
  atom) to `1` (every atom moving). In a glass this is half the physics — it is
  how propagons, diffusons and locons are separated, and how the boson-peak
  region is identified.
- **Reduced DOS** `g(ν)/ν²` as a CSV column and plot panel. Debye predicts
  `g ~ ν²`, so this is flat for a crystal and shows a peak in a glass — the boson
  peak.
- A **coordination census**: how many oxygens are 1-, 2-, 3-coordinated, and the
  mean Si–O–Si angle. In a quenched glass some oxygens are non-bridging, and the
  census says how much of the structure the decomposition actually covers. Worth
  reading as a glass-quality check in its own right.

Extra outputs when this is on:

- `<date>_vdos_dynmat_modes.csv` — one row per mode: frequency in all four units,
  participation ratio, `frac_stretch/bend/rock`, and per-element fractions. This
  is what makes "which modes are in this peak" answerable.
- `<date>_vdos_dynmat_character.png` — three panels: character-resolved DOS,
  participation ratio vs frequency, and the reduced DOS.

It needs `dynmat_ref.lammpstrj`, the coordinates LAMMPS writes immediately after
`minimize`. The last frame of `dump.lammpstrj` will not do: it predates the
minimization, and the relaxation rotates exactly the bond directions this
analysis projects onto.

### Reading the diagnostics

Three zero modes are the acoustic translations and are expected. The script snaps
modes within one bin width of zero to exactly zero (finite-difference noise
scatters them across zero, and without this some fall out of range and quietly
cost the spectrum weight), then reports anything still negative:

```
  Lowest 6 signed frequencies (meV): 0.0000, 0.0000, 0.0000, 18.5958, 32.1654, 43.4362
  Modes snapped to zero (|nu| < 1 meV, one bin width): 3  (3 acoustic translations are expected)
  Imaginary modes: 0
```

Imaginary modes beyond the acoustic ones mean the minimization did not reach a
true minimum. Tighten `DYNMAT_MIN_FTOL` or reduce `DYNMAT_DISPLACEMENT` before
trusting the spectrum. Two other checks worth doing once on a new system: the
reported `max|D - D.T| / max|D|` should be small, and the spectrum should be
stable across `DYNMAT_DISPLACEMENT` of 1e-5 to 1e-3 (too small is numerical
noise, too large samples anharmonicity).
