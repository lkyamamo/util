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
| `rdf_freud.py` | Radial distribution function g(r) and coordination number n(r) for all element pairs | `rdfs.csv`, `rdfs.png`, `nrs.csv`, `nrs.png` |
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
| `DUMP_FILE` | Trajectory path | |
| `R_MAX` | Max r in Å for g(r) | Must be < half the shortest box dimension |
| `COL_ELEMENT` | Column index of element symbol in `ITEM: ATOMS` | 0-indexed; default layout: `id element x y z` → set to `1` |
| `COL_X`, `COL_Y`, `COL_Z` | Column indices of x, y, z coordinates | Default: `2, 3, 4` |

Optional:

| Variable | What it controls | Default |
|---|---|---|
| `BINS` | Number of r-bins | `200` |
| `OUTPUT_CSV` | g(r) CSV path; `None` to skip | `"rdfs.csv"` |
| `OUTPUT_NR_CSV` | n(r) CSV path; `None` to skip | `"nrs.csv"` |

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
| `TIME_UNIT` | Time between consecutive dumped frames in **femtoseconds** | Falls back to `DT` if unset — same quantity dsf.py calls `DT` |

Optional (see the module docstring for the full METHOD/WINDOW/NORMALIZATION
explanation):

| Variable | What it controls | Default |
|---|---|---|
| `METHOD` | `'vacf_cosine_transform'` (matches `analysis/dynamics/src/msd.cpp`) or `'fft_periodogram'` (faster, this script's own approach) | `'vacf_cosine_transform'` |
| `CORR_LENGTH` | VACF max lag / Welch-segment length, in fs | 75% of total trajectory duration |
| `CORR_INTERVAL` | Spacing between VACF reference frames / segment starts, in fs | 10% of `CORR_LENGTH` |
| `MAX_FREQUENCY_EV` | Upper frequency limit of the output grid, in eV (`vacf_cosine_transform` only) | `0.1` |
| `NUM_GRIDS` | Frequency grid points (`vacf_cosine_transform` only) | `5000` |
| `WINDOW` | `'cosine_lag'`/`'none'` under `vacf_cosine_transform`; `'hann'`/`'none'` under `fft_periodogram` | Matches `METHOD` |
| `NORMALIZATION` | `'phonon'` (mole-fraction-weighted, matches msd.cpp) or `'unit_area'` | `'phonon'` |
| `N_FRAMES`, `STRIDE` | Max frames to read / read every Nth frame | `0` (all), `1` |
| `VDOS_THREADS` | scipy FFT thread count (`fft_periodogram` only); `0` = all cores | `0` |

---

### `msd.py` — required changes

| Variable | What it controls | Notes |
|---|---|---|
| `DUMP_FILE` | Dynamics trajectory path (env var `DYNAMICS_TRAJ`, default `dynamics.lammpstrj`) | Requires `element x y z` columns (wrapped, not `xu yu zu`; velocities not read) |
| `TIME_UNIT` | Time between consecutive dumped frames in **femtoseconds** | Falls back to `DT` if unset — same quantity dsf.py/vdos.py use |

Optional (see the module docstring for the full METHOD explanation — no
unwrapped coordinates needed, minimum-image correction on each
reference-to-current displacement instead):

| Variable | What it controls | Default |
|---|---|---|
| `CORR_LENGTH` | Max time lag / reference-frame spacing window, in fs — same meaning as vdos.py's `CORR_LENGTH` | 75% of total trajectory duration |
| `CORR_INTERVAL` | Spacing between reference frames, in fs — same meaning as vdos.py's `CORR_INTERVAL` | 10% of `CORR_LENGTH` |
| `MSD_FIT_FRACTION` | Fraction of the tail of the correlation window used for the diffusion-coefficient linear fit (excludes the early ballistic regime) | `0.5` |
| `N_FRAMES`, `STRIDE` | Max frames to read / read every Nth frame | `0` (all), `1` |

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
| `rdf_freud.py` | `element` (symbol) or any string | `x y z` (real, Å) | Column indices set manually via `COL_*` |
| `bad_freud.py` | `element` (symbol) or any string | `x y z` (real, Å) | Column indices set manually via `COL_*` |
| `dsf.py` | `element` (symbol) — **required** | `x y z` or `xs ys zs` or `xu yu zu` | Reads `dynamics.lammpstrj`; column layout auto-detected from `ITEM: ATOMS` header |
| `vdos.py` | `element` (symbol) — **required** | `vx vy vz` (velocities) — **required**; positions not read | Reads `dynamics.lammpstrj`; column layout auto-detected from `ITEM: ATOMS` header |
| `msd.py` | `element` (symbol) — **required** | `x y z` (wrapped, real, Å) — **required**; velocities not read | Reads `dynamics.lammpstrj`; column layout auto-detected from `ITEM: ATOMS` header |
