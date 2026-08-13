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
| `rdf_freud.py` | Radial distribution function g(r) and coordination number n(r) for all element pairs | `rdfs.csv`, `rdfs.png`, `nrs.csv`, `nrs.png` |
| `bad_freud.py` | Bond angle distribution P(θ) for all A-B-C triplets | `bads.csv`, `bads.png` |
| `dsf.py` | Static structure factor S(q) and dynamic structure factor S(q,ω) | `sq.csv`, `sq.png`, `dsf.csv`, `dsf.png` |
| `vdos.py` | Vibrational density of states (aligned with `analysis/dynamics/src/msd.cpp` by default) | `vdos.csv`, `vdos.png` |
| `msd.py` | Mean square displacement and self-diffusion coefficient (10⁻⁵ cm²/s) per element | `msd.csv`, `msd.png` |
| `vdos_dynmat.py` | Vibrational density of states from the LAMMPS dynamical matrix — harmonic, 0 K, no trajectory | `vdos_dynmat.csv`, `vdos_dynmat.png` |

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
| `NORMALIZATION` (`VDOS_NORMALIZATION`) | `'phonon'` (mole-fraction-weighted, matches msd.cpp) or `'unit_area'` | `'phonon'` |
| `N_FRAMES`, `STRIDE` (`VDOS_N_FRAMES`, `VDOS_STRIDE`) | Max frames to read / read every Nth frame | `0` (all), `1` |
| `VDOS_THREADS` | scipy FFT thread count (`fft_periodogram` only); `0` = all cores | `0` |

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
| `rdf_freud.py` | `element` (symbol) or any string | `x y z` (real, Å) | Column indices set manually via `COL_*` |
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
| `VDOS_DYNMAT_NORMALIZATION` | `phonon` | `phonon` (∫ = 3 per atom) or `unit_area` |
| `VDOS_DYNMAT_PARTIAL` | `yes` | `no` skips per-element curves, uses `eigvalsh`, halves memory and runtime |
| `VDOS_DYNMAT_ASR` | `none` | `simple` imposes the acoustic sum rule |
| `VDOS_DYNMAT_THREADS` | unset | BLAS threads; the pipeline's `OMP_NUM_THREADS` already covers this |
| `VDOS_DYNMAT_OUTPUT` | `vdos_dynmat` | output basename |

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
