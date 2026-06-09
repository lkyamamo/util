# Structural Analysis Pipeline

Scripts for computing structural and dynamical properties from a LAMMPS MD trajectory.
All scripts read from the same dump file and produce independent outputs — they can be
run in any order or simultaneously.

---

## Scripts

| Script | What it computes | Output files |
|---|---|---|
| `rdf_freud.py` | Radial distribution function g(r) and coordination number n(r) for all element pairs | `rdfs.csv`, `rdfs.png`, `nrs.csv`, `nrs.png` |
| `bad_freud.py` | Bond angle distribution P(θ) for all A-B-C triplets | `bads.csv`, `bads.png` |
| `dsf.py` | Static structure factor S(q) and dynamic structure factor S(q,ω) | `sq.csv`, `sq.png`, `dsf.csv`, `dsf.png` |

---

## How to Run

### On the cluster (SLURM)

```bash
sbatch dsf_submit.slurm
```

The slurm script runs all three Python scripts in sequence under the same job.
Toggle which scripts run at the top of `dsf_submit.slurm`:

```bash
RUN_DSF=1   # 1 = run, 0 = skip
RUN_RDF=1
RUN_BAD=1
```

### Locally (all cores)

```bash
./run_dsf.sh          # all cores
./run_dsf.sh 8        # 8 threads
```

`run_dsf.sh` only runs `dsf.py`. To run the freud scripts locally, call them directly:

```bash
python rdf_freud.py
python bad_freud.py
```

---

## What to Modify Before Running

### All scripts — trajectory input

Every script reads from the same dump file. Set `DUMP_FILE` consistently in each:

| Script | Variable | Default |
|---|---|---|
| `rdf_freud.py` | `DUMP_FILE` | `"dump.lammpstrj"` |
| `bad_freud.py` | `DUMP_FILE` | `"dump.lammpstrj"` |
| `dsf.py` | `DUMP_FILE` | `"dump.lammpstrj"` |

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
| `DUMP_FILE` | Trajectory path | Requires `element` column in dump (not numeric type) |
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

### `dsf_submit.slurm` — required changes

| Item | Location | Notes |
|---|---|---|
| `RUN_DSF` / `RUN_RDF` / `RUN_BAD` | Lines 23–25 (run flags) | Set `1` to run, `0` to skip |
| `--time` | SBATCH header | Increase for large trajectories |
| `--cpus-per-task` | SBATCH header | Sets thread count; passed to `OMP_NUM_THREADS` |
| venv path | `source /home1/lkyamamo/venv/struc_analysis/bin/activate` | Update if environment moves |

---

### `run_dsf.sh` — required changes

No changes required for HPC use — activates the same venv as `dsf_submit.slurm`:
`/home1/lkyamamo/venv/struc_analysis/bin/activate`

---

## Environment

### HPC (cluster)

The shared venv on the HPC contains all dependencies:

```bash
source /home1/lkyamamo/venv/struc_analysis/bin/activate
```

This is what `dsf_submit.slurm` activates automatically.

### Local install

```bash
pip install dynasor matplotlib freud
pip install icc_rt    # optional: 5–10× speedup for dynasor's numba backend
```

---

## Dump Format Requirements

| Script | Element column | Coordinate columns | Notes |
|---|---|---|---|
| `rdf_freud.py` | `element` (symbol) or any string | `x y z` (real, Å) | Column indices set manually via `COL_*` |
| `bad_freud.py` | `element` (symbol) or any string | `x y z` (real, Å) | Column indices set manually via `COL_*` |
| `dsf.py` | `element` (symbol) — **required** | `x y z` or `xs ys zs` or `xu yu zu` | Column layout auto-detected from `ITEM: ATOMS` header |
