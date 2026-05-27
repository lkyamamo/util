# Dynamics / MSD workflow

## Directory layout

**Top level** (scripts and docs you edit/run):

| File | Purpose |
|------|---------|
| `msd_submit.slurm` | Slurm **job array** — one MSD run per temperature |
| `lammps_custom_to_input_xyz.py` | LAMMPS custom dump → MSD input |
| `1.robust_lammps_to_input_xyz.py` | Legacy dump converter (hardcoded box) |
| `2.xyz_split.sh` | Split `calculation.xyz` for `qmd_corr.exe` |
| `run.sh` | Local full pipeline (convert → split → corr + MSD) |
| `README.md` | This file |

**`src/`** (C++ sources and build; not needed for routine use):

- `msd.cpp`, `corr.cpp`, `histogram.cpp`, `nbr.cpp`, `Makefile`
- `bkup/`, `stat/` — older sources and Fortran utilities

Build executables: `make -C src msd` (outputs `src/qmd_msd.exe`).  
`make -C src corr` requires the Boost module (`BOOST_ROOT`) on the cluster.

---

## `msd_submit.slurm`

Submits **one Slurm job array** so each temperature runs as its own parallel task. This is the reliable way to run several `qmd_msd.exe` jobs at once on the cluster (each array task is a separate allocation step).

### Slurm resources

| Directive | Value | Meaning |
|-----------|-------|---------|
| `--array=0-7` | 8 tasks | One task per temperature (must match `TEMPERATURES` length) |
| `--ntasks=8` | 8 | Tasks per array job (see cluster docs; MSD itself uses one core) |
| `--output=msd_%a.log` | per task | Log file per array index (`%a` = array task id) |
| `--exclusive` | node | Full node per array task |
| `--time=24:00:00` | 24 h | Wall time per array task |

### What each array task does

1. Picks a temperature: `T="${TEMPERATURES[$SLURM_ARRAY_TASK_ID]}"` (e.g. index `3` → `60` → input `60C.custom`).
2. Reads **`${INPUT_DIR}/${T}C.custom`** (default `../input_files/60C.custom` relative to the submit directory).
3. Builds **`src/qmd_msd.exe`** (`make msd` in `DYNAMICS_DIR/src`).
4. Converts the dump → **`${RUN_DIR}/msd_${T}C/${T}C_calculation.xyz`** via `lammps_custom_to_input_xyz.py`.
5. Runs **`qmd_msd.exe`** in **`msd_${T}C/`** (writes `msd.dat` and `dos.dat` there).

Array tasks do not share output directories, so nothing is overwritten between temperatures.

### What to edit before `sbatch`

At the top of `msd_submit.slurm`:

```bash
TEMPERATURES=(15 30 45 60 75 90 120 150)   # must match --array=0-7 (N temps → --array=0-(N-1))
INPUT_DIR="../input_files"
DYNAMICS_DIR="/path/to/util/analysis/dynamics"
```

Also set MSD arguments passed to `qmd_msd.exe` (if not set, the script may pass empty values — add these lines before the `qmd_msd.exe` call):

```bash
CORR_LENGTH=5000
CORR_INTERVAL=10
TIME_UNIT=2.5
MAX_FREQ=0.5
```

Keep **`#SBATCH --array=0-7`** in sync with **`TEMPERATURES`**: 8 temperatures → `--array=0-7`; 5 temperatures → `--array=0-4`, etc.

### Submit

```bash
cd /path/to/your/run_directory          # where outputs (msd_*C/) should go
sbatch /path/to/util/analysis/dynamics/msd_submit.slurm
```

- **`SLURM_SUBMIT_DIR`** becomes `RUN_DIR` (output folders `msd_15C/`, `msd_30C/`, … are created there).
- Put **`{T}C.custom`** files in `INPUT_DIR` (or adjust `INPUT_DIR` to point at your dumps).
- Set **`DYNAMICS_DIR`** to the `analysis/dynamics` folder in this repo (absolute path recommended).

### Outputs (per temperature)

```
msd_30C/
  30C_calculation.xyz
  msd.dat
  dos.dat
msd_30.log          # from msd_%a.log (array index, not temperature)
```

Check **`msd_<array_index>.log`** in the submit directory for Slurm stdout/stderr for each task.

---

## Python scripts

### `lammps_custom_to_input_xyz.py`

Converts a **standard LAMMPS custom dump** (with `ITEM: TIMESTEP`, box bounds, and `element x y z vx vy vz`) into the concatenated XYZ format expected by `qmd_msd.exe`.

- Parses each frame’s header and reads lattice parameters from `ITEM: BOX BOUNDS` (supports triclinic tilt).
- Writes one output file with frames back-to-back: atom count, `a b c α β γ`, then `element x y z vx vy vz` per atom.
- Optional `type_map` if the dump uses numeric `type` instead of `element`.
- Config at top of file: `input_filename`, `output_filename`, `start_frame`, `end_frame`, `type_map`.

Used by `msd_submit.slurm` (paths passed at runtime).

### `1.robust_lammps_to_input_xyz.py`

Older converter for a **non-standard concatenated dump** with a fixed 9-line header per frame (not full `ITEM:` blocks). Uses **hardcoded** lattice constants (`a`, `b`, `c`, angles) instead of reading the box from the file. Maps the first character of each atom line through `type_map` (e.g. `1` → `O`).

Used by `run.sh` when input is the legacy `*.custom` layout. Prefer `lammps_custom_to_input_xyz.py` for standard LAMMPS custom dumps.

Fortran/plot utilities live under `src/stat/` (not part of the MSD Slurm pipeline).

---

## `qmd_msd.exe` usage

Build: `make -C src msd` → produces `src/qmd_msd.exe` (and `src/nnqmd_msd.exe` with a different atom-line format).

### Command

```bash
./src/qmd_msd.exe <trajectory.xyz> [corr_length] [corr_interval] [time_unit] [max_frequency_ev]
```

### Arguments

| Arg | Name | Units | Description |
|-----|------|-------|-------------|
| 1 | trajectory file | — | Single concatenated XYZ file (`calculation.xyz` / `{T}C_calculation.xyz`). All frames in one file. |
| 2 | `corr_length` | fs | Maximum time lag for MSD/VAC correlation. |
| 3 | `corr_interval` | fs | Spacing between new reference frames (multiple time origins). |
| 4 | `time_unit` | fs/frame | Elapsed time between consecutive MD frames in the input file. |
| 5 | `max_frequency_ev` | eV | Upper frequency limit for the DOS grid written to `dos.dat`. |

If only the trajectory file is given, defaults are derived from the number of frames (e.g. `corr_length` ≈ 75% of total trajectory time, `time_unit` = 1 fs).

### Input trajectory format (per frame)

```
<num_atoms>
<a> <b> <c> <alpha> <beta> <gamma>
<element> <x> <y> <z> <vx> <vy> <vz>
...
```

### Outputs (written to current working directory)

| File | Contents |
|------|----------|
| `msd.dat` | Time [fs], MSD and VAC per element, sample counts, truncated VAC, current correlation `Jt` |
| `dos.dat` | Frequency [meV], DOS per element, total DOS, `JDoS` |

`msd_submit.slurm` runs each temperature in its own `msd_{T}C/` subdirectory so these fixed filenames do not collide.

### Example

```bash
./src/qmd_msd.exe calculation.xyz 5000 10 2.5 0.5
```

Correlation length 5000 fs, new reference every 10 fs, 2.5 fs between frames, DOS computed up to 0.5 eV.
