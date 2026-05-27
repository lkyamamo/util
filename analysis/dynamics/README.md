# Dynamics / MSD workflow

## Directory layout

**Top level** (scripts and docs you edit/run):

| File | Purpose |
|------|---------|
| `msd_submit.slurm` | Batch MSD job on the cluster |
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

Sanity check after layout changes:

```bash
cd analysis/dynamics
bash verify_setup.sh
```

---

## `msd_submit.slurm`

Batch Slurm job that runs MSD analysis for multiple temperatures in one submission.

1. Reads settings from the **edit block** at the top (paths, temperature list, MSD parameters).
2. Resolves `INPUT_DIR` and `DYNAMICS_DIR` to absolute paths (relative paths are taken from the submit directory).
3. Builds `src/qmd_msd.exe` via `make msd` in `src/`.
4. Runs each temperature **in parallel** (up to `MAX_PARALLEL`, default `SLURM_NTASKS`):
   - Reads `{T}C.custom` from `INPUT_DIR`
   - Converts it to MSD input with `lammps_custom_to_input_xyz.py`
   - Runs `qmd_msd.exe` in a dedicated output folder `msd_{T}C/` under the submit directory
   - Writes `msd.dat`, `dos.dat`, `{T}C_calculation.xyz`, and `run.log` there

Submit from any directory; edit `INPUT_DIR`, `DYNAMICS_DIR`, and the other variables before `sbatch`.

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

`msd_submit.slurm` runs each temperature in its own subdirectory so these fixed filenames do not collide.

### Example

```bash
./src/qmd_msd.exe calculation.xyz 5000 10 2.5 0.5
```

Correlation length 5000 fs, new reference every 10 fs, 2.5 fs between frames, DOS computed up to 0.5 eV.
