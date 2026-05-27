# util repo summary

This repository is a personal utilities + scientific workflow workspace (LAMMPS/VASP/QXMD analysis, plotting, visualization, HPC job scripts, and grading automation). Most code is organized by workflow area rather than as a single installable package.

## Top-level directories

- **`calculations/`**: Science/MD/DFT analysis and helper codes.
  - **`calculations/codes/`**: C++/Fortran analysis executables + build scripts.
    - **`msd.cpp`**: Reads concatenated `all.xyz` (custom XYZ w/ lattice + velocities) and writes `msd.dat` and `dos.dat`.
    - **`corr.cpp`**: Correlation-function tool (paired with `msd.cpp`).
    - **`Makefile`**: `make msd` / `make corr` builds `qmd_*.exe` (and `nnqmd_*.exe` variants).
    - **`1.robust_lammps_to_input_xyz.py`, `2.xyz_split.sh`**: Convert trajectories into the formats expected by the analysis executables.
    - **`stat/`**: Fortran/statistics utilities (`stat.f90`, `voxels.f90`, etc.) with its own `Makefile`.
  - **`calculations/lammps/`**: LAMMPS-focused workflows (conversion, dielectric, MSD prep, OVITO helpers, source mods).
    - **`msd_code/`**: Converters that generate the `all.xyz` format consumed by `calculations/codes/msd.cpp`.
    - **`new-dielectric-constant/`**: MPI/batch scripts and utilities for dipole/dielectric post-processing and dump-format conversions.
    - **`lammps_src_mod/`**: Local LAMMPS source modifications / prototype force code.
  - **`calculations/md_analysis/`**: Python notebooks/scripts and sample dumps/data for trajectory/bond/structure analysis.
  - **`calculations/VASP/`**: Small VASP helper scripts (INCAR/structure utilities, parameter helpers).
  - **`calculations/qxmd_utils/`**: Fortran/Python utilities for QXMD-related workflows (e.g., participation numbers).
  - **`calculations/general_calc/`**: General-purpose calculation scripts and small datasets (`box_size.py`, parameter tables).
  - **`calculations/kenichi_codes/`**: Additional/legacy analysis code variants (often mirroring tools in `calculations/codes/`).

- **`lammps/`**: General LAMMPS utilities.
  - **`input_description.py`**: Parses a LAMMPS input script and writes a human-readable description of phases/commands.
  - **`lammps_archive.sh`**: Helper script for archiving LAMMPS-related artifacts.

- **`plotting/`**: Plotting utilities.
  - **`plotting/lammps/`**: Scripts to parse/combine thermodynamics CSV/log outputs.
  - **`plotting/vasp/`**: VASP output parsing helpers.

- **`visualization/`**: OVITO-based visualization helpers + batch scripts.
  - `ovito_visual_pipeline.py`, `ovito_render_rotation.py`, configs, and Slurm scripts for rendering (`pyproject.toml` present).

- **`md_setup/`**: Setup utilities for building MD inputs (e.g., bubble/geometry generation).
  - `bubble_setup.py`, `create_bubble.in`, plus a small Python project (`pyproject.toml`, `uv.lock`).

- **`NAS/`**: Scripts and docs for syncing/transferring data from HPC → local NAS with manifests/rsync/tarballs.
  - Entry points: `create_manifest.sh`, `upload_directory.sh`, `rsync_list.sh`, `process_large_directories.sh`.
  - Docs: `README.md`, `TUTORIAL.md`, `FILE_STRUCTURE.md`.

- **`grading/`**: Coursework/grading automation (Python project with `pyproject.toml` / `uv.lock`).
  - **`grading/ai/`**: AI-assisted grading / report tooling (CSV/JSON artifacts and scripts).
  - **`grading/instructor_in_loop/`**: Human-in-the-loop checklist/evidence extraction helpers.

- **`submission_files/`**: HPC scheduler submission scripts (`*.slurm`, `*.pbs`) for LAMMPS/VASP/etc.
  - Note: additional workflow-specific Slurm scripts also live under `calculations/kenichi_codes/` and some `calculations/lammps/*` subtrees.

- **`carc/`**: Small cluster convenience scripts (e.g., `sinfo2_command.sh`).

## Common workflows (high level)

- **Generate MSD inputs from LAMMPS**:
  - LAMMPS dump → `calculations/lammps/msd_code/*` converter → `all.xyz` → compile/run `calculations/codes/msd.cpp`.

- **Build analysis executables**:
  - `cd calculations/codes && make msd && make corr`

- **Transfer data from HPC to NAS**:
  - Use scripts in `NAS/` (manifest + rsync/tarball workflow; see `NAS/TUTORIAL.md`).

## Notes / conventions

- Many scripts assume they are run from a particular working directory (paths are often relative).
- Several folders contain “workflow snapshots” (Slurm scripts, example dumps, notebooks) rather than reusable libraries.
- Notebooks (e.g. under `calculations/md_analysis/` and `calculations/lammps/msd_code/`) often serve as the most up-to-date “how-to” for a workflow.

