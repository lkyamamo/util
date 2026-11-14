# Working in util

Quick reference for AI assistants and for me. Read this before editing.

Research workflows for LAMMPS/VASP/QXMD molecular dynamics of silica and water on the USC CARC cluster
(endeavour): building systems, running them, tracking runs, analyzing and visualizing results, and moving
data to the NAS. Organized by workflow stage, not as one package. Most code runs on the cluster from
`/home1/lkyamamo/util`; the Mac checkout is for editing, plotting and visualization.

## Repo map

- `md_setup/` — build MD inputs: bubbles (`md_setup/create_bubble/`), ring and water permeation paths
  (`md_setup/ring_path/`, `md_setup/water_path/`), hand placement, mirroring/rotation, broken-bond repair.
  Shared I/O: `md_setup/lammps_data.py`, `md_setup/structure_io.py`. uv project.
- `starting-structures/` — reference crystal structures (`.data`), the one place `.data` files are tracked.
- `simulation/` — run folders and helpers per code: `simulation/lammps/` (dated run folders such as the
  dielectric runs, LAMMPS source mods, MSD input converters), `simulation/vasp/` (INCAR/POSCAR helpers).
- `jobs/` — running on the cluster.
  - `jobs/slurm/` — base SLURM/PBS submit scripts for LAMMPS and VASP.
  - `jobs/pipeline/` — multi-stage sbatch pipelines: `jobs/pipeline/lammps_to_analysis/` (run → distribution
    analysis), `jobs/pipeline/lammps_to_temperature_sweep/`, `jobs/pipeline/plotting/`. Each reads a
    gitignored `.conf` copied from its committed `.conf.example`.
  - `jobs/ledger/` — `simledger`, the read-only run ledger (SQLite + Markdown cards). Stdlib only, has tests.
    See `jobs/ledger/README.md`.
  - `jobs/carc/` — cluster convenience scripts.
- `analysis/` — post-processing, by what is measured: `analysis/distributions/` (RDF, BAD, S(q), VDOS, MSD;
  `analysis/distributions/20260608_GrNrBaSqw/` is the current Python set the pipeline copies),
  `analysis/dynamics/` (C++/Fortran MSD and correlation codes), `analysis/bonds/`, `analysis/structure/`,
  `analysis/energetics/`, `analysis/electronic/` (QXMD participation numbers), `analysis/shock/` (voxel
  analysis), `analysis/general/` (box size / density calculators), `analysis/kenichi_codes/` (inherited
  reference codes). uv project.
- `plotting/` — thermo and OSZICAR parsers (`plotting/lammps/thermodynamics.py` is called by `jobs/slurm/lammps_submit.slurm`).
- `visualization/` — OVITO renders, rotations, reaction-path GIFs, bounding-box subsets. uv project.
- `NAS/` — HPC → NAS transfer. `NAS/pipeline/nas.sh` (`sync` / `status` / `manifest` / `relink`) is current;
  `NAS/legacy/` is superseded.
- `grading/` — coursework grading automation, unrelated to the research code. uv project.
- `tools/` — the commit gate: `tools/repo_check.py`, its suites in `tools/verification.json`, `tools/setup.sh`.
  See `tools/README.md`.
- `plans/` — plans for work in progress, one file per scope, from `plans/_TEMPLATE.md`. Deleted when the work
  ships.

## Development workflow

Every code change goes through these steps, in this order. The gate enforces each one (`tools/README.md`).

1. **Plan.** Write `plans/<slug>.md` from `plans/_TEMPLATE.md`. Its `## Scope` lists every file the work touches;
   its `## Verification` is the verification plan: known answer, invariants, independent implementation,
   regression, tolerances and sources, slow checks, done when. "n/a" needs a reason.
2. **Approve, then commit the plan.** Only the owner sets `**State:** approved`. The plan is committed before
   any code it covers; Claude Code cannot write uncovered code at all (a local `PreToolUse` hook).
3. **Code with its verification.** Each commit carries `Plan: plans/<slug>.md`, changes the tests of any code it
   changes in `analysis/`, `jobs/ledger/` or `tools/` (or names one with `Verified-by: <test> — <reason>`), and
   updates the folder README's `## Verification`: every file, every test of it, method, what is checked,
   tolerance. Every analysis script is run by a test.
4. **Commit.** Pre-commit runs the checks and the touched suites on the staged tree, then records the passing run
   as the README's `### Last verified` block, replacing the previous one.
5. **Push.** Pre-push runs the verification run (every pushed commit, every fast suite) and refuses the push on
   failure; CI repeats it, and runs the slow checks nightly.
6. **Ship.** Delete the plan once its work is done.

Run `tools/setup.sh` once per checkout to turn the hooks on.

## Commit gate

`tools/repo_check.py` also refuses a change that:

- adds a generated file (anything `.gitignore` matches, force-added), a compiled binary, a file over 1 MB, or a
  backup copy (`bkup/`, `*.bak`, `*_old.py`, `* copy*`)
- adds a line hardcoding `/home1/lkyamamo`, `/scratch*/lkyamamo` or `/Users/loganyamamoto` in code, other
  than as an overridable default (`${VAR:-/path}`, argparse `default=`), in a comment, or on a line marked
  `allow-path`
- adds a file identical to one already tracked
- has a `.py` that does not compile or a shell/SLURM script `bash -n` rejects
- leaves this map out of date: a top-level directory missing, or a path named here or in any `.md` that the
  change deleted or moved
- imports a package `analysis/pyproject.toml` does not declare (a warning in other projects)

It warns on a source file sharing its name with one elsewhere. `python3 tools/repo_check.py --all` lists every
existing violation; `tools/tests/baselines/debt.json` pins that list, so it changes only deliberately.

## Conventions

- **One copy of each tool.** Never copy a script or Fortran/C++ source into a new analysis folder to change
  it. Parameterize the shared one (flags, a config file) and call it from the run folder. If a fork is truly
  needed, say in its header which file it forked and why.
- **No backups in the tree.** No `bkup/`, `_old`, `_v2` copies; git history keeps every version. Delete code
  instead of commenting it out.
- **Paths come from config.** Cluster and laptop paths go in a gitignored `.conf`/`.config` with a committed
  `.example`, or a flag with an overridable default. Scripts locate siblings relative to their own location
  (`Path(__file__).parent`, `$(dirname "$0")`).
- **Inputs you edit per use are not tracked.** If a config is rewritten for every system (box sizes,
  compositions), commit a `.example` and gitignore the real one, or take the values as arguments.
- **Generated output never goes in git.** Dumps, logs, CSVs, figures, `.h5`, compiled executables and Fortran
  `.mod` files are regenerated by whatever made them. Build from source with the folder's `Makefile`.
- **Python dependencies are declared.** A third-party import goes in that project's `pyproject.toml`
  (`uv add <pkg>`) and `uv.lock` is committed. Code that runs on the cluster names the venv it needs.
- **Tests hit real code.** Build a synthetic input whose answer is known by construction (`analysis/testkit/`)
  and run the real script on it. Analysis tests are pytest files named `test_*.py` beside the code; long
  convergence and sampling checks are marked `@pytest.mark.slow`. A new suite area goes in
  `tools/verification.json`.
- **Docs say what is true now.** Update this map and the folder's README in the same commit that adds, moves
  or deletes something. Replace text that stops being true; do not append dated logs. A doc carries only the last
  run that verified the current code (its `### Last verified` block); older runs live in git history, run results
  in the ledger, the story in the commit message. Add a new `.md` only for a new tool or workflow.
- **Commit messages** are one imperative sentence saying what changed (`Track SLURM jobs per run: ...`); body
  for the why. Work on a branch and merge.

## Environment

- Local Python: the pyenv 3.11.13 build is linked against a removed MacPorts `libsqlite3` and cannot import
  `sqlite3`. The hook picks the first `python3` that can; set `UTIL_PYTHON` to override.
- Analysis tests use the analysis project's `.venv` (override with `ANALYSIS_PYTHON`).
- On the cluster, hooks and cron call `simledger` with the spack Python 3.11.9 by full path.
