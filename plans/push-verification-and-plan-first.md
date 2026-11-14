# Development workflow: plan first, verification with code, verify on push

**State:** approved (2026-11-13)

## Goal

Make these rules hold for every change, enforced by git hooks and CI rather than habit:

1. **A plan exists before code.** A commit that changes code names an approved plan in `plans/`. The plan was
   committed before the code commit, is complete, and its Scope covers every code file the commit changes. Every
   plan has a verification plan.
2. **Analysis code ships with its verification.** Every Python script in `analysis/` is run by at least one test.
   A commit that changes analysis code also changes its tests. Verification uses all four methods: known answer,
   invariants, independent implementation and regression.
3. **Documentation includes the tests.** Each code folder's README has a `## Verification` section. It names every
   test that exercises each file there, with the method, what is checked and the tolerance. For code without
   automated tests, it records the manual checks done. Documentation of the current state carries **only the last
   run that verified it**: a `### Last verified` block, replaced on every run and never appended.
4. **Every push runs a verification run.** All fast suites and every gate check run over the pushed range, and a
   failure blocks the push. Slow checks run on demand and nightly in CI.

This is plan 1 of 3. Plan 2, `plans/consolidate-analysis-code.md`, moves analysis Python into `analysis/`. Plan 3
(written after plan 1 lands) backfills verification for `analysis/distributions/20260608_GrNrBaSqw/`.

## Scope

- `tools/repo_check.py`: the new checks (below), a `--commit SHA` mode (per-commit checks against that commit's
  parent) and a `--push` mode (pre-push and CI). The suite registry moves out to `tools/verification.json`.
- `tools/verification.json` (new): the suites. Each entry has a name, area prefixes, test globs, interpreter,
  fast command, slow command, and whether it is enforced (pairing and run-by-test).
- `tools/tests/test_repo_check.py` (new): verification for the gate itself (see Verification).
- `tools/README.md` (new): what the gate checks, how to run it, and its `## Verification` section.
- `.githooks/_python.sh` (new): the interpreter lookup the hooks share.
- `.githooks/pre-commit`: use `_python.sh`.
- `.githooks/commit-msg` (new): the per-commit plan and pairing checks. They need the message, which pre-commit
  does not have.
- `.githooks/pre-push` (new): the verification run.
- `.github/workflows/ci.yml`: the verification run on push and PR, plus a nightly `schedule` job for the slow
  checks. Python 3.11, matching the cluster.
- `analysis/pyproject.toml`: declare every third-party import (numpy, scipy, freud-analysis, h5py, ase, ovito,
  pandas, matplotlib, dynasor, pyvista, psutil), add pytest as a dev dependency, and configure pytest
  (`testpaths`, a `slow` marker, `-m "not slow"` by default).
- `analysis/uv.lock` (new).
- `analysis/conftest.py` (new): puts `analysis/testkit/` on the path.
- `analysis/testkit/` (new): shared test helpers.
  - Synthetic-system generators: ideal gas, a perfect lattice (from `starting-structures/`), a random walk with a
    set D, harmonic oscillators at a set ω, and a dipole series with a set variance.
  - Writers for LAMMPS dump, custom and data formats.
  - `run_script(path, args)` and `load_script(path)`, for scripts in folders whose names (e.g. `20260608_GrNrBaSqw`)
    are not importable.
- `analysis/bonds/test_bond_events.py` and `analysis/shock/etc/test_pipeline.py`: convert to pytest and move their
  synthesizers into `analysis/testkit/`. The shock test writes to `tmp_path` instead of `etc/test_output/`.
- `analysis/bonds/README.md`, `analysis/shock/README.md` (new): `## Verification` sections for the two existing
  suites.
- `jobs/ledger/README.md`: a `## Verification` section for `jobs/ledger/tests/test_simledger.py`.
- `plans/_TEMPLATE.md` (new): the required plan sections.
- `CLAUDE.md`: the workflow (plan → approve → code + tests + README → push), `plans/` and `tools/` in the repo
  map, and the new rules in the commit-gate section.
- `.gitignore`: re-include `analysis/testkit/data/`, for small committed reference inputs (each under the 1 MB
  limit).
- `.claude/settings.json` (local; `.claude/` is gitignored): a `PreToolUse` hook that blocks Claude's Edit/Write
  on a code file no approved plan's Scope covers.

It also covers the first gate already in the working tree (`tools/repo_check.py`, `tools/setup.sh`,
`.githooks/pre-commit`, CI, `CLAUDE.md`, `.gitignore`). That code was written before this rule existed and lands
with this plan.

## Approach

### Definitions

- **Code:** `.py .sh .cpp .c .h .f90 .f`. SLURM/PBS scripts and LAMMPS inputs get the syntax check only.
- **Analysis code:** a non-test `.py` under `analysis/`, excluding `analysis/testkit/`.
- **Enforced areas:** `analysis/`, `jobs/ledger/` and `tools/` plus `.githooks/`. Pairing, run-by-test and
  declared dependencies apply here. The plan and README rules apply to all code.
- **A test of file F:** a `test_*.py` that imports F's module or names F's filename (to run it). For
  `jobs/ledger/` and `tools/`, the suite's test globs.

### Plans (`plans/<slug>.md`, from `plans/_TEMPLATE.md`)

- A `**State:**` line: `draft`, `approved` or `in progress`. Only you set `approved`.
- Required sections, each non-empty: `## Goal`, `## Scope` (backticked paths or globs), `## Approach`,
  `## Verification`, `## Out of scope`.
- `## Verification` is the verification plan. Its required subsections:
  - `### Known answer`: synthetic inputs whose result is known by construction, and the expected values.
  - `### Invariants`: sum rules, limits, symmetries.
  - `### Independent implementation`: a cross-check against another code (e.g. Kenichi's Fortran vs freud).
  - `### Regression`: pinned outputs on a small real input, and what a change in them means.
  - `### Tolerances and sources`: each tolerance with its justification, and where each expected value comes from
    (analytic, literature, the other code).
  - `### Slow checks`: convergence and sampling studies, marked slow.
  - `### Done when`

  Each must hold content. "n/a" passes only with a reason (e.g. "n/a: no second implementation exists").
- A plan is deleted when its work ships. Git history keeps it. At most about five plans are open at once.

### New gate checks

Per commit (`commit-msg` locally, and per commit in pre-push and CI). They apply when the commit changes code. A
pure rename (git similarity 100%) does not count as changing code.

- **plan-trailer:** the message has a `Plan: plans/<slug>.md` trailer.
- **plan-before:** that plan is tracked in the commit's parent, so it was committed first.
- **plan-approved:** its State is `approved` or `in progress`.
- **plan-complete:** every required section and verification subsection is present and holds content.
- **plan-scope:** every changed code file matches a path or glob in the plan's Scope.
- **paired** (enforced areas): for each suite whose area the commit changes, the commit also changes that suite's
  tests. The escape is a `Verified-by: <test file> — <reason>` trailer, for a pure refactor or a fix the existing
  tests already catch. The named test must exist and be a test of a changed file.

On the staged or ranged change (pre-commit, pre-push, CI):

- **run-by-test** (analysis code): every changed or added analysis script has at least one test of it.
- **readme-verification** (all code): the README in each changed code file's folder has `## Verification`, naming
  the file and every test of it.
- **deps:** in `analysis/`, an undeclared third-party import is an error. Elsewhere it stays a warning.

- **last-verified** (all code): that README has exactly one `### Last verified` block under `## Verification`, and
  the commit updates it.
  - **Folders with a suite:** pre-commit writes the block itself after the suite passes, then stages the README.
    The run that gated the commit is the run recorded. The block holds:
    - the date
    - the parent commit the change sits on
    - each command run, with its pass count and runtime
  - **Slow tests:** they get their own line, rewritten only by `tools/repo_check.py --record --slow`.
  - **Folders without a suite:** you edit the block by hand. It says what was checked manually and when.
  - **History:** older runs are not kept in the doc. Git history has them.

Untested existing code is strict on first edit: the commit that first changes it also adds its test and its README
section. `python3 tools/repo_check.py --all` lists the rest as debt.

### Verification run (`.githooks/pre-push` → `tools/repo_check.py --push`)

- It reads the refs being pushed from stdin. For each range: the gate checks over the range, the per-commit checks
  on every non-merge commit, and **every** suite's fast command (not only the touched ones): `repo_check`,
  `simledger`, and `analysis` (`pytest` with slow tests excluded).
- Any failure blocks the push. It must stay under a minute; a test that would exceed that is marked slow.
- CI runs the same on push and PR, so a push with `--no-verify` is still caught. The nightly job runs
  `pytest -m slow` in `analysis/`.

### Claude Code hook

- A `PreToolUse` hook on Edit/Write/NotebookEdit runs `tools/repo_check.py --plan-covers <path>`.
- When the path is code and no approved plan's Scope covers it, the hook exits 2 with the reason, which blocks the
  write and tells Claude to write or update a plan first.
- Plans, docs and configs are always writable.

## Verification

### Known answer

`tools/tests/test_repo_check.py` (stdlib `unittest`) builds a temporary git repo with a fake `analysis/` tree and
drives the real hooks. Each case asserts refused or allowed and the error text:

- **Commit-time plan and pairing:**
  - code with no `Plan:` trailer: refused
  - a plan added in the same commit as the code: refused (not committed first)
  - a draft plan: refused
  - a plan missing `### Invariants`: refused
  - a plan whose `### Regression` is just "n/a": refused
  - a file outside the plan's Scope: refused
  - a complete approved plan with code and its test: allowed
- **Test pairing and run-by-test:**
  - analysis code with no test change: refused
  - the same with a valid `Verified-by:`: allowed
  - a `Verified-by:` naming a test that does not exercise the file: refused
  - a new analysis script no test runs: refused
- **Documentation and dependencies:**
  - a README with no `## Verification`, or one that omits the file: refused
  - code changed without its README's `### Last verified` block changing: refused
  - two `### Last verified` blocks: refused
  - a pre-commit on a suite folder rewrites the block in place, with the parent commit and pass counts, and stages
    it; a second run replaces it rather than appending
  - an undeclared import in `analysis/`: refused; elsewhere: a warning
- **Exemptions:**
  - a docs/plan/config-only commit with no trailer: allowed
  - a pure rename: allowed
- **Pre-push:**
  - a range with one bad commit: refused, naming it
  - a range where a suite fails: refused
- **Claude Code hook:** `--plan-covers` on a covered path: exit 0; on an uncovered path: exit 2.
- **Existing checks:** one refused and one allowed case each for artifacts, user paths, duplicates, syntax and docs
  (the throwaway scenario script from the last session, made permanent).

### Invariants

- A commit that changes no code is never blocked by the plan or pairing checks.
- Running a check twice gives the same result.
- `--all` on the repo after this lands reports no plan, pairing or README errors for files this plan touched.

### Independent implementation

- The same cases run two ways, through the real git hooks and by calling `repo_check.py` directly. Both must
  agree.
- CI on GitHub and pre-push on the Mac must give the same verdict on the same range. This is checked once by
  pushing a deliberately failing branch.

### Regression

- `--all` output on the repo is pinned in `tools/tests/baselines/debt.json`. A change in the debt list (files fixed
  or new findings) fails the test until the baseline is updated in the same commit. Debt therefore only moves
  deliberately.

### Tolerances and sources

- Pre-push runtime is under 60 s on this Mac; it is currently ~10 s (simledger 1.3 s, bonds 6 s, shock 2 s).
- The expected verdicts come from the rules in this plan.

### Slow checks

- n/a: every gate test is fast, since the temporary repos are tiny.

### Done when

- `tools/tests/test_repo_check.py`, the converted bonds and shock tests, and simledger all pass.
- `pytest` in `analysis/` collects both converted tests.
- The pinned debt baseline matches.
- In a scratch clone, the real hooks refuse and allow the listed commits.
- A real `git push` runs the verification run, and a deliberately failing push is refused.
- CI passes on GitHub, and the nightly job runs once.
- The Claude Code hook blocks an uncovered write in this session.

## Out of scope

- Moving analysis Python into `analysis/`: plan 2.
- Backfilling tests for `analysis/distributions/20260608_GrNrBaSqw/`: plan 3.
- Fixing the other existing debt: tracked binaries, hardcoded paths, duplicate Fortran/C++ sources.
- Physics validation beyond what each analysis plan's verification section specifies.
- GitHub branch protection requiring CI. It's a repo setting to turn on by hand.

## Decisions taken as defaults (change before approving if wrong)

- **`Verified-by:` escape allowed** for refactors and fixes the existing tests catch. It must name a real test of a
  changed file. Run-by-test still applies.
- **Claude Code hook on**, on this machine only.
- **Code set:** source files only. SLURM/PBS scripts and LAMMPS inputs get syntax checks, not tests.
