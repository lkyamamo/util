# tools: the commit gate

`repo_check.py` enforces the development workflow in `CLAUDE.md`. The git
hooks in `.githooks/` run it locally; `.github/workflows/ci.yml` runs it on
GitHub. `setup.sh` turns the hooks on for a checkout (once per clone).

| File | What it is |
|---|---|
| `repo_check.py` | All the checks, and the per-commit, range, push, debt and record modes (its docstring lists them) |
| `verification.json` | The verification suites: area, test globs, interpreter, fast and slow commands, the README for each area |
| `setup.sh` | `git config core.hooksPath .githooks` |
| `../.githooks/_python.sh` | Finds a `python3` that can import `sqlite3` (`UTIL_PYTHON` overrides) |
| `../.githooks/pre-commit` | Staged change: checks, the touched suites run on the staged tree, then each affected README's `### Last verified` block rewritten and staged |
| `../.githooks/commit-msg` | `Plan:` and `Verified-by:` trailers |
| `../.githooks/pre-push` | The verification run: every pushed commit, the range, every fast suite on the pushed tree |
| `tests/test_repo_check.py` | This gate's verification (below) |
| `tests/baselines/debt.json` | The existing debt, pinned (`--all --json`) |

## What a code commit needs

1. A plan in `plans/` (from `plans/_TEMPLATE.md`), committed **before** the code, `**State:** approved`, with every
   section and verification subsection filled in, and the changed files in its `## Scope`.
2. `Plan: plans/<slug>.md` in the commit message.
3. For code in a verified area (`analysis/`, `jobs/ledger/`, `tools/`): its tests changed in the same commit, or
   `Verified-by: <test file> — <reason>` naming a test that exercises it. Every analysis script must be run by a test.
4. Its README's `## Verification` naming the file and every test of it, with one `### Last verified` block holding
   only the latest run. Pre-commit writes that block for suite areas; elsewhere you write it by hand.

Docs, plans, configs, SLURM scripts, deletions and pure renames need none of this.

## Commands

```bash
tools/setup.sh                                   # once per clone
python3 tools/repo_check.py --all                # the existing debt
python3 tools/repo_check.py --record --slow      # run slow suites, record them in the READMEs
python3 tools/repo_check.py --base origin/master --verify   # what CI runs
```

`git commit --no-verify` and `git push --no-verify` skip the hooks; CI still runs the verification run on the push.

## Verification

`tests/test_repo_check.py` tests `repo_check.py`, `verification.json` and the hooks (`_python.sh`, `pre-commit`,
`commit-msg`, `pre-push`, `setup.sh`'s `core.hooksPath`) by building small git repositories with the real hooks
installed and a fake analysis suite, then making the commits and pushes each rule exists for.

| Tests | Method | Checked |
|---|---|---|
| `test_code_without_plan_trailer_is_refused`, `test_plan_in_same_commit_as_code_is_refused`, `test_draft_plan_is_refused`, `test_plan_missing_a_verification_subsection_is_refused`, `test_bare_na_is_refused`, `test_file_outside_plan_scope_is_refused` | known answer | Each plan rule refuses the commit built to break it, naming the rule |
| `test_complete_plan_with_code_and_test_is_allowed_and_recorded`, `test_second_run_replaces_the_block` | known answer, invariant | A valid commit passes; its README block names the parent commit and pass count; a second commit replaces the block rather than appending |
| `test_analysis_change_without_test_change_is_refused`, `test_verified_by_an_exercising_test_is_allowed`, `test_verified_by_a_test_that_does_not_exercise_the_file_is_refused`, `test_new_script_no_test_runs_is_refused` | known answer | Pairing, the `Verified-by:` escape, run-by-test |
| `test_readme_without_verification_section_is_refused`, `test_readme_that_omits_the_file_is_refused`, `test_manual_code_without_block_update_is_refused`, `test_manual_code_with_block_update_is_allowed`, `test_two_last_verified_blocks_are_refused` | known answer | README verification and the single last-run block |
| `test_undeclared_import_is_an_error_in_analysis`, `test_undeclared_import_is_a_warning_elsewhere` | known answer | Dependency strictness by area |
| `test_docs_and_plan_only_commit_needs_no_plan`, `test_pure_rename_needs_no_plan` | invariant | A commit that changes no code is never blocked by the plan rules |
| `test_push_with_a_bad_commit_is_refused`, `test_push_with_a_failing_suite_is_refused`, `test_push_of_good_commits_is_allowed` | known answer | The pre-push verification run |
| `test_ci_mode_agrees_with_pre_push` | independent implementation | CI's `--base --verify` and the pre-push hook reach the same verdict on the same range |
| `test_plan_covers`, `test_draft_plan_does_not_cover` | known answer | The Claude Code write-time check (`--plan-covers`): exit 2 for uncovered code, 0 otherwise |
| `test_force_added_ignored_file_is_refused`, `test_hardcoded_user_path_is_refused`, `test_overridable_user_path_is_allowed`, `test_identical_file_is_refused`, `test_syntax_error_is_refused`, `test_unmapped_top_level_directory_is_refused` | known answer | The original checks: artifacts, user paths, duplicates, syntax, the CLAUDE.md map |
| `test_all_mode_is_deterministic` | invariant | `--all` gives the same findings twice |
| `test_debt_matches_baseline` | regression | `--all` on this repo equals `tests/baselines/debt.json`; debt moves only in a commit that regenerates it |

Mutation check (done once, by hand): disabling the pairing rule fails the two pairing tests, and disabling the
approval rule fails `test_draft_plan_is_refused`.

### Last verified

<!-- Rewritten by tools/repo_check.py on every run that verifies this folder; never appended. -->
- Fast: 2026-11-13 on top of `00e7f05`: `python -m unittest discover -s tools/tests -q` in `./`: 35 tests OK (39.0 s)
- Slow: none defined
