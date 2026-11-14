"""Verification for the commit gate (tools/repo_check.py and .githooks/).

Each test builds a small git repository with the real hooks installed and a
fake analysis suite, makes a commit (or push) the rules should refuse or allow,
and checks the verdict and the error. The hooks are driven for real: these are
the commits a person or Claude would make.

Run from the repo root:  python3 -m unittest discover -s tools/tests -v
"""

import json
import os
import shutil
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent
TREE = HERE.parent.parent  # the tree under test: a checkout, or the gate's export of one
GATE = TREE / "tools" / "repo_check.py"
HOOKS = TREE / ".githooks"
REAL_REPO = Path(os.environ.get("REPO_CHECK_ROOT") or TREE)
BASELINE = HERE / "baselines" / "debt.json"


def clean_env(**extra):
    """No GIT_* variables: inside a hook they point at the real repo's index."""
    env = {k: v for k, v in os.environ.items()
           if not k.startswith("GIT_") and k not in ("REPO_CHECK_ROOT", "REPO_CHECK_REF", "ANALYSIS_PYTHON")}
    env.update(UTIL_PYTHON=sys.executable, GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@t",
               GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@t", **extra)
    return env


FAKE_SUITES = {
    "suites": [{
        "name": "fake",
        "areas": ["analysis/"],
        "tests": ["analysis/*test_*.py"],
        "python": "stdlib",
        "cwd": "analysis",
        "fast": ["-c", "import glob, runpy\n"
                       "for p in sorted(glob.glob('**/test_*.py', recursive=True)):\n"
                       "    runpy.run_path(p, run_name='__main__')\n"
                       "print('3 passed')"],
        "slow": None,
        "per_file": True,
        "doc": None,
    }]
}

PLAN = """# {title}

**State:** {state}

## Goal

Test plan.

## Scope

{scope}

## Approach

Build it.

## Verification

### Known answer

double(2) == 4.

### Invariants

{invariants}

### Independent implementation

n/a: arithmetic needs no second implementation.

### Regression

{regression}

### Tolerances and sources

Exact integers.

### Slow checks

n/a: everything is fast.

### Done when

Tests pass.

## Out of scope

Nothing else.
"""


def plan(title="Work", state="approved", scope="`analysis/foo/`, `analysis/bar/`, `md_setup/`",
         invariants="double is linear.", regression="n/a: no real input exists for this."):
    text = PLAN.format(title=title, state=state, scope=scope, invariants=invariants, regression=regression)
    if invariants is None:
        text = text.replace("### Invariants\n\nNone\n\n", "")
    return text


CALC = "def double(x):\n    return 2 * x\n"
TEST_CALC = ("import os, sys\nsys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))\n"
             "from calc import double\nassert double(2) == 4\n")
FOO_README = """# foo

## Verification

`test_calc.py` checks `calc.py` (known answer: double(2) == 4).

### Last verified

<!-- Rewritten by tools/repo_check.py on every run that verifies this folder; never appended. -->
- Fast: not run yet
- Slow: none defined
"""
MD_README = """# md_setup

## Verification

`build.py`: checked by hand.

### Last verified

- 2026-01-01: ran build.py, output looked right.
"""
CLAUDE = "# map\n\n- `analysis/`\n- `tools/`\n- `plans/`\n- `md_setup/`\n"


class GateTest(unittest.TestCase):
    """Each test gets a copy of one template repository with the gate installed."""

    @classmethod
    def setUpClass(cls):
        cls.tmp = Path(tempfile.mkdtemp(prefix="gate_test_"))
        template = cls.tmp / "template"
        template.mkdir()
        files = {
            "CLAUDE.md": CLAUDE,
            ".gitignore": "*.log\n",
            "tools/verification.json": json.dumps(FAKE_SUITES, indent=1),
            "analysis/pyproject.toml": "[project]\nname = \"a\"\ndependencies = []\n",
            "analysis/foo/calc.py": CALC,
            "analysis/foo/test_calc.py": TEST_CALC,
            "analysis/foo/README.md": FOO_README,
            "md_setup/pyproject.toml": "[project]\nname = \"m\"\ndependencies = []\n",
            "md_setup/build.py": "print('built')\n",
            "md_setup/README.md": MD_README,
            "plans/work.md": plan(),
            "plans/draft.md": plan("Draft", state="draft"),
            "plans/incomplete.md": plan("Incomplete", invariants=None),
            "plans/na.md": plan("NA", regression="n/a"),
        }
        for rel, text in files.items():
            p = template / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(text)
        # The incomplete plan drops its Invariants section entirely.
        inc = template / "plans/incomplete.md"
        inc.write_text(inc.read_text().replace("### Invariants\n\nNone\n\n", ""))
        shutil.copy(GATE, template / "tools/repo_check.py")
        shutil.copytree(HOOKS, template / ".githooks")
        run(["git", "init", "-q"], template)
        run(["git", "config", "core.hooksPath", ".githooks"], template)
        run(["git", "add", "-A"], template)
        run(["git", "commit", "-q", "--no-verify", "-m", "initial"], template)
        cls.template = template

    @classmethod
    def tearDownClass(cls):
        shutil.rmtree(cls.tmp, ignore_errors=True)

    def setUp(self):
        self.root = Path(tempfile.mkdtemp(prefix="repo_", dir=self.tmp))
        shutil.rmtree(self.root)
        shutil.copytree(self.template, self.root, symlinks=True)

    # helpers ---------------------------------------------------------------

    def write(self, rel, text):
        p = self.root / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text)

    def git(self, *args, check=True):
        return run(["git", *args], self.root, check=check)

    def commit(self, message, files=None, verify=True, add_all=True):
        for rel, text in (files or {}).items():
            self.write(rel, text)
        if add_all:
            self.git("add", "-A")
        msg = self.root.parent / f"msg_{self.root.name}.txt"
        msg.write_text(message)
        args = ["commit", "-F", str(msg)] + ([] if verify else ["--no-verify"])
        return self.git(*args, check=False)

    def gate(self, *args, stdin=None):
        return run([sys.executable, str(self.root / "tools/repo_check.py"), *args], self.root,
                   check=False, stdin=stdin)

    def assertRefused(self, result, *needles):
        out = result.stdout + result.stderr
        self.assertNotEqual(result.returncode, 0, f"expected a refusal:\n{out}")
        for needle in needles:
            self.assertIn(needle, out)

    def assertAllowed(self, result):
        self.assertEqual(result.returncode, 0, f"expected success:\n{result.stdout}{result.stderr}")

    def head(self):
        return self.git("rev-parse", "--short", "HEAD").stdout.strip()

    def foo_change(self, n=3):
        return {"analysis/foo/calc.py": CALC + f"\n# revision {n}\n",
                "analysis/foo/test_calc.py": TEST_CALC + f"assert double({n}) == {2 * n}\n"}

    # plan before code --------------------------------------------------------

    def test_code_without_plan_trailer_is_refused(self):
        self.assertRefused(self.commit("Change calc", self.foo_change()), "plan-trailer")

    def test_plan_in_same_commit_as_code_is_refused(self):
        files = dict(self.foo_change(), **{"plans/new.md": plan("New", scope="`analysis/foo/`")})
        self.assertRefused(self.commit("Change\n\nPlan: plans/new.md", files), "plan-before")

    def test_draft_plan_is_refused(self):
        self.assertRefused(self.commit("Change\n\nPlan: plans/draft.md", self.foo_change()), "plan-approved")

    def test_plan_missing_a_verification_subsection_is_refused(self):
        self.assertRefused(self.commit("Change\n\nPlan: plans/incomplete.md", self.foo_change()),
                           "plan-complete", "Invariants")

    def test_bare_na_is_refused(self):
        self.assertRefused(self.commit("Change\n\nPlan: plans/na.md", self.foo_change()),
                           "Regression", "n/a without a reason")

    def test_file_outside_plan_scope_is_refused(self):
        # Valid in every other way, so pre-commit passes and commit-msg is what refuses it.
        files = {"analysis/baz/x.py": "X = 1\n",
                 "analysis/baz/test_x.py": "import os, sys\nsys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))\nimport x\n",
                 "analysis/baz/README.md": FOO_README.replace("`test_calc.py` checks `calc.py`", "`test_x.py` checks `x.py`")}
        self.assertRefused(self.commit("Add baz\n\nPlan: plans/work.md", files), "plan-scope", "analysis/baz/x.py")

    def test_complete_plan_with_code_and_test_is_allowed_and_recorded(self):
        parent = self.head()
        self.assertAllowed(self.commit("Change calc\n\nPlan: plans/work.md", self.foo_change()))
        readme = self.git("show", "HEAD:analysis/foo/README.md").stdout
        self.assertEqual(readme.count("### Last verified"), 1)
        self.assertIn(f"on top of `{parent}`", readme)
        self.assertIn("3 passed", readme)
        self.assertEqual(self.git("status", "--porcelain").stdout.strip(), "")

    def test_second_run_replaces_the_block(self):
        self.assertAllowed(self.commit("One\n\nPlan: plans/work.md", self.foo_change(3)))
        first = self.head()
        self.assertAllowed(self.commit("Two\n\nPlan: plans/work.md", self.foo_change(4)))
        readme = (self.root / "analysis/foo/README.md").read_text()
        self.assertEqual(readme.count("### Last verified"), 1)
        self.assertEqual(readme.count("- Fast:"), 1)
        self.assertIn(f"on top of `{first}`", readme)

    # verification with code ------------------------------------------------

    def test_analysis_change_without_test_change_is_refused(self):
        files = {"analysis/foo/calc.py": CALC + "# tweak\n"}
        self.assertRefused(self.commit("Tweak\n\nPlan: plans/work.md", files), "paired", "analysis/foo/calc.py")

    def test_verified_by_an_exercising_test_is_allowed(self):
        files = {"analysis/foo/calc.py": CALC + "# tweak\n"}
        msg = "Tweak\n\nPlan: plans/work.md\nVerified-by: analysis/foo/test_calc.py — comment only"
        self.assertAllowed(self.commit(msg, files))

    def test_verified_by_a_test_that_does_not_exercise_the_file_is_refused(self):
        files = {"analysis/foo/calc.py": CALC + "# tweak\n"}
        msg = "Tweak\n\nPlan: plans/work.md\nVerified-by: analysis/foo/test_nothing.py — trust me"
        self.assertRefused(self.commit(msg, files), "paired", "does not exercise it")

    def test_new_script_no_test_runs_is_refused(self):
        files = {"analysis/foo/extra.py": "Y = 2\n", "analysis/foo/README.md": FOO_README.replace(
            "`calc.py`", "`calc.py` and `extra.py`")}
        self.assertRefused(self.commit("Add extra\n\nPlan: plans/work.md", files), "run-by-test", "extra.py")

    # documentation ---------------------------------------------------------

    def test_readme_without_verification_section_is_refused(self):
        files = {"analysis/bar/tool.py": "Z = 3\n",
                 "analysis/bar/test_tool.py": "import os, sys\nsys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))\nimport tool\n",
                 "analysis/bar/README.md": "# bar\n\nA tool.\n"}
        self.assertRefused(self.commit("Add bar\n\nPlan: plans/work.md", files),
                           "readme-verification", "no ## Verification")

    def test_readme_that_omits_the_file_is_refused(self):
        files = {"analysis/bar/tool.py": "Z = 3\n",
                 "analysis/bar/test_tool.py": "import os, sys\nsys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))\nimport tool\n",
                 "analysis/bar/README.md": "# bar\n\n## Verification\n\nSomething.\n\n### Last verified\n\n- Fast: not run yet\n"}
        self.assertRefused(self.commit("Add bar\n\nPlan: plans/work.md", files),
                           "readme-verification", "tool.py")

    def test_manual_code_without_block_update_is_refused(self):
        files = {"md_setup/build.py": "print('built twice')\n"}
        self.assertRefused(self.commit("Build\n\nPlan: plans/work.md", files), "last-verified", "did not")

    def test_manual_code_with_block_update_is_allowed(self):
        files = {"md_setup/build.py": "print('built twice')\n",
                 "md_setup/README.md": MD_README.replace("2026-01-01", "2026-02-02")}
        self.assertAllowed(self.commit("Build\n\nPlan: plans/work.md", files))

    def test_two_last_verified_blocks_are_refused(self):
        files = {"md_setup/build.py": "print('x')\n",
                 "md_setup/README.md": MD_README + "\n### Last verified\n\n- 2025: older run\n"}
        self.assertRefused(self.commit("Build\n\nPlan: plans/work.md", files), "last-verified", "found 2")

    def test_undeclared_import_is_an_error_in_analysis(self):
        files = self.foo_change()
        files["analysis/foo/calc.py"] = "import requests\n" + CALC
        self.assertRefused(self.commit("Use requests\n\nPlan: plans/work.md", files), "deps", "requests")

    def test_undeclared_import_is_a_warning_elsewhere(self):
        files = {"md_setup/build.py": "import requests\n",
                 "md_setup/README.md": MD_README.replace("2026-01-01", "2026-03-03")}
        result = self.commit("Use requests\n\nPlan: plans/work.md", files)
        self.assertAllowed(result)
        self.assertIn("warning [deps]", result.stdout + result.stderr)  # git sends hook output to stderr

    # exemptions --------------------------------------------------------------

    def test_docs_and_plan_only_commit_needs_no_plan(self):
        files = {"CLAUDE.md": CLAUDE + "\nMore words.\n", "plans/work.md": plan(title="Work, reworded")}
        self.assertAllowed(self.commit("Reword docs", files))

    def test_pure_rename_needs_no_plan(self):
        self.git("mv", "md_setup/build.py", "md_setup/build_it.py")
        self.write("md_setup/README.md", MD_README.replace("build.py", "build_it.py"))
        self.assertAllowed(self.commit("Rename build"))

    # verification run on push ----------------------------------------------

    def remote(self):
        bare = self.root.parent / f"{self.root.name}.git"
        run(["git", "init", "-q", "--bare", str(bare)], self.root.parent)
        self.git("remote", "add", "origin", str(bare))
        self.git("push", "-q", "--no-verify", "origin", "HEAD:refs/heads/master")

    def test_push_with_a_bad_commit_is_refused(self):
        self.remote()
        self.assertAllowed(self.commit("Sneak a change", self.foo_change(), verify=False))
        result = self.git("push", "origin", "HEAD:refs/heads/master", check=False)
        self.assertRefused(result, "plan-trailer", "Sneak a change")

    def test_push_with_a_failing_suite_is_refused(self):
        self.remote()
        files = {"analysis/foo/test_calc.py": TEST_CALC + "assert double(2) == 5\n",
                 "analysis/foo/README.md": FOO_README.replace("not run yet", "claimed")}
        self.assertAllowed(self.commit("Break test\n\nPlan: plans/work.md", files, verify=False))
        result = self.git("push", "origin", "HEAD:refs/heads/master", check=False)
        self.assertRefused(result, "[tests]", "fake failed")

    def test_push_of_good_commits_is_allowed(self):
        self.remote()
        self.assertAllowed(self.commit("Change calc\n\nPlan: plans/work.md", self.foo_change()))
        result = self.git("push", "origin", "HEAD:refs/heads/master", check=False)
        self.assertAllowed(result)
        self.assertIn("verification run (1 commit(s))", result.stdout + result.stderr)

    def test_ci_mode_agrees_with_pre_push(self):
        """Independent paths to one verdict: the pre-push hook and CI's --base --verify."""
        self.remote()
        base = self.git("rev-parse", "HEAD").stdout.strip()
        self.assertAllowed(self.commit("Sneak a change", self.foo_change(), verify=False))
        push = self.git("push", "origin", "HEAD:refs/heads/master", check=False)
        ci = self.gate("--base", base, "--verify")
        self.assertNotEqual(push.returncode, 0)
        self.assertNotEqual(ci.returncode, 0)
        self.assertIn("plan-trailer", ci.stdout)

    # Claude Code write-time check ------------------------------------------

    def test_plan_covers(self):
        self.assertEqual(self.gate("--plan-covers", "analysis/foo/calc.py").returncode, 0)
        uncovered = self.gate("--plan-covers", "analysis/zzz/new.py")
        self.assertEqual(uncovered.returncode, 2)
        self.assertIn("no approved plan", uncovered.stderr)
        self.assertEqual(self.gate("--plan-covers", "analysis/zzz/notes.md").returncode, 0)
        self.assertEqual(self.gate("--plan-covers", "/tmp/elsewhere.py").returncode, 0)
        payload = json.dumps({"tool_input": {"file_path": str(self.root / "analysis/zzz/new.py")}})
        self.assertEqual(self.gate("--plan-covers", "-", stdin=payload).returncode, 2)
        payload = json.dumps({"tool_input": {"file_path": str(self.root / "analysis/foo/calc.py")}})
        self.assertEqual(self.gate("--plan-covers", "-", stdin=payload).returncode, 0)

    def test_draft_plan_does_not_cover(self):
        self.write("plans/only_draft.md", plan("Only draft", state="draft", scope="`analysis/qq/`"))
        self.assertEqual(self.gate("--plan-covers", "analysis/qq/a.py").returncode, 2)

    # the earlier checks ----------------------------------------------------

    def test_force_added_ignored_file_is_refused(self):
        self.write("analysis/foo/run.log", "output\n")
        self.git("add", "-f", "analysis/foo/run.log")
        self.assertRefused(self.commit("Add log", add_all=False), "artifacts", "run.log")

    def test_hardcoded_user_path_is_refused(self):
        files = {"md_setup/build.py": "DATA = '/scratch1/lkyamamo/run'\n",
                 "md_setup/README.md": MD_README.replace("2026-01-01", "2026-04-04")}
        self.assertRefused(self.commit("Path\n\nPlan: plans/work.md", files), "user-paths")

    def test_overridable_user_path_is_allowed(self):
        self.write("md_setup/run.slurm", '#!/bin/bash\nlmp="${LMP:-/home1/lkyamamo/lmp}"\n')
        self.git("add", "-A")
        self.assertNotIn("user-paths", self.gate("--no-tests").stdout)

    def test_identical_file_is_refused(self):
        blob = "x" * 300 + "\n"
        self.assertRefused(self.commit("Data", {"analysis/foo/a.json": blob, "analysis/foo/b.json": blob}),
                           "duplicates")

    def test_syntax_error_is_refused(self):
        files = self.foo_change()
        files["analysis/foo/calc.py"] = "def double(x:\n"
        self.assertRefused(self.commit("Break\n\nPlan: plans/work.md", files), "syntax")

    def test_unmapped_top_level_directory_is_refused(self):
        self.assertRefused(self.commit("New dir", {"newtop/notes.md": "hi\n"}), "docs", "newtop/")

    # invariants ------------------------------------------------------------

    def test_all_mode_is_deterministic(self):
        first = self.gate("--all", "--json", "--no-tests").stdout
        second = self.gate("--all", "--json", "--no-tests").stdout
        self.assertEqual(first, second)
        self.assertTrue(json.loads(first) is not None)


class DebtBaseline(unittest.TestCase):
    """Regression: the repo's existing debt only changes deliberately.

    When this fails, the change either fixed debt or added some. Regenerate the
    baseline in the same commit and say which in the commit message:
        python3 tools/repo_check.py --all --json --no-tests > tools/tests/baselines/debt.json
    """

    def test_debt_matches_baseline(self):
        args = [sys.executable, str(GATE), "--all", "--json", "--no-tests"]
        ref = os.environ.get("REPO_CHECK_REF")
        if ref and ref != ":":
            args += ["--ref", ref]
        env = {k: v for k, v in os.environ.items() if k != "REPO_CHECK_ROOT"}
        env["REPO_CHECK_ROOT"] = str(REAL_REPO)
        result = subprocess.run(args, capture_output=True, text=True, env=env, cwd=REAL_REPO)
        self.assertEqual(result.returncode, 0, result.stderr)
        current = json.loads(result.stdout)
        baseline = json.loads(BASELINE.read_text())
        fixed = [f for f in baseline if f not in current]
        added = [f for f in current if f not in baseline]
        self.assertFalse(fixed or added, (
            f"debt changed: {len(fixed)} fixed, {len(added)} new. Regenerate the baseline in this commit.\n"
            + "\n".join(f"  new:   {f}" for f in added[:20])
            + "\n" + "\n".join(f"  fixed: {f}" for f in fixed[:20])
        ))


def run(cmd, cwd, check=True, stdin=None):
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True, env=clean_env(), input=stdin)
    if check and result.returncode:
        raise AssertionError(f"{' '.join(cmd)} failed:\n{result.stdout}{result.stderr}")
    return result


if __name__ == "__main__":
    unittest.main()
