#!/usr/bin/env python3
"""The checks every change to this repo must pass. See tools/README.md.

Run by the git hooks in .githooks/ and by CI:

    pre-commit   repo_check.py                    staged change: checks, the touched suites,
                                                  then rewrites and stages each affected
                                                  README's "### Last verified" block
    commit-msg   repo_check.py --commit-msg FILE  plan and pairing checks on the message
    pre-push     repo_check.py --push             the verification run: every commit being
                                                  pushed, the range, and every fast suite
    CI           repo_check.py --base REF --verify

Other modes:
    --all [--ref REF]       every tracked file: the existing debt
    --commit SHA            the per-commit checks for one commit
    --record [--slow]       run suites on the index and rewrite README blocks
    --plan-covers PATH|-    exit 2 unless an approved plan's Scope covers PATH
                            ("-" reads a Claude Code hook payload from stdin)
    --no-tests              skip running suites

Checks:
  artifacts            force-added ignored file, compiled binary, file over 1 MB, backup copy
  user-paths           new line hardcoding /home1/lkyamamo, /scratch*/lkyamamo, /Users/loganyamamoto
  duplicates           file identical to a tracked one (error); same source name elsewhere (warning)
  syntax               .py that does not compile; shell or hook script bash -n rejects
  docs                 CLAUDE.md map out of date; a doc naming a path this change removes
  deps                 import missing from pyproject.toml (error in analysis/, else warning)
  run-by-test          analysis script no test runs
  readme-verification  code whose README lacks "## Verification" naming it and its tests
  last-verified        README without exactly one "### Last verified" block, or not updated
  plan-*               code commit without an approved, complete, earlier-committed plan
                       whose Scope covers it
  paired               code in a verified area changed without its tests
  tests                a suite failing
"""

import argparse
import ast
import datetime
import fnmatch
import json
import os
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
import time
from collections import defaultdict
from pathlib import Path, PurePosixPath

ROOT = Path(os.environ.get("REPO_CHECK_ROOT") or Path(__file__).resolve().parent.parent)
CONFIG = Path(__file__).resolve().parent / "verification.json"

CODE_EXTS = {".py", ".sh", ".cpp", ".c", ".h", ".f90", ".f"}
HOOK_DIR = ".githooks/"
SHELL_EXTS = {".sh", ".slurm", ".pbs", ".sbatch"}

MAX_BYTES = 1_000_000
LARGE_FILE_ALLOW = ["starting-structures/*", "*/uv.lock", "uv.lock"]

BACKUP_PATTERNS = [
    "*/bkup/*", "*/backup/*", "*/old/*",
    "*.bak", "*.orig", "*_old.*", "*_bak.*", "*_backup.*", "* copy*", "*_copy.*",
]

USER_PATH = re.compile(r"/home1/lkyamamo|/scratch\d*/lkyamamo|/Users/loganyamamoto")
USER_PATH_EXTS = {
    ".py", ".sh", ".slurm", ".pbs", ".sbatch", ".f90", ".f", ".cpp", ".h",
    ".input", ".in", ".toml", ".json", ".ipynb",
}  # not .conf/.config: those are where paths belong
USER_PATH_MARKER = "allow-path"
# A user path is fine as an overridable default: ${VAR:-/path} in shell, or an
# argparse default=. Comments and tests (synthetic fixtures) are exempt too.
OVERRIDABLE = re.compile(
    r"\$\{\w+:-[^}]*(" + USER_PATH.pattern + r")|default\s*=\s*[\"'](" + USER_PATH.pattern + r")"
)
USER_PATH_EXEMPT = ["*/tests/*", "*/test_*.py", "*.example", "tools/repo_check.py"]

DUP_MIN_BYTES = 256
SOURCE_EXTS = {".py", ".f90", ".f", ".cpp", ".c", ".h", ".sh", ".slurm", ".pbs"}
COMMON_NAMES = {
    "main.py", "__init__.py", "__main__.py", "cli.py", "config.py", "utils.py",
    "test.py", "run.sh", "setup.py", "plot.py", "main.f90", "main.cpp", "conftest.py",
}

# import name -> distribution name, where they differ
IMPORT_TO_DIST = {
    "sklearn": "scikit-learn", "yaml": "pyyaml", "PIL": "pillow",
    "freud": "freud-analysis", "MDAnalysis": "mdanalysis", "cv2": "opencv-python",
    "skimage": "scikit-image", "dotenv": "python-dotenv", "PySide6": "pyside6",
    "mpl_toolkits": "matplotlib", "log_lammps_reader": "log-lammps-reader",
    "docx": "python-docx",
}
STRICT_DEPS_AREA = "analysis/"

PLAN_DIR = "plans/"
PLAN_TEMPLATE = "plans/_TEMPLATE.md"
PLAN_SECTIONS = ["goal", "scope", "approach", "verification", "out of scope"]
VERIFICATION_SUBSECTIONS = [
    "known answer", "invariants", "independent implementation", "regression",
    "tolerances and sources", "slow checks", "done when",
]
PLAN_STATES = ("approved", "in progress")
PLAN_TRAILER = re.compile(r"^Plan:\s*(\S+)\s*$", re.M)
VERIFIED_TRAILER = re.compile(r"^Verified-by:\s*(\S+)\s+(?:—|--?)\s+(\S.*)$", re.M)
NA_ONLY = re.compile(r"^[-*\s]*n/?a[\s.:;,-]*$", re.I)

LAST_VERIFIED = "### Last verified"
BLOCK_NOTE = "<!-- Rewritten by tools/repo_check.py on every run that verifies this folder; never appended. -->"

ZERO_SHA = "0" * 40


# --------------------------------------------------------------------------
# git and config


def git(*args, input=None, ok=(0,), cwd=None):
    r = subprocess.run(["git", *args], cwd=cwd or ROOT, input=input, capture_output=True)
    if r.returncode not in ok:
        raise subprocess.CalledProcessError(r.returncode, ["git", *args], r.stdout, r.stderr)
    return r.stdout


def git_text(*args, **kw):
    return git(*args, **kw).decode().strip()


def rev_exists(rev):
    return subprocess.run(
        ["git", "cat-file", "-e", f"{rev}^{{commit}}"], cwd=ROOT, capture_output=True
    ).returncode == 0


def load_suites():
    return json.loads(CONFIG.read_text())["suites"]


SUITES = load_suites()


def matches(path, patterns):
    return any(fnmatch.fnmatch(path, p) or fnmatch.fnmatch("/" + path, p) for p in patterns)


def suite_for(path):
    return next((s for s in SUITES if any(path.startswith(a) for a in s["areas"])), None)


def is_code(path):
    return Path(path).suffix in CODE_EXTS or path.startswith(HOOK_DIR)


def is_verification(path):
    return any(matches(path, s["tests"]) for s in SUITES)


def is_test_file(path):
    return is_verification(path) and PurePosixPath(path).name.startswith("test_") and path.endswith(".py")


def is_analysis_code(path):
    return path.startswith(STRICT_DEPS_AREA) and path.endswith(".py") and not is_verification(path)


# --------------------------------------------------------------------------
# the file set under check


class Change:
    """What changed, the tree it changed into, and how to read both sides.

    staged  index vs HEAD          range   merge-base(base, head) vs head
    commit  SHA^ vs SHA            all     every file of the index (or of ref)
    """

    def __init__(self, mode, base=None, head="HEAD", ref=None):
        self.mode = mode
        self._cache = {}
        has_head = rev_exists("HEAD")
        if mode == "staged":
            self.new, self.old = ":", ("HEAD" if has_head else None)
        elif mode == "range":
            self.new = head
            self.old = git_text("merge-base", base, head) if base else None
        elif mode == "commit":
            self.new = head
            parents = git_text("rev-list", "--parents", "-n", "1", head).split()[1:]
            self.old = parents[0] if parents else None
        else:
            self.new, self.old = (ref or ":"), None
        self.tracked = self._tree(self.new)
        self.old_tracked = self._tree(self.old) if self.old else {}
        self.entries = self._entries()

    def _tree(self, ref):
        tracked = {}
        if ref == ":":
            for rec in filter(None, git("ls-files", "-s", "-z").decode().split("\0")):
                meta, path = rec.split("\t", 1)
                tracked[path] = meta.split()[1]
        else:
            for rec in filter(None, git("ls-tree", "-r", "-z", ref).decode().split("\0")):
                meta, path = rec.split("\t", 1)
                tracked[path] = meta.split()[2]
        return tracked

    def _diff_args(self):
        if self.mode == "staged":
            return ["--cached"] + ([] if self.old else [git_text("hash-object", "-t", "tree", "/dev/null")])
        return [self.old or git_text("hash-object", "-t", "tree", "/dev/null"), self.new]

    def _entries(self):
        """List of (status, path, old_path, similarity); status is A, M, D or R."""
        if self.mode == "all":
            return [("A", p, None, 0) for p in sorted(self.tracked)]
        out = git("diff", *self._diff_args(), "--name-status", "-M", "-z").decode().split("\0")
        entries, i = [], 0
        while i < len(out) and out[i]:
            code = out[i]
            status = code[0]
            if status in "RC":
                sim = int(code[1:] or 0)
                entries.append(("R" if status == "R" else "A", out[i + 2], out[i + 1], sim))
                i += 3
            else:
                entries.append((status if status in "AMD" else "M", out[i + 1], None, 0))
                i += 2
        return entries

    def read(self, path, old=False):
        ref = self.old if old else self.new
        key = (ref, path)
        if key not in self._cache:
            spec = f":{path}" if ref == ":" else f"{ref}:{path}"
            self._cache[key] = git("show", spec)
        return self._cache[key]

    def text(self, path, old=False):
        return self.read(path, old).decode(errors="replace")

    def present(self):
        return [p for s, p, _, _ in self.entries if s != "D" and p in self.tracked]

    def new_paths(self):
        return [p for s, p, _, _ in self.entries if s in "AR" and p in self.tracked]

    def gone(self):
        return [(old if s == "R" else p) for s, p, old, _ in self.entries if s in "DR"]

    def touched(self):
        return [p for _, p, _, _ in self.entries] + [o for _, _, o, _ in self.entries if o]

    def structural(self):
        return any(s in "ADR" for s, _, _, _ in self.entries)

    def code_changes(self):
        """Code files added or edited. Deletions and pure (100%) renames change no code."""
        return [
            p for s, p, _, sim in self.entries
            if p in self.tracked and is_code(p) and (s in "AM" or (s == "R" and sim < 100))
        ]

    def added_lines(self, path):
        if self.mode == "all":
            return self.text(path).splitlines()
        diff = git("diff", *self._diff_args(), "-U0", "--", path).decode(errors="replace")
        return [
            line[1:] for line in diff.splitlines()
            if line.startswith("+") and not line.startswith("+++")
        ]


class Report:
    def __init__(self):
        self.errors = []
        self.warnings = []

    def error(self, check, path, msg):
        self.errors.append((check, path, msg))

    def warn(self, check, path, msg):
        self.warnings.append((check, path, msg))

    def extend(self, other, prefix=""):
        self.errors += [(c, prefix + p, m) for c, p, m in other.errors]
        self.warnings += [(c, prefix + p, m) for c, p, m in other.warnings]

    def lines(self):
        out = []
        for level, items in (("warning", self.warnings), ("error", self.errors)):
            for check, path, msg in sorted(set(items)):
                out.append(f"{level} [{check}] {path}: {msg}")
        return out

    def print(self, label="repo_check"):
        for line in self.lines():
            print(line)
        print(f"{label}: {len(set(self.errors))} error(s), {len(set(self.warnings))} warning(s)")


# --------------------------------------------------------------------------
# which tests exercise which code, and which README documents it


def referenced_names(source):
    """Basenames and module names a test file uses: strings (not docstrings) and imports."""
    try:
        tree = ast.parse(source)
    except SyntaxError:
        return set()
    docstrings = {
        id(node.value) for node in ast.walk(tree)
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Constant)
    }
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str) and id(node) not in docstrings:
            names.update(re.findall(r"[\w.-]+\.py\b", node.value))
        elif isinstance(node, ast.Import):
            names.update(a.name.split(".")[-1] for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module and node.level == 0:
            names.add(node.module.split(".")[-1])
            names.update(a.name for a in node.names)
    return names


class Coverage:
    """Maps code files to the tests that exercise them, on one tree."""

    def __init__(self, change):
        self.change = change
        self._refs = {}

    def refs(self, test):
        if test not in self._refs:
            self._refs[test] = referenced_names(self.change.read(test))
        return self._refs[test]

    def tests_of(self, path):
        suite = suite_for(path)
        if not suite:
            return []
        tests = [p for p in self.change.tracked if is_test_file(p) and suite_for(p) is suite]
        if not suite["per_file"]:
            return sorted(tests)
        name, stem = PurePosixPath(path).name, PurePosixPath(path).stem
        return sorted(t for t in tests if name in self.refs(t) or stem in self.refs(t))

    def code_tested_by(self, test):
        """Analysis code files a per-file test exercises; same-topic files win a name clash."""
        refs = self.refs(test)
        by_name = defaultdict(list)
        for p in self.change.tracked:
            if is_code(p) and not is_verification(p) and suite_for(p) is suite_for(test):
                by_name[PurePosixPath(p).name].append(p)
                by_name[PurePosixPath(p).stem].append(p)
        topic = "/".join(test.split("/")[:2]) + "/"
        found = set()
        for ref in refs:
            candidates = by_name.get(ref, [])
            local = [c for c in candidates if c.startswith(topic)]
            found.update(local or candidates)
        return found


def doc_for(path, tracked):
    """The README documenting a code file: nearest one up to its area root, else its suite's doc."""
    parts = PurePosixPath(path).parts
    stop = 2 if parts[0] == "analysis" else 1
    d = PurePosixPath(path).parent
    while len(d.parts) >= stop:
        candidate = str(d / "README.md")
        if candidate in tracked:
            return candidate
        d = d.parent
    suite = suite_for(path)
    return suite.get("doc") if suite else None


def section(text, heading, level=2):
    """Body of the first '#'*level heading named heading (case-insensitive), or None."""
    marks = "#" * level
    lines = text.splitlines()
    body, inside = [], False
    for line in lines:
        m = re.match(r"^(#{1,6})\s+(.*?)\s*$", line)
        if m:
            if inside and len(m.group(1)) <= level:
                break
            if len(m.group(1)) == level and m.group(2).lower() == heading.lower():
                inside = True
                continue
        if inside:
            body.append(line)
    return "\n".join(body) if inside else None


def last_verified_block(text):
    verification = section(text, "verification")
    if verification is None:
        return None
    return section(verification, "last verified", level=3)


# --------------------------------------------------------------------------
# checks on a change


def check_artifacts(change, report):
    paths = change.new_paths()
    if not paths:
        return
    # -v so a path re-included by a "!pattern" can be told apart; exit 1 means no match
    fields = git(
        "check-ignore", "--no-index", "--stdin", "-z", "-v",
        input="\0".join(paths).encode(), ok=(0, 1),
    ).decode().split("\0")
    ignored = {
        fields[i + 3] for i in range(0, len(fields) - 3, 4)
        if not fields[i + 2].startswith("!")
    }
    for path in paths:
        if path in ignored:
            report.error("artifacts", path, "matches .gitignore (force-added?); generated output stays out of git")
        if matches(path, BACKUP_PATTERNS):
            report.error("artifacts", path, "backup copy; delete it, git history keeps old versions")
        data = change.read(path)
        if data[:4] == b"\x7fELF" or data[:4] in (
            b"\xfe\xed\xfa\xce", b"\xfe\xed\xfa\xcf", b"\xce\xfa\xed\xfe",
            b"\xcf\xfa\xed\xfe", b"\xca\xfe\xba\xbe",
        ):
            report.error("artifacts", path, "compiled binary; commit the source and build it")
        if len(data) > MAX_BYTES and not matches(path, LARGE_FILE_ALLOW):
            report.error("artifacts", path, f"{len(data) / 1e6:.1f} MB exceeds the {MAX_BYTES / 1e6:.0f} MB limit")


def check_user_paths(change, report):
    for path in change.present():
        if Path(path).suffix not in USER_PATH_EXTS or matches(path, USER_PATH_EXEMPT):
            continue
        for line in change.added_lines(path):
            m = USER_PATH.search(line)
            stripped = line.lstrip()
            comment = stripped.startswith("#") and not stripped.startswith("#SBATCH")
            if m and USER_PATH_MARKER not in line and not comment and not OVERRIDABLE.search(line):
                report.error(
                    "user-paths", path,
                    f"hardcoded {m.group(0)}; read it from a .conf/.config or flag "
                    f"(or mark the line '{USER_PATH_MARKER}'): {line.strip()[:100]}",
                )


def check_duplicates(change, report):
    by_blob = defaultdict(list)
    by_name = defaultdict(list)
    for path, blob in change.tracked.items():
        by_blob[blob].append(path)
        by_name[PurePosixPath(path).name].append(path)
    for path in change.present():
        name = PurePosixPath(path).name
        twins = [p for p in by_blob[change.tracked[path]] if p != path]
        if twins and name not in COMMON_NAMES and len(change.read(path)) >= DUP_MIN_BYTES:
            report.error("duplicates", path, f"identical to {', '.join(twins)}; share one copy instead")
            continue
        if Path(path).suffix in SOURCE_EXTS and name not in COMMON_NAMES:
            others = [p for p in by_name[name] if p != path]
            if others:
                report.warn("duplicates", path, f"same name as {', '.join(others)}; one source or say which is current")


def check_syntax(change, report):
    for path in change.present():
        suffix = Path(path).suffix
        if suffix == ".py":
            try:
                compile(change.read(path), path, "exec")
            except SyntaxError as e:
                report.error("syntax", path, f"line {e.lineno}: {e.msg}")
        elif suffix in SHELL_EXTS or path.startswith(HOOK_DIR):
            data = change.read(path)
            first = data.split(b"\n", 1)[0]
            if first.startswith(b"#!") and b"sh" not in first:
                continue
            if any(sh in first for sh in (b"csh", b"tcsh", b"zsh")):
                continue
            r = subprocess.run(["bash", "-n"], input=data, capture_output=True)
            if r.returncode:
                report.error("syntax", path, r.stderr.decode(errors="replace").strip().replace("bash: ", ""))


DOC_PATH = re.compile(r"`([A-Za-z0-9_.-]+(?:/[A-Za-z0-9_.*-]+)*/?)`")


def check_docs(change, report):
    claude_md = "CLAUDE.md"
    if change.mode == "all" or change.structural() or claude_md in change.touched():
        if claude_md not in change.tracked:
            report.error("docs", claude_md, "missing; the repo map lives here")
        else:
            text = change.text(claude_md)
            top_dirs = {p.split("/")[0] for p in change.tracked if "/" in p}
            listed = set(DOC_PATH.findall(text))
            for d in sorted(top_dirs):
                if not d.startswith(".") and f"{d}/" not in listed:
                    report.error("docs", claude_md, f"top-level `{d}/` is not in the repo map")
            tracked_dirs = {str(parent) for p in change.tracked for parent in PurePosixPath(p).parents}
            for ref in sorted(listed):
                if "/" not in ref or ref.split("/")[0] not in top_dirs:
                    continue
                target = ref.rstrip("/")
                if "*" in ref:
                    exists = any(fnmatch.fnmatch(p, ref) for p in change.tracked)
                else:
                    exists = target in change.tracked or target in tracked_dirs
                if not exists:
                    report.error("docs", claude_md, f"names `{ref}`, which is not tracked")

    gone = change.gone()
    if not gone:
        return
    # plans describe moves and deletions, so they name paths on purpose
    docs = [p for p in change.tracked if p.endswith(".md") and not p.startswith(PLAN_DIR)]
    for doc in docs:
        text = change.text(doc)
        for old in gone:
            if old in text:
                report.error("docs", doc, f"still names `{old}`, which this change removes or moves")


def pyproject_deps(pyproject_text):
    m = re.search(r"^dependencies\s*=\s*\[(.*?)\]", pyproject_text, re.S | re.M)
    if not m:
        return set()
    names = re.findall(r"[\"']([A-Za-z0-9_.\-]+)", m.group(1))
    return {re.sub(r"[_.]", "-", n).lower() for n in names}


def check_deps(change, report):
    stdlib = getattr(sys, "stdlib_module_names", None)
    if stdlib is None:
        return
    projects = {str(PurePosixPath(p).parent) for p in change.tracked if PurePosixPath(p).name == "pyproject.toml"}
    for path in change.present():
        if not path.endswith(".py"):
            continue
        project = next((str(d) for d in PurePosixPath(path).parents if str(d) in projects), None)
        if project is None:
            continue
        try:
            tree = ast.parse(change.read(path))
        except SyntaxError:
            continue
        imported = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                imported |= {a.name.split(".")[0] for a in node.names}
            elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
                imported.add(node.module.split(".")[0])
        prefix = "" if project == "." else project + "/"
        local = set()
        for p in change.tracked:
            if p.startswith(prefix):
                local |= {PurePosixPath(part).stem for part in PurePosixPath(p[len(prefix):]).parts}
        declared = pyproject_deps(change.text(f"{prefix}pyproject.toml"))
        # pytest is the test runner, declared as a dev dependency
        dev = {"pytest"} if is_verification(path) else set()
        for mod in sorted(imported - set(stdlib) - local - {"__future__"} - dev):
            dist = IMPORT_TO_DIST.get(mod, mod).replace("_", "-").lower()
            if dist not in declared:
                msg = f"imports {mod}, not declared in {prefix}pyproject.toml (uv add {dist})"
                if path.startswith(STRICT_DEPS_AREA):
                    report.error("deps", path, msg)
                else:
                    report.warn("deps", path, msg)


def check_run_by_test(change, report, coverage):
    targets = change.present() if change.mode == "all" else change.code_changes()
    for path in targets:
        if is_analysis_code(path) and not coverage.tests_of(path):
            report.error("run-by-test", path, "no test runs this script; add one that runs or imports it")


def documented_targets(change, coverage):
    """(code files whose README must name them and their tests, READMEs whose block must update)."""
    changed = change.present() if change.mode == "all" else change.code_changes()
    targets, docs = set(), set()
    for path in changed:
        doc = doc_for(path, change.tracked)
        if is_verification(path):
            if doc and doc in change.tracked and LAST_VERIFIED in change.text(doc):
                docs.add(doc)
            if is_test_file(path) and suite_for(path) and suite_for(path)["per_file"]:
                targets |= coverage.code_tested_by(path)
            continue
        targets.add(path)
    for path in targets:
        doc = doc_for(path, change.tracked)
        if doc:
            docs.add(doc)
    return sorted(targets), sorted(docs)


def check_readme_verification(change, report, coverage, targets):
    for path in targets:
        doc = doc_for(path, change.tracked)
        if not doc or doc not in change.tracked:
            here = str(PurePosixPath(path).parent / "README.md")
            report.error("readme-verification", path, f"no README documents it; add {here} with a ## Verification section")
            continue
        verification = section(change.text(doc), "verification")
        if verification is None:
            report.error("readme-verification", doc, f"no ## Verification section (documents {path})")
            continue
        missing = [n for n in [PurePosixPath(path).name] + [PurePosixPath(t).name for t in coverage.tests_of(path)]
                   if n not in verification]
        if missing:
            report.error("readme-verification", doc, f"## Verification does not name {', '.join(missing)} (for {path})")


def check_last_verified(change, report, docs):
    for doc in docs:
        if doc not in change.tracked:
            continue
        text = change.text(doc)
        count = len(re.findall(r"^### Last verified\s*$", text, re.M))
        if count != 1 or last_verified_block(text) is None:
            report.error("last-verified", doc,
                         f"needs exactly one '{LAST_VERIFIED}' block under ## Verification (found {count})")
            continue
        if change.mode == "all":
            continue
        old = last_verified_block(change.text(doc, old=True)) if doc in change.old_tracked else None
        if old == last_verified_block(text):
            report.error("last-verified", doc,
                         "code it documents changed but its Last verified block did not; record the run that verified this change")


def parse_plan(text):
    """Problems with a plan's completeness, and its Scope entries."""
    problems = []
    state = re.search(r"^\*\*State:\*\*\s*(.*)$", text, re.M | re.I)
    state_value = state.group(1).strip().lower() if state else ""
    for name in PLAN_SECTIONS:
        body = section(text, name)
        if body is None or not body.strip():
            problems.append(f"missing or empty ## {name.title()}")
    verification = section(text, "verification") or ""
    for name in VERIFICATION_SUBSECTIONS:
        body = section(verification, name, level=3)
        if body is None or not body.strip():
            problems.append(f"verification plan missing or empty ### {name.title()}")
        elif NA_ONLY.match(body.strip()):
            problems.append(f"### {name.title()} says n/a without a reason")
    scope = [s for s in DOC_PATH.findall(section(text, "scope") or "") if "/" in s or "." in s]
    return state_value, problems, scope


def in_scope(path, scope):
    for entry in scope:
        if path == entry or fnmatch.fnmatch(path, entry):
            return True
        if path.startswith(entry.rstrip("/") + "/"):
            return True
    return False


def check_plan(change, report, message):
    code = change.code_changes()
    if not code:
        return
    m = PLAN_TRAILER.search(message or "")
    if not m:
        report.error("plan-trailer", code[0],
                     f"this commit changes code ({len(code)} file(s)) but its message has no 'Plan: plans/<slug>.md' trailer")
        return
    plan = m.group(1)
    if not plan.startswith(PLAN_DIR) or plan == PLAN_TEMPLATE:
        report.error("plan-trailer", plan, "the Plan trailer must name a plan in plans/")
        return
    if plan not in change.old_tracked:
        report.error("plan-before", plan, "the plan must be committed before the code it covers (it is not in the parent commit)")
        return
    state, problems, scope = parse_plan(change.text(plan, old=True))
    if not state.startswith(PLAN_STATES):
        report.error("plan-approved", plan, f"State is '{state or 'missing'}'; it must be approved (or in progress) before code")
    for p in problems:
        report.error("plan-complete", plan, p)
    for path in code:
        if not in_scope(path, scope):
            report.error("plan-scope", path, f"not in the Scope of {plan}")


def check_paired(change, report, coverage, message):
    changed = set(change.code_changes())
    verified_by = {m.group(1): m.group(2) for m in VERIFIED_TRAILER.finditer(message or "")}
    for path in sorted(changed):
        if is_verification(path):
            continue
        suite = suite_for(path)
        if not suite:
            continue
        tests = coverage.tests_of(path)
        if any(t in changed for t in tests):
            continue
        named = [t for t in verified_by if t in tests]
        if named:
            continue
        bad = [t for t in verified_by if t not in tests]
        hint = f" ('Verified-by: {bad[0]}' does not exercise it)" if bad else ""
        report.error("paired", path,
                     f"changed without a change to its tests ({', '.join(tests) or 'none exist'}); update them, "
                     f"or add 'Verified-by: <test> — <reason>' naming a test that covers this change{hint}")


# --------------------------------------------------------------------------
# suites


def suite_python(suite):
    if suite["python"] == "stdlib":
        return sys.executable
    env = os.environ.get("ANALYSIS_PYTHON")
    if env:
        return env
    venv = ROOT / "analysis/.venv/bin/python"
    return str(venv) if venv.exists() else None


def export_tree(treeish, dest):
    """Write a tree (commit, tree id, or ':' for the index) into dest."""
    if treeish == ":":
        treeish = git_text("write-tree")
    data = git("archive", "--format=tar", treeish)
    with tempfile.TemporaryFile() as f:
        f.write(data)
        f.seek(0)
        with tarfile.open(fileobj=f) as tar:
            tar.extractall(dest)


def summarize(output):
    m = re.findall(r"(\d+ passed[^\n]*?)(?: in [\d.]+s)?\s*=*\s*$", output, re.M)
    if m:
        return m[-1].strip()
    m = re.search(r"^Ran (\d+) tests? in", output, re.M)
    if m:
        return f"{m.group(1)} tests OK"
    return "passed"


def run_suites(suites, tree, report, kind="fast", ref=":"):
    """Run each suite's fast (or slow) command in tree, exported from ref.

    Returns {name: result line, "none defined", or None when it failed or could not run}.
    """
    results = {}
    for suite in suites:
        args = suite.get(kind)
        if not args:
            continue
        python = suite_python(suite)
        if python is None:
            report.warn("tests", suite["cwd"], f"{suite['name']}: no interpreter (set ANALYSIS_PYTHON); not run")
            results[suite["name"]] = None
            continue
        print(f"running {suite['name']} {kind} suite...", flush=True)
        # No GIT_* variables: inside a hook they point at this repo's index, and a
        # suite that runs git (the gate's own tests) must not touch it.
        env = {k: v for k, v in os.environ.items() if not k.startswith("GIT_")}
        env.update(REPO_CHECK_ROOT=str(ROOT), REPO_CHECK_REF=ref, PYTHONDONTWRITEBYTECODE="1")
        start = time.time()
        r = subprocess.run([python, *args], cwd=Path(tree) / suite["cwd"], capture_output=True, env=env)
        secs = time.time() - start
        output = (r.stdout + r.stderr).decode(errors="replace")
        if kind == "slow" and r.returncode == 5:  # pytest: nothing collected
            results[suite["name"]] = "none defined"
            continue
        if r.returncode:
            tail = output.strip().splitlines()[-20:]
            report.error("tests", suite["cwd"], f"{suite['name']} failed:\n    " + "\n    ".join(tail))
            results[suite["name"]] = None
            continue
        cmd = "python " + " ".join(args)
        results[suite["name"]] = f"`{cmd}` in `{suite['cwd'].rstrip('/') or '.'}/`: {summarize(output)} ({secs:.1f} s)"
    return results


def suites_touching(paths):
    return [s for s in SUITES if any(p.startswith(a) for p in paths for a in s["areas"])]


def write_block(path, kind, line):
    """Rewrite the Fast or Slow line of a README's Last verified block, creating the block if needed."""
    lines = Path(path).read_text().splitlines()
    today = datetime.date.today().isoformat()
    parent = git_text("rev-parse", "--short", "HEAD") if rev_exists("HEAD") else "the root commit"
    entry = f"- {kind.title()}: none defined" if line == "none defined" else \
        f"- {kind.title()}: {today} on top of `{parent}`: {line}"

    heading = next((i for i, l in enumerate(lines) if l.strip() == LAST_VERIFIED), None)
    if heading is None:
        start = next((i for i, l in enumerate(lines) if l.strip() == "## Verification"), None)
        if start is None:
            return False
        end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith("## ")), len(lines))
        while end > start + 1 and not lines[end - 1].strip():
            end -= 1
        lines[end:end] = ["", LAST_VERIFIED]
        heading = end + 1
    end = next((i for i in range(heading + 1, len(lines)) if re.match(r"^#{1,3}\s", lines[i])), len(lines))
    old = [l for l in lines[heading + 1:end] if l.startswith("- ")]
    fast = next((l for l in old if l.startswith("- Fast:")), "- Fast: not run yet")
    slow = next((l for l in old if l.startswith("- Slow:")), "- Slow: none defined")
    if kind == "fast":
        fast = entry
    else:
        slow = entry
    block = ["", BLOCK_NOTE, fast, slow]
    if end < len(lines):
        block.append("")
    lines[heading + 1:end] = block
    Path(path).write_text("\n".join(lines) + "\n")
    return True


# --------------------------------------------------------------------------
# modes


def range_checks(change, report, coverage):
    for check in (check_artifacts, check_user_paths, check_duplicates, check_syntax, check_docs, check_deps):
        check(change, report)
    check_run_by_test(change, report, coverage)


def doc_checks(change, report, coverage):
    targets, docs = documented_targets(change, coverage)
    check_readme_verification(change, report, coverage, targets)
    check_last_verified(change, report, docs)


def mode_pre_commit(args):
    change = Change("staged")
    if not change.entries:
        return 0
    report = Report()
    coverage = Coverage(change)
    range_checks(change, report, coverage)

    targets, docs = documented_targets(change, coverage)
    results = {}
    if not args.no_tests:
        suites = suites_touching(change.touched())
        if suites:
            with tempfile.TemporaryDirectory(prefix="repo_check_") as tree:
                export_tree(":", tree)
                results = run_suites(suites, tree, report)

    recorded = set()
    for doc in docs:
        suite = suite_for(doc)
        line = results.get(suite["name"]) if suite else None
        if not line or doc not in change.tracked:
            continue
        unstaged = git("diff", "--name-only", "--", doc).decode().strip()
        if unstaged:
            report.error("last-verified", doc, "has unstaged edits; stage or stash them so the run can be recorded")
            continue
        if write_block(ROOT / doc, "fast", line):
            git("add", "--", doc)
            recorded.add(doc)

    change = Change("staged")
    coverage = Coverage(change)
    check_readme_verification(change, report, coverage, targets)
    failed_suites = {s["name"] for s in SUITES if results.get(s["name"], "") is None}
    pending = [d for d in docs if suite_for(d) and suite_for(d)["name"] in failed_suites]
    check_last_verified(change, report, [d for d in docs if d not in pending])
    report.print()
    return 1 if report.errors else 0


def mode_commit_msg(args):
    message = Path(args.commit_msg).read_text()
    message = "\n".join(l for l in message.splitlines() if not l.startswith("#"))
    change = Change("staged")
    report = Report()
    check_plan(change, report, message)
    check_paired(change, report, Coverage(change), message)
    if report.errors or report.warnings:
        report.print("repo_check commit-msg")
    return 1 if report.errors else 0


def commit_report(sha):
    change = Change("commit", head=sha)
    message = git("log", "-1", "--format=%B", sha).decode()
    report = Report()
    coverage = Coverage(change)
    check_plan(change, report, message)
    check_paired(change, report, coverage, message)
    doc_checks(change, report, coverage)
    return report


def verify_range(base, head, report, run_tests=True):
    """The verification run for base..head: range checks, per-commit checks, every fast suite."""
    if base is None:
        commits = git_text("rev-list", "--no-merges", head, "--not", "--remotes").split()
        oldest = commits[-1] if commits else head
        parents = git_text("rev-list", "--parents", "-n", "1", oldest).split()[1:]
        base = parents[0] if parents else None
    else:
        commits = git_text("rev-list", "--no-merges", f"{base}..{head}").split()
    change = Change("range", base=base, head=head)
    range_checks(change, report, Coverage(change))
    for sha in reversed(commits):
        subject = git_text("log", "-1", "--format=%h %s", sha)
        report.extend(commit_report(sha), prefix=f"[{subject[:50]}] ")
    if run_tests:
        with tempfile.TemporaryDirectory(prefix="repo_check_verify_") as tree:
            export_tree(head, tree)
            run_suites(SUITES, tree, report, ref=head)
    return commits


def mode_base(args):
    report = Report()
    if args.verify:
        commits = verify_range(args.base, args.head, report, run_tests=not args.no_tests)
        print(f"verified {len(commits)} commit(s) in {args.base[:12]}..{args.head}")
    else:
        change = Change("range", base=args.base, head=args.head)
        range_checks(change, report, Coverage(change))
        if not args.no_tests:
            with tempfile.TemporaryDirectory(prefix="repo_check_") as tree:
                export_tree(args.head, tree)
                run_suites(suites_touching(change.touched()), tree, report, ref=args.head)
    report.print()
    return 1 if report.errors else 0


def mode_push(args):
    report = Report()
    total = 0
    for line in sys.stdin.read().splitlines():
        parts = line.split()
        if len(parts) != 4:
            continue
        local_ref, local_sha, remote_ref, remote_sha = parts
        if local_sha == ZERO_SHA:
            continue
        base = remote_sha if remote_sha != ZERO_SHA and rev_exists(remote_sha) else None
        print(f"verification run: {local_ref} -> {remote_ref}", flush=True)
        total += len(verify_range(base, local_sha, report, run_tests=not args.no_tests))
    report.print(f"verification run ({total} commit(s))")
    if report.errors:
        print("push refused: fix the errors above (git push --no-verify skips this; CI still runs it)")
    return 1 if report.errors else 0


def mode_all(args):
    change = Change("all", ref=args.ref)
    report = Report()
    coverage = Coverage(change)
    range_checks(change, report, coverage)
    targets = [p for p in change.present() if is_code(p) and not is_verification(p)]
    check_readme_verification(change, report, coverage, targets)
    check_last_verified(change, report, sorted({d for d in (doc_for(t, change.tracked) for t in targets) if d}))
    if args.json:
        print(json.dumps(sorted(set(report.errors)) + sorted(set(report.warnings)), indent=1))
        return 0
    report.print()
    return 1 if report.errors else 0


def mode_record(args):
    report = Report()
    kind = "slow" if args.slow else "fast"
    with tempfile.TemporaryDirectory(prefix="repo_check_") as tree:
        export_tree(":", tree)
        results = run_suites(SUITES, tree, report, kind=kind)
    tracked = set(git("ls-files").decode().split("\n"))
    for doc in sorted(p for p in tracked if p.endswith("README.md")):
        suite = suite_for(doc)
        if suite and results.get(suite["name"]) and LAST_VERIFIED in (ROOT / doc).read_text():
            write_block(ROOT / doc, kind, results[suite["name"]])
            print(f"recorded {kind} run in {doc}")
    report.print()
    return 1 if report.errors else 0


def mode_plan_covers(args):
    target = args.plan_covers
    if target == "-":
        try:
            payload = json.load(sys.stdin)
        except ValueError:
            return 0
        tool_input = payload.get("tool_input") or {}
        target = tool_input.get("file_path") or tool_input.get("notebook_path")
        if not target:
            return 0
    path = Path(target)
    if not path.is_absolute():
        path = Path.cwd() / path
    try:
        rel = path.resolve().relative_to(ROOT.resolve()).as_posix()
    except ValueError:
        return 0
    if not is_code(rel):
        return 0
    plans = sorted((ROOT / PLAN_DIR).glob("*.md"))
    for plan in plans:
        if plan.name.startswith("_"):
            continue
        state, _, scope = parse_plan(plan.read_text())
        if state.startswith(PLAN_STATES) and in_scope(rel, scope):
            return 0
    sys.stderr.write(
        f"{rel} is code, and no approved plan in plans/ covers it. Write or extend a plan "
        f"(plans/_TEMPLATE.md) with this path in its ## Scope and a verification plan, and get it "
        f"approved before writing this code.\n"
    )
    return 2


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    mode = ap.add_mutually_exclusive_group()
    mode.add_argument("--commit-msg", metavar="FILE", help="plan and pairing checks for a commit message")
    mode.add_argument("--commit", metavar="SHA", help="the per-commit checks for one commit")
    mode.add_argument("--base", metavar="REF", help="check REF...HEAD (CI)")
    mode.add_argument("--push", action="store_true", help="pre-push verification run (refs on stdin)")
    mode.add_argument("--all", action="store_true", help="every tracked file: existing debt")
    mode.add_argument("--record", action="store_true", help="run suites and rewrite README Last verified blocks")
    mode.add_argument("--plan-covers", metavar="PATH", help="exit 2 unless an approved plan covers PATH")
    ap.add_argument("--head", default="HEAD", help="with --base: the head of the range")
    ap.add_argument("--verify", action="store_true", help="with --base: per-commit checks and every fast suite")
    ap.add_argument("--ref", help="with --all: a commit instead of the index")
    ap.add_argument("--json", action="store_true", help="with --all: findings as JSON")
    ap.add_argument("--slow", action="store_true", help="with --record: the slow suites")
    ap.add_argument("--no-tests", action="store_true", help="skip the test suites")
    args = ap.parse_args()

    if args.plan_covers:
        return mode_plan_covers(args)
    if args.commit_msg:
        return mode_commit_msg(args)
    if args.commit:
        report = commit_report(args.commit)
        report.print()
        return 1 if report.errors else 0
    if args.base:
        return mode_base(args)
    if args.push:
        return mode_push(args)
    if args.all:
        return mode_all(args)
    if args.record:
        return mode_record(args)
    return mode_pre_commit(args)


if __name__ == "__main__":
    sys.exit(main())
