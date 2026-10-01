"""Find projects, runs, sub-runs and analysis directories. Read-only.

Parallel filesystems are slow on metadata, so everything goes through
os.scandir with symlinks never followed, and each run tree is walked once.
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

from .classify import classify
from .config import Config

PROJECT_DIRS = ("runs", "analysis", "potentials", "structures")

# A directory holding any of these is a calculation (a "calc dir").
CALC_MARKERS = {"log.lammps", "OUTCAR", "OSZICAR", "vasprun.xml", "STREAM_OUTPUT", "INCAR"}
# Directories whose contents are inputs/staging, never calculations.
NON_CALC_DIRS = {"input_files", "setup"}


@dataclass
class FileEntry:
    relpath: str        # relative to the owning run/analysis dir, '/'-separated
    size: int
    mtime: float
    category: str
    symlink: Optional[str] = None   # link target, when the entry is a symlink
    broken: bool = False            # symlink whose target does not exist


@dataclass
class Project:
    name: str
    path: Path
    present: Dict[str, bool]


@dataclass
class RunDir:
    project: str
    run_id: str
    label: str           # directory name, e.g. "vasp-interface-0252"
    group_path: str      # dirs between runs/ and the run dir, e.g. "dft_surface_coverage"
    path: Path
    files: List[FileEntry] = field(default_factory=list)
    truncated: bool = False  # walk hit max depth somewhere

    @property
    def run_key(self) -> str:
        return f"{self.project}/{self.run_id}"


@dataclass
class AnalysisDir:
    project: str
    name: str
    path: Path
    run_ids: List[str]
    hint: str            # remainder of the name after the id part, e.g. "T60C_msd"
    files: List[FileEntry] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Walking
# ---------------------------------------------------------------------------

def walk_files(root: Path, cfg: Config) -> Tuple[List[FileEntry], bool]:
    """All files under root (symlinks recorded, never followed)."""
    out: List[FileEntry] = []
    truncated = False
    stack: List[Tuple[str, int]] = [(str(root), 0)]
    root_s = str(root)
    while stack:
        d, depth = stack.pop()
        try:
            it = os.scandir(d)
        except OSError:
            continue
        with it:
            for e in it:
                if e.name in cfg.exclude:
                    continue
                rel = os.path.relpath(e.path, root_s).replace(os.sep, "/")
                try:
                    if e.is_symlink():
                        target = os.readlink(e.path)
                        broken = not os.path.exists(e.path)
                        st = e.stat(follow_symlinks=False)
                        out.append(FileEntry(rel, 0, st.st_mtime, classify(e.name), target, broken))
                    elif e.is_dir(follow_symlinks=False):
                        if depth + 1 >= cfg.max_walk_depth:
                            truncated = True
                        else:
                            stack.append((e.path, depth + 1))
                    else:
                        st = e.stat(follow_symlinks=False)
                        out.append(FileEntry(rel, st.st_size, st.st_mtime, classify(e.name)))
                except OSError:
                    continue
    out.sort(key=lambda f: f.relpath)
    return out, truncated


def _subdirs(path: Path, cfg: Config) -> Iterator[os.DirEntry]:
    try:
        with os.scandir(path) as it:
            entries = [e for e in it if e.name not in cfg.exclude]
    except OSError:
        return
    for e in sorted(entries, key=lambda e: e.name):
        try:
            if e.is_dir(follow_symlinks=False):
                yield e
        except OSError:
            continue


# ---------------------------------------------------------------------------
# Projects and runs
# ---------------------------------------------------------------------------

def find_projects(root: Path, cfg: Config, max_depth: int = 3) -> List[Project]:
    """Directories containing runs/, searched from root downward (root included)."""
    projects: List[Project] = []

    def visit(path: Path, depth: int) -> None:
        if (path / "runs").is_dir():
            present = {d: (path / d).is_dir() for d in PROJECT_DIRS}
            projects.append(Project(path.name, path, present))
            return
        if depth >= max_depth:
            return
        for e in _subdirs(path, cfg):
            visit(Path(e.path), depth + 1)

    visit(root, 0)
    return projects


def find_runs(project: Project, cfg: Config, warnings: List[Tuple[str, str]],
              max_group_depth: int = 4) -> List[RunDir]:
    """Run dirs under <project>/runs; intermediate dirs become group_path."""
    runs: List[RunDir] = []

    def visit(path: Path, groups: List[str]) -> int:
        found = 0
        for e in _subdirs(path, cfg):
            sub = Path(e.path)
            run_id = cfg.match_id(e.name)
            if run_id is not None:
                runs.append(RunDir(project.name, run_id, e.name, "/".join(groups), sub))
                found += 1
            elif len(groups) < max_group_depth:
                n = visit(sub, groups + [e.name])
                if n == 0:
                    warnings.append((str(sub), "directory under runs/ contains no run ids"))
                found += n
        return found

    visit(project.path / "runs", [])
    return runs


# ---------------------------------------------------------------------------
# Calc dirs (sub-runs) inside a run
# ---------------------------------------------------------------------------

def calc_dirs(files: List[FileEntry]) -> Dict[str, List[FileEntry]]:
    """Group files by directory and keep the dirs that look like calculations.

    Keys are directory relpaths ('' = run dir itself)."""
    by_dir: Dict[str, List[FileEntry]] = {}
    for f in files:
        d = f.relpath.rsplit("/", 1)[0] if "/" in f.relpath else ""
        by_dir.setdefault(d, []).append(f)
    calcs = {}
    for d, fs in by_dir.items():
        if any(part in NON_CALC_DIRS for part in d.split("/") if part):
            continue
        names = {f.relpath.rsplit("/", 1)[-1] for f in fs}
        if names & CALC_MARKERS or any(f.category == "lammps_log" for f in fs):
            calcs[d] = fs
    return calcs


def group_frames(calcs) -> Tuple[List[str], Dict[str, List[str]]]:
    """Split calc dirs into single calcs and frame groups.

    Numbered sibling calc dirs (e.g. NEB images or path frames
    `size_03__ring_00008__roll_060/1 .. /29`) form one group keyed by their
    parent dir, ordered by number."""
    by_parent: Dict[str, List[str]] = {}
    for c in calcs:
        parent, _, last = c.rpartition("/")
        if last.isdigit():
            by_parent.setdefault(parent, []).append(c)
    groups = {p: sorted(cs, key=lambda c: int(c.rpartition("/")[2]))
              for p, cs in by_parent.items() if len(cs) >= 2}
    grouped = {c for cs in groups.values() for c in cs}
    singles = sorted(c for c in calcs if c not in grouped)
    return singles, groups


def subrun_label(calc_dir: str) -> str:
    """'T45C/run' -> 'T45C', 'run/output_atom1' -> 'output_atom1', 'run' -> ''."""
    return "/".join(p for p in calc_dir.split("/") if p and p != "run")


# ---------------------------------------------------------------------------
# Analysis dirs
# ---------------------------------------------------------------------------

_ID_PREFIX_RE = re.compile(r"^\d+(?:-\d+)?(?:_\d+(?:-\d+)?)*")


def parse_analysis_name(name: str, width: int) -> Tuple[List[str], str, List[str]]:
    """Run ids named by an analysis dir, plus the descriptive remainder.

    '0239_T60C_msd'           -> ['0239'], 'T60C_msd'
    '0281-ring-selection'     -> ['0281'], 'ring-selection'
    '0242_0243_compare'       -> ['0242', '0243'], 'compare'
    '0242-0245_rdf'           -> ['0242' .. '0245'], 'rdf'
    '025_6-8_reaction_energy' -> ['0256', '0257', '0258'], 'reaction_energy'
    """
    warnings: List[str] = []
    m = _ID_PREFIX_RE.match(name)
    if not m:
        return [], name, ["analysis dir name does not start with a run id"]
    rest = name[m.end():].lstrip("_-")
    tokens = m.group(0).split("_")
    ids: List[str] = []
    prefix: Optional[str] = None
    for i, tok in enumerate(tokens):
        a, _, b = tok.partition("-")
        if len(a) < width:
            if prefix is None and not b and i + 1 < len(tokens):
                prefix = a          # '025' in '025_6-8'
                continue
            if prefix is None or len(prefix) + len(a) != width:
                warnings.append(f"cannot expand id token {tok!r}")
                continue
            a = prefix + a
        elif len(a) > width:
            warnings.append(f"id token {tok!r} longer than {width} digits")
            continue
        if b:
            if len(b) < width:
                b = a[: width - len(b)] + b
            lo, hi = int(a), int(b)
            if hi < lo or hi - lo > 200:
                warnings.append(f"implausible id range {tok!r}")
                continue
            ids.extend(str(n).zfill(width) for n in range(lo, hi + 1))
        else:
            ids.append(a)
    return ids, rest, warnings


def find_analysis(project: Project, cfg: Config) -> List[AnalysisDir]:
    adir = project.path / "analysis"
    if not adir.is_dir():
        return []
    out = []
    for e in _subdirs(adir, cfg):
        ids, hint, warns = parse_analysis_name(e.name, cfg.id_width)
        out.append(AnalysisDir(project.name, e.name, Path(e.path), ids, hint, warnings=warns))
    return out


# ---------------------------------------------------------------------------
# Resolving an arbitrary path (e.g. from a hook) to project/run/subrun
# ---------------------------------------------------------------------------

def locate(path: str, cfg: Config) -> Optional[Dict[str, str]]:
    """Find project, run id and sub-run label from a path's components.

    Works on paths recorded before data was moved (e.g. /scratch1 -> /scratch2),
    since it only looks at '<project>/runs/.../<id>/...'."""
    parts = [p for p in str(path).replace(os.sep, "/").split("/") if p]
    for i in range(len(parts) - 1, -1, -1):
        if parts[i] != "runs" or i == 0:
            continue
        for j in range(i + 1, len(parts)):
            run_id = cfg.match_id(parts[j])
            if run_id is not None:
                return {
                    "project": parts[i - 1],
                    "run_id": run_id,
                    "group_path": "/".join(parts[i + 1: j]),
                    "subrun": subrun_label("/".join(parts[j + 1:])),
                }
        break
    return None
