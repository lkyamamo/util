"""`simledger survey`: describe a tree without parsing files or touching the DB.

Used to check the layout assumptions (run-id names, grouping depth, sub-runs,
analysis names) against real data before trusting a scan."""

from __future__ import annotations

import re
from collections import Counter
from pathlib import Path
from typing import List, Tuple

from .config import Config
from .discover import (calc_dirs, find_analysis, find_projects, find_runs, group_frames,
                       subrun_label, walk_files)
from .export import human_bytes


def survey(cfg: Config, root: Path, out=print) -> None:
    root = root.resolve()
    projects = find_projects(root, cfg)
    out(f"# Survey of {root}\n")
    if not projects:
        out("No projects found (no directory containing runs/).")
        return
    for p in projects:
        warnings: List[Tuple[str, str]] = []
        runs = find_runs(p, cfg, warnings)
        out(f"## {p.name}  ({p.path})")
        out("   dirs present: " + ", ".join(f"{k}={'yes' if v else 'NO'}" for k, v in p.present.items()))
        shapes = Counter(re.sub(r"\d", "#", r.label) for r in runs)
        groups = Counter(r.group_path or "(none)" for r in runs)
        ids = Counter(r.run_id for r in runs)
        out(f"   runs: {len(runs)}  ids {min(ids) if ids else '-'}..{max(ids) if ids else '-'}")
        out("   name shapes: " + ", ".join(f"{s} x{n}" for s, n in shapes.most_common()))
        out("   groups: " + ", ".join(f"{g} x{n}" for g, n in groups.most_common()))
        dups = [i for i, n in ids.items() if n > 1]
        if dups:
            out(f"   DUPLICATE ids: {', '.join(dups)}")
        cats: Counter = Counter()
        sizes: Counter = Counter()
        biggest = []
        sub_counts = Counter()
        sub_examples = {}
        broken = 0
        for r in runs:
            files, truncated = walk_files(r.path, cfg)
            if truncated:
                warnings.append((str(r.path), "walk truncated at max depth"))
            for f in files:
                cats[f.category] += 1
                sizes[f.category] += f.size
                broken += f.broken
                biggest.append((f.size, f"{r.label}/{f.relpath}"))
            singles, frames = group_frames(calc_dirs(files))
            units = [subrun_label(c) or "(top)" for c in singles] + \
                    [f"{subrun_label(p) or p}[{len(fs)} frames]" for p, fs in sorted(frames.items())]
            sub_counts[len(units)] += 1
            if len(units) > 1 or frames:
                sub_examples[r.label] = units[:4] + (["..."] if len(units) > 4 else [])
        out("   sub-runs per run: " + ", ".join(f"{k} x{n}" for k, n in sorted(sub_counts.items())))
        for label, ex in list(sub_examples.items())[:6]:
            out(f"      e.g. {label}: {', '.join(ex)}")
        out("   files: " + ", ".join(f"{c} {n} ({human_bytes(sizes[c])})" for c, n in cats.most_common()))
        if broken:
            out(f"   broken symlinks: {broken}")
        biggest.sort(reverse=True)
        out("   largest: " + ", ".join(f"{n} ({human_bytes(s)})" for s, n in biggest[:5]))
        analysis = find_analysis(p, cfg)
        unparsed = [a.name for a in analysis if not a.run_ids or a.warnings]
        multi = [f"{a.name}->{','.join(a.run_ids)}" for a in analysis if len(a.run_ids) > 1]
        unknown = sorted({i for a in analysis for i in a.run_ids} - set(ids))
        out(f"   analysis dirs: {len(analysis)}")
        if multi:
            out("      multi-run: " + ", ".join(multi))
        if unparsed:
            out("      unparsed names: " + ", ".join(unparsed))
        if unknown:
            out("      ids with no run in this project: " + ", ".join(unknown))
        if warnings:
            out("   warnings:")
            for path, msg in warnings[:30]:
                out(f"      {path}: {msg}")
            if len(warnings) > 30:
                out(f"      ... {len(warnings) - 30} more")
        out("")
