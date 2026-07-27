#!/usr/bin/env python3
"""Ledger for completed MD/DFT SLURM runs, identified by a run id.

A run's id is the name of the directory two levels above its working/
output directory (the directory SLURM writes logs/dumps into):

    <id>/<subdir>/<run_dir>/...   ->   run_id = <id>

Ledger entries live as one JSON file per run id under runs/, plus one
JSON file per run id tracking analyses performed against it under
analysis/. Both directories sit alongside this script.

Subcommands:
    record        Record a finished run, reading target_dir + its slurm script.
    check         Read-only: report whether target_dir is already ledgered.
    resolve       Merge a conflicting run-id entry back into the ledger.
    search        Filter/list recorded runs.
    add-analysis  Append an analysis record for a run.
    show          Pretty-print a run's ledger entry and its analyses.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime
from pathlib import Path
from typing import Optional

LEDGER_ROOT = Path(__file__).resolve().parent
RUNS_DIR = LEDGER_ROOT / "runs"
ANALYSIS_DIR = LEDGER_ROOT / "analysis"

NESTED_KEYS = ("slurm_directives", "specifications")

SBATCH_LINE_RE = re.compile(r"^#SBATCH\s+(\S.*)$")


# ---------------------------------------------------------------------------
# Small helpers
# ---------------------------------------------------------------------------

def _now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _load_json(path: Path) -> dict:
    with path.open() as f:
        return json.load(f)


def _save_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        json.dump(data, f, indent=2, sort_keys=False)
        f.write("\n")


def _run_file(run_id: str) -> Path:
    return RUNS_DIR / f"{run_id}.json"


def _analysis_file(run_id: str) -> Path:
    return ANALYSIS_DIR / f"{run_id}_analysis.json"


def resolve_run_id(target_dir: Path) -> str:
    """Run id = the directory name two levels above target_dir."""
    resolved = target_dir.resolve()
    parents = resolved.parents
    if len(parents) < 2:
        raise ValueError(
            f"{resolved} is too shallow to have a grandparent directory; "
            "pass --run-id explicitly."
        )
    return parents[1].name


def find_slurm_script(target_dir: Path) -> Optional[Path]:
    candidates = sorted(target_dir.glob("*.slurm")) + sorted(target_dir.glob("*.pbs"))
    if len(candidates) == 1:
        return candidates[0]
    return None


def parse_slurm_directives(script_path: Path) -> dict:
    directives: dict = {}
    text = script_path.read_text(errors="replace")
    for line in text.splitlines():
        m = SBATCH_LINE_RE.match(line.strip())
        if not m:
            continue
        rest = m.group(1).strip()
        if rest.startswith("--"):
            token = rest[2:]
            if "=" in token:
                key, _, val = token.partition("=")
            else:
                parts = token.split(None, 1)
                key, val = parts[0], (parts[1] if len(parts) > 1 else "")
        elif rest.startswith("-"):
            parts = rest[1:].split(None, 1)
            key, val = parts[0], (parts[1] if len(parts) > 1 else "")
        else:
            continue
        directives[key] = val
    return directives


def parse_fields(field_args: list[str]) -> dict:
    fields: dict = {}
    for item in field_args or []:
        if "=" not in item:
            raise ValueError(f"--field expects key=value, got: {item!r}")
        key, _, val = item.partition("=")
        fields[key] = val
    return fields


def build_entry(
    run_id: str,
    target_dir: Path,
    slurm_script: Optional[Path],
    status: str,
    notes: str,
    run_type: Optional[str],
    fields: dict,
) -> dict:
    return {
        "run_id": run_id,
        "recorded_at": _now(),
        "run_directory": str(target_dir),
        "run_type": run_type,
        "slurm_script": str(slurm_script) if slurm_script else None,
        "slurm_directives": parse_slurm_directives(slurm_script) if slurm_script else {},
        "status": status,
        "notes": notes or "",
        "specifications": fields,
    }


# ---------------------------------------------------------------------------
# record
# ---------------------------------------------------------------------------

def cmd_record(args: argparse.Namespace) -> int:
    target_dir = Path(args.target_dir).resolve()
    if not target_dir.is_dir():
        print(f"ERROR: {target_dir} is not a directory", file=sys.stderr)
        return 1

    run_id = args.run_id or resolve_run_id(target_dir)

    if args.slurm_script:
        slurm_script = Path(args.slurm_script).resolve()
    else:
        slurm_script = find_slurm_script(target_dir)
        if slurm_script is None:
            candidates = sorted(target_dir.glob("*.slurm")) + sorted(target_dir.glob("*.pbs"))
            if candidates:
                print(
                    f"WARNING: multiple slurm scripts found in {target_dir}; "
                    "pass --slurm-script explicitly. Recording with no slurm_directives.",
                    file=sys.stderr,
                )
            else:
                print(
                    f"WARNING: no *.slurm/*.pbs found in {target_dir}; "
                    "recording with no slurm_directives.",
                    file=sys.stderr,
                )

    try:
        fields = parse_fields(args.field)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    entry = build_entry(
        run_id, target_dir, slurm_script, args.status, args.notes, args.run_type, fields
    )

    run_file = _run_file(run_id)
    if not run_file.exists():
        _save_json(run_file, entry)
        print(f"Recorded run {run_id} -> {run_file}")
        return 0

    existing = _load_json(run_file)
    if existing.get("run_directory") == str(target_dir):
        _save_json(run_file, entry)
        print(f"Run {run_id} already recorded for this directory — updated existing entry.")
        return 0

    # Genuine id collision: two different directories claim the same run id.
    print(
        f"WARNING: run id {run_id} is already ledgered for a different directory "
        f"({existing.get('run_directory')}). Not overwriting {run_file}.",
        file=sys.stderr,
    )
    conflict_file = target_dir / f"ledger_entry.{run_id}.conflict.json"
    _save_json(conflict_file, entry)
    print(f"Wrote conflicting entry to {conflict_file}", file=sys.stderr)

    existing["conflict"] = {"detected_at": _now(), "incoming_file": str(conflict_file)}
    _save_json(run_file, existing)

    print(f"Run: ledger.py resolve {run_id} --incoming {conflict_file}", file=sys.stderr)
    return 1


# ---------------------------------------------------------------------------
# check
# ---------------------------------------------------------------------------

def cmd_check(args: argparse.Namespace) -> int:
    target_dir = Path(args.target_dir).resolve()
    run_id = args.run_id or resolve_run_id(target_dir)
    run_file = _run_file(run_id)

    if not run_file.exists():
        print(f"run_id={run_id}: no ledger entry yet.")
        return 0

    existing = _load_json(run_file)
    existing_dir = existing.get("run_directory")
    if existing_dir == str(target_dir):
        print(f"run_id={run_id}: exact match — already recorded for this directory.")
        return 0

    print(
        f"run_id={run_id}: CONFLICT — already recorded for a different directory:\n"
        f"  existing: {existing_dir}\n"
        f"  this dir: {target_dir}"
    )
    return 1


# ---------------------------------------------------------------------------
# resolve
# ---------------------------------------------------------------------------

def _prompt_choice(field_name: str, existing_val, incoming_val):
    print(f"\nField '{field_name}' differs:")
    print(f"  [1] existing: {existing_val!r}")
    print(f"  [2] incoming: {incoming_val!r}")
    print("  [3] enter a custom value")
    while True:
        choice = input("Choice [1/2/3]: ").strip()
        if choice == "1":
            return existing_val
        if choice == "2":
            return incoming_val
        if choice == "3":
            return input("Custom value: ")
        print("Please enter 1, 2, or 3.")


def _merge_dict(existing: dict, incoming: dict, prefix: str = "") -> dict:
    merged = {}
    for key in sorted(set(existing) | set(incoming)):
        ev = existing.get(key)
        iv = incoming.get(key)
        label = f"{prefix}{key}"
        if key in NESTED_KEYS and isinstance(ev, dict) and isinstance(iv, dict):
            merged[key] = _merge_dict(ev, iv, prefix=f"{key}.")
        elif ev == iv:
            merged[key] = ev
        else:
            merged[key] = _prompt_choice(label, ev, iv)
    return merged


def cmd_resolve(args: argparse.Namespace) -> int:
    run_id = args.run_id
    run_file = _run_file(run_id)
    if not run_file.exists():
        print(f"ERROR: no ledger entry for run {run_id}; nothing to resolve.", file=sys.stderr)
        return 1

    incoming_path = Path(args.incoming).resolve()
    if not incoming_path.exists():
        print(f"ERROR: incoming file {incoming_path} does not exist.", file=sys.stderr)
        return 1

    existing = _load_json(run_file)
    incoming = _load_json(incoming_path)

    backup_file = run_file.with_suffix(run_file.suffix + ".bak")
    _save_json(backup_file, existing)

    existing_data = {k: v for k, v in existing.items() if k != "conflict"}
    incoming_data = {k: v for k, v in incoming.items() if k != "conflict"}
    merged = _merge_dict(existing_data, incoming_data)

    _save_json(run_file, merged)

    resolved_path = incoming_path.with_name(incoming_path.name + ".resolved")
    incoming_path.rename(resolved_path)

    print(f"Merged entry written to {run_file} (backup at {backup_file}).")
    print(f"Incoming conflict file archived as {resolved_path}.")
    return 0


# ---------------------------------------------------------------------------
# search
# ---------------------------------------------------------------------------

def _lookup(entry: dict, key: str):
    if key in entry:
        return entry[key]
    for nested_key in NESTED_KEYS:
        nested = entry.get(nested_key) or {}
        if key in nested:
            return nested[key]
    return None


def cmd_search(args: argparse.Namespace) -> int:
    try:
        wanted_fields = parse_fields(args.field)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    if not RUNS_DIR.is_dir():
        print("No runs recorded yet.")
        return 0

    matches = []
    for run_file in sorted(RUNS_DIR.glob("*.json")):
        entry = _load_json(run_file)

        if args.status and entry.get("status") != args.status:
            continue
        if args.conflicts_only and "conflict" not in entry:
            continue
        if wanted_fields and any(str(_lookup(entry, k)) != v for k, v in wanted_fields.items()):
            continue
        if args.contains:
            haystack = f"{entry.get('run_id', '')} {entry.get('notes', '')}".lower()
            if args.contains.lower() not in haystack:
                continue

        matches.append(entry)

    if not matches:
        print("No matching runs.")
        return 0

    for entry in matches:
        flag = " [CONFLICT]" if "conflict" in entry else ""
        print(
            f"{entry.get('run_id')}{flag}  status={entry.get('status')}  "
            f"run_type={entry.get('run_type')}  dir={entry.get('run_directory')}"
        )
    return 0


# ---------------------------------------------------------------------------
# add-analysis
# ---------------------------------------------------------------------------

def cmd_add_analysis(args: argparse.Namespace) -> int:
    if args.dir:
        target_dir = Path(args.dir).resolve()
        run_id = args.run_id or resolve_run_id(target_dir)
    else:
        run_id = args.run_id

    run_file = _run_file(run_id)
    if not run_file.exists():
        print(
            f"WARNING: run {run_id} has no ledger entry yet (recording analysis anyway).",
            file=sys.stderr,
        )

    try:
        fields = parse_fields(args.field)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    analysis_file = _analysis_file(run_id)
    data = _load_json(analysis_file) if analysis_file.exists() else {"run_id": run_id, "analyses": []}

    data["analyses"].append(
        {
            "type": args.analysis_type,
            "recorded_at": _now(),
            "output": args.output,
            "notes": args.notes or "",
            "fields": fields,
        }
    )
    _save_json(analysis_file, data)
    print(f"Recorded '{args.analysis_type}' analysis for run {run_id} -> {analysis_file}")
    return 0


# ---------------------------------------------------------------------------
# show
# ---------------------------------------------------------------------------

def cmd_show(args: argparse.Namespace) -> int:
    run_file = _run_file(args.run_id)
    if run_file.exists():
        print(json.dumps(_load_json(run_file), indent=2))
    else:
        print(f"No ledger entry for run {args.run_id}.")

    analysis_file = _analysis_file(args.run_id)
    if analysis_file.exists():
        print()
        print(json.dumps(_load_json(analysis_file), indent=2))
    return 0


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Ledger for completed MD/DFT SLURM runs.")
    sub = p.add_subparsers(dest="command", required=True)

    p_record = sub.add_parser("record", help="Record a finished run.")
    p_record.add_argument("target_dir", help="Run's working/output directory.")
    p_record.add_argument("--run-id", help="Override the auto-detected run id.")
    p_record.add_argument("--slurm-script", help="Explicit path to the submitted slurm script.")
    p_record.add_argument("--status", default="completed")
    p_record.add_argument("--notes", default="")
    p_record.add_argument("--run-type", choices=["job", "setup"], default=None)
    p_record.add_argument("--field", action="append", default=[], help="key=value, repeatable.")
    p_record.set_defaults(func=cmd_record)

    p_check = sub.add_parser("check", help="Report whether target_dir is already ledgered.")
    p_check.add_argument("target_dir")
    p_check.add_argument("--run-id", help="Override the auto-detected run id.")
    p_check.set_defaults(func=cmd_check)

    p_resolve = sub.add_parser("resolve", help="Merge a conflicting entry back into the ledger.")
    p_resolve.add_argument("run_id")
    p_resolve.add_argument("--incoming", required=True, help="Path to the *.conflict.json file.")
    p_resolve.set_defaults(func=cmd_resolve)

    p_search = sub.add_parser("search", help="Filter/list recorded runs.")
    p_search.add_argument("--field", action="append", default=[], help="key=value, repeatable.")
    p_search.add_argument("--contains", help="Substring match against run_id/notes.")
    p_search.add_argument("--status")
    p_search.add_argument("--conflicts-only", action="store_true")
    p_search.set_defaults(func=cmd_search)

    p_analysis = sub.add_parser("add-analysis", help="Append an analysis record for a run.")
    p_analysis.add_argument("run_id", nargs="?", help="Run id (omit if using --dir).")
    p_analysis.add_argument("analysis_type")
    p_analysis.add_argument("--dir", help="Analyzed run's data directory (resolves the run id).")
    p_analysis.add_argument("--output", help="Path to the analysis output.")
    p_analysis.add_argument("--notes", default="")
    p_analysis.add_argument("--field", action="append", default=[], help="key=value, repeatable.")
    p_analysis.set_defaults(func=cmd_add_analysis)

    p_show = sub.add_parser("show", help="Pretty-print a run's ledger entry and analyses.")
    p_show.add_argument("run_id")
    p_show.set_defaults(func=cmd_show)

    return p


def main(argv: Optional[list[str]] = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if args.command == "add-analysis" and bool(args.run_id) == bool(args.dir):
        parser.error("add-analysis requires exactly one of a positional run_id or --dir")

    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
