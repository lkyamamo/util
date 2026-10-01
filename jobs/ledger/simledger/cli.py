"""Command-line interface. Run `simledger -h` or `simledger <command> -h`."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

from . import __version__, config, export, inbox, search, store
from .discover import locate


def _fields(items: List[str]) -> dict:
    out = {}
    for item in items or []:
        key, sep, val = item.partition("=")
        if not sep:
            raise SystemExit(f"--field expects key=value, got {item!r}")
        out[key.strip()] = val.strip()
    return out


def _conn(cfg, readonly=True):
    if not cfg.db_path.exists():
        raise SystemExit(f"No ledger at {cfg.db_path}. Run `simledger scan <root>` first "
                         f"(or set SIMLEDGER_HOME / --home).")
    return store.connect(cfg.db_path, readonly=readonly)


def resolve_run(conn, cfg, ref: str) -> str:
    """Run key from 'project/0247', '0247', or a path inside a run."""
    if conn.execute("SELECT 1 FROM runs WHERE run_key=?", (ref,)).fetchone():
        return ref
    if "/" in ref or Path(ref).exists():
        loc = locate(str(Path(ref).absolute()) if Path(ref).exists() else ref, cfg)
        if loc:
            key = f"{loc['project']}/{loc['run_id']}"
            if conn.execute("SELECT 1 FROM runs WHERE run_key=?", (key,)).fetchone():
                return key
            ref = loc["run_id"]
    rid = cfg.match_id(ref) or ref
    rows = conn.execute("SELECT run_key FROM runs WHERE run_id=? OR label=?", (rid, ref)).fetchall()
    if len(rows) == 1:
        return rows[0][0]
    if not rows:
        raise SystemExit(f"No run matches {ref!r}")
    raise SystemExit(f"{ref!r} is ambiguous: " + ", ".join(r[0] for r in rows))


# ---------------------------------------------------------------------------
# Commands
# ---------------------------------------------------------------------------

def cmd_scan(args, cfg) -> int:
    from .scan import scan
    from .survey import survey
    if args.dry_run:
        survey(cfg, Path(args.root))
        return 0
    result = scan(cfg, Path(args.root), full=args.full, log=lambda m: print(m, file=sys.stderr))
    conn = result["conn"]
    if not args.no_export:
        export.export_all(conn, cfg, result)
    c = result["counts"]
    print(f"scan {result['scan_id']}: {c['projects']} projects, {c['runs']} runs "
          f"({c['new']} new, {c['changed']} changed, {c['unchanged']} unchanged, {c['missing']} missing), "
          f"{c['analysis']} analysis dirs, {c['events']} hook events, {c['warnings']} warnings")
    print(f"ledger: {cfg.db_path}\ncards:  {cfg.cards}" + ("" if args.no_export else f"  (report: {cfg.cards / 'last_scan.md'})"))
    return 0


def cmd_survey(args, cfg) -> int:
    from .survey import survey
    survey(cfg, Path(args.root))
    return 0


def _print_runs(rows, analysis_hits=()) -> None:
    if not rows and not analysis_hits:
        print("No matching runs.")
        return
    for r in rows:
        print(f"{r['run_key']:<34} {r['status']:<12} {r['code']:<7} {(r['calc_type'] or ''):<11} "
              f"{(r['group_path'] or ''):<22} {r['label']}" + ("  [missing]" if r["missing"] else ""))
    if rows:
        print(f"({len(rows)} runs)")
    for a in analysis_hits:
        print(f"analysis: {a}")


def cmd_search(args, cfg) -> int:
    conn = _conn(cfg)
    rows, hits = search.run(conn, " ".join(args.query), args.text or "", args.all)
    _print_runs(rows, hits)
    return 0


def cmd_list(args, cfg) -> int:
    q = []
    if args.project:
        q.append(f"project={args.project}")
    if args.status:
        q.append(f"status={args.status}")
    if args.code:
        q.append(f"code={args.code}")
    if args.group:
        q.append(f"group~{args.group}")
    rows, _ = search.run(_conn(cfg), " ".join(q), include_missing=args.all)
    _print_runs(rows)
    return 0


def cmd_show(args, cfg) -> int:
    conn = _conn(cfg)
    key = resolve_run(conn, cfg, args.run)
    print(export.render_card(conn, key), end="")
    if args.files:
        for f in conn.execute("SELECT relpath, category, size, symlink, broken FROM files WHERE run_key=? "
                              "ORDER BY relpath", (key,)):
            link = f"  -> {f['symlink']}{' (broken)' if f['broken'] else ''}" if f["symlink"] else ""
            print(f"{f['category']:<16} {export.human_bytes(f['size']):>10}  {f['relpath']}{link}")
    return 0


def cmd_note(args, cfg) -> int:
    conn = _conn(cfg, readonly=False)
    key = resolve_run(conn, cfg, args.run)
    conn.execute("INSERT INTO user_notes (run_key, note, tags, created) VALUES (?,?,?,?)",
                 (key, args.text, args.tag, store.now()))
    store.add_fts(conn, key, "note", f"{args.text} {args.tag or ''}")
    conn.commit()
    path = export.card_path(cfg, key)
    if path.parent.is_dir():
        export._write(path, export.render_card(conn, key))
    print(f"note added to {key}")
    return 0


def cmd_export(args, cfg) -> int:
    conn = _conn(cfg)
    if args.csv:
        export.write_csv(conn, Path(args.csv))
        print(f"wrote {args.csv}")
    else:
        export.export_all(conn, cfg)
        print(f"wrote cards to {cfg.cards}")
    return 0


def cmd_sql(args, cfg) -> int:
    conn = _conn(cfg)
    cur = conn.execute(args.query)
    cols = [d[0] for d in cur.description or []]
    if cols:
        print("\t".join(cols))
    for row in cur:
        print("\t".join("" if v is None else str(v) for v in row))
    return 0


def cmd_hook(args, cfg) -> int:
    """Never fails the job: problems are reported on stderr and exit status is 0."""
    try:
        fields = _fields(args.field)
        if args.hook_kind == "run":
            path = inbox.hook_run(cfg, args.dir, args.exit_code, args.script, fields, args.note)
        else:
            if not args.dir and not args.run:
                print("simledger hook analysis: give --dir or --run", file=sys.stderr)
                return 0
            path = inbox.hook_analysis(cfg, args.dir, args.run, args.type, args.output, fields, args.note)
        loc = locate(args.dir, cfg) if getattr(args, "dir", None) else None
        where = f"{loc['project']}/{loc['run_id']}" + (f" [{loc['subrun']}]" if loc and loc["subrun"] else "") \
            if loc else "(run not recognised from path; will retry at scan)"
        print(f"simledger: recorded {args.hook_kind} event for {where} -> {path}", file=sys.stderr)
    except Exception as exc:  # noqa: BLE001 - a ledger problem must never fail a job
        print(f"simledger: WARNING: hook failed: {exc}", file=sys.stderr)
    return 0


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="simledger", description="Read-only ledger of LAMMPS/VASP runs.")
    p.add_argument("--home", help="Ledger directory (default $SIMLEDGER_HOME or ~/ledger).")
    p.add_argument("--version", action="version", version=__version__)
    sub = p.add_subparsers(dest="command", required=True)

    s = sub.add_parser("scan", help="Scan a tree and update the ledger + cards.")
    s.add_argument("root")
    s.add_argument("--full", action="store_true", help="Reparse every run, not just changed ones.")
    s.add_argument("--dry-run", action="store_true", help="Survey only; write nothing.")
    s.add_argument("--no-export", action="store_true", help="Update the DB but not the Markdown cards.")
    s.set_defaults(func=cmd_scan)

    s = sub.add_parser("survey", help="Describe a tree's layout; no parsing, no DB writes.")
    s.add_argument("root")
    s.set_defaults(func=cmd_survey)

    s = sub.add_parser("search", help='Filter runs, e.g. "code=vasp ENCUT>=520" or free text.')
    s.add_argument("query", nargs="*", help="key=value, key>=n, key~sub, or bare words (full text).")
    s.add_argument("--text", help="Full-text terms.")
    s.add_argument("--all", action="store_true", help="Include runs missing from the latest scan.")
    s.set_defaults(func=cmd_search)

    s = sub.add_parser("list", help="List runs.")
    s.add_argument("--project")
    s.add_argument("--status")
    s.add_argument("--code")
    s.add_argument("--group")
    s.add_argument("--all", action="store_true")
    s.set_defaults(func=cmd_list)

    s = sub.add_parser("show", help="Print a run card (run id, project/id, or a path inside a run).")
    s.add_argument("run")
    s.add_argument("--files", action="store_true", help="Also list every file.")
    s.set_defaults(func=cmd_show)

    s = sub.add_parser("note", help="Attach a note to a run (kept across rescans).")
    s.add_argument("run")
    s.add_argument("text")
    s.add_argument("--tag", help="Comma-separated tags.")
    s.set_defaults(func=cmd_note)

    s = sub.add_parser("export", help="Rewrite the cards, or write runs CSV to a path.")
    s.add_argument("--csv", help="Write runs.csv to this path instead.")
    s.set_defaults(func=cmd_export)

    s = sub.add_parser("sql", help="Run a read-only SQL query against the ledger.")
    s.add_argument("query")
    s.set_defaults(func=cmd_sql)

    s = sub.add_parser("hook", help="Record a job-time event from a SLURM script (never fails the job).")
    hk = s.add_subparsers(dest="hook_kind", required=True)
    h = hk.add_parser("run", help="A simulation job finished.")
    h.add_argument("--dir", required=True, help="The job's run directory.")
    h.add_argument("--exit-code", type=int, help="Exit status of the main program.")
    h.add_argument("--script", help="Path to the job script, if `scontrol` cannot provide it.")
    h.add_argument("--field", action="append", default=[], help="key=value, repeatable.")
    h.add_argument("--note", default="")
    h = hk.add_parser("analysis", help="An analysis job finished.")
    h.add_argument("--dir", help="The analysed run's data directory.")
    h.add_argument("--run", help="Run id or project/run id, instead of --dir.")
    h.add_argument("--type", required=True, help="Analysis type, e.g. msd, dielectric.")
    h.add_argument("--output", help="Output file or directory.")
    h.add_argument("--field", action="append", default=[], help="key=value, repeatable.")
    h.add_argument("--note", default="")
    s.set_defaults(func=cmd_hook)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    cfg = config.load(Path(args.home) if args.home else None)
    return args.func(args, cfg)
