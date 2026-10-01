"""Human-readable output regenerated from the database: run cards, project
summaries, INDEX.md, runs.csv and last_scan.md under <home>/cards."""

from __future__ import annotations

import csv
import json
import os
from collections import Counter, OrderedDict
from pathlib import Path
from typing import Dict, List, Optional

from .config import Config
from .describe import fmt_time

CSV_PARAMS = ["ensemble", "T_target", "P_target", "box", "density", "composition", "simulated_time",
              "units", "pair_style", "timestep", "total_steps",
              "ENCUT", "PREC", "ISMEAR", "IBRION", "ISIF", "NSW", "GGA", "METAGGA", "IVDW",
              "partition", "job_id", "exit_code"]
CSV_RESULTS = ["energy_sigma0", "energy_per_atom", "pressure", "ionic_steps", "steps_completed",
               "performance", "elapsed_s"]


def human_bytes(n: Optional[float]) -> str:
    n = float(n or 0)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1024 or unit == "TB":
            return f"{n:.0f} {unit}" if unit == "B" else f"{n:.1f} {unit}"
        n /= 1024
    return ""


def human_time(s: Optional[float]) -> str:
    if not s:
        return ""
    s = int(s)
    return f"{s // 3600}:{s % 3600 // 60:02d}:{s % 60:02d}"


def _cell(v) -> str:
    if v is None:
        return ""
    if isinstance(v, float):
        v = f"{v:.6g}"
    return str(v).replace("|", "\\|").replace("\n", " ")


def _table(headers: List[str], rows: List[List]) -> str:
    out = ["| " + " | ".join(headers) + " |", "|" + "---|" * len(headers)]
    out += ["| " + " | ".join(_cell(c) for c in r) + " |" for r in rows]
    return "\n".join(out)


def _yaml_val(v) -> str:
    if v is None:
        return "null"
    if isinstance(v, (int, float)):
        return str(v)
    s = str(v)
    return json.dumps(s) if (not s or any(c in s for c in ":#'\"[]{}\n,") or s != s.strip()) else s


def card_path(cfg: Config, run_key: str) -> Path:
    project, _, run_id = run_key.partition("/")
    return cfg.cards / project / f"{run_id.replace('/', '_')}.md"


# ---------------------------------------------------------------------------
# Run card
# ---------------------------------------------------------------------------

def render_card(conn, run_key: str) -> str:
    r = conn.execute("SELECT * FROM runs WHERE run_key=?", (run_key,)).fetchone()
    if r is None:
        return f"No run {run_key}\n"
    subs = conn.execute("SELECT * FROM subruns WHERE run_key=? ORDER BY relpath", (run_key,)).fetchall()
    params = conn.execute("SELECT * FROM params WHERE run_key=? ORDER BY source, key", (run_key,)).fetchall()
    results = conn.execute("SELECT * FROM results WHERE run_key=? ORDER BY subrun, key", (run_key,)).fetchall()
    notes = conn.execute("SELECT * FROM user_notes WHERE run_key=? ORDER BY note_id", (run_key,)).fetchall()
    events = conn.execute("SELECT * FROM events WHERE run_key=? ORDER BY created", (run_key,)).fetchall()
    analyses = conn.execute(
        "SELECT ar.*, a.path, a.n_files, a.readme FROM analysis_runs ar "
        "LEFT JOIN analysis a ON a.analysis_key = ar.analysis_key WHERE ar.run_key=? "
        "ORDER BY ar.analysis_key", (run_key,)).fetchall()

    front = OrderedDict((k, r[k]) for k in (
        "run_key", "project", "run_id", "label", "group_path", "code", "calc_type", "status",
        "n_atoms", "formula", "cores", "nodes", "job_name", "job_id", "start_time", "end_time",
        "wall_time_s", "n_subruns", "path", "last_changed", "missing"))
    system = json.loads(r["system"]) if r["system"] else {}
    for k in ("ensemble", "calc_type_text", "functional", "T_target", "P_target", "box", "density",
              "composition", "simulated_time"):
        if system.get(k) is not None:
            front[k] = system[k]
    if r["summary"]:
        front["summary"] = r["summary"]
    tags = sorted({n["tags"] for n in notes if n["tags"]})
    if tags:
        front["tags"] = ",".join(tags)
    out = ["---"] + [f"{k}: {_yaml_val(v)}" for k, v in front.items()] + ["---", ""]
    out.append(f"# {r['run_key']} — {r['label']}")
    out.append("")
    out.append(f"**{r['status']}** ({r['status_evidence'] or 'no evidence'}) · {r['code']}"
               + (f" · {r['calc_type']}" if r["calc_type"] else "")
               + (f" · group `{r['group_path']}`" if r["group_path"] else ""))
    out.append("")
    out.append(f"`{r['path']}`")
    if r["summary"]:
        out += ["", f"> {r['summary']}"]
    if r["missing"]:
        out += ["", "> **Missing:** this directory was not found in the latest scan."]
    out.append("")

    if notes:
        out.append("## Notes")
        for n in notes:
            out.append(f"- {n['created'][:10]}: {n['note']}" + (f" _(tags: {n['tags']})_" if n["tags"] else ""))
        out.append("")
    if r["readme"]:
        out += ["## README", "", r["readme"].strip(), ""]

    multi = subs and (len(subs) > 1 or subs[0]["relpath"] != "run")
    if multi:
        out.append(f"## Sub-runs ({len(subs)})")
        out.append("")
        out.append(_table(["label", "dir", "conditions", "status", "atoms", "wall", "evidence"],
                          [[s["label"], s["relpath"], s["conditions"], s["status"], s["n_atoms"],
                            human_time(s["wall_time_s"]), s["status_evidence"]] for s in subs]))
        out.append("")

    if system:
        ref = next((s for s in subs if s["system"] == r["system"]), None)
        title = "## System" + (f" (from sub-run {ref['label']})" if multi and ref and ref["label"] else "")
        out += [title, ""] + _system_table(system) + [""]
    protos = [(s["label"], json.loads(s["protocol"])) for s in subs if s["protocol"]]
    if protos:
        out += ["## Protocol", ""]
        unit = (system.get("time_unit") or "") if system else ""
        for label, proto in protos[:12]:
            if len(protos) > 1 or (label and label != "(top)"):
                out += [f"### {label or '(main)'}", ""]
            out += [_protocol_table(proto, unit), ""]
        if len(protos) > 12:
            out += [f"_{len(protos) - 12} more sub-runs; see `simledger sql` on subruns.protocol._", ""]

    out += _params_section([p for p in params if p["source"] != "system"])
    if results:
        out.append("## Results")
        out.append("")
        rows = [[x["subrun"] or "(main)", x["key"], x["value_text"] if x["value_num"] is None else x["value_num"],
                 x["unit"]] for x in results]
        out.append(_table(["sub-run", "key", "value", "unit"], rows[:300]))
        if len(rows) > 300:
            out.append(f"\n_{len(rows) - 300} more rows in the database._")
        out.append("")

    if events:
        out.append("## Jobs (from SLURM hooks)")
        out.append("")
        rows = []
        for e in events:
            d = json.loads(e["data"])
            sl = d.get("slurm") or {}
            rows.append([e["created"], e["kind"], e["subrun"] or "", sl.get("SLURM_JOB_ID"),
                         d.get("exit_code"), sl.get("SLURM_JOB_NODELIST"),
                         d.get("analysis_type") or d.get("note") or ""])
        out.append(_table(["time", "kind", "sub-run", "job", "exit", "nodes", "detail"], rows))
        last_script = next((json.loads(e["data"]).get("script_text") for e in reversed(events)
                            if json.loads(e["data"]).get("script_text")), None)
        if last_script:
            out += ["", "<details><summary>Last submitted script</summary>", "", "```bash",
                    last_script.rstrip(), "```", "", "</details>"]
        out.append("")

    if analyses:
        out.append("## Analysis")
        out.append("")
        rows = [[a["analysis_key"], a["subrun"] or "", a["source"], a["n_files"], a["path"]] for a in analyses]
        out.append(_table(["analysis", "sub-run", "linked by", "files", "path"], rows))
        figs = conn.execute(
            "SELECT analysis_key, relpath FROM analysis_files WHERE category='figure' AND analysis_key IN "
            f"({','.join('?' * len(analyses))}) ORDER BY analysis_key, relpath LIMIT 40",
            [a["analysis_key"] for a in analyses]).fetchall()
        if figs:
            out += ["", "Figures: " + ", ".join(f"`{f['analysis_key'].split('/', 1)[1]}/{f['relpath']}`" for f in figs)]
        out.append("")

    issues = conn.execute("SELECT path, message FROM run_warnings WHERE run_key=?", (run_key,)).fetchall()
    if issues:
        out += ["## Issues", ""] + [f"- `{i['path']}`: {i['message']}" for i in issues] + [""]
    out += _files_section(conn, run_key)
    return "\n".join(out) + "\n"


LAMMPS_ROWS = [
    ("Atoms", lambda s: f"{s['n_atoms']}" + (f" ({s['n_types']} types)" if s.get("n_types") else "")
     if s.get("n_atoms") else None),
    ("Composition", lambda s: s.get("composition")),
    ("Elements by type", lambda s: s.get("elements")),
    ("Box", lambda s: s.get("box") and (s["box"] + (f", tilt {s['box_tilt']}" if s.get("box_tilt") else ""))),
    ("Volume", lambda s: s.get("volume") and f"{s['volume']:.6g} {s.get('length_unit', '')}³"),
    ("Density (start)", lambda s: s.get("density") and f"{s['density']:.4g} g/cm³"),
    ("Units / atom style", lambda s: " / ".join(x for x in (s.get("units"), s.get("atom_style")) if x) or None),
    ("Potential", lambda s: s.get("pair_style") and (s["pair_style"] + (f" ({s['potential_files']})"
                                                                         if s.get("potential_files") else ""))),
    ("Starting structure", lambda s: s.get("data_file")),
    ("Structure origin", lambda s: s.get("structure_origin") and f"{s['structure_origin']} (link target, not readable)"),
    ("Initial velocities", lambda s: s.get("velocity_T") is not None and f"{s['velocity_T']:g} K"),
    ("Timestep", lambda s: s.get("timestep") and f"{s['timestep']:g} {s.get('time_unit', '')}"),
    ("Length", lambda s: s.get("total_steps") and
     f"{s['total_steps']:,} steps = {fmt_time(s.get('simulated_time'), s.get('time_unit', ''))}"),
    ("Main ensemble", lambda s: s.get("ensemble")),
    ("Target T", lambda s: s.get("T_range")),
    ("Target P", lambda s: s.get("P_target") is not None and f"{s['P_target']:g} {s.get('pressure_unit', '')}"),
    ("Measured T (2nd half)", lambda s: s.get("T_measured") is not None and f"{s['T_measured']:.5g} K"),
    ("Measured P (2nd half)", lambda s: s.get("P_measured") is not None and
     f"{s['P_measured']:.5g} {s.get('pressure_unit', '')}"),
    ("Measured density", lambda s: s.get("density_measured") and f"{s['density_measured']:.4g} g/cm³"),
    ("Described from", lambda s: s.get("description_source")),
    ("Unresolved variables", lambda s: s.get("unresolved")),
]
VASP_ROWS = [
    ("Calculation", lambda s: s.get("calc_type_text")),
    ("Atoms", lambda s: s.get("n_atoms") and f"{s['n_atoms']} ({s.get('formula')})"),
    ("Composition", lambda s: s.get("composition")),
    ("Cell", lambda s: s.get("box")),
    ("Volume", lambda s: s.get("volume") and (f"{s['volume']:.6g} Å³" + (
        f" → {s['volume_final']:.6g} Å³ final" if s.get("volume_final") and s["volume_final"] != s["volume"] else ""))),
    ("Density", lambda s: s.get("density") and f"{s['density']:.4g} g/cm³"),
    ("Functional", lambda s: s.get("functional")),
    ("ENCUT", lambda s: s.get("encut") and f"{s['encut']:g} eV"),
    ("k-points", lambda s: s.get("kpoints")),
    ("Spin", lambda s: s.get("ispin") == 2 and "polarized (ISPIN=2)"),
    ("EDIFF / EDIFFG", lambda s: (s.get("ediff") or s.get("ediffg")) and f"{s.get('ediff') or '-'} / {s.get('ediffg') or '-'}"),
    ("ISIF / NSW", lambda s: (s.get("isif") is not None or s.get("nsw")) and
     f"{s.get('isif') if s.get('isif') is not None else '-'} / {int(s['nsw']) if s.get('nsw') else '-'}"),
    ("Selective dynamics", lambda s: s.get("selective_dynamics") and "yes"),
    ("MD temperature", lambda s: s.get("T_range")),
    ("MD length", lambda s: s.get("simulated_time") and
     f"{int(s['nsw'])} × {s['timestep']:g} fs = {fmt_time(s['simulated_time'], 'ps')}"),
    ("Measured T (2nd half)", lambda s: s.get("T_measured") and f"{s['T_measured']:.5g} K"),
    ("POTCARs", lambda s: s.get("potcars")),
]


def _system_table(s: Dict) -> List[str]:
    rows = VASP_ROWS if s.get("code") == "vasp" else LAMMPS_ROWS
    out = []
    for label, fn in rows:
        try:
            v = fn(s)
        except (TypeError, ValueError, KeyError):
            v = None
        if v not in (None, False, ""):
            out.append([label, v])
    lines = [_table(["", ""], out)] if out else []
    if s.get("vasp"):
        lines += ["", "VASP in the same directory:", ""] + _system_table(s["vasp"])
    return lines


def _protocol_table(proto: List[Dict], time_unit: str) -> str:
    rows = []
    for i, p in enumerate(proto, 1):
        T = p.get("T")
        rows.append([i, p["ensemble"], (f"{T[0]:g}" if T[0] == T[1] else f"{T[0]:g}→{T[1]:g}") if T else "",
                     (f"{p['P'][0]:g}" if p["P"][0] == p["P"][1] else f"{p['P'][0]:g}→{p['P'][1]:g}") if p.get("P") else "",
                     (f"{p['steps']:,}" if p.get("steps") is not None else "?") + (f" × {p['count']}" if p["count"] > 1 else ""),
                     f"{p['timestep']:g}" if p.get("timestep") else "",
                     fmt_time(p["time"] * p["count"], time_unit) if p.get("time") else "",
                     f"{p['T_mean']:.5g} ± {p['T_std']:.2g}" if p.get("T_mean") is not None and p.get("T_std") is not None else "",
                     f"{p['P_mean']:.4g}" if p.get("P_mean") is not None else "",
                     f"{p['density_mean']:.4g}" if p.get("density_mean") else ""])
    return _table(["#", "ensemble", "T target (K)", "P target", "steps", "dt", "time", "⟨T⟩ (K)", "⟨P⟩", "⟨ρ⟩ (g/cm³)"],
                  rows)


def _params_section(params) -> List[str]:
    """Parameters shared by all sub-runs once; varying ones as a sub-run x key table."""
    if not params:
        return []
    by_key: Dict[tuple, Dict[str, str]] = OrderedDict()
    subs = []
    for p in params:
        if p["subrun"] not in subs:
            subs.append(p["subrun"])
        by_key.setdefault((p["source"], p["key"]), {})[p["subrun"]] = p["value_text"]
    common, varying = [], []
    for (source, key), vals in by_key.items():
        if len(vals) == len(subs) and len(set(vals.values())) == 1:
            common.append((source, key, next(iter(vals.values()))))
        elif len(subs) > 1:
            varying.append(((source, key), vals))
        else:
            common.append((source, key, next(iter(vals.values()))))
    out = ["## Parameters", ""]
    out.append(_table(["source", "key", "value"], [list(c) for c in common]))
    out.append("")
    if varying:
        out.append("### Varying across sub-runs")
        out.append("")
        if len(varying) <= 8:
            headers = ["sub-run"] + [k for (_, k), _ in varying]
            rows = [[s or "(main)"] + [vals.get(s) for _, vals in varying] for s in subs]
            out.append(_table(headers, rows))
        else:
            rows = [[s or "(main)", src, k, vals.get(s)] for (src, k), vals in varying for s in subs if s in vals]
            out.append(_table(["sub-run", "source", "key", "value"], rows[:400]))
        out.append("")
    return out


def _files_section(conn, run_key: str) -> List[str]:
    rows = conn.execute("SELECT category, COUNT(*) n, SUM(size) b, SUM(broken) broken FROM files "
                        "WHERE run_key=? GROUP BY category ORDER BY b DESC", (run_key,)).fetchall()
    if not rows:
        return []
    out = ["## Files", ""]
    out.append(_table(["category", "files", "size"], [[r["category"], r["n"], human_bytes(r["b"])] for r in rows]))
    broken = sum(r["broken"] or 0 for r in rows)
    if broken:
        out.append(f"\n{broken} broken symlinks (targets moved or deleted).")
    big = conn.execute("SELECT relpath, size FROM files WHERE run_key=? ORDER BY size DESC LIMIT 8",
                       (run_key,)).fetchall()
    out.append("\nLargest: " + ", ".join(f"`{b['relpath']}` ({human_bytes(b['size'])})" for b in big))
    out.append("")
    return out


# ---------------------------------------------------------------------------
# Summaries
# ---------------------------------------------------------------------------

def _key_result(conn, run_key: str) -> str:
    row = conn.execute(
        "SELECT key, value_num, unit FROM results WHERE run_key=? AND subrun='' "
        "AND key IN ('energy_per_atom','energy_sigma0') ORDER BY key='energy_per_atom' DESC LIMIT 1",
        (run_key,)).fetchone()
    if row and row["value_num"] is not None:
        return f"{row['key']}={row['value_num']:.6g} {row['unit'] or ''}".strip()
    return ""


def project_readme(conn, project: str) -> str:
    runs = conn.execute("SELECT * FROM runs WHERE project=? ORDER BY run_id", (project,)).fetchall()
    out = [f"# {project}", ""]
    p = conn.execute("SELECT * FROM projects WHERE name=?", (project,)).fetchone()
    if p:
        out += [f"`{p['path']}` · last scanned {p['last_scanned']}", ""]
    counts = Counter(r["status"] for r in runs if not r["missing"])
    out += [", ".join(f"{n} {s}" for s, n in counts.most_common()), ""]
    rows = []
    for r in runs:
        rows.append([f"[{r['run_id']}]({r['run_id'].replace('/', '_')}.md)", r["label"], r["group_path"],
                     r["status"] + (" (missing)" if r["missing"] else ""), r["summary"] or r["code"],
                     human_time(r["wall_time_s"]), _key_result(conn, r["run_key"])])
    out.append(_table(["id", "dir", "group", "status", "summary", "wall", "key result"], rows))
    for kind in ("potentials", "structures"):
        files = conn.execute("SELECT * FROM project_files WHERE project=? AND kind=? ORDER BY name",
                             (project, kind)).fetchall()
        if files:
            out += ["", f"## {kind.title()}", ""]
            out.append(_table(["file", "size", "description", "mentions runs"],
                              [[f["name"], human_bytes(f["size"]), f["description"], f["mentioned_ids"]]
                               for f in files]))
    ana = conn.execute(
        "SELECT a.name, a.n_files, a.missing, GROUP_CONCAT(DISTINCT ar.run_key) runs FROM analysis a "
        "LEFT JOIN analysis_runs ar ON ar.analysis_key = a.analysis_key WHERE a.project=? "
        "GROUP BY a.analysis_key ORDER BY a.name", (project,)).fetchall()
    if ana:
        out += ["", "## Analysis", ""]
        out.append(_table(["analysis", "files", "runs"],
                          [[a["name"] + (" (missing)" if a["missing"] else ""), a["n_files"], a["runs"]] for a in ana]))
    return "\n".join(out) + "\n"


def write_csv(conn, path: Path) -> None:
    runs = conn.execute("SELECT * FROM runs ORDER BY project, run_id").fetchall()
    base = ["run_key", "project", "run_id", "label", "group_path", "code", "calc_type", "status", "summary",
            "n_atoms", "formula", "cores", "nodes", "job_id", "start_time", "end_time", "wall_time_s",
            "n_subruns", "n_files", "total_bytes", "missing", "path"]
    with path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(base + CSV_PARAMS + CSV_RESULTS)
        for r in runs:
            pv = {p["key"]: p["value_text"] for p in conn.execute(
                "SELECT key, value_text FROM params WHERE run_key=? AND subrun=''", (r["run_key"],))}
            rv = {x["key"]: x["value_text"] for x in conn.execute(
                "SELECT key, value_text FROM results WHERE run_key=? AND subrun=''", (r["run_key"],))}
            w.writerow([r[c] for c in base] + [pv.get(k) for k in CSV_PARAMS] + [rv.get(k) for k in CSV_RESULTS])


def scan_report(result: Dict) -> str:
    c = result["counts"]
    out = [f"# Scan {result['scan_id']}", "",
           f"{c['projects']} projects, {c['runs']} runs ({c['new']} new, {c['changed']} changed, "
           f"{c['unchanged']} unchanged, {c['missing']} missing), {c['analysis']} analysis dirs, "
           f"{c['events']} hook events ingested, {c['warnings']} warnings.", ""]
    if result["new"]:
        out += ["## New runs", ""] + [f"- {k}" for k in result["new"]] + [""]
    if result["changed"]:
        out += ["## Changed runs", ""] + [f"- {k}" for k in result["changed"]] + [""]
    if result["warnings"]:
        out += ["## Warnings from this scan", ""] + [f"- `{p}`: {m}" for p, m in result["warnings"]] + [""]
    if result.get("open_issues"):
        out += ["## Open issues in recorded runs", "",
                "Found when each run was last parsed; they clear when the run changes and parses cleanly.", ""]
        out += [f"- {k}: `{p}`: {m}" for k, p, m in result["open_issues"]] + [""]
    return "\n".join(out)


def _write(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp")
    tmp.write_text(text)
    os.replace(tmp, path)


def export_all(conn, cfg: Config, scan_result: Optional[Dict] = None) -> None:
    cfg.cards.mkdir(parents=True, exist_ok=True)
    projects = [p["name"] for p in conn.execute("SELECT name FROM projects ORDER BY name")]
    index = ["# Simulation ledger", ""]
    for project in projects:
        _write(cfg.cards / project / "README.md", project_readme(conn, project))
        index += [f"## [{project}]({project}/README.md)", ""]
        for r in conn.execute("SELECT run_key, run_id, label, status, code, calc_type, missing FROM runs "
                              "WHERE project=? ORDER BY run_id", (project,)):
            _write(card_path(cfg, r["run_key"]), render_card(conn, r["run_key"]))
            index.append(f"- [{r['run_id']}]({project}/{r['run_id'].replace('/', '_')}.md) {r['label']} — "
                         f"{r['code']} {r['calc_type'] or ''} — {r['status']}" + (" (missing)" if r["missing"] else ""))
        index.append("")
    _write(cfg.cards / "INDEX.md", "\n".join(index) + "\n")
    write_csv(conn, cfg.cards / "runs.csv")
    if scan_result is not None:
        _write(cfg.cards / "last_scan.md", scan_report(scan_result))
