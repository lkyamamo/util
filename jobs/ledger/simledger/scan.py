"""`simledger scan`: discover runs, reparse what changed, apply hook events."""

from __future__ import annotations

import hashlib
import json
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from . import PARSER_VERSION, inbox, store
from .classify import code_from_categories
from .config import Config
from .discover import (AnalysisDir, FileEntry, Project, RunDir, calc_dirs, find_analysis,
                       find_projects, find_runs, group_frames, subrun_label, walk_files)
from .parsers import iter_lines, read_head, read_tail
from .parsers import lammps, readme, scheduler, vasp

Warnings = List[Tuple[str, str]]

FINAL_STATES = {"completed", "failed", "unconverged"}


def fingerprint(files: List[FileEntry]) -> str:
    h = hashlib.sha1()
    for f in files:
        h.update(f"{f.relpath}\0{f.size}\0{int(f.mtime)}\0{f.symlink or ''}\n".encode())
    return h.hexdigest()


def _name(f: FileEntry) -> str:
    return f.relpath.rsplit("/", 1)[-1]


def _dir(relpath: str) -> str:
    return relpath.rsplit("/", 1)[0] if "/" in relpath else ""


def _join(d: str, name: str) -> str:
    return f"{d}/{name}" if d else name


# ---------------------------------------------------------------------------
# One calc dir (sub-run)
# ---------------------------------------------------------------------------

def _input_dir(calc: str, by_path: Dict[str, FileEntry]) -> Optional[str]:
    """input_files/ next to the calc dir, or at the run root."""
    parent = _dir(calc)
    for cand in (_join(parent, "input_files"), "input_files"):
        if any(p.startswith(cand + "/") for p in by_path):
            return cand
    return None


def _pick(by_path: Dict[str, FileEntry], dirs: List[Optional[str]], names: List[str]) -> Optional[str]:
    """First readable file (regular, or a symlink that resolves) among dirs x names."""
    for d in dirs:
        if d is None:
            continue
        for n in names:
            p = _join(d, n)
            f = by_path.get(p)
            if f and not f.broken:
                return p
    return None


def _lammps_input(by_path, dirs) -> Optional[str]:
    for d in dirs:
        if d is None:
            continue
        cands = [p for p, f in by_path.items()
                 if _dir(p) == d and f.category == "lammps_input" and not f.broken]
        cands.sort(key=lambda p: (not _name(by_path[p]).startswith("in."), p))
        if cands:
            return cands[0]
    return None


def _job_script(by_path, dirs) -> Optional[str]:
    for d in dirs:
        if d is None:
            continue
        cands = sorted(p for p, f in by_path.items()
                       if _dir(p) == d and f.category == "job_script" and not f.broken
                       and _name(f).endswith((".slurm", ".pbs", ".sbatch")))
        if cands:
            return cands[0]
    return None


def parse_calc(root: Path, calc: str, by_path: Dict[str, FileEntry], cfg: Config,
               warnings: Warnings, driver: bool = False) -> Dict:
    """Status, params and results for one calc dir. Paths are relative to the run dir.

    A driver is the top of a run whose calculations live in sub-runs (e.g. a
    temperature sweep's submit script + STREAM_OUTPUT); its template inputs are
    not attributed to it."""
    inp = None if driver else _input_dir(calc, by_path)
    # The calc dir's own (non-symlink) copy is what actually ran; input_files/ is the template.
    in_dirs = [None if driver else calc, inp]
    names = {_name(f) for p, f in by_path.items() if _dir(p) == calc}
    params: Dict[str, Dict] = {}
    results: Dict = {}
    units: Dict[str, str] = {}
    info: Dict = {"relpath": calc, "code": "unknown", "calc_type": None, "n_atoms": None,
                  "start_time": None, "end_time": None, "wall_time_s": None, "fts": []}
    evidence: List[str] = []
    status = None
    now = time.time()

    # --- scheduler -------------------------------------------------------
    script = _job_script(by_path, [calc, _dir(calc), ""])
    if script:
        text = read_head(root / script, cfg.max_text_bytes)
        sched = scheduler.parse_script(text)
        params["slurm"] = dict(sched["directives"])
        if sched["modules"]:
            params["slurm"]["modules"] = " ".join(sched["modules"])
        for k, v in sched["executables"].items():
            params["slurm"][f"exe:{k}"] = v
        info["script"] = script
        info["job_name"] = sched["directives"].get("job-name")
    stdout_failures: List[str] = []
    if "STREAM_OUTPUT" in names:
        p = root / _join(calc, "STREAM_OUTPUT")
        so = scheduler.parse_stdout(read_head(p, cfg.head_bytes), read_tail(p, cfg.tail_bytes))
        info["start_time"], info["end_time"] = so["start"], so["end"]
        stdout_failures = so["failures"]
        if so["failures"]:
            evidence.append("STREAM_OUTPUT: " + ",".join(so["failures"]))

    # --- LAMMPS ----------------------------------------------------------
    in_path = _lammps_input(by_path, in_dirs)
    if in_path:
        text = read_head(root / in_path, cfg.max_text_bytes)
        li = lammps.parse_input(text)
        info["code"] = "lammps"
        lp = dict(li["settings"])
        for fx in li["fixes"]:
            lp[f"fix:{fx['id']}"] = f"{fx['style']} {fx['args']}".strip()
            if "temp" in fx:
                lp["T_start"], _, rest = fx["temp"].partition(" ")
                lp["T_stop"] = rest.split(" ")[0] if rest else None
            for kw in ("iso", "aniso"):
                if kw in fx:
                    lp["P_start"] = fx[kw].split(" ")[0]
        if li["pair_coeff"]:
            lp["pair_coeff"] = " | ".join(li["pair_coeff"])
        lp["total_steps"] = li["total_steps"]
        lp["n_run_commands"] = len(li["runs"])
        lp["minimize"] = li["minimize"] or None
        lp["simulated_time"] = li["simulated_time"]
        if li["unresolved"]:
            lp["unresolved_vars"] = ",".join(li["unresolved"])
        for d in li["dumps"]:
            lp[f"dump:{d['id']}"] = f"{d['style']} every {d['every']} -> {d['file']} {d['columns']}".strip()
        params["lammps"] = lp
        styles = {fx["style"] for fx in li["fixes"]}
        info["calc_type"] = ("md_npt" if styles & {"npt", "nph"} else
                             "md_nvt" if styles & {"nvt", "langevin", "temp/berendsen", "temp/rescale", "temp/csvr"} else
                             "md_nve" if "nve" in styles else
                             "minimize" if li["minimize"] else None)
        info["input"] = in_path
        info["fts"].append(text[:20000])
        if li["time_unit"]:
            units["simulated_time"] = li["time_unit"]
    log = _join(calc, "log.lammps")
    if log in by_path and not by_path[log].broken:
        lg = lammps.parse_log(read_head(root / log, cfg.head_bytes), read_tail(root / log, cfg.tail_bytes))
        info["code"] = "lammps"
        info["n_atoms"] = lg["n_atoms"]
        info["wall_time_s"] = lg["wall_time_s"]
        results.update({"lammps_version": lg["version"], "wall_time_s": lg["wall_time_s"],
                        "performance": lg["performance"], "procs": lg.get("procs"),
                        "steps_completed": sum(s["steps"] for s in lg["segments"]) or None,
                        "log_warnings": lg["warnings"] or None})
        if lg["finished"]:
            status = "completed"
            evidence.append("log.lammps: Total wall time")
        elif lg["errors"]:
            status = "failed"
            evidence.append("log.lammps: " + lg["errors"][-1][:200])
        elif now - by_path[log].mtime < cfg.running_window_s:
            status = "running"
            evidence.append("log.lammps modified recently, no end marker")
        else:
            status = "incomplete"
            evidence.append("log.lammps has no Total wall time")

    # --- VASP ------------------------------------------------------------
    incar = _pick(by_path, in_dirs, ["INCAR"])
    tags: Dict[str, str] = {}
    if incar:
        text = read_head(root / incar, cfg.max_text_bytes)
        tags = vasp.parse_incar(text)
        params["incar"] = tags
        info["code"] = "both" if info["code"] == "lammps" else "vasp"
        info["calc_type"] = vasp.calc_type(tags)
        info["fts"].append(text)
    kp = _pick(by_path, in_dirs, ["KPOINTS"])
    if kp:
        k = vasp.parse_kpoints(read_head(root / kp, 8192))
        params["kpoints"] = {key: v for key, v in k.items() if key != "comment"}
    pos = _pick(by_path, in_dirs, ["POSCAR"])
    if pos:
        ps = vasp.parse_poscar(read_head(root / pos, cfg.head_bytes))
        if ps:
            info["n_atoms"] = info["n_atoms"] or ps.get("n_atoms")
            info["formula"] = ps.get("formula")
            params["poscar"] = {"formula": ps.get("formula"), "n_atoms": ps.get("n_atoms"),
                                "volume": ps.get("volume"), "selective_dynamics": ps.get("selective_dynamics"),
                                "lattice": " ".join(map(str, ps.get("lattice_lengths", [])))}
    contcar = _pick(by_path, [calc], ["CONTCAR"])
    if contcar and pos:
        cs = vasp.parse_poscar(read_head(root / contcar, cfg.head_bytes))
        ps = params.get("poscar", {})
        if cs.get("volume") and ps.get("volume"):
            results["final_volume"] = cs["volume"]
            results["volume_change_pct"] = round(100 * (cs["volume"] / ps["volume"] - 1), 3)
    potcar = _pick(by_path, in_dirs, ["POTCAR"])
    if potcar:
        titles = vasp.potcar_titles(iter_lines(root / potcar, 50 * 1024 * 1024))
        if titles:
            params["potcar"] = {"titles": "; ".join(t["titel"] for t in titles),
                                "enmax_max": max((t.get("enmax", 0) for t in titles), default=None)}
    outcar = _join(calc, "OUTCAR")
    if outcar in by_path and not by_path[outcar].broken:
        oc = vasp.parse_outcar(read_head(root / outcar, cfg.head_bytes), read_tail(root / outcar, cfg.tail_bytes))
        if info["code"] == "unknown":
            info["code"] = "vasp"
        results.update({"vasp_version": oc.get("version"), "energy_sigma0": oc["energy_sigma0"],
                        "toten": oc["toten"], "e_fermi": oc["e_fermi"], "pressure": oc["pressure_kB"],
                        "elapsed_s": oc["elapsed_s"], "cores": oc.get("cores")})
        units.update({"energy_sigma0": "eV", "toten": "eV", "e_fermi": "eV", "pressure": "kB",
                      "energy_per_atom": "eV/atom", "elapsed_s": "s"})
        if oc["energy_sigma0"] is not None and info["n_atoms"]:
            results["energy_per_atom"] = oc["energy_sigma0"] / info["n_atoms"]
        info["wall_time_s"] = info["wall_time_s"] or oc["elapsed_s"]
        relax = (vasp._tag_num(tags, "NSW") or 0) > 0 and vasp._tag_num(tags, "IBRION") in (1, 2, 3)
        so_tail = read_tail(root / _join(calc, "STREAM_OUTPUT"), 64 * 1024) if "STREAM_OUTPUT" in names else ""
        reached = oc["reached_accuracy"] or "reached required accuracy" in so_tail
        if oc["finished"]:
            status = "unconverged" if relax and not reached else "completed"
            evidence.append("OUTCAR: timing block" + ("" if not relax else
                            (", reached required accuracy" if reached else ", relaxation did not converge")))
        elif now - by_path[outcar].mtime < cfg.running_window_s:
            status = "running"
            evidence.append("OUTCAR modified recently, no timing block")
        else:
            status = status or "incomplete"
            evidence.append("OUTCAR has no timing block")
    osz = _join(calc, "OSZICAR")
    if osz in by_path and not by_path[osz].broken:
        oz = vasp.parse_oszicar(read_tail(root / osz, 64 * 1024))
        if oz:
            results["ionic_steps"] = oz["ionic_steps"]
            if "T" in oz:
                results["final_T"] = oz["T"]

    if stdout_failures and status not in ("completed",):
        status = "failed"
    if status is None:
        status = "unknown" if "STREAM_OUTPUT" in names else "not_started"
        evidence.append("no LAMMPS log or OUTCAR")

    # inputs edited after being copied into the calc dir (only against the
    # calc's own input_files/, not a run-level template shared by many frames)
    if inp and inp == _join(_dir(calc), "input_files"):
        for p, f in by_path.items():
            if _dir(p) == inp and not f.symlink:
                twin = by_path.get(_join(calc, _name(f)))
                if twin and not twin.symlink and twin.size != f.size:
                    warnings.append((str(root / twin.relpath), f"differs from {inp}/{_name(f)}"))

    info.update(status=status, status_evidence="; ".join(evidence), params=params,
                results=results, units=units)
    return info


# ---------------------------------------------------------------------------
# One run
# ---------------------------------------------------------------------------

def aggregate_status(statuses: List[str], unit: str = "calcs") -> Tuple[str, str]:
    if not statuses:
        return "not_started", f"no {unit}"
    counts = ", ".join(f"{statuses.count(s)} {s}" for s in sorted(set(statuses)))
    if all(s == "completed" for s in statuses):
        return "completed", f"{len(statuses)}/{len(statuses)} {unit} completed"
    if any(s == "running" for s in statuses):
        return "running", f"{unit}: {counts}"
    if all(s in ("failed", "unconverged") for s in statuses):
        return (statuses[0] if len(set(statuses)) == 1 else "failed"), f"{unit}: {counts}"
    if any(s == "completed" for s in statuses):
        return "partial", f"{unit}: {counts}"
    return statuses[0], f"{unit}: {counts}"


def parse_frames(root: Path, parent: str, frames: List[str], by_path, cfg: Config,
                 warnings: Warnings) -> Dict:
    """One sub-run for numbered frame dirs (NEB images / path frames)."""
    infos = [parse_calc(root, f, by_path, cfg, warnings) for f in frames]
    info = dict(infos[0], relpath=parent)
    info["status"], info["status_evidence"] = aggregate_status([i["status"] for i in infos], "frames")
    energies = [i["results"].get("energy_sigma0") for i in infos]
    valid = [e for e in energies if e is not None]
    results = {"n_frames": len(infos),
               "frames_completed": sum(i["status"] == "completed" for i in infos),
               "frame_energies": json.dumps(energies)}
    if valid:
        results.update(energy_first=energies[0], energy_last=energies[-1], energy_min=min(valid),
                       energy_max=max(valid), energy_span=max(valid) - min(valid),
                       energy_max_frame=int(frames[energies.index(max(valid))].rpartition("/")[2]))
    info["results"] = results
    info["units"] = {k: "eV" for k in ("energy_first", "energy_last", "energy_min", "energy_max", "energy_span")}
    info["wall_time_s"] = sum(i["wall_time_s"] or 0 for i in infos) or None
    return info


def parse_run(run: RunDir, cfg: Config, warnings: Warnings) -> Dict:
    by_path = {f.relpath: f for f in run.files}
    calcs = calc_dirs(run.files)
    singles, frame_groups = group_frames(calcs)
    n_units = len(singles) + len(frame_groups)
    subruns = []
    for c in singles:
        info = parse_calc(run.path, c, by_path, cfg, warnings, driver=(c == "" and n_units > 1))
        subruns.append(info)
    for parent, frames in sorted(frame_groups.items()):
        subruns.append(parse_frames(run.path, parent, frames, by_path, cfg, warnings))
    for s in subruns:
        s["label"] = (subrun_label(s["relpath"]) or "") if s["relpath"] else "(top)"
    if not subruns and not any(f.relpath.startswith("input_files/") or "/input_files/" in f.relpath
                               for f in run.files):
        warnings.append((str(run.path), "no calculations or input_files found"))

    # The main calc is run/ when present, else the only calc. Its params/results
    # are stored with subrun '' so plain `key=value` searches hit them.
    main = next((s for s in subruns if s["relpath"] == "run"), subruns[0] if len(subruns) == 1 else None)
    for s in subruns:
        s["is_main"] = s is main
        if s is not main and s["label"] == "":
            s["label"] = s["relpath"]
    statuses = [s["status"] for s in subruns if s["label"] != "(top)"] or [s["status"] for s in subruns]
    status, evidence = aggregate_status(statuses)
    if len(subruns) == 1:
        evidence = subruns[0]["status_evidence"]

    # README text at the run's top level
    readme_text = ""
    for f in run.files:
        if "/" not in f.relpath and f.category == "readme" and not f.symlink:
            readme_text += read_head(run.path / f.relpath, cfg.max_text_bytes) + "\n"

    ref = main or next((s for s in subruns if s["label"] != "(top)"), subruns[0] if subruns else {})
    slurm = next((s["params"]["slurm"] for s in ([ref] + subruns) if s and s.get("params", {}).get("slurm")), {})
    starts = [s["start_time"] for s in subruns if s.get("start_time")]
    ends = [s["end_time"] for s in subruns if s.get("end_time")]
    codes = {s["code"] for s in subruns} - {"unknown"}
    if len(codes) == 1:
        code = codes.pop()
    elif codes:
        code = "both"
    else:
        code = code_from_categories(f.category for f in run.files)

    if run.truncated:
        warnings.append((str(run.path), f"walk stopped at depth {cfg.max_walk_depth}"))
    broken = sum(1 for f in run.files if f.broken)

    row = {
        "run_key": run.run_key, "project": run.project, "run_id": run.run_id, "label": run.label,
        "group_path": run.group_path, "path": str(run.path), "code": code,
        "calc_type": ref.get("calc_type") or next((s["calc_type"] for s in subruns if s.get("calc_type")), None),
        "status": status, "status_evidence": evidence,
        "n_atoms": ref.get("n_atoms") or next((s["n_atoms"] for s in subruns if s.get("n_atoms")), None),
        "formula": ref.get("formula") or next((s["formula"] for s in subruns if s.get("formula")), None),
        "start_time": min(starts) if starts else None, "end_time": max(ends) if ends else None,
        "wall_time_s": sum(s["wall_time_s"] or 0 for s in subruns) or None,
        "cores": _int(slurm.get("ntasks")), "nodes": _int(slurm.get("nodes")),
        "job_name": ref.get("job_name") or slurm.get("job-name"), "n_subruns": len(subruns),
        "n_files": len(run.files), "total_bytes": sum(f.size for f in run.files),
        "readme": readme_text.strip() or None,
    }
    fts = " ".join(filter(None, [run.label, run.group_path, run.project, readme_text,
                                 row["job_name"], f"broken_symlinks={broken}" if broken else None]))
    fts += " " + " ".join(t for s in subruns for t in s["fts"])
    return {"row": row, "subruns": subruns, "fts": fts}


def _int(v) -> Optional[int]:
    try:
        return int(str(v).strip())
    except (TypeError, ValueError):
        return None


def write_run(conn, run: RunDir, parsed: Dict, fp: str, first_seen: Optional[str]) -> None:
    ts = store.now()
    store.clear_run(conn, run.run_key)
    row = dict(parsed["row"], fingerprint=fp, parser_version=PARSER_VERSION,
               first_seen=first_seen or ts, last_scanned=ts, last_changed=ts, missing=0)
    store.upsert_row(conn, "runs", ["run_key"], row)
    for s in parsed["subruns"]:
        conn.execute("INSERT OR REPLACE INTO subruns VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                     (run.run_key, s["label"], s["relpath"], s["code"], s["calc_type"], s["status"],
                      s["status_evidence"], s["n_atoms"], s["start_time"], s["end_time"], s["wall_time_s"]))
        sub = "" if s["is_main"] else s["label"]
        for source, items in s["params"].items():
            store.add_params(conn, run.run_key, sub, source, items)
        store.add_results(conn, run.run_key, sub, s["results"], s["units"])
    conn.executemany("INSERT INTO files VALUES (?,?,?,?,?,?,?)",
                     [(run.run_key, f.relpath, f.category, f.size, f.mtime, f.symlink, int(f.broken))
                      for f in run.files])
    store.add_fts(conn, run.run_key, "run", parsed["fts"])


# ---------------------------------------------------------------------------
# Projects: potentials/structures, analysis
# ---------------------------------------------------------------------------

def scan_project_files(conn, project: Project, cfg: Config) -> None:
    conn.execute("DELETE FROM project_files WHERE project = ?", (project.name,))
    for kind in ("potentials", "structures"):
        d = project.path / kind
        if not d.is_dir():
            continue
        files, _ = walk_files(d, cfg)
        sections: Dict[str, str] = {}
        for f in files:
            if f.category == "readme" and "/" not in f.relpath:
                sections.update(readme.file_sections(read_head(d / f.relpath, cfg.max_text_bytes)))
        for f in files:
            if f.category == "readme":
                continue
            desc = sections.get(f.relpath) or sections.get(_name(f))
            conn.execute("INSERT OR REPLACE INTO project_files VALUES (?,?,?,?,?,?,?)",
                         (project.name, kind, f.relpath, f.size, f.mtime, desc,
                          ",".join(readme.mentioned_ids(desc or "", cfg.id_width)) or None))


def scan_analysis(conn, a: AnalysisDir, cfg: Config, run_index: Dict[str, List[str]],
                  subrun_index: Dict[str, List[str]], warnings: Warnings, full: bool) -> str:
    key = f"{a.project}/{a.name}"
    a.files, _ = walk_files(a.path, cfg)
    fp = fingerprint(a.files)
    prev = conn.execute("SELECT fingerprint, first_seen FROM analysis WHERE analysis_key=?", (key,)).fetchone()
    changed = full or prev is None or prev["fingerprint"] != fp
    for w in a.warnings:
        warnings.append((str(a.path), w))
    if changed:
        text = ""
        for f in a.files:
            if f.category == "readme" and f.relpath.count("/") <= 1:
                text += read_head(a.path / f.relpath, cfg.max_text_bytes) + "\n"
        ts = store.now()
        store.upsert_row(conn, "analysis", ["analysis_key"], {
            "analysis_key": key, "project": a.project, "name": a.name, "path": str(a.path),
            "hint": a.hint, "run_ids": ",".join(a.run_ids), "readme": text.strip() or None,
            "n_files": len(a.files), "total_bytes": sum(f.size for f in a.files),
            "fingerprint": fp, "first_seen": prev["first_seen"] if prev else ts, "last_scanned": ts,
            "missing": 0})
        conn.execute("DELETE FROM analysis_files WHERE analysis_key=?", (key,))
        conn.executemany("INSERT INTO analysis_files VALUES (?,?,?,?,?)",
                         [(key, f.relpath, f.category, f.size, f.mtime) for f in a.files])
        if store.has_fts(conn):
            conn.execute("DELETE FROM fts WHERE run_key=? AND kind='analysis'", (key,))
        figs = " ".join(_name(f) for f in a.files if f.category in ("figure", "data", "analysis_script"))
        store.add_fts(conn, key, "analysis", f"{a.name} {a.hint} {text} {figs}")
    else:
        conn.execute("UPDATE analysis SET last_scanned=?, missing=0 WHERE analysis_key=?", (store.now(), key))

    # Name links are cheap and depend on which runs exist, so redo them every scan.
    conn.execute("DELETE FROM analysis_runs WHERE analysis_key=? AND source='name'", (key,))
    hint_tokens = {t.lower() for t in a.hint.replace("-", "_").split("_") if t}
    for rid in a.run_ids:
        keys = run_index.get(rid, [])
        own = [k for k in keys if k.startswith(a.project + "/")]
        targets = own or keys
        if not targets:
            warnings.append((str(a.path), f"run id {rid} not found in any project"))
            continue
        if not own:
            warnings.append((str(a.path), f"run id {rid} found only in another project: {', '.join(keys)}"))
        for rk in targets:
            sub = next((s for s in subrun_index.get(rk, []) if s.lower() in hint_tokens), None)
            conn.execute("INSERT INTO analysis_runs VALUES (?,?,?,'name')", (key, rk, sub))
    return "changed" if changed else "unchanged"


# ---------------------------------------------------------------------------
# Hook events -> runs
# ---------------------------------------------------------------------------

def apply_events(conn) -> None:
    """Re-derive everything that comes from hook events (idempotent)."""
    inbox.link_events(conn)
    conn.execute("DELETE FROM params WHERE source IN ('user', 'hook')")
    conn.execute("DELETE FROM analysis_runs WHERE source = 'event'")
    for ev in conn.execute("SELECT * FROM events WHERE run_key IS NOT NULL ORDER BY created"):
        data = json.loads(ev["data"])
        rk, sub = ev["run_key"], ev["subrun"] or ""
        if data.get("fields"):
            store.add_params(conn, rk, sub, "user", data["fields"])
        slurm = data.get("slurm") or {}
        if ev["kind"] == "run":
            store.add_params(conn, rk, sub, "hook", {
                "job_id": slurm.get("SLURM_JOB_ID"), "nodelist": slurm.get("SLURM_JOB_NODELIST"),
                "partition": slurm.get("SLURM_JOB_PARTITION"), "exit_code": data.get("exit_code")})
            if slurm.get("SLURM_JOB_ID"):
                conn.execute("UPDATE runs SET job_id=? WHERE run_key=?", (slurm["SLURM_JOB_ID"], rk))
            code = data.get("exit_code")
            if code not in (None, 0, "0"):
                ev_note = f"hook: job {slurm.get('SLURM_JOB_ID', '?')} exited {code}"
                if sub:
                    conn.execute("UPDATE subruns SET status='failed', status_evidence=? WHERE run_key=? AND label=?",
                                 (ev_note, rk, sub))
                if not sub or conn.execute("SELECT n_subruns FROM runs WHERE run_key=?", (rk,)).fetchone()[0] <= 1:
                    conn.execute("UPDATE runs SET status='failed', status_evidence=? WHERE run_key=?", (ev_note, rk))
        elif ev["kind"] == "analysis":
            akey = _analysis_key_for(conn, data.get("output")) or f"event:{data.get('analysis_type')}:{ev['event_id']}"
            conn.execute("INSERT INTO analysis_runs VALUES (?,?,?,'event')", (akey, rk, sub or None))


def _analysis_key_for(conn, output: Optional[str]) -> Optional[str]:
    if not output:
        return None
    parts = [p for p in output.split("/") if p]
    if "analysis" in parts:
        i = len(parts) - 1 - parts[::-1].index("analysis")
        if 0 < i < len(parts) - 1:
            key = f"{parts[i - 1]}/{parts[i + 1]}"
            if conn.execute("SELECT 1 FROM analysis WHERE analysis_key=?", (key,)).fetchone():
                return key
    return None


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def scan(cfg: Config, root: Path, full: bool = False, log=print) -> Dict:
    root = root.resolve()
    home = cfg.home.resolve()
    if home == root or root in home.parents:
        raise SystemExit(f"ledger home {home} is inside the scanned tree {root}; refusing to write there")

    conn = store.connect(cfg.db_path)
    warnings: Warnings = []
    started = store.now()
    scan_id = conn.execute("INSERT INTO scans (root, started) VALUES (?,?)", (str(root), started)).lastrowid
    counts = {"events": inbox.ingest(conn, cfg, warnings), "projects": 0, "runs": 0, "new": 0,
              "changed": 0, "unchanged": 0, "missing": 0, "analysis": 0, "analysis_changed": 0}
    new_runs: List[str] = []
    changed_runs: List[str] = []

    projects = find_projects(root, cfg)
    if not projects:
        warnings.append((str(root), "no project (directory containing runs/) found"))
    run_index: Dict[str, List[str]] = {}
    for r in conn.execute("SELECT run_key, run_id FROM runs WHERE missing=0"):
        run_index.setdefault(r["run_id"], []).append(r["run_key"])

    for project in projects:
        counts["projects"] += 1
        missing = [d for d, ok in project.present.items() if not ok]
        if missing:
            warnings.append((str(project.path), "missing project dirs: " + ", ".join(missing)))
        ts = store.now()
        prev = conn.execute("SELECT first_seen FROM projects WHERE name=?", (project.name,)).fetchone()
        store.upsert_row(conn, "projects", ["name"], {
            "name": project.name, "path": str(project.path), "present": json.dumps(project.present),
            "first_seen": prev["first_seen"] if prev else ts, "last_scanned": ts})
        scan_project_files(conn, project, cfg)

        seen = set()
        for run in find_runs(project, cfg, warnings):
            if run.run_key in seen:
                warnings.append((str(run.path), f"duplicate run id {run.run_id} in project {project.name}; "
                                                 f"recorded as {run.run_key}~{run.label}"))
                run = RunDir(run.project, f"{run.run_id}~{run.label}", run.label, run.group_path, run.path)
            seen.add(run.run_key)
            others = [k for k in run_index.get(run.run_id, []) if not k.startswith(project.name + "/")]
            if others:
                warnings.append((str(run.path), f"run id {run.run_id} also used in {', '.join(others)}"))
            run_index.setdefault(run.run_id, [])
            if run.run_key not in run_index[run.run_id]:
                run_index[run.run_id].append(run.run_key)

            run.files, run.truncated = walk_files(run.path, cfg)
            fp = fingerprint(run.files)
            counts["runs"] += 1
            prev = conn.execute("SELECT fingerprint, parser_version, first_seen, status, path FROM runs "
                                "WHERE run_key=?", (run.run_key,)).fetchone()
            if prev and prev["path"] != str(run.path) and Path(prev["path"]).exists():
                warnings.append((str(run.path), f"{run.run_key} previously recorded at {prev['path']}"))
            stale = (full or prev is None or prev["fingerprint"] != fp
                     or prev["parser_version"] != PARSER_VERSION or prev["status"] == "running")
            if not stale:
                conn.execute("UPDATE runs SET last_scanned=?, missing=0 WHERE run_key=?",
                             (store.now(), run.run_key))
                counts["unchanged"] += 1
                continue
            run_warns: Warnings = []
            parsed = parse_run(run, cfg, run_warns)
            write_run(conn, run, parsed, fp, prev["first_seen"] if prev else None)
            conn.executemany("INSERT INTO run_warnings VALUES (?,?,?)", [(run.run_key, p, m) for p, m in run_warns])
            warnings.extend(run_warns)
            if prev is None:
                counts["new"] += 1
                new_runs.append(run.run_key)
            else:
                counts["changed"] += 1
                changed_runs.append(run.run_key)
            conn.commit()
            if (counts["new"] + counts["changed"]) % 20 == 0:
                log(f"  ... {counts['runs']} runs scanned")

        stale_rows = conn.execute("SELECT run_key FROM runs WHERE project=? AND missing=0", (project.name,)).fetchall()
        gone = [r["run_key"] for r in stale_rows if r["run_key"] not in seen]
        for rk in gone:
            conn.execute("UPDATE runs SET missing=1 WHERE run_key=?", (rk,))
        counts["missing"] += len(gone)

        subrun_index: Dict[str, List[str]] = {}
        for s in conn.execute("SELECT run_key, label FROM subruns"):
            subrun_index.setdefault(s["run_key"], []).append(s["label"])
        seen_a = set()
        for a in find_analysis(project, cfg):
            counts["analysis"] += 1
            seen_a.add(f"{a.project}/{a.name}")
            if scan_analysis(conn, a, cfg, run_index, subrun_index, warnings, full) == "changed":
                counts["analysis_changed"] += 1
        for r in conn.execute("SELECT analysis_key FROM analysis WHERE project=? AND missing=0",
                              (project.name,)).fetchall():
            if r["analysis_key"] not in seen_a:
                conn.execute("UPDATE analysis SET missing=1 WHERE analysis_key=?", (r["analysis_key"],))
        conn.commit()

    apply_events(conn)
    unlinked = conn.execute("SELECT COUNT(*) FROM events WHERE run_key IS NULL").fetchone()[0]
    if unlinked:
        warnings.append(("inbox", f"{unlinked} hook events not matched to any run yet"))
    conn.executemany("INSERT INTO scan_warnings VALUES (?,?,?)", [(scan_id, p, m) for p, m in warnings])
    counts["warnings"] = len(warnings)
    conn.execute("UPDATE scans SET finished=?, counts=? WHERE scan_id=?",
                 (store.now(), json.dumps(counts), scan_id))
    conn.commit()
    open_issues = [(r["run_key"], r["path"], r["message"]) for r in conn.execute(
        "SELECT w.* FROM run_warnings w JOIN runs r USING (run_key) WHERE r.missing = 0 ORDER BY run_key")]
    return {"scan_id": scan_id, "counts": counts, "new": new_runs, "changed": changed_runs,
            "warnings": warnings, "open_issues": open_issues, "conn": conn}
