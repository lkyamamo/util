"""`simledger scan`: discover runs, reparse what changed, apply hook events."""

from __future__ import annotations

import hashlib
import itertools
import json
import os
import re
import sqlite3
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from . import PARSER_VERSION, describe, inbox, results as rx, sacct, store
from .classify import code_from_categories
from .config import Config
from .discover import (AnalysisDir, FileEntry, Project, RunDir, calc_dirs, find_analysis,
                       find_projects, find_runs, group_frames, locate, subrun_label, walk_files)
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


def _job_scripts(by_path, dirs) -> List[str]:
    """Every *.slurm/*.pbs/*.sbatch in dirs, nearest dir first (duplicates by name dropped)."""
    out: List[str] = []
    seen = set()
    for d in dict.fromkeys(dirs):
        if d is None:
            continue
        for p in sorted(p for p, f in by_path.items()
                        if _dir(p) == d and f.category == "job_script" and not f.broken
                        and _name(f).endswith((".slurm", ".pbs", ".sbatch"))):
            if _name(by_path[p]) not in seen:
                seen.add(_name(by_path[p]))
                out.append(p)
    return out


def _lammps_logs(by_path: Dict[str, FileEntry], calc: str) -> List[str]:
    """log.lammps then log.lammps.resumeN in order (dielectric restarts write these)."""
    logs = [p for p, f in by_path.items() if _dir(p) == calc and f.category == "lammps_log" and not f.broken]

    def order(p):
        n = _name(by_path[p])
        m = re.search(r"(\d+)(?!.*\d)", n)   # last number in the name: resume3, log_050
        return (0 if n == "log.lammps" else 1, int(m.group(1)) if m else 0, n)
    return sorted(logs, key=order)


def _candidates(name: str, dirs: List[Optional[str]], by_path: Dict[str, FileEntry]) -> List[FileEntry]:
    """Inventory entries a command's file name can refer to, nearest dir first, input_files/ last."""
    base = name.rsplit("/", 1)[-1]
    out = []
    for d in dirs + ["input_files"]:
        if d is None:
            continue
        for cand in (_join(d, name), _join(d, base)):
            f = by_path.get(os.path.normpath(cand).replace(os.sep, "/"))
            if f and f not in out:
                out.append(f)
    return out


def _data_file(name: str, dirs: List[Optional[str]], by_path: Dict[str, FileEntry]) -> Tuple[Optional[str], Optional[str]]:
    """(readable relpath, target of a broken link with that name) for a read_data/write_data file.

    The target records where a starting structure came from (often another run)."""
    broken_target = None
    for f in _candidates(name, dirs, by_path):
        if not f.broken:
            return f.relpath, broken_target
        # keep the last one: run/x -> input_files/x -> <real origin>, so input_files/ wins
        broken_target = f.symlink
    return None, broken_target


def _link_origin(name: str, dirs: List[Optional[str]], by_path: Dict[str, FileEntry]) -> Optional[str]:
    """Where a linked-in file (e.g. OH.usc -> potentials/20260802_OH_B_mod_v11.usc) really comes from:
    the last link in the chain run/x -> input_files/x -> <source>."""
    targets = [f.symlink for f in _candidates(name, dirs, by_path) if f.symlink]
    return targets[-1] if targets else None


def _remap(target: str, tree: Optional[Path]) -> Optional[Path]:
    """Where a link target recorded before a move now lives under the scanned tree.

    /scratch1/lkyamamo/<project>/runs/0177/run/final.data -> <tree>/<project>/runs/0177/run/final.data:
    leading components are dropped until the remainder exists under tree."""
    if tree is None:
        return None
    parts = [p for p in target.split("/") if p]
    for i in range(len(parts) - 1):
        cand = tree.joinpath(*parts[i:])
        if os.path.lexists(cand):
            return cand
    return None


def follow_moved_link(target: str, tree: Optional[Path], hops: int = 5) -> Tuple[Optional[Path], str]:
    """(readable file, deepest target) following a chain of links that may point at moved paths."""
    for _ in range(hops):
        cand = _remap(target, tree) if not os.path.isfile(target) else Path(target)
        if cand is None:
            return None, target
        if cand.is_file():
            return cand, str(cand)
        if cand.is_symlink():
            target = os.readlink(cand)
            continue
        return None, target
    return None, target


def _lammps_description(root: Path, calc: str, logs: List[str], in_path: Optional[str],
                        in_dirs: List[Optional[str]], by_path: Dict[str, FileEntry], cfg: Config,
                        warnings: Warnings, tree: Optional[Path] = None) -> Dict:
    """System/protocol from the log (what actually ran); data file for per-element
    counts and masses; the input script only when there is no usable log."""
    log_info: Dict = {}
    total = sum(by_path[l].size for l in logs)
    if logs and total <= cfg.max_log_bytes:
        log_info = lammps.parse_log_stream(itertools.chain.from_iterable(
            iter_lines(root / l, cfg.max_log_bytes) for l in logs))
    elif logs:
        warnings.append((str(root / logs[0]), f"logs total {total} bytes > max_log_bytes; "
                                              "protocol taken from the input script"))
    source = "log"
    if not log_info.get("segments") and in_path:   # log without echoed run/fix commands
        inp = lammps.parse_log_stream(iter_lines(root / in_path, cfg.max_text_bytes), substituted=False)
        for key in ("units", "atom_style", "pair_style", "pair_coeff", "masses", "data_files",
                    "timestep", "segments", "unresolved", "velocity_T"):
            if not log_info.get(key):
                log_info[key] = inp.get(key)
        log_info.setdefault("box", None)
        source = "input script" + (" (log has no echoed commands)" if logs else " (no log yet)")
    data_rel, origin, data_note = None, None, None
    data: Dict = {}
    for name in log_info.get("data_files") or []:
        data_rel, target = _data_file(name, [calc] + in_dirs, by_path)
        origin = origin or target
        if data_rel:
            break
    data_path: Optional[Path] = root / data_rel if data_rel else None
    if not data_rel and origin:   # broken link: the file may have moved with the data
        found, origin = follow_moved_link(origin, tree)
        if found:
            data_path, data_rel = found, str(found)
            data_note = "followed a broken link to its moved location"
    if not data_path:   # starting structure unreadable: the written end state has the same atoms
        for name in log_info.get("write_data") or []:
            data_rel, _ = _data_file(name, [calc], by_path)
            if data_rel:
                data_path = root / data_rel
                data_note = "end state; starting data file unavailable"
                break
    if data_path:
        try:
            count = data_path.stat().st_size <= cfg.max_data_bytes
        except OSError:
            count = False
        data = lammps.parse_data_file(iter_lines(data_path, cfg.max_data_bytes if count else cfg.head_bytes),
                                      log_info.get("atom_style"), count_types=count)
    desc = describe.lammps_system(log_info, data, data_rel)
    desc["system"]["description_source"] = source
    pot_sources = []
    for pot in (desc["system"].get("potential_files") or "").split(", "):
        target = _link_origin(pot, [calc] + in_dirs, by_path) if pot else None
        if target:
            pot_sources.append(target)
    if pot_sources:
        desc["system"]["potential_source"] = ", ".join(pot_sources)
    if data_note:
        desc["system"]["data_file"] = f"{data_rel} ({data_note})"
    if origin:
        desc["system"]["structure_origin"] = origin
        loc = locate(origin, cfg)
        here = locate(str(root), cfg)
        if loc and not (here and (loc["project"], loc["run_id"]) == (here["project"], here["run_id"])):
            desc["system"]["structure_origin_run"] = f"{loc['project']}/{loc['run_id']}"
    return desc


def parse_calc(root: Path, calc: str, by_path: Dict[str, FileEntry], cfg: Config,
               warnings: Warnings, driver: bool = False, tree: Optional[Path] = None) -> Dict:
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
    scripts = _job_scripts(by_path, [calc, _dir(calc), ""])
    if scripts:
        launches: List[str] = []
        for i, script in enumerate(scripts):
            sched = scheduler.parse_script(read_head(root / script, cfg.max_text_bytes))
            launches += sched["launch"]
            if i == 0:   # the calc dir's own script describes this calc
                params["slurm"] = dict(sched["directives"])
                if sched["modules"]:
                    params["slurm"]["modules"] = " ".join(sched["modules"])
                for k, v in sched["executables"].items():
                    params["slurm"][f"exe:{k}"] = v
                if sched["launch"]:
                    params["slurm"]["launch"] = sched["launch"][0][:300]
                info["script"] = script
                info["job_name"] = sched["directives"].get("job-name")
        if len(scripts) > 1:
            params["slurm"]["scripts"] = ", ".join(scripts)
        lv = scheduler.launch_vars(launches)
        if lv:
            params["lmp_var"] = lv
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
    logs = _lammps_logs(by_path, calc)
    if logs or in_path:
        info.update(_lammps_description(root, calc, logs, in_path, in_dirs, by_path, cfg, warnings, tree))
    if logs:
        log = logs[-1]   # a resumed run's newest log decides its status
        lg = lammps.parse_log(read_head(root / logs[0], cfg.head_bytes), read_tail(root / log, cfg.tail_bytes))
        info["code"] = "lammps"
        info["n_atoms"] = info.get("n_atoms") or lg["n_atoms"]
        info["wall_time_s"] = lg["wall_time_s"]
        results.update({"lammps_version": lg["version"], "wall_time_s": lg["wall_time_s"],
                        "performance": lg["performance"], "procs": lg.get("procs"),
                        "steps_completed": sum(s["steps"] for s in lg["segments"]) or None,
                        "log_warnings": lg["warnings"] or None})
        lname = _name(by_path[log])
        if lg["finished"]:
            status = "completed"
            evidence.append(f"{lname}: Total wall time")
        elif lg["errors"]:
            status = "failed"
            evidence.append(f"{lname}: " + lg["errors"][-1][:200])
        elif now - by_path[log].mtime < cfg.running_window_s:
            status = "running"
            evidence.append(f"{lname} modified recently, no end marker")
        else:
            status = "incomplete"
            evidence.append(f"{lname} has no Total wall time")

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
    k: Dict = {}
    if kp:
        k = vasp.parse_kpoints(read_head(root / kp, 8192))
        params["kpoints"] = {key: v for key, v in k.items() if key != "comment"}
    pos = _pick(by_path, in_dirs, ["POSCAR"])
    ps: Dict = {}
    cs: Dict = {}
    titles: List[Dict] = []
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
        if cs.get("volume") and ps.get("volume"):
            results["final_volume"] = cs["volume"]
            results["volume_change_pct"] = round(100 * (cs["volume"] / ps["volume"] - 1), 3)
    potcar = _pick(by_path, in_dirs, ["POTCAR"])
    if potcar:
        titles = vasp.potcar_titles(iter_lines(root / potcar, 50 * 1024 * 1024)) or []
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
    osz_md: Dict = {}
    if osz in by_path and not by_path[osz].broken:
        oz = vasp.parse_oszicar(read_tail(root / osz, 64 * 1024))
        if oz:
            results["ionic_steps"] = oz["ionic_steps"]
            if "T" in oz:
                results["final_T"] = oz["T"]
                osz_md = vasp.oszicar_md(iter_lines(root / osz, cfg.max_log_bytes))
    if incar or pos:
        desc = describe.vasp_system(tags, k, ps, cs, titles, osz_md)
        if info.get("system"):   # LAMMPS + VASP in one dir: keep the LAMMPS one, note VASP
            info["system"]["vasp"] = desc["system"]
        else:
            info.update(desc)

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

    sysd = info.get("system")
    if sysd:
        params["system"] = {k: v for k, v in sysd.items()
                            if isinstance(v, (int, float, str)) and not isinstance(v, bool) and k != "code"}
        info["n_atoms"] = info["n_atoms"] or sysd.get("n_atoms")
        info["formula"] = info.get("formula") or sysd.get("formula")
        ens = (sysd.get("ensemble") or "").split(" ")[0]
        if info["code"] == "lammps" and ens:
            info["calc_type"] = {"NVT": "md_nvt", "NPT": "md_npt", "NVE": "md_nve", "NPH": "md_nph",
                                 "minimize": "minimize"}.get(ens, info["calc_type"])
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
                 warnings: Warnings, tree: Optional[Path] = None) -> Dict:
    """One sub-run for numbered frame dirs (NEB images / path frames)."""
    infos = [parse_calc(root, f, by_path, cfg, warnings, tree=tree) for f in frames]
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
    info["summary"] = describe.frames_summary(infos[0], len(infos), results)
    info["conditions"] = f"{len(infos)} frames" + (f", span {results['energy_span']:.3g} eV"
                                                    if results.get("energy_span") is not None else "")
    info["protocol"] = []
    info["wall_time_s"] = sum(i["wall_time_s"] or 0 for i in infos) or None
    return info


def _main_label(subruns: List[Dict]) -> Optional[str]:
    """Label of the sub-run whose values are stored as the run's own ('' sub-run)."""
    real = [s for s in subruns if s["label"] != "(top)"]
    main = next((s for s in subruns if s["relpath"] == "run"), real[0] if len(real) == 1 else None)
    return main["label"] if main else None


def target_subrun(v: "rx.Value", file_dir: Optional[str], subs: List[Tuple[str, str, Optional[float]]],
                  main: Optional[str], default: str = "") -> Tuple[str, Optional[str]]:
    """(sub-run label, note) a result value belongs to.

    subs: (label, relpath, T_target). Explicit labels win, then temperature
    (sweeps), then the directory the file sits in."""
    real = [s for s in subs if s[0] != "(top)"]
    label, note = None, None
    if v.subrun:
        label = next((s[0] for s in subs if s[0] == v.subrun or s[0].endswith("/" + v.subrun)), v.subrun)
    elif v.T is not None:
        hits = [s for s in real if s[2] is not None and abs(s[2] - v.T) < 0.6]
        if hits:
            # several sub-runs at that T (e.g. T303 and T303/run2-files): the plainest label
            label = min(hits, key=lambda s: (len(s[0]), s[0]))[0]
        elif len(real) == 1:
            label = real[0][0]
        else:
            note = f"T = {v.T:g} K (no single sub-run at that temperature)"
    elif file_dir is not None:
        best = max((s for s in subs if s[1] and (file_dir == s[1] or file_dir.startswith(s[1] + "/"))),
                   key=lambda s: len(s[1]), default=None)
        label = best[0] if best else None
    if label is None:
        label = default
    if main is not None and label == main:
        label = ""
    return label, note


def run_file_results(run: RunDir, subruns: List[Dict], warnings: Warnings, limit: int = 300):
    """Results from analysis outputs written inside the run dir (sweep tables, SUMMARY.txt...)."""
    subs = [(s["label"], s["relpath"], (s.get("system") or {}).get("T_target")) for s in subruns]
    main = _main_label(subruns)
    found, texts = [], []
    cands = [f for f in run.files if not f.broken and rx.is_candidate(_name(f))
             and not any(part in ("input_files", "setup") for part in f.relpath.split("/")[:-1])]
    for f in cands[:limit]:
        path = run.path / f.relpath
        size = f.size if not f.symlink else (path.stat().st_size if path.exists() else 0)
        got = rx.extract(path, size)
        if not got:
            continue
        warnings.extend((str(path), w) for w in got.warnings)
        if got.text:
            texts.append(got.text)
        for v in got.values:
            label, note = target_subrun(v, _dir(f.relpath), subs, main)
            found.append((label, v, f"file:{f.relpath}", "; ".join(x for x in (v.note, note) if x) or None))
    return found, texts


NON_SIM_ROLES = ("calc", "aggregat", "analysis", "distribution", "collect", "plot")


def run_pipeline(run: RunDir, subruns: List[Dict], cfg: Config) -> Dict:
    """Pipeline settings (.conf), pipeline logs and every job id the run's files mention.

    Jobs are mapped to sub-runs by the temperature of their log section, by a
    `cascade` role, or by the directory of a *_<jobid>.out file."""
    subs = [(s["label"], s["relpath"], (s.get("system") or {}).get("T_target")) for s in subruns]
    main = _main_label(subruns)
    params: Dict[str, Dict] = {}
    jobs: List[Tuple[str, str, str, str, Optional[str]]] = []   # (subrun, job_id, role, source, detail)
    texts: List[str] = []
    for f in run.files:
        if f.broken:
            continue
        depth, name, path = f.relpath.count("/"), _name(f), run.path / f.relpath
        if name.endswith(".conf") and depth <= 1 and f.size <= 1024 * 1024:
            conf = scheduler.parse_conf(read_head(path, 1024 * 1024))
            if conf:
                params[f"conf:{f.relpath}"] = conf
        elif f.category == "log" and depth <= 1 and f.size <= cfg.max_log_bytes:
            pl = scheduler.parse_pipeline_log(iter_lines(path, cfg.max_log_bytes))
            if pl["script"]:
                params.setdefault("pipeline", {}).update(
                    {"script": pl["script"], "started": pl["started"], "id": pl["sweep_id"]})
            if pl["settings"]:
                params.setdefault("pipeline", {}).update(pl["settings"])
            for j in pl["jobs"]:
                if "aggregat" in j["role"]:   # combines every temperature: belongs to the whole run
                    label = ""
                elif j["T"] is not None:
                    label, _ = target_subrun(rx.Value("", None, T=j["T"]), None, subs, main)
                elif j["role"] == "cascade":
                    label = next((s[0] for s in subs if s[0] == "cascade"), "")
                    label = "" if label == main else label
                else:
                    label = ""
                jobs.append((label, j["job_id"], j["role"], f"log:{f.relpath}", j["section"]))
        jid = scheduler.job_id_from_filename(name)
        if jid:
            role = re.sub(r"[_.-]?\d{6,9}\.(out|err|log)$", "", name) or "slurm"
            label, _ = target_subrun(rx.Value("", None), _dir(f.relpath), subs, main)
            jobs.append((label, jid, role, f"file:{f.relpath}", None))
    for conf in params.values():
        for k in ("INPUT_SCRIPT", "STARTING_STRUCTURE", "POTENTIAL_FILE"):
            if conf.get(k):
                texts.append(conf[k])
    return {"params": params, "jobs": jobs, "texts": texts}


def parse_run(run: RunDir, cfg: Config, warnings: Warnings, tree: Optional[Path] = None) -> Dict:
    by_path = {f.relpath: f for f in run.files}
    calcs = calc_dirs(run.files)
    singles, frame_groups = group_frames(calcs)
    n_units = len(singles) + len(frame_groups)
    subruns = []
    for c in singles:
        info = parse_calc(run.path, c, by_path, cfg, warnings, driver=(c == "" and n_units > 1), tree=tree)
        subruns.append(info)
    for parent, frames in sorted(frame_groups.items()):
        subruns.append(parse_frames(run.path, parent, frames, by_path, cfg, warnings, tree=tree))
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
    summary = (main or {}).get("summary") or describe.multi_summary(subruns)
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
        "summary": summary,
    }
    file_results, result_texts = run_file_results(run, subruns, warnings)
    pipe = run_pipeline(run, subruns, cfg)
    result_texts += pipe["texts"]
    # the pipeline conf names the starting structure and potential version, even when links are gone
    sysd = ref.get("system")
    for conf in pipe["params"].values():
        if sysd is None:
            break
        if conf.get("STARTING_STRUCTURE") and not sysd.get("structure_origin"):
            sysd["structure_origin"] = conf["STARTING_STRUCTURE"]
            loc = locate(conf["STARTING_STRUCTURE"], cfg)
            if loc and (loc["project"], loc["run_id"]) != (run.project, run.run_id):
                sysd["structure_origin_run"] = f"{loc['project']}/{loc['run_id']}"
        if conf.get("POTENTIAL_FILE") and not sysd.get("potential_source"):
            sysd["potential_source"] = conf["POTENTIAL_FILE"]
    fts = " ".join(filter(None, [run.label, run.group_path, run.project, readme_text, summary, *result_texts,
                                 row["job_name"], f"broken_symlinks={broken}" if broken else None]))
    fts += " " + " ".join(t for s in subruns for t in s["fts"])
    row["system"] = json.dumps(ref.get("system")) if ref.get("system") else None
    return {"row": row, "subruns": subruns, "fts": fts, "file_results": file_results,
            "pipeline_params": pipe["params"], "jobs": pipe["jobs"]}


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
        store.upsert_row(conn, "subruns", ["run_key", "label"], {
            "run_key": run.run_key, "label": s["label"], "relpath": s["relpath"], "code": s["code"],
            "calc_type": s["calc_type"], "status": s["status"], "status_evidence": s["status_evidence"],
            "n_atoms": s["n_atoms"], "start_time": s["start_time"], "end_time": s["end_time"],
            "wall_time_s": s["wall_time_s"], "summary": s.get("summary"), "conditions": s.get("conditions"),
            "system": json.dumps(s["system"]) if s.get("system") else None,
            "protocol": json.dumps(s["protocol"]) if s.get("protocol") else None})
        sub = "" if s["is_main"] else s["label"]
        for source, items in s["params"].items():
            store.add_params(conn, run.run_key, sub, source, items)
        store.add_results(conn, run.run_key, sub, s["results"], s["units"])
    for source, items in parsed.get("pipeline_params", {}).items():
        store.add_params(conn, run.run_key, "", source, items)
    conn.executemany("INSERT INTO jobs (run_key, subrun, analysis_key, job_id, role, source, detail) "
                     "VALUES (?,?,NULL,?,?,?,?)",
                     [(run.run_key, sub, jid, role, src, det) for sub, jid, role, src, det in parsed.get("jobs", [])])
    for label, v, source, note in parsed.get("file_results", []):
        store.add_results(conn, run.run_key, label, {v.key: v.value}, {v.key: v.unit}, source, {v.key: note})
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


def describe_analysis(conn, key: str, a: AnalysisDir, cfg: Config, warnings: Warnings,
                      limit: int = 300) -> Tuple[Optional[str], str]:
    """Extract result values into analysis_results; return (description, searchable text).

    The description is each top-level script's docstring summary; job scripts
    and small command files (e.g. silanol-run.txt) add searchable text."""
    conn.execute("DELETE FROM analysis_results WHERE analysis_key=?", (key,))
    conn.execute("DELETE FROM jobs WHERE analysis_key=? AND source != 'hook'", (key,))
    descs, texts = [], []
    for f in a.files[:2000]:
        if f.broken:
            continue
        name, depth, path = _name(f), f.relpath.count("/"), a.path / f.relpath
        jid = scheduler.job_id_from_filename(name)
        if jid:
            role = re.sub(r"[_.-]?\d{6,9}\.(out|err|log)$", "", name) or "slurm"
            conn.execute("INSERT INTO jobs (run_key, subrun, analysis_key, job_id, role, source) "
                         "VALUES (NULL, NULL, ?, ?, ?, ?)", (key, jid, role, f"file:{f.relpath}"))
        if name.endswith(".py") and depth <= 1 and f.size <= cfg.max_text_bytes:
            d = rx.script_description(read_head(path, 16384))
            if d:
                descs.append(f"{f.relpath}: {d}")
        elif f.category == "job_script" and depth <= 1:
            sched = scheduler.parse_script(read_head(path, 65536))
            texts.append(" ".join(filter(None, [sched["directives"].get("job-name"), *sched["launch"][:2]])))
        elif name.endswith(".txt") and f.size <= 4096 and name not in rx.TEXT_RESULT_NAMES and depth <= 2:
            texts.append(read_head(path, 4096))
        if limit and rx.is_candidate(name):
            limit -= 1
            got = rx.extract(path, f.size)
            if not got:
                continue
            for w in got.warnings:
                warnings.append((str(path), w))
                conn.execute("INSERT INTO analysis_results (analysis_key, file, kind, key, value_text) "
                             "VALUES (?,?,?,?,?)", (key, f.relpath, got.kind, "__warning__", w))
            if got.text:
                texts.append(got.text)
            conn.executemany(
                "INSERT INTO analysis_results VALUES (?,?,?,?,?,?,?,?,?,?)",
                [(key, f.relpath, got.kind, v.key, store._num(v.value), store._text(v.value), v.unit,
                  v.subrun, v.T, v.note) for v in got.values])
    return (" | ".join(descs) or None), " ".join(t for t in texts if t)


def materialize_analysis_results(conn) -> None:
    """Copy analysis values onto every run/sub-run the analysis is linked to."""
    conn.execute("DELETE FROM results WHERE source LIKE 'analysis:%'")
    subs: Dict[str, List[Tuple[str, str, Optional[float]]]] = {}
    for s in conn.execute("SELECT run_key, label, relpath, system FROM subruns"):
        T = (json.loads(s["system"]) or {}).get("T_target") if s["system"] else None
        subs.setdefault(s["run_key"], []).append((s["label"], s["relpath"], T))
    links: Dict[Tuple[str, str], Optional[str]] = {}
    for l in conn.execute("SELECT analysis_key, run_key, subrun FROM analysis_runs"):
        k = (l["analysis_key"], l["run_key"])
        links[k] = links.get(k) or l["subrun"]          # a sub-run link beats a whole-run link
    rows: Dict[str, List] = {}
    for r in conn.execute("SELECT * FROM analysis_results WHERE key != '__warning__'"):
        rows.setdefault(r["analysis_key"], []).append(r)
    for (akey, run_key), link_sub in links.items():
        run_subs = subs.get(run_key, [])
        real = [s for s in run_subs if s[0] != "(top)"]
        main = next((s[0] for s in run_subs if s[1] == "run"), real[0][0] if len(real) == 1 else None)
        for r in rows.get(akey, []):
            v = rx.Value(r["key"], None, subrun=r["subrun"], T=r["T"])
            label, note = target_subrun(v, None, run_subs, main, default=link_sub or "")
            conn.execute(
                "INSERT INTO results (run_key, subrun, key, value_num, value_text, unit, source, note) "
                "VALUES (?,?,?,?,?,?,?,?)",
                (run_key, label, r["key"], r["value_num"], r["value_text"], r["unit"],
                 f"analysis:{akey.split('/', 1)[1]}/{r['file']}", "; ".join(x for x in (r["note"], note) if x) or None))


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
        description, extra = describe_analysis(conn, key, a, cfg, warnings)
        ts = store.now()
        store.upsert_row(conn, "analysis", ["analysis_key"], {
            "analysis_key": key, "project": a.project, "name": a.name, "path": str(a.path),
            "hint": a.hint, "run_ids": ",".join(a.run_ids), "readme": text.strip() or None,
            "n_files": len(a.files), "total_bytes": sum(f.size for f in a.files),
            "fingerprint": fp, "first_seen": prev["first_seen"] if prev else ts, "last_scanned": ts,
            "missing": 0, "description": description})
        conn.execute("DELETE FROM analysis_files WHERE analysis_key=?", (key,))
        conn.executemany("INSERT INTO analysis_files VALUES (?,?,?,?,?)",
                         [(key, f.relpath, f.category, f.size, f.mtime) for f in a.files])
        if store.has_fts(conn):
            conn.execute("DELETE FROM fts WHERE run_key=? AND kind='analysis'", (key,))
        figs = " ".join(_name(f) for f in a.files if f.category in ("figure", "data", "analysis_script"))
        store.add_fts(conn, key, "analysis", f"{a.name} {a.hint} {text} {figs} {description or ''} {extra}")
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
    conn.execute("DELETE FROM jobs WHERE source = 'hook'")
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
                conn.execute("INSERT INTO jobs (run_key, subrun, analysis_key, job_id, role, source, detail) "
                             "VALUES (?,?,NULL,?,?,'hook',?)",
                             (rk, sub, slurm["SLURM_JOB_ID"], "run", f"exit {data.get('exit_code')}"))
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
            if slurm.get("SLURM_JOB_ID"):
                conn.execute("INSERT INTO jobs (run_key, subrun, analysis_key, job_id, role, source) "
                             "VALUES (?,?,?,?,?,'hook')",
                             (rk, sub, akey, slurm["SLURM_JOB_ID"], data.get("analysis_type") or "analysis"))


def update_sacct(conn, log=print) -> int:
    """Ask sacct about job ids not cached yet or not finished; returns how many were updated."""
    if not sacct.available():
        return 0
    cached = {r["job_id"]: r["state"] for r in conn.execute("SELECT job_id, state FROM sacct")}
    ids = [r[0] for r in conn.execute("SELECT DISTINCT job_id FROM jobs")
           if sacct.needs_query(cached.get(r[0]))]
    if not ids:
        return 0
    got = sacct.query(ids)
    for jid, rec in got.items():
        store.upsert_row(conn, "sacct", ["job_id"], {**{k: rec.get(k) for k in (
            "job_id", "job_name", "state", "state_detail", "elapsed_s", "start", "end", "nodes", "ncpus",
            "exit_code", "partition", "nodelist", "timelimit", "max_rss_kb")}, "queried": store.now()})
    return len(got)


SACCT_FAILED = {"FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED", "BOOT_FAIL", "DEADLINE"}


def apply_job_states(conn) -> None:
    """Where the files cannot tell (no end marker), let the scheduler's verdict decide.

    Uses the newest simulation job (not analysis jobs) of each run/sub-run."""
    rows = conn.execute(
        "SELECT j.run_key, COALESCE(j.subrun, '') sub, j.job_id, j.role, s.state, s.state_detail "
        "FROM jobs j JOIN sacct s USING (job_id) WHERE j.run_key IS NOT NULL").fetchall()
    newest: Dict[Tuple[str, str], sqlite3.Row] = {}
    for r in rows:
        if any(x in (r["role"] or "") for x in NON_SIM_ROLES):
            continue
        k = (r["run_key"], r["sub"])
        if k not in newest or int(r["job_id"]) > int(newest[k]["job_id"]):
            newest[k] = r
    touched = set()
    for (rk, sub), r in newest.items():
        run_subs = conn.execute("SELECT label, relpath, status FROM subruns WHERE run_key=?", (rk,)).fetchall()
        real = [s for s in run_subs if s["label"] != "(top)"]
        main = next((s for s in run_subs if s["relpath"] == "run"), real[0] if len(real) == 1 else None)
        target = main if sub == "" else next((s for s in run_subs if s["label"] == sub), None)
        status = target["status"] if target else conn.execute(
            "SELECT status FROM runs WHERE run_key=?", (rk,)).fetchone()[0]
        if status not in ("incomplete", "unknown", "not_started", "running"):
            continue
        state = r["state"]
        new = ("failed" if state in SACCT_FAILED else "running" if state == "RUNNING" else
               "pending" if state == "PENDING" else None)
        if new is None or new == status:
            continue
        evidence = f"sacct: job {r['job_id']} {r['state_detail']}"
        if target:
            conn.execute("UPDATE subruns SET status=?, status_evidence=? WHERE run_key=? AND label=?",
                         (new, evidence, rk, target["label"]))
        if not target or target is main or len(real) <= 1:
            conn.execute("UPDATE runs SET status=?, status_evidence=? WHERE run_key=?", (new, evidence, rk))
        else:
            touched.add(rk)
    for rk in touched:   # re-aggregate multi-sub-run runs
        sts = [s[0] for s in conn.execute("SELECT status FROM subruns WHERE run_key=? AND label != '(top)'", (rk,))]
        status, evidence = aggregate_status(sts)
        conn.execute("UPDATE runs SET status=?, status_evidence=? WHERE run_key=?", (status, evidence, rk))


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

def scan(cfg: Config, root: Path, full: bool = False, log=print, use_sacct: bool = True) -> Dict:
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

    # analysis dirs have no per-row parser version; reparse them all when parsing changed
    prev_pv = conn.execute("SELECT value FROM meta WHERE key='parser_version'").fetchone()
    reparse_analysis = full or prev_pv is None or prev_pv[0] != str(PARSER_VERSION)
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
            parsed = parse_run(run, cfg, run_warns, tree=root)
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
            if scan_analysis(conn, a, cfg, run_index, subrun_index, warnings, reparse_analysis) == "changed":
                counts["analysis_changed"] += 1
        for r in conn.execute("SELECT analysis_key FROM analysis WHERE project=? AND missing=0",
                              (project.name,)).fetchall():
            if r["analysis_key"] not in seen_a:
                conn.execute("UPDATE analysis SET missing=1 WHERE analysis_key=?", (r["analysis_key"],))
        conn.commit()

    apply_events(conn)
    materialize_analysis_results(conn)
    if use_sacct:
        counts["sacct"] = update_sacct(conn, log)
    apply_job_states(conn)
    unlinked = conn.execute("SELECT COUNT(*) FROM events WHERE run_key IS NULL").fetchone()[0]
    if unlinked:
        warnings.append(("inbox", f"{unlinked} hook events not matched to any run yet"))
    conn.execute("INSERT OR REPLACE INTO meta VALUES ('parser_version', ?)", (str(PARSER_VERSION),))
    conn.executemany("INSERT INTO scan_warnings VALUES (?,?,?)", [(scan_id, p, m) for p, m in warnings])
    counts["warnings"] = len(warnings)
    conn.execute("UPDATE scans SET finished=?, counts=? WHERE scan_id=?",
                 (store.now(), json.dumps(counts), scan_id))
    conn.commit()
    open_issues = [(r["run_key"], r["path"], r["message"]) for r in conn.execute(
        "SELECT w.* FROM run_warnings w JOIN runs r USING (run_key) WHERE r.missing = 0 ORDER BY run_key")]
    open_issues += [(r["analysis_key"], r["file"], r["value_text"]) for r in conn.execute(
        "SELECT x.* FROM analysis_results x JOIN analysis a USING (analysis_key) "
        "WHERE x.key = '__warning__' AND a.missing = 0 ORDER BY analysis_key")]
    return {"scan_id": scan_id, "counts": counts, "new": new_runs, "changed": changed_runs,
            "warnings": warnings, "open_issues": open_issues, "conn": conn}
