"""Job-time events: written by `simledger hook` from SLURM scripts, ingested by `scan`.

Each event is one JSON file, written to a temp name and renamed into
<home>/inbox so concurrent jobs never collide and a scan never sees a
half-written file. Hooks only append files; they never open the database.
"""

from __future__ import annotations

import json
import os
import shutil
import socket
import subprocess
import uuid
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from . import store
from .config import Config
from .discover import locate

SLURM_VARS = ("SLURM_JOB_ID", "SLURM_JOB_NAME", "SLURM_JOB_PARTITION", "SLURM_JOB_NODELIST",
              "SLURM_JOB_NUM_NODES", "SLURM_NTASKS", "SLURM_CPUS_PER_TASK", "SLURM_SUBMIT_DIR",
              "SLURM_JOB_ACCOUNT", "SLURM_ARRAY_JOB_ID", "SLURM_ARRAY_TASK_ID",
              "SLURM_JOB_START_TIME", "SLURM_JOB_END_TIME")


def _batch_script(job_id: Optional[str], script: Optional[str], max_bytes: int) -> Tuple[Optional[str], Optional[str]]:
    """Text of the submitted script.

    Under sbatch, $0 is slurmd's spool copy, so prefer asking SLURM for the
    script; fall back to an explicitly given path."""
    if job_id:
        try:
            res = subprocess.run(["scontrol", "write", "batch_script", job_id, "-"],
                                 capture_output=True, text=True, timeout=15)
            if res.returncode == 0 and res.stdout.strip():
                return res.stdout[:max_bytes], "scontrol"
        except (OSError, subprocess.SubprocessError):
            pass
    if script and Path(script).is_file():
        try:
            return Path(script).read_text(errors="replace")[:max_bytes], str(Path(script).resolve())
        except OSError:
            pass
    return None, None


def write_event(cfg: Config, kind: str, data: Dict) -> Path:
    cfg.inbox.mkdir(parents=True, exist_ok=True)
    created = datetime.now()
    job = os.environ.get("SLURM_JOB_ID", str(os.getpid()))
    name = f"{created:%Y%m%dT%H%M%S}_{kind}_{job}_{uuid.uuid4().hex[:8]}.json"
    event = {"kind": kind, "created": created.isoformat(timespec="seconds"),
             "host": socket.gethostname(), "user": os.environ.get("USER"),
             "slurm": {k: os.environ[k] for k in SLURM_VARS if k in os.environ}, **data}
    tmp = cfg.inbox / f".{name}.tmp"
    with tmp.open("w") as f:
        json.dump(event, f, indent=1)
        f.write("\n")
    os.replace(tmp, cfg.inbox / name)
    return cfg.inbox / name


def hook_run(cfg: Config, run_dir: str, exit_code: Optional[int], script: Optional[str],
             fields: Dict[str, str], note: str) -> Path:
    job_id = os.environ.get("SLURM_JOB_ID")
    text, source = _batch_script(job_id, script, cfg.max_text_bytes)
    return write_event(cfg, "run", {
        "dir": str(Path(run_dir).absolute()),
        "exit_code": exit_code,
        "script_text": text,
        "script_source": source,
        "fields": fields,
        "note": note,
    })


def hook_analysis(cfg: Config, data_dir: Optional[str], run: Optional[str], analysis_type: str,
                  output: Optional[str], fields: Dict[str, str], note: str) -> Path:
    return write_event(cfg, "analysis", {
        "dir": str(Path(data_dir).absolute()) if data_dir else None,
        "run": run,
        "analysis_type": analysis_type,
        "output": str(Path(output).absolute()) if output else None,
        "fields": fields,
        "note": note,
    })


def ingest(conn, cfg: Config, warnings: List[Tuple[str, str]]) -> int:
    """Move inbox/*.json into the events table (archived under inbox/processed/)."""
    if not cfg.inbox.is_dir():
        return 0
    processed = cfg.inbox / "processed" / datetime.now().strftime("%Y%m")
    n = 0
    for path in sorted(cfg.inbox.glob("*.json")):
        try:
            event = json.loads(path.read_text())
        except (OSError, ValueError) as exc:
            bad = cfg.inbox / "bad"
            bad.mkdir(exist_ok=True)
            shutil.move(str(path), str(bad / path.name))
            warnings.append((str(path), f"unreadable inbox event moved to bad/: {exc}"))
            continue
        loc = locate(event["dir"], cfg) if event.get("dir") else None
        if loc is None and event.get("run"):
            proj, _, rid = event["run"].rpartition("/")
            loc = {"project": proj or None, "run_id": cfg.match_id(rid) or rid, "subrun": ""}
        conn.execute(
            "INSERT OR IGNORE INTO events VALUES (?,?,?,?,?,?,?,?,NULL)",
            (path.name, event.get("kind"), event.get("created"), store.now(), json.dumps(event),
             loc and loc["project"], loc and loc["run_id"], loc["subrun"] if loc else None),
        )
        if loc is None:
            warnings.append((event.get("dir") or str(path), "inbox event does not point inside <project>/runs/<id>"))
        processed.mkdir(parents=True, exist_ok=True)
        shutil.move(str(path), str(processed / path.name))
        n += 1
    return n


def link_events(conn) -> None:
    """Attach events to runs by project + run id, or by run id alone if unique."""
    conn.execute("""
        UPDATE events SET run_key = project || '/' || run_id
        WHERE run_key IS NULL AND project IS NOT NULL
          AND EXISTS (SELECT 1 FROM runs r WHERE r.run_key = events.project || '/' || events.run_id)""")
    conn.execute("""
        UPDATE events SET run_key = (SELECT r.run_key FROM runs r WHERE r.run_id = events.run_id)
        WHERE run_key IS NULL AND run_id IS NOT NULL
          AND (SELECT COUNT(*) FROM runs r WHERE r.run_id = events.run_id) = 1""")
