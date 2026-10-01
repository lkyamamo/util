"""Read-only `sacct` lookups for job ids the ledger has found.

One batched call per scan for ids not yet known or not yet finished;
answers are cached in the `sacct` table. Without `sacct` (e.g. on a laptop)
nothing is queried and cards show what the files say.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from typing import Dict, Iterable, List

FIELDS = ["JobIDRaw", "JobName", "State", "Elapsed", "Start", "End", "NNodes", "NCPUS", "ExitCode",
          "Partition", "MaxRSS", "NodeList", "Timelimit"]
FINAL_STATES = {"COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "OUT_OF_MEMORY", "NODE_FAIL", "PREEMPTED",
                "BOOT_FAIL", "DEADLINE"}


def available() -> bool:
    return shutil.which("sacct") is not None


def elapsed_seconds(s: str):
    """'D-HH:MM:SS', 'HH:MM:SS' or 'MM:SS' -> seconds."""
    if not s:
        return None
    days = 0
    if "-" in s:
        d, s = s.split("-", 1)
        days = int(d)
    parts = [int(float(p)) for p in s.split(":")]
    while len(parts) < 3:
        parts.insert(0, 0)
    h, m, sec = parts[-3:]
    return days * 86400 + h * 3600 + m * 60 + sec


def rss_kb(s: str):
    m = re.match(r"^([\d.]+)([KMGT]?)$", s or "")
    if not m:
        return None
    return float(m.group(1)) * {"": 1 / 1024, "K": 1, "M": 1024, "G": 1024 ** 2, "T": 1024 ** 3}[m.group(2)]


def parse(output: str) -> Dict[str, Dict]:
    jobs: Dict[str, Dict] = {}
    rss: Dict[str, float] = {}
    for line in output.splitlines():
        cols = line.split("|")
        if len(cols) != len(FIELDS):
            continue
        rec = dict(zip(FIELDS, cols))
        base, _, step = rec["JobIDRaw"].partition(".")
        kb = rss_kb(rec["MaxRSS"])
        if kb is not None:
            rss[base] = max(rss.get(base, 0), kb)
        if step:
            continue
        jobs[base] = {
            "job_id": base, "job_name": rec["JobName"], "state": rec["State"].split()[0] if rec["State"] else None,
            "state_detail": rec["State"], "elapsed_s": elapsed_seconds(rec["Elapsed"]),
            "start": None if rec["Start"] in ("Unknown", "None", "") else rec["Start"],
            "end": None if rec["End"] in ("Unknown", "None", "") else rec["End"],
            "nodes": int(rec["NNodes"]) if rec["NNodes"].isdigit() else None,
            "ncpus": int(rec["NCPUS"]) if rec["NCPUS"].isdigit() else None,
            "exit_code": rec["ExitCode"], "partition": rec["Partition"], "nodelist": rec["NodeList"],
            "timelimit": rec["Timelimit"],
        }
    for base, kb in rss.items():
        if base in jobs:
            jobs[base]["max_rss_kb"] = kb
    return jobs


def query(job_ids: Iterable[str], chunk: int = 200, timeout: int = 120) -> Dict[str, Dict]:
    ids = sorted({j for j in job_ids if j and j.isdigit()})
    out: Dict[str, Dict] = {}
    for i in range(0, len(ids), chunk):
        try:
            res = subprocess.run(["sacct", "-j", ",".join(ids[i:i + chunk]), "-P", "-n",
                                  "--format=" + ",".join(FIELDS)],
                                 capture_output=True, text=True, timeout=timeout)
        except (OSError, subprocess.SubprocessError):
            break
        if res.returncode == 0:
            out.update(parse(res.stdout))
    return out


def needs_query(cached_state) -> bool:
    return cached_state not in FINAL_STATES
