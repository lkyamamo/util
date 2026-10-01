"""Result extraction from analysis outputs (dielectric, diffusion, EOS, ring barriers, ...).

Each extractor reads one small file and returns values with enough context
(a sub-run label or a temperature) for the scanner to attach them to the
right run or sub-run. Values that the ledger computes rather than reads
(bulk modulus fits, RDF peaks) say so in their note.
"""

from __future__ import annotations

import csv
import io
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Dict, List, Optional, Tuple

from .parsers import read_head, to_num

MAX_RESULT_FILE = 2 * 1024 * 1024
EV_PER_A3_TO_GPA = 160.21766


@dataclass
class Value:
    key: str
    value: object
    unit: Optional[str] = None
    subrun: Optional[str] = None      # explicit sub-run label (e.g. a ring case)
    T: Optional[float] = None         # temperature (K) the value belongs to, for sweeps
    note: Optional[str] = None


@dataclass
class Extracted:
    kind: str
    values: List[Value] = field(default_factory=list)
    warnings: List[str] = field(default_factory=list)
    text: str = ""                    # extra searchable text (columns, commands, docstrings)


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _rows(text: str) -> Tuple[List[str], List[List[str]]]:
    reader = csv.reader(io.StringIO(text))
    rows = [r for r in reader if r]
    return (rows[0], rows[1:]) if rows else ([], [])


def _fit(xs: List[float], ys: List[float], degree: int) -> Optional[List[float]]:
    """Least-squares polynomial coefficients [c0, c1, (c2)] via normal equations."""
    n = degree + 1
    if len(xs) < n + 1:
        return None
    a = [[sum(x ** (i + j) for x in xs) for j in range(n)] for i in range(n)]
    b = [sum(y * x ** i for x, y in zip(xs, ys)) for i in range(n)]
    for i in range(n):              # Gaussian elimination with partial pivoting
        p = max(range(i, n), key=lambda r: abs(a[r][i]))
        a[i], a[p], b[i], b[p] = a[p], a[i], b[p], b[i]
        if abs(a[i][i]) < 1e-300:
            return None
        for r in range(i + 1, n):
            f = a[r][i] / a[i][i]
            a[r] = [x - f * y for x, y in zip(a[r], a[i])]
            b[r] -= f * b[i]
    c = [0.0] * n
    for i in reversed(range(n)):
        c[i] = (b[i] - sum(a[i][j] * c[j] for j in range(i + 1, n))) / a[i][i]
    return c


def traceback_warning(text: str, name: str) -> Optional[str]:
    if "Traceback (most recent call last)" not in text:
        return None
    last = [l.strip() for l in text.splitlines() if l.strip()]
    err = next((l for l in reversed(last) if re.match(r"^\w+(Error|Exception)\b", l)), last[-1] if last else "")
    return f"{name} contains a Python traceback (the script that wrote it failed): {err[:200]}"


# ---------------------------------------------------------------------------
# extractors
# ---------------------------------------------------------------------------

EPS_RE = re.compile(r"eps_x\s*=\s*([-\d.eE+]+),\s*eps_y\s*=\s*([-\d.eE+]+),\s*eps_z\s*=\s*([-\d.eE+]+),\s*"
                    r"eps_total\s*=\s*([-\d.eE+]+)")
DEV_RE = re.compile(r"deviation\s*=\s*([-\d.eE+]+)")
FRAMES_RE = re.compile(r"Processing complete:\s*([\d,]+)\s*frames")


def dielectric_summary(text: str, name: str) -> Optional[Extracted]:
    m = None
    for m in EPS_RE.finditer(text):
        pass
    if not m:
        return None
    out = Extracted("dielectric")
    for k, v in zip(("eps_x", "eps_y", "eps_z", "eps_total"), m.groups()):
        out.values.append(Value(k, to_num(v)))
    d = DEV_RE.findall(text)
    if d:
        out.values.append(Value("dipole_deviation", to_num(d[-1]), "e·Å"))
    f = FRAMES_RE.findall(text)
    if f:
        out.values.append(Value("dielectric_frames", int(f[-1].replace(",", ""))))
    return out


def vs_temperature(text: str, name: str) -> Optional[Extracted]:
    header, rows = _rows(text)
    tcol = next((c for c in ("temperature_K", "temperature", "T", "T_K") if c in header), None)
    ccol = "temperature_C" if "temperature_C" in header else None
    if not rows or not (tcol or ccol):
        return None
    out = Extracted(name.rsplit("_vs_", 1)[0], text=" ".join(header))
    unit = "1e-5 cm2/s" if name.startswith("diffusion") else None
    for row in rows:
        rec = dict(zip(header, row))
        T = to_num(rec.get(tcol)) if tcol else None
        if T is None and ccol:
            c = to_num(rec.get(ccol))
            T = c + 273.15 if c is not None else None
        for col, val in rec.items():
            if col in (tcol, ccol):
                continue
            v = to_num(val)
            if v is not None:
                out.values.append(Value(col, v, unit if col.startswith("D_") else None, T=T))
    return out


def diffusion_csv(text: str, name: str) -> Optional[Extracted]:
    header, rows = _rows(text)
    if len(header) != 2 or header[0] != "label" or not header[1].startswith("D"):
        return None
    m = re.match(r"D_(.+)$", header[1])
    unit = m.group(1).replace("_", " ").replace("cm2 s", "cm2/s") if m else None
    out = Extracted("diffusion")
    for label, val in rows:
        v = to_num(val)
        if v is not None:
            out.values.append(Value(f"D_{label}", v, unit))
    return out


def eos_summary(text: str, name: str) -> Optional[Extracted]:
    """Bulk modulus as plot_eos.py defines it, computed here from its input table."""
    header, rows = _rows(text)
    recs = [dict(zip(header, r)) for r in rows]
    out = Extracted("equation_of_state", text=" ".join(header))
    vcol = "nve_volume_A3" if "nve_volume_A3" in header else "volume_A3"
    pts = [(to_num(r.get(vcol)), to_num(r.get("nve_Press_mean_GPa"))) for r in recs]
    pts = [(v, p) for v, p in pts if v is not None and p is not None]
    if len(pts) >= 3:
        c = _fit([v for v, _ in pts], [p for _, p in pts], 1)
        if c and c[1] != 0:
            v0 = -c[0] / c[1]
            out.values += [Value("bulk_modulus_P", -v0 * c[1], "GPa",
                                 note="derived: B = -V0·dP/dV, linear fit of nve_Press_mean_GPa vs volume"),
                           Value("V0_P", v0, "Å³", note="derived: volume where the P(V) fit crosses 0")]
    ept = [(to_num(r.get("volume_A3")), to_num(r.get("cg_TotEng_eV"))) for r in recs]
    ept = [(v, e) for v, e in ept if v is not None and e is not None]
    if len(ept) >= 4:
        c = _fit([v for v, _ in ept], [e for _, e in ept], 2)
        if c and c[2] > 0:
            v0 = -c[1] / (2 * c[2])
            out.values += [Value("bulk_modulus_E0K", 2 * c[2] * v0 * EV_PER_A3_TO_GPA, "GPa",
                                 note="derived: B = V0·d²E/dV², parabola fit of cg_TotEng_eV vs volume"),
                           Value("V0_E0K", v0, "Å³", note="derived: minimum of the E(V) parabola")]
    out.values.append(Value("eos_points", len(recs)))
    return out if out.values else None


SURVEY_ROW_RE = re.compile(r"^\s*(\d+)\s+(\d+)\s*\|\s*([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)\s*\|"
                           r"\s*([-\d.]+)\s*eV\s*\|\s*([-\d.]+)\s*A")
CASE_ROW_RE = re.compile(r"^\s*(\S+)\s+(\d+)\s+([-\d.]+)\s+(\d+)\s+([-+\d.]+)\s+(\d+)\s*$")


def ring_summary(text: str, name: str) -> Optional[Extracted]:
    out = Extracted("ring_barrier")
    m = re.search(r"^(\d+) cases, (\d+) of them with energies", text, re.M)
    if m:
        out.values += [Value("cases", int(m.group(1))), Value("cases_with_energy", int(m.group(2)))]
    for line in text.splitlines():
        m = SURVEY_ROW_RE.match(line)
        if m:
            n = m.group(1)
            out.values += [Value(f"rings_n{n}", int(m.group(2))),
                           Value(f"barrier_p10_n{n}", to_num(m.group(3)), "eV"),
                           Value(f"barrier_median_n{n}", to_num(m.group(4)), "eV"),
                           Value(f"barrier_p90_n{n}", to_num(m.group(5)), "eV"),
                           Value(f"barrier_min_n{n}", to_num(m.group(6)), "eV")]
            continue
        m = CASE_ROW_RE.match(line)
        if m and m.group(1) != "case":
            case = m.group(1)
            out.values += [Value("barrier", to_num(m.group(3)), "eV", subrun=case),
                           Value("barrier_frame", int(m.group(4)), subrun=case),
                           Value("tail_dE", to_num(m.group(5)), "eV", subrun=case),
                           Value("frames_unconverged", int(m.group(6)), subrun=case)]
    for m in re.finditer(r"^\s*(ring size n|aperture)\s*:\s*r\s*=\s*([-+\d.]+)", text, re.M):
        out.values.append(Value(f"barrier_corr_{m.group(1).split()[-1]}", to_num(m.group(2))))
    return out if out.values else None


def polars_table(text: str, name: str) -> Optional[Extracted]:
    """A one-row polars DataFrame printout (thermodynamics.py output.txt)."""
    rows = [l for l in text.splitlines() if "│" in l]
    if len(rows) < 4:
        return None
    cells = [[c.strip() for c in r.strip().strip("│").split("┆")] for r in rows]
    header, values = cells[0], cells[-1]
    if len(header) != len(values) or len(re.findall(r"shape: \((\d+),", text)) != 1 or \
            re.search(r"shape: \(1,", text) is None:
        return None
    out = Extracted("thermo_average")
    for k, v in zip(header, values):
        n = to_num(v)
        if n is not None:
            out.values.append(Value(f"avg_{k}", n))
    return out if out.values else None


def rdf_csv(text: str, name: str) -> Optional[Extracted]:
    header, rows = _rows(text)
    if len(header) < 2 or not header[0].lower().startswith("r"):
        return None
    data = [[to_num(x) for x in r] for r in rows]
    data = [r for r in data if len(r) == len(header) and all(v is not None for v in r)]
    if len(data) < 10:
        return None
    out = Extracted("rdf", text=" ".join(header))
    for j, col in enumerate(header[1:], 1):
        if "-" not in col:          # element pairs only (skip total/neutron weights)
            continue
        i = max(range(len(data)), key=lambda k: data[k][j])
        out.values += [Value(f"rdf_peak_r_{col}", data[i][0], "Å", note="derived: position of the g(r) maximum"),
                       Value(f"rdf_peak_g_{col}", data[i][j], note="derived: g(r) maximum")]
    return out if out.values else None


def generic_csv(text: str, name: str) -> Optional[Extracted]:
    """Any other small CSV: columns for search; values when it has a single row."""
    header, rows = _rows(text)
    if not header:
        return None
    out = Extracted("table", text=f"{name}: " + " ".join(header))
    if len(rows) == 1:
        for k, v in zip(header, rows[0]):
            n = to_num(v)
            if n is not None:
                out.values.append(Value(f"{Path(name).stem}.{k}", n))
    return out


# (kind, filename predicate, extractor), first match wins
EXTRACTORS: List[Tuple[str, Callable[[str], bool], Callable[[str, str], Optional[Extracted]]]] = [
    ("dielectric", lambda n: n == "summary.txt", dielectric_summary),
    ("vs_temperature", lambda n: n.endswith("_vs_temperature.csv"), vs_temperature),
    ("diffusion", lambda n: n.endswith("diffusion.csv"), diffusion_csv),
    ("eos", lambda n: n == "eos_summary.csv", eos_summary),
    ("ring_barrier", lambda n: n == "SUMMARY.txt", ring_summary),
    ("thermo_average", lambda n: n == "output.txt", polars_table),
    ("rdf", lambda n: n.endswith(".csv") and "rdf" in n.lower(), rdf_csv),
    ("table", lambda n: n.endswith(".csv"), generic_csv),
]

TEXT_RESULT_NAMES = ("summary.txt", "SUMMARY.txt", "output.txt")


def is_candidate(name: str) -> bool:
    return any(pred(name) for _, pred, _ in EXTRACTORS)


def extract(path: Path, size: int) -> Optional[Extracted]:
    name = path.name
    if size > MAX_RESULT_FILE:
        return None
    for kind, pred, fn in EXTRACTORS:
        if not pred(name):
            continue
        text = read_head(path, MAX_RESULT_FILE)
        warn = traceback_warning(text, name) if name in TEXT_RESULT_NAMES else None
        try:
            got = fn(text, name)
        except (ValueError, IndexError, ZeroDivisionError):
            got = None
        if warn:
            got = got or Extracted(kind)
            got.warnings.append(warn)
        if got is not None:
            return got
    return None


# ---------------------------------------------------------------------------
# describing analysis directories
# ---------------------------------------------------------------------------

def script_description(text: str) -> Optional[str]:
    """First sentence/line of a Python module docstring or leading comment block."""
    body = text.lstrip()
    if body.startswith("#!"):
        body = body.split("\n", 1)[1].lstrip() if "\n" in body else ""
    m = re.match(r'^(?:[rubRUB]{0,2})("""|\'\'\')(.*?)\1', body, re.S)
    if m:
        para = m.group(2).strip().split("\n\n", 1)[0]
        return " ".join(para.split())[:300] or None
    comments = []
    for line in body.splitlines():
        if line.startswith("#"):
            comments.append(line.lstrip("#").strip())
        elif comments or line.strip():
            break
    return " ".join(c for c in comments if c)[:300] or None
