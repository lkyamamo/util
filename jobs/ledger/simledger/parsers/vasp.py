"""VASP inputs (INCAR, KPOINTS, POSCAR, POTCAR titles) and outputs (OUTCAR/OSZICAR tails).

POTCAR content is never stored; only TITEL/VRHFIN/ZVAL/ENMAX values are kept.
"""

from __future__ import annotations

import math
import re
from typing import Dict, Iterable, List, Optional

from . import to_num

KEY_TAGS = ("ENCUT", "PREC", "ISMEAR", "SIGMA", "EDIFF", "EDIFFG", "IBRION", "ISIF", "NSW",
            "POTIM", "ISPIN", "MAGMOM", "GGA", "METAGGA", "LDAU", "IVDW", "LHFCALC", "HFSCREEN",
            "MDALGO", "TEBEG", "TEEND", "SMASS", "NPAR", "NCORE", "KPAR", "ALGO", "LREAL",
            "ICHAIN", "IMAGES", "LEPSILON", "LCALCEPS", "LORBIT", "ICHARG", "ISTART")


def parse_incar(text: str) -> Dict[str, str]:
    tags: Dict[str, str] = {}
    for raw in text.splitlines():
        line = re.split(r"[!#]", raw, 1)[0]
        for stmt in line.split(";"):
            if "=" not in stmt:
                continue
            key, _, val = stmt.partition("=")
            key = key.strip().upper()
            if key:
                tags[key] = val.strip()
    return tags


def _tag_num(tags: Dict[str, str], key: str) -> Optional[float]:
    return to_num(tags[key].split()[0]) if tags.get(key) else None


def calc_type(tags: Dict[str, str]) -> str:
    ibrion = _tag_num(tags, "IBRION")
    nsw = _tag_num(tags, "NSW") or 0
    isif = _tag_num(tags, "ISIF")
    if tags.get("IMAGES") or tags.get("ICHAIN"):
        return "neb"
    if ibrion is not None and ibrion in (5, 6, 7, 8):
        return "phonon"
    if ibrion == 0:
        mdalgo = _tag_num(tags, "MDALGO")
        smass = _tag_num(tags, "SMASS")
        if isif is not None and isif >= 3:
            return "md_npt"
        if mdalgo in (1, 2, 3, 5) or (mdalgo is None and smass is not None and smass >= 0):
            return "md_nvt"
        return "md_nve"   # default SMASS = -3 is microcanonical
    if tags.get("LEPSILON", "").upper().startswith((".T", "T")) or \
            tags.get("LCALCEPS", "").upper().startswith((".T", "T")):
        return "dielectric"
    if nsw > 0 and ibrion in (1, 2, 3):
        return "relax_cell" if (isif or 2) >= 3 else "relax_ions"
    if _tag_num(tags, "ICHARG") == 11:
        return "nscf"
    return "static"


def parse_kpoints(text: str) -> Dict:
    lines = [l.strip() for l in text.splitlines()]
    out: Dict = {"comment": lines[0] if lines else None}
    if len(lines) < 3:
        return out
    n = to_num(lines[1].split()[0]) if lines[1] else None
    mode = lines[2][:1].lower() if lines[2] else ""
    if n == 0 and mode in ("g", "m"):
        out["scheme"] = "Gamma" if mode == "g" else "Monkhorst-Pack"
        if len(lines) > 3:
            out["mesh"] = lines[3]
        if len(lines) > 4 and lines[4]:
            out["shift"] = lines[4]
    elif mode == "l":
        out["scheme"] = "line"
        out["n_per_segment"] = n
    elif n == 0 and mode == "a":
        out["scheme"] = "auto"
        out["mesh"] = lines[3] if len(lines) > 3 else None
    else:
        out["scheme"] = "explicit"
        out["n_kpoints"] = n
    return out


def parse_poscar(text: str) -> Dict:
    lines = text.splitlines()
    if len(lines) < 7:
        return {}
    try:
        scale = float(lines[1].split()[0])
        vecs = [[float(x) for x in lines[i].split()[:3]] for i in (2, 3, 4)]
    except (ValueError, IndexError):
        return {}
    species_line = lines[5].split()
    if all(to_num(x) is not None for x in species_line):
        species, counts = [], [int(x) for x in species_line]   # VASP 4 format
    else:
        species = species_line
        try:
            counts = [int(x) for x in lines[6].split()]
        except ValueError:
            return {"species": species}
    if scale < 0:   # negative scale = target volume
        raw_a, raw_b, raw_c = vecs
        raw_vol = abs(raw_a[0] * (raw_b[1] * raw_c[2] - raw_b[2] * raw_c[1]) - raw_a[1] * (raw_b[0] * raw_c[2] - raw_b[2] * raw_c[0])
                      + raw_a[2] * (raw_b[0] * raw_c[1] - raw_b[1] * raw_c[0]))
        scale = (-scale / raw_vol) ** (1 / 3)
    vecs = [[scale * x for x in v] for v in vecs]
    a, b, c = vecs
    volume = abs(a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
                 + a[2] * (b[0] * c[1] - b[1] * c[0]))
    lengths = [sum(x * x for x in v) ** 0.5 for v in vecs]

    def angle(u, v, lu, lv):
        cosv = sum(x * y for x, y in zip(u, v)) / (lu * lv)
        return math.degrees(math.acos(max(-1.0, min(1.0, cosv))))
    angles = [angle(b, c, lengths[1], lengths[2]), angle(a, c, lengths[0], lengths[2]),
              angle(a, b, lengths[0], lengths[1])]
    rest = lines[7].strip().lower() if len(lines) > 7 else ""
    formula = "".join(f"{s}{n}" for s, n in zip(species, counts)) if species else None
    return {"species": species, "counts": counts, "n_atoms": sum(counts), "formula": formula,
            "lattice_lengths": [round(x, 4) for x in lengths], "lattice_angles": [round(x, 2) for x in angles],
            "volume": round(volume, 4),
            "selective_dynamics": rest.startswith("s")}


def potcar_titles(lines: Iterable[str]) -> List[Dict]:
    out: List[Dict] = []
    cur: Dict = {}
    for line in lines:
        s = line.strip()
        # VRHFIN precedes TITEL within each element's header; either may open an entry.
        if s.startswith("VRHFIN"):
            cur = {"vrhfin": s.split("=", 1)[1].split(":")[0].strip()}
            out.append(cur)
        elif s.startswith("TITEL"):
            if not cur or "titel" in cur:
                cur = {}
                out.append(cur)
            cur["titel"] = s.split("=", 1)[1].strip()
        elif s.startswith("POMASS") and cur:
            m = re.search(r"ZVAL\s*=\s*([\d.]+)", s)
            if m:
                cur["zval"] = float(m.group(1))
            m = re.search(r"POMASS\s*=\s*([\d.]+)", s)
            if m:
                cur["pomass"] = float(m.group(1))
        elif s.startswith("ENMAX") and cur:
            m = re.search(r"ENMAX\s*=\s*([\d.]+)", s)
            if m:
                cur["enmax"] = float(m.group(1))
    return out


E0_RE = re.compile(r"energy\(sigma->0\)\s*=\s*([-\d.Ee+]+)")
TOTEN_RE = re.compile(r"free\s+energy\s+TOTEN\s*=\s*([-\d.Ee+]+)")
EFERMI_RE = re.compile(r"E-fermi\s*:\s*([-\d.Ee+]+)")
PRESSURE_RE = re.compile(r"external pressure =\s*([-\d.Ee+]+)")
ELAPSED_RE = re.compile(r"Elapsed time \(sec\):\s*([\d.]+)")
VERSION_RE = re.compile(r"^\s*(vasp\.\S+)")
CORES_RE = re.compile(r"running (?:on\s+)?(\d+) total cores|running\s+(\d+)\s+mpi-ranks")


def _last(rx, text) -> Optional[float]:
    vals = rx.findall(text)
    if not vals:
        return None
    v = vals[-1]
    if isinstance(v, tuple):
        v = next(x for x in v if x)
    return to_num(v)


def parse_outcar(head: str, tail: str) -> Dict:
    out: Dict = {}
    for line in head.splitlines()[:20]:
        m = VERSION_RE.match(line)
        if m:
            out["version"] = m.group(1)
            break
    m = CORES_RE.search(head)
    if m:
        out["cores"] = int(m.group(1) or m.group(2))
    out["finished"] = "General timing and accounting" in tail
    out["reached_accuracy"] = "reached required accuracy" in tail
    out["energy_sigma0"] = _last(E0_RE, tail)
    out["toten"] = _last(TOTEN_RE, tail)
    out["e_fermi"] = _last(EFERMI_RE, tail)
    out["pressure_kB"] = _last(PRESSURE_RE, tail)
    out["elapsed_s"] = _last(ELAPSED_RE, tail)
    return out


OSZ_IONIC_RE = re.compile(r"^\s*(\d+)\s+(?:T=\s*([\d.]+)\s+)?.*?F=\s*([-\d.Ee+]+)\s+E0=\s*([-\d.Ee+]+)")


def parse_oszicar(tail: str) -> Dict:
    steps = [m for m in map(OSZ_IONIC_RE.match, tail.splitlines()) if m]
    if not steps:
        return {}
    last = steps[-1]
    out = {"ionic_steps": int(last.group(1)), "F": to_num(last.group(3)), "E0": to_num(last.group(4))}
    if last.group(2):
        out["T"] = to_num(last.group(2))
    return out


GGA_NAMES = {"PE": "PBE", "PS": "PBEsol", "RP": "RPBE", "91": "PW91", "AM": "AM05", "RE": "revPBE",
             "B3": "B3LYP", "BO": "optB86b", "MK": "optB88", "OR": "optPBE", "ML": "vdW-DF2"}
IVDW_NAMES = {"1": "D2", "10": "D2", "11": "D3", "12": "D3(BJ)", "13": "D4", "2": "TS", "20": "TS",
              "21": "TS-SCS", "202": "MBD", "4": "dDsC"}


def functional(tags: Dict[str, str], potcar: List[Dict]) -> Optional[str]:
    meta = tags.get("METAGGA", "").strip().upper()
    gga = tags.get("GGA", "").strip().upper()
    if tags.get("LHFCALC", "").upper().startswith((".T", "T")):
        hf = to_num(tags.get("HFSCREEN", "0").split()[0]) or 0
        name = "HSE06" if abs(hf - 0.2) < 1e-3 else ("HSE03" if abs(hf - 0.3) < 1e-3 else
                                                     ("PBE0" if hf == 0 else f"hybrid(HFSCREEN={hf})"))
    elif meta:
        name = {"R2SCAN": "r2SCAN", "SCAN": "SCAN", "RTPSS": "revTPSS", "TPSS": "TPSS"}.get(meta, meta)
    elif gga:
        name = GGA_NAMES.get(gga, gga)
    elif potcar and all("PBE" in t.get("titel", "") for t in potcar):
        name = "PBE"
    elif potcar:
        name = "LDA"
    else:
        return None
    ivdw = tags.get("IVDW", "").split()[0] if tags.get("IVDW") else ""
    if ivdw and ivdw != "0":
        name += f"+{IVDW_NAMES.get(ivdw, 'IVDW' + ivdw)}"
    if tags.get("LDAU", "").upper().startswith((".T", "T")):
        name += "+U"
    return name


def md_settings(tags: Dict[str, str]) -> Dict:
    """Temperatures, timestep and length of an IBRION=0 run."""
    if _tag_num(tags, "IBRION") != 0:
        return {}
    potim = _tag_num(tags, "POTIM")
    nsw = _tag_num(tags, "NSW")
    out = {"T_start": _tag_num(tags, "TEBEG"), "T_end": _tag_num(tags, "TEEND") or _tag_num(tags, "TEBEG"),
           "potim_fs": potim, "nsw": int(nsw) if nsw else None}
    if potim and nsw:
        out["time_ps"] = potim * nsw / 1000
    return out


def oszicar_md(lines: Iterable[str]) -> Dict:
    """Mean temperature over the second half of an MD OSZICAR."""
    temps = []
    for line in lines:
        m = re.search(r"\bT=\s*([\d.]+)", line)
        if m:
            temps.append(float(m.group(1)))
    if not temps:
        return {}
    half = temps[len(temps) // 2:]
    mean = sum(half) / len(half)
    return {"md_steps": len(temps), "T_mean": mean,
            "T_std": (sum((t - mean) ** 2 for t in half) / len(half)) ** 0.5}
