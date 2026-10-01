"""Turn parsed LAMMPS/VASP data into a description of the simulated system.

Each builder returns a `system` dict (stored as JSON and as searchable
params), a `protocol` list (LAMMPS run segments), a one-line `summary` and a
short `conditions` string for sub-run tables.
"""

from __future__ import annotations

import re
from typing import Dict, List, Optional, Tuple

from .parsers import lammps, vasp

AMU_PER_A3_TO_G_PER_CM3 = 1.66053907
ELEMENT_RE = re.compile(r"^[A-Z][a-z]?$")

TIME_UNIT = {"metal": "ps", "real": "fs", "lj": "tau", "si": "s", "cgs": "s", "electron": "fs",
             "micro": "us", "nano": "ns"}
PRESSURE_UNIT = {"metal": "bar", "real": "atm", "lj": "", "si": "Pa", "cgs": "dyn/cm2",
                 "electron": "Pa", "micro": "pg/(um ns2)", "nano": "ag/(nm ns2)"}
LENGTH_UNIT = {"metal": "Å", "real": "Å", "lj": "σ", "si": "m", "cgs": "cm", "electron": "bohr",
               "micro": "um", "nano": "nm"}
VASP_TYPES = {"static": "single point", "relax_ions": "ionic relaxation", "relax_cell": "cell relaxation",
              "md_nvt": "NVT MD", "md_npt": "NPT MD", "md_nve": "NVE MD", "neb": "NEB", "phonon": "phonons",
              "dielectric": "dielectric (DFPT)", "nscf": "non-SCF"}


def _int(x: Optional[float]) -> Optional[int]:
    return int(x) if x is not None else None


def _g(x: Optional[float], digits: int = 4) -> str:
    if x is None:
        return "?"
    return f"{x:.{digits}g}"


def fmt_time(value: Optional[float], unit: str) -> str:
    if value is None:
        return "?"
    if unit == "fs" and value >= 1000:
        value, unit = value / 1000, "ps"
    if unit == "ps" and value >= 1000:
        value, unit = value / 1000, "ns"
    return f"{value:.4g} {unit}"


def fmt_box(lengths: List[float], unit: str = "Å", angles: Optional[List[float]] = None) -> str:
    if not lengths:
        return "?"
    a, b, c = lengths
    ortho = not angles or all(abs(x - 90) < 0.01 for x in angles)
    if ortho and max(lengths) - min(lengths) < 1e-3:
        return f"{a:.4g} {unit} cubic"
    s = f"{a:.4g} × {b:.4g} × {c:.4g} {unit}"
    if not ortho:
        s += " (α,β,γ = " + ", ".join(f"{x:.4g}" for x in angles) + "°)"
    return s


def fmt_range(pair, unit: str) -> Optional[str]:
    if not pair:
        return None
    a, b = pair
    return f"{a:g} {unit}" if a == b else f"{a:g}→{b:g} {unit}"


def composition_text(comp: Dict[str, int]) -> str:
    return ", ".join(f"{el} {n}" for el, n in comp.items())


# ---------------------------------------------------------------------------
# LAMMPS
# ---------------------------------------------------------------------------

def _pair_elements(pair_coeff: List[str], n_types: Optional[int]) -> Optional[List[str]]:
    """Elements from `pair_coeff * * file E1 E2 ...` (manybody styles map types in order)."""
    for pc in reversed(pair_coeff):
        toks = pc.split()
        if len(toks) > 3 and toks[0] == "*" and toks[1] == "*":
            els = toks[3:]
            if all(ELEMENT_RE.match(e) or e == "NULL" for e in els) and (n_types is None or len(els) == n_types):
                return els
    return None


def type_elements(log: Dict, data: Dict) -> Tuple[Dict[str, str], Dict[str, float]]:
    """Element and mass for each atom type, from the best source available."""
    masses: Dict[str, float] = {}
    labels: Dict[str, str] = {}
    for src in (data.get("masses") or {}, log.get("masses") or {}):   # log `mass` commands win
        for t, m in src.items():
            masses[t] = m["mass"]
            lab = (m.get("label") or "").split()
            if lab and ELEMENT_RE.match(lab[0]):
                labels[t] = lab[0]
    n_types = data.get("n_types") or (len(masses) or None)
    pair_els = _pair_elements(log.get("pair_coeff") or [], n_types)
    types = sorted(set(masses) | set((data.get("type_counts") or {})) |
                   ({str(i + 1) for i in range(len(pair_els))} if pair_els else set()), key=lambda t: int(t) if t.isdigit() else 0)
    out = {}
    for t in types:
        el = labels.get(t)
        if not el and pair_els and t.isdigit() and int(t) <= len(pair_els) and pair_els[int(t) - 1] != "NULL":
            el = pair_els[int(t) - 1]
        if not el and t in masses:
            el = lammps.guess_element(masses[t])
        out[t] = el or f"type{t}"
        if t not in masses and el in lammps.ELEMENT_MASSES:
            masses[t] = lammps.ELEMENT_MASSES[el]
    return out, masses


def lammps_system(log: Dict, data: Dict, data_file: Optional[str]) -> Dict:
    units = log.get("units") or "metal"
    tu, pu, lu = TIME_UNIT.get(units, ""), PRESSURE_UNIT.get(units, ""), LENGTH_UNIT.get(units, "")
    elements, masses = type_elements(log, data)
    counts = dict(data.get("type_counts") or {})
    structures = list(dict.fromkeys(f.rsplit("/", 1)[-1] if "/" not in f else f for f in log.get("data_files") or []))
    multi = len(structures) > 1
    # several structures read one after another: describe the first (the data file read)
    n_atoms = (data.get("n_atoms") or log.get("n_atoms")) if multi else (log.get("n_atoms") or data.get("n_atoms"))
    # `replicate` multiplies the data file's atoms; scale per-type counts to what ran
    if not multi and counts and data.get("n_atoms") and n_atoms and n_atoms != data["n_atoms"] and not log.get("atoms_changed"):
        factor = n_atoms / data["n_atoms"]
        if abs(factor - round(factor)) < 1e-9 and (not log.get("replicate") or log["replicate"] == round(factor)):
            counts = {t: n * int(round(factor)) for t, n in counts.items()}
        else:
            counts = {}
    comp: Dict[str, int] = {}
    for t, n in sorted(counts.items(), key=lambda kv: int(kv[0]) if kv[0].isdigit() else 0):
        comp[elements.get(t, f"type{t}")] = comp.get(elements.get(t, f"type{t}"), 0) + n

    box = log.get("box_at_run") or log.get("box") or ({"lengths": [data["box"].get(k) for k in "xyz"], "tilt": data.get("tilt")}
                             if len(data.get("box") or {}) == 3 else None)
    lengths = box["lengths"] if box else None
    volume = lengths[0] * lengths[1] * lengths[2] if lengths else None
    total_mass = sum(masses.get(t, 0) * n for t, n in counts.items()) if counts and all(t in masses for t in counts) else None
    density = (total_mass * AMU_PER_A3_TO_G_PER_CM3 / volume) if (total_mass and volume and lu == "Å") else None

    protocol = []
    total_steps, total_time = 0, 0.0
    for seg in lammps.collapse_segments(log.get("segments") or []):
        steps = seg.get("steps") or 0
        dt = seg.get("timestep")
        th = seg.get("thermo") or {}
        row = {"count": seg["count"], "kind": seg["kind"], "ensemble": seg["ensemble"],
               "T": seg.get("T"), "P": seg.get("P"), "timestep": dt, "steps": seg.get("steps"),
               "time": steps * dt if (seg["kind"] == "run" and dt) else None,
               "T_mean": th.get("Temp", {}).get("mean"), "T_std": th.get("Temp", {}).get("std"),
               "P_mean": th.get("Press", {}).get("mean"),
               "volume_mean": th.get("Volume", {}).get("mean"),
               "density_mean": th.get("Density", {}).get("mean"),
               "fixes": seg.get("fixes")}
        if row["density_mean"] is None and row["volume_mean"] and total_mass and lu == "Å":
            row["density_mean"] = total_mass * AMU_PER_A3_TO_G_PER_CM3 / row["volume_mean"]
        protocol.append(row)
        if seg["kind"] == "run":
            total_steps += steps * seg["count"]
            total_time += (row["time"] or 0) * seg["count"]
    runs = [r for r in protocol if r["kind"] == "run"]
    # the stage with the most steps; ties go to the later one (production usually comes last)
    main = max(reversed(runs), key=lambda r: (r["steps"] or 0) * r["count"], default=None)

    pot_files = sorted({pc.split()[2] for pc in log.get("pair_coeff") or [] if len(pc.split()) > 2
                        and pc.split()[0] == "*" and "." in pc.split()[2]})
    system = {
        "code": "lammps", "units": units, "atom_style": log.get("atom_style") or data.get("atoms_style"),
        "n_atoms": n_atoms, "n_types": data.get("n_types") or (len(elements) or None),
        "elements": " ".join(elements[t] for t in elements) or None,
        "composition": composition_text(comp) or None,
        "formula": "".join(f"{el}{n}" for el, n in comp.items()) or None,
        "box": fmt_box(lengths, lu) if lengths else None,
        "box_lengths": [round(x, 4) for x in lengths] if lengths else None,
        "box_tilt": box.get("tilt") if box else None,
        "volume": round(volume, 3) if volume else None,
        "density": round(density, 4) if density else None,
        "pair_style": log.get("pair_style"),
        "potential_files": ", ".join(pot_files) or None,
        "data_file": data_file,
        "velocity_T": log.get("velocity_T"),
        "timestep": main["timestep"] if main else log.get("timestep"),
        "time_unit": tu, "pressure_unit": pu, "length_unit": lu,
        "total_steps": total_steps or None,
        "simulated_time": round(total_time, 6) if total_time else None,
        "ensemble": main["ensemble"] if main else ("minimize" if protocol else None),
        "T_target": (main["T"][1] if main and main["T"] else None),
        "T_range": fmt_range(main["T"], "K") if main and main["T"] else None,
        "P_target": (main["P"][1] if main and main["P"] else None),
        "T_measured": round(main["T_mean"], 3) if main and main["T_mean"] is not None else None,
        "P_measured": round(main["P_mean"], 3) if main and main["P_mean"] is not None else None,
        "density_measured": round(main["density_mean"], 4) if main and main["density_mean"] else None,
        "n_stages": len(protocol) or None,
        "n_structures": len(structures) if multi else None,
        "structures": (", ".join(structures[:8]) + (f", … ({len(structures)} total)" if len(structures) > 8 else ""))
        if multi else None,
        "unresolved": ", ".join(log.get("unresolved") or []) or None,
        "lammps_version": log.get("version"),
    }
    return {"system": system, "protocol": protocol,
            "summary": lammps_summary(system, protocol), "conditions": lammps_conditions(system)}


def _stage_text(r: Dict, tu: str, pu: str) -> str:
    if r["kind"] == "minimize":
        return "minimize"
    parts = [r["ensemble"]]
    if r["T"]:
        parts.append(fmt_range(r["T"], "K"))
    if r["P"]:
        parts.append(fmt_range(r["P"], pu))
    t = fmt_time((r["time"] or 0) * r["count"], tu) if r["time"] else f"{r['steps']} steps"
    return " ".join(parts) + f" {t}" + (f" ({r['count']}×)" if r["count"] > 1 else "")


def lammps_summary(s: Dict, protocol: List[Dict]) -> str:
    kind = "MD" if s.get("total_steps") else ("Energy minimization" if protocol else "LAMMPS run")
    head = f"{s['ensemble']} MD" if (kind == "MD" and len({r['ensemble'] for r in protocol if r['kind'] == 'run'}) == 1) \
        else ("MD" if kind == "MD" else kind)
    if s.get("n_structures"):
        out = f"{head} of {s['n_structures']} structures; first"
        if s.get("data_file"):
            out += f" ({s['data_file'].split(' (')[0].rsplit('/', 1)[-1]})"
        out += f": {s['n_atoms']} atoms" if s.get("n_atoms") else ""
    else:
        out = f"{head} of {s['n_atoms']} atoms" if s.get("n_atoms") else head
    if s.get("composition"):
        out += f" ({s['composition']})"
    if s.get("box"):
        out += f" in a {s['box']} box"
    if s.get("density"):
        out += f", {s['density']:.3g} g/cm³"
    stages = [r for r in protocol if r["kind"] == "minimize" or (r["kind"] == "run" and r.get("steps"))]
    if len({(r["ensemble"], r["T"], r["P"]) for r in stages}) > 1:
        out += "; stages: " + " → ".join(_stage_text(r, s["time_unit"], s["pressure_unit"]) for r in stages[:6])
        if len(stages) > 6:
            out += f" → … ({len(stages)} stages)"
        if s.get("simulated_time"):
            out += f"; total {fmt_time(s['simulated_time'], s['time_unit'])}"
    else:
        if s.get("T_range"):
            out += f" at {s['T_range']}"
        if s.get("P_target") is not None:
            out += f", {s['P_target']:g} {s['pressure_unit']}"
        if s.get("simulated_time"):
            out += f"; {fmt_time(s['simulated_time'], s['time_unit'])}"
            if s.get("timestep"):
                out += f" ({s['total_steps']:,} steps × {s['timestep']:g} {s['time_unit']})"
    if s.get("pair_style"):
        out += f"; pair {s['pair_style'].split()[0]}"
        if s.get("potential_files"):
            out += f" ({s['potential_files']})"
    return out + "."


def lammps_conditions(s: Dict) -> str:
    parts = [s.get("ensemble") or ""]
    if s.get("T_range"):
        parts.append(s["T_range"])
    if s.get("P_target") is not None:
        parts.append(f"{s['P_target']:g} {s['pressure_unit']}")
    if s.get("simulated_time"):
        parts.append(fmt_time(s["simulated_time"], s["time_unit"]))
    return ", ".join(p for p in parts if p)


# ---------------------------------------------------------------------------
# VASP
# ---------------------------------------------------------------------------

def vasp_system(tags: Dict[str, str], kpoints: Dict, poscar: Dict, contcar: Dict,
                potcar: List[Dict], osz_md: Dict) -> Dict:
    species = poscar.get("species") or []
    counts = poscar.get("counts") or []
    pomass = [t.get("pomass") for t in potcar] if len(potcar) == len(species) else []
    masses = [pm or lammps.ELEMENT_MASSES.get(el) for el, pm in zip(species, pomass or [None] * len(species))]
    total_mass = sum(m * n for m, n in zip(masses, counts)) if masses and all(masses) else None
    vol = poscar.get("volume")
    density = total_mass * AMU_PER_A3_TO_G_PER_CM3 / vol if (total_mass and vol) else None
    md = vasp.md_settings(tags)
    ctype = vasp.calc_type(tags) if tags else None
    kmesh = kpoints.get("mesh") or (f"{kpoints.get('n_kpoints')} explicit" if kpoints.get("n_kpoints") else None)
    system = {
        "code": "vasp", "calc_type": ctype, "calc_type_text": VASP_TYPES.get(ctype, ctype),
        "n_atoms": poscar.get("n_atoms"), "formula": poscar.get("formula"),
        "composition": composition_text(dict(zip(species, counts))) or None,
        "box": fmt_box(poscar.get("lattice_lengths"), "Å", poscar.get("lattice_angles")) if poscar.get("lattice_lengths") else None,
        "box_lengths": poscar.get("lattice_lengths"), "box_angles": poscar.get("lattice_angles"),
        "volume": vol, "density": round(density, 4) if density else None,
        "volume_final": contcar.get("volume") if contcar else None,
        "selective_dynamics": poscar.get("selective_dynamics"),
        "functional": vasp.functional(tags, potcar),
        "encut": vasp._tag_num(tags, "ENCUT"), "kpoints": (" ".join(f"{kpoints.get('scheme', '')} {kmesh}".split())
                                                           if kmesh else kpoints.get("scheme")),
        "ispin": _int(vasp._tag_num(tags, "ISPIN")), "ediff": tags.get("EDIFF"), "ediffg": tags.get("EDIFFG"),
        "isif": _int(vasp._tag_num(tags, "ISIF")), "nsw": _int(vasp._tag_num(tags, "NSW")),
        "T_target": md.get("T_end"), "T_range": fmt_range((md["T_start"], md["T_end"]), "K") if md.get("T_start") else None,
        "timestep": md.get("potim_fs"), "time_unit": "fs" if md else None,
        "simulated_time": md.get("time_ps"),
        "T_measured": round(osz_md["T_mean"], 2) if osz_md.get("T_mean") else None,
        "md_steps_done": osz_md.get("md_steps"),
        "potcars": "; ".join(t["titel"] for t in potcar) or None,
    }
    return {"system": system, "protocol": [], "summary": vasp_summary(system),
            "conditions": vasp_conditions(system)}


def vasp_summary(s: Dict) -> str:
    out = f"VASP {s.get('calc_type_text') or 'calculation'} of {s.get('formula') or '?'} ({s.get('n_atoms') or '?'} atoms)"
    if s.get("box"):
        out += f", cell {s['box']}"
    if s.get("density"):
        out += f", {s['density']:.3g} g/cm³"
    tech = [x for x in (s.get("functional"), f"ENCUT {s['encut']:g} eV" if s.get("encut") else None,
                        f"k {s['kpoints']}" if s.get("kpoints") else None,
                        "spin-polarized" if s.get("ispin") == 2 else None,
                        "selective dynamics" if s.get("selective_dynamics") else None) if x]
    if tech:
        out += "; " + ", ".join(tech)
    if s.get("T_range"):
        out += f"; {s['T_range']}"
        if s.get("simulated_time"):
            out += f", {fmt_time(s['simulated_time'], 'ps')} ({int(s['nsw'])} × {s['timestep']:g} fs)"
    if s.get("volume_final") and s.get("volume") and s.get("isif") and s["isif"] >= 3:
        out += f"; volume {s['volume']:.4g} → {s['volume_final']:.4g} Å³"
    return out + "."


def vasp_conditions(s: Dict) -> str:
    parts = [s.get("calc_type_text") or "", s.get("functional") or ""]
    if s.get("T_range"):
        parts.append(s["T_range"])
    if s.get("simulated_time"):
        parts.append(fmt_time(s["simulated_time"], "ps"))
    return ", ".join(p for p in parts if p)


# ---------------------------------------------------------------------------
# Runs made of several sub-runs
# ---------------------------------------------------------------------------

def frames_summary(first: Dict, n_frames: int, results: Dict) -> str:
    s = first.get("system") or {}
    out = f"{n_frames} frames of {s.get('formula') or '?'} ({s.get('n_atoms') or '?'} atoms)"
    tech = [x for x in (s.get("calc_type_text"), s.get("functional"),
                        f"ENCUT {s['encut']:g} eV" if s.get("encut") else None) if x]
    if tech:
        out += ", " + ", ".join(tech)
    if results.get("energy_span") is not None:
        out += f"; energy span {results['energy_span']:.4g} eV, highest at frame {results.get('energy_max_frame')}"
    return out + "."


def multi_summary(subruns: List[Dict]) -> Optional[str]:
    units = [s for s in subruns if s.get("label") != "(top)" and s.get("summary")]
    if not units:
        return None
    frames = [s for s in units if (s.get("results") or {}).get("n_frames")]
    if frames and len(frames) == len(units):
        n = sorted({s["results"]["n_frames"] for s in frames})
        systems = [s.get("system") or {} for s in frames]
        atoms = sorted({x.get("n_atoms") for x in systems} - {None})
        formulas = sorted({x.get("formula") for x in systems} - {None})
        first = systems[0]
        out = f"{len(frames)} frame sets ({'/'.join(map(str, n))} frames each) of "
        out += (f"{formulas[0]} ({atoms[0]} atoms)" if len(formulas) == 1 else
                f"{len(formulas)} structures ({atoms[0]}–{atoms[-1]} atoms)" if atoms else "?")
        tech = [x for x in (first.get("calc_type_text"), first.get("functional"),
                            f"ENCUT {first['encut']:g} eV" if first.get("encut") else None) if x]
        if tech:
            out += ", " + ", ".join(tech)
        spans = [s["results"]["energy_span"] for s in frames if s["results"].get("energy_span") is not None]
        if spans:
            out += (f"; energy span {spans[0]:.4g} eV" if len(spans) == 1 else
                    f"; energy span {min(spans):.3g}–{max(spans):.3g} eV across sets")
        return out + "."
    temps = sorted({(s.get("system") or {}).get("T_target") for s in units} - {None})
    head = f"{len(units)} sub-run" + ("s" if len(units) != 1 else "")
    if len(temps) > 1:
        head += f" over T = {temps[0]:g}–{temps[-1]:g} K ({len(temps)} temperatures)"
    return f"{head}. {units[0]['label']}: {units[0]['summary']}"
