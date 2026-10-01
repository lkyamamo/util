"""Query syntax -> SQL.

    code=vasp ENCUT>=520 status=completed      field filters, ANDed
    pair_style~reaxff  group~ring             ~ = case-insensitive substring
    incar.ENCUT>=520   slurm.partition=priya   source-qualified parameter keys
    T_start=300                                any param/result key (any sub-run)

Bare words (no operator) are full-text terms (FTS5 over READMEs, inputs, notes).
"""

from __future__ import annotations

import re
from typing import List, Tuple

from .parsers import to_num
from .store import has_fts

RUN_FIELDS = {
    "run_key": "run_key", "key": "run_key", "project": "project", "id": "run_id", "run_id": "run_id",
    "label": "label", "dir": "label", "group": "group_path", "group_path": "group_path", "path": "path",
    "code": "code", "type": "calc_type", "calc_type": "calc_type", "status": "status",
    "atoms": "n_atoms", "n_atoms": "n_atoms", "formula": "formula", "cores": "cores", "nodes": "nodes",
    "job": "job_id", "job_id": "job_id", "job_name": "job_name", "wall": "wall_time_s",
    "wall_time_s": "wall_time_s", "start": "start_time", "end": "end_time", "missing": "missing",
    "subruns": "n_subruns",
}

COND_RE = re.compile(r"^([A-Za-z_][\w.:/-]*?)(>=|<=|!=|=|>|<|~)(.*)$")


def _cond(col: str, op: str, val: str, numeric_col: str = None) -> Tuple[str, List]:
    if op == "~":
        return f"{col} LIKE ?", [f"%{val}%"]
    num = to_num(val)
    if op in (">", "<", ">=", "<="):
        if num is None:
            return f"{col} {op} ?", [val]
        return f"{numeric_col or col} {op} ?", [num]
    sql_op = "=" if op == "=" else "!="
    if num is not None and numeric_col:
        if op == "=":
            return f"({numeric_col} = ? OR {col} = ? COLLATE NOCASE)", [num, val]
        return f"({numeric_col} IS NULL OR {numeric_col} != ?)", [num]
    return f"{col} {sql_op} ? COLLATE NOCASE", [val]


def build(query: str, include_missing: bool = False) -> Tuple[str, List, List[str]]:
    """Return (sql, args, text_terms) selecting rows from runs r."""
    where: List[str] = []
    args: List = []
    terms: List[str] = []
    for tok in _split(query):
        m = COND_RE.match(tok)
        if not m:
            terms.append(tok)
            continue
        key, op, val = m.group(1), m.group(2), m.group(3).strip("\"'")
        if key.lower() in RUN_FIELDS:
            col = f"r.{RUN_FIELDS[key.lower()]}"
            numeric = col if col.endswith(("n_atoms", "cores", "nodes", "wall_time_s", "missing", "n_subruns")) else None
            if numeric and op in ("=", "!="):
                numeric = None
            sql, a = _cond(col, op, val, numeric)
            where.append(sql)
            args += a
            if key.lower() == "missing":
                include_missing = True
            continue
        source = None
        if "." in key:
            source, key = key.split(".", 1)
        p_sql, p_args = _cond("p.value_text", op, val, "p.value_num")
        x_sql, x_args = _cond("x.value_text", op, val, "x.value_num")
        src_sql = " AND p.source = ? COLLATE NOCASE" if source else ""
        where.append(
            f"(EXISTS (SELECT 1 FROM params p WHERE p.run_key = r.run_key AND p.key = ? COLLATE NOCASE"
            f"{src_sql} AND {p_sql})"
            + ("" if source else
               f" OR EXISTS (SELECT 1 FROM results x WHERE x.run_key = r.run_key AND x.key = ? COLLATE NOCASE"
               f" AND {x_sql})")
            + ")")
        args += [key] + ([source] if source else []) + p_args + ([] if source else [key] + x_args)
    if not include_missing:
        where.append("r.missing = 0")
    sql = "SELECT r.* FROM runs r" + (" WHERE " + " AND ".join(where) if where else "") + \
          " ORDER BY r.project, r.run_id"
    return sql, args, terms


def _split(query: str) -> List[str]:
    # whitespace split, keeping "quoted values" together
    return [t for t in re.findall(r'(?:[^\s"]+(?:"[^"]*")?)+|"[^"]*"', query) if t]


def fts_query(terms: List[str]) -> str:
    return " ".join('"' + t.strip('"').replace('"', '""') + '"' for t in terms)


def run(conn, query: str, text: str = "", include_missing: bool = False):
    """Matching runs, plus matching analysis dirs when there are text terms."""
    sql, args, terms = build(query, include_missing)
    if text:
        terms += _split(text)
    analysis_hits = []
    if terms:
        if has_fts(conn):
            q = fts_query(terms)
            keys = {row[0]: row[1] for row in conn.execute(
                "SELECT run_key, kind FROM fts WHERE fts MATCH ? ORDER BY rank", (q,))}
        else:
            keys = {}
            like = " AND ".join("body LIKE ?" for _ in terms)
            for row in conn.execute(f"SELECT run_key, kind FROM fts WHERE {like}", [f"%{t}%" for t in terms]):
                keys[row[0]] = row[1]
        run_keys = [k for k, kind in keys.items() if kind in ("run", "note")]
        analysis_hits = [k for k, kind in keys.items() if kind == "analysis"]
        sql = sql.replace(" ORDER BY", f" {'AND' if ' WHERE ' in sql else 'WHERE'} r.run_key IN "
                                        f"({','.join('?' * len(run_keys)) or 'NULL'}) ORDER BY")
        args += run_keys
    return conn.execute(sql, args).fetchall(), analysis_hits
