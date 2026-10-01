"""SQLite schema and writes. The database is the ledger's source of truth.

Everything a scan derives from files (runs, subruns, params, results, files,
analysis) is replaced wholesale when a run is reparsed. Hook events and
user notes are never touched by a rescan.
"""

from __future__ import annotations

import json
import sqlite3
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, Optional

SCHEMA_VERSION = 2

SCHEMA = """
CREATE TABLE IF NOT EXISTS meta (key TEXT PRIMARY KEY, value TEXT);

CREATE TABLE IF NOT EXISTS projects (
    name TEXT PRIMARY KEY, path TEXT, present TEXT, first_seen TEXT, last_scanned TEXT);

CREATE TABLE IF NOT EXISTS runs (
    run_key TEXT PRIMARY KEY, project TEXT, run_id TEXT, label TEXT, group_path TEXT,
    path TEXT, code TEXT, calc_type TEXT, status TEXT, status_evidence TEXT,
    n_atoms INTEGER, formula TEXT, start_time TEXT, end_time TEXT, wall_time_s REAL,
    cores INTEGER, nodes INTEGER, job_name TEXT, job_id TEXT, n_subruns INTEGER,
    n_files INTEGER, total_bytes INTEGER, readme TEXT, fingerprint TEXT,
    parser_version INTEGER, first_seen TEXT, last_scanned TEXT, last_changed TEXT,
    missing INTEGER DEFAULT 0, summary TEXT, system TEXT);
CREATE INDEX IF NOT EXISTS runs_id ON runs(run_id);

CREATE TABLE IF NOT EXISTS subruns (
    run_key TEXT, label TEXT, relpath TEXT, code TEXT, calc_type TEXT, status TEXT,
    status_evidence TEXT, n_atoms INTEGER, start_time TEXT, end_time TEXT, wall_time_s REAL,
    summary TEXT, conditions TEXT, system TEXT, protocol TEXT,
    PRIMARY KEY (run_key, label));

CREATE TABLE IF NOT EXISTS params (
    run_key TEXT, subrun TEXT, source TEXT, key TEXT, value_text TEXT, value_num REAL);
CREATE INDEX IF NOT EXISTS params_key ON params(key, run_key);

CREATE TABLE IF NOT EXISTS results (
    run_key TEXT, subrun TEXT, key TEXT, value_num REAL, value_text TEXT, unit TEXT);
CREATE INDEX IF NOT EXISTS results_key ON results(key, run_key);

CREATE TABLE IF NOT EXISTS files (
    run_key TEXT, relpath TEXT, category TEXT, size INTEGER, mtime REAL,
    symlink TEXT, broken INTEGER);
CREATE INDEX IF NOT EXISTS files_run ON files(run_key);

CREATE TABLE IF NOT EXISTS analysis (
    analysis_key TEXT PRIMARY KEY, project TEXT, name TEXT, path TEXT, hint TEXT,
    run_ids TEXT, readme TEXT, n_files INTEGER, total_bytes INTEGER, fingerprint TEXT,
    first_seen TEXT, last_scanned TEXT, missing INTEGER DEFAULT 0);

CREATE TABLE IF NOT EXISTS analysis_runs (
    analysis_key TEXT, run_key TEXT, subrun TEXT, source TEXT);
CREATE INDEX IF NOT EXISTS analysis_runs_run ON analysis_runs(run_key);

CREATE TABLE IF NOT EXISTS analysis_files (
    analysis_key TEXT, relpath TEXT, category TEXT, size INTEGER, mtime REAL);

CREATE TABLE IF NOT EXISTS project_files (
    project TEXT, kind TEXT, name TEXT, size INTEGER, mtime REAL, description TEXT,
    mentioned_ids TEXT, PRIMARY KEY (project, kind, name));

CREATE TABLE IF NOT EXISTS events (
    event_id TEXT PRIMARY KEY, kind TEXT, created TEXT, ingested TEXT, data TEXT,
    project TEXT, run_id TEXT, subrun TEXT, run_key TEXT);
CREATE INDEX IF NOT EXISTS events_run ON events(run_key);

CREATE TABLE IF NOT EXISTS user_notes (
    note_id INTEGER PRIMARY KEY AUTOINCREMENT, run_key TEXT, note TEXT, tags TEXT, created TEXT);

CREATE TABLE IF NOT EXISTS scans (
    scan_id INTEGER PRIMARY KEY AUTOINCREMENT, root TEXT, started TEXT, finished TEXT, counts TEXT);

CREATE TABLE IF NOT EXISTS scan_warnings (scan_id INTEGER, path TEXT, message TEXT);

-- Problems found while parsing a run; kept until the run is reparsed.
CREATE TABLE IF NOT EXISTS run_warnings (run_key TEXT, path TEXT, message TEXT);
"""

FTS_SCHEMA = "CREATE VIRTUAL TABLE IF NOT EXISTS fts USING fts5(run_key UNINDEXED, kind UNINDEXED, body)"


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def connect(db_path: Path, readonly: bool = False) -> sqlite3.Connection:
    if readonly:
        conn = sqlite3.connect(f"file:{db_path}?mode=ro", uri=True)
    else:
        db_path.parent.mkdir(parents=True, exist_ok=True)
        conn = sqlite3.connect(str(db_path), timeout=60)
        conn.executescript(SCHEMA)
        try:
            conn.execute(FTS_SCHEMA)
        except sqlite3.OperationalError:  # no FTS5: text search falls back to LIKE
            pass
        migrate(conn)
        conn.execute("INSERT OR REPLACE INTO meta VALUES ('schema_version', ?)", (str(SCHEMA_VERSION),))
    conn.row_factory = sqlite3.Row
    return conn


# Columns added after a table was first created: (table, column, type).
ADDED_COLUMNS = [
    ("runs", "summary", "TEXT"), ("runs", "system", "TEXT"),
    ("subruns", "summary", "TEXT"), ("subruns", "conditions", "TEXT"),
    ("subruns", "system", "TEXT"), ("subruns", "protocol", "TEXT"),
]


def migrate(conn: sqlite3.Connection) -> None:
    for table, col, typ in ADDED_COLUMNS:
        have = {r[1] for r in conn.execute(f"PRAGMA table_info({table})")}
        if col not in have:
            conn.execute(f"ALTER TABLE {table} ADD COLUMN {col} {typ}")


def has_fts(conn: sqlite3.Connection) -> bool:
    return conn.execute("SELECT 1 FROM sqlite_master WHERE name='fts'").fetchone() is not None


def _num(v):
    if isinstance(v, bool):
        return None
    if isinstance(v, (int, float)):
        return float(v)
    from .parsers import to_num
    return to_num(v) if isinstance(v, str) else None


def _text(v) -> Optional[str]:
    if v is None:
        return None
    if isinstance(v, (dict, list)):
        return json.dumps(v, sort_keys=True)
    return str(v)


def clear_run(conn: sqlite3.Connection, run_key: str) -> None:
    for table in ("subruns", "params", "results", "files", "run_warnings"):
        conn.execute(f"DELETE FROM {table} WHERE run_key = ?", (run_key,))
    if has_fts(conn):
        conn.execute("DELETE FROM fts WHERE run_key = ? AND kind != 'note'", (run_key,))


def upsert_row(conn: sqlite3.Connection, table: str, key_cols: Iterable[str], row: Dict) -> None:
    cols = list(row)
    placeholders = ",".join("?" for _ in cols)
    updates = ",".join(f"{c}=excluded.{c}" for c in cols if c not in key_cols)
    conn.execute(
        f"INSERT INTO {table} ({','.join(cols)}) VALUES ({placeholders}) "
        f"ON CONFLICT ({','.join(key_cols)}) DO UPDATE SET {updates}",
        [row[c] for c in cols],
    )


def add_params(conn, run_key: str, subrun: str, source: str, items: Dict) -> None:
    conn.executemany(
        "INSERT INTO params VALUES (?,?,?,?,?,?)",
        [(run_key, subrun, source, k, _text(v), _num(v)) for k, v in items.items() if v is not None],
    )


def add_results(conn, run_key: str, subrun: str, items: Dict, units: Optional[Dict] = None) -> None:
    units = units or {}
    conn.executemany(
        "INSERT INTO results VALUES (?,?,?,?,?,?)",
        [(run_key, subrun, k, _num(v), _text(v), units.get(k)) for k, v in items.items() if v is not None],
    )


def add_fts(conn, run_key: str, kind: str, body: str) -> None:
    if body and has_fts(conn):
        conn.execute("INSERT INTO fts VALUES (?,?,?)", (run_key, kind, body))
