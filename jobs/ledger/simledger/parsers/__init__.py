"""Parsers: pure functions from file text to dicts, tolerant of truncated input.

Files are only ever opened for reading, and large files only partially
(``read_head`` / ``read_tail``).
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Iterator


def read_head(path: Path, nbytes: int) -> str:
    try:
        with open(path, "rb") as f:
            return f.read(nbytes).decode("utf-8", errors="replace")
    except OSError:
        return ""


def read_tail(path: Path, nbytes: int) -> str:
    """Last nbytes of a file, starting at a line boundary when truncated."""
    try:
        with open(path, "rb") as f:
            f.seek(0, os.SEEK_END)
            size = f.tell()
            f.seek(max(0, size - nbytes))
            data = f.read()
    except OSError:
        return ""
    text = data.decode("utf-8", errors="replace")
    if size > nbytes and "\n" in text:
        text = text.split("\n", 1)[1]
    return text


def iter_lines(path: Path, max_bytes: int) -> Iterator[str]:
    """Stream lines, stopping after max_bytes."""
    seen = 0
    try:
        with open(path, "rb") as f:
            for raw in f:
                seen += len(raw)
                if seen > max_bytes:
                    return
                yield raw.decode("utf-8", errors="replace")
    except OSError:
        return


def to_num(val: str):
    """float(val) for plain numbers (Fortran 1.0d-5 included), else None."""
    s = str(val).strip().lower().replace("d", "e", 1) if val is not None else ""
    try:
        return float(s)
    except ValueError:
        return None
