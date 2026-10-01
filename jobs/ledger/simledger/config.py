"""Defaults, overridable from ``<ledger home>/simledger.toml``.

Example simledger.toml:

    run_id_regex = '(?:^|[-_])(\\d{4})$'
    exclude = ["__pycache__", ".git"]
    max_walk_depth = 8
"""

from __future__ import annotations

import os
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

try:  # Python 3.11+
    import tomllib
except ImportError:  # pragma: no cover - older interpreters just skip the file
    tomllib = None


def default_home() -> Path:
    return Path(os.environ.get("SIMLEDGER_HOME", "~/ledger")).expanduser()


@dataclass
class Config:
    home: Path
    # A directory under runs/ is a run when its name matches this; group(1) is the id.
    run_id_regex: str = r"(?:^|[-_])(\d{4})$"
    exclude: List[str] = field(default_factory=lambda: ["__pycache__", ".git", ".ipynb_checkpoints"])
    max_walk_depth: int = 8          # below a run dir
    head_bytes: int = 64 * 1024      # how much of a file's start to read
    tail_bytes: int = 256 * 1024     # how much of a log's end to read
    max_text_bytes: int = 256 * 1024  # largest README/script stored as text
    running_window_s: int = 3600     # log touched this recently + no end marker -> running

    @property
    def db_path(self) -> Path:
        return self.home / "ledger.db"

    @property
    def inbox(self) -> Path:
        return self.home / "inbox"

    @property
    def cards(self) -> Path:
        return self.home / "cards"

    @property
    def id_re(self) -> "re.Pattern[str]":
        return re.compile(self.run_id_regex)

    @property
    def id_width(self) -> int:
        """Digits in a full run id (4 for the default regex)."""
        m = re.search(r"\\d\{(\d+)\}", self.run_id_regex)
        return int(m.group(1)) if m else 4

    def match_id(self, name: str) -> Optional[str]:
        m = self.id_re.search(name)
        return m.group(1) if m else None


def load(home: Optional[Path] = None) -> Config:
    home = Path(home).expanduser() if home else default_home()
    cfg = Config(home=home.resolve())
    toml_path = cfg.home / "simledger.toml"
    if tomllib and toml_path.is_file():
        with toml_path.open("rb") as f:
            overrides = tomllib.load(f)
        for key, val in overrides.items():
            if hasattr(cfg, key) and key != "home":
                setattr(cfg, key, val)
    return cfg
