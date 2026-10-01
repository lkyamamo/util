"""README files: whole text for search, plus `## <filename>` sections.

potentials/README.md and structures/README.md describe each file under a
heading named after it:

    ## 0039.data
    0232 end (finalized mod)
"""

from __future__ import annotations

import re
from typing import Dict, List

HEADING_RE = re.compile(r"^#{1,3}\s+(\S.*?)\s*$")


def file_sections(text: str) -> Dict[str, str]:
    sections: Dict[str, List[str]] = {}
    current = None
    for line in text.splitlines():
        m = HEADING_RE.match(line)
        if m:
            current = m.group(1).strip("` ")
            sections.setdefault(current, [])
            continue
        if current is not None:
            sections[current].append(line)
    return {k: "\n".join(v).strip() for k, v in sections.items()}


def mentioned_ids(text: str, width: int = 4) -> List[str]:
    return sorted(set(re.findall(rf"(?<!\d)(\d{{{width}}})(?!\d)", text)))


def title(text: str) -> str:
    for line in text.splitlines():
        m = HEADING_RE.match(line)
        if m:
            return m.group(1)
        if line.strip():
            return line.strip()[:120]
    return ""
