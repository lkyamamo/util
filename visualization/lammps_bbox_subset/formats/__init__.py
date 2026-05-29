"""Format registry and auto-detection for trajectory adapters.

Adding a new format: implement the TrajectoryAdapter protocol in a new module
under formats/, then register it in ADAPTERS below.
"""

from __future__ import annotations

from formats import lammps_custom, lammps_xyz
from formats.base import FrameSpec

# Order matters for auto-detection: more-specific detectors first.
# lammps_custom must come before lammps_xyz because some .xyz files are
# lammps_custom dumps (they start with 'ITEM: TIMESTEP').
_REGISTRY: dict[str, object] = {
    "lammps_custom": lammps_custom,
    "lammps_xyz": lammps_xyz,
}


def get_format(name: str) -> object:
    """Return the adapter module for the given format name."""
    if name not in _REGISTRY:
        raise ValueError(
            f"Unknown format {name!r}. Available: {sorted(_REGISTRY)}"
        )
    return _REGISTRY[name]


def auto_detect(path: str) -> object:
    """Probe the file and return the matching adapter module.

    Detection order: lammps_custom → lammps_xyz.
    Raises ValueError if no adapter matches.
    """
    for name, adapter in _REGISTRY.items():
        if adapter.detect(path):  # type: ignore[attr-defined]
            return adapter
    raise ValueError(
        f"Could not auto-detect format of {path!r}. "
        f"Specify --format explicitly. Available: {sorted(_REGISTRY)}"
    )


def resolve_adapter(fmt: str, path: str) -> object:
    """Return adapter for fmt='auto' (detect) or a named format."""
    if fmt == "auto":
        return auto_detect(path)
    return get_format(fmt)
