"""Bounding-box definition shared across all format adapters."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class BBox:
    xmin: float
    xmax: float
    ymin: float
    ymax: float
    zmin: float
    zmax: float

    def validate(self) -> None:
        for lo, hi, axis in (
            (self.xmin, self.xmax, "x"),
            (self.ymin, self.ymax, "y"),
            (self.zmin, self.zmax, "z"),
        ):
            if lo >= hi:
                raise ValueError(
                    f"bbox {axis}min ({lo}) must be strictly less than {axis}max ({hi})"
                )
