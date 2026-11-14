"""Soft assertions: record every check in a test, fail once at the end.

A test of a script's output checks dozens of values; reporting all the wrong
ones at once is more useful than stopping at the first.
"""

import math


class Checker:
    def __init__(self):
        self.failures = []

    def chk(self, label, got, expected, tol=1e-3):
        """Compare got to expected: within tol for floats (NaN equals NaN), exactly otherwise."""
        if isinstance(expected, float) or isinstance(got, float):
            if got is None:
                ok = False
            elif math.isnan(float(expected)):
                ok = math.isnan(float(got))
            else:
                ok = abs(float(got) - float(expected)) <= tol
        else:
            ok = got == expected
        if not ok:
            self.failures.append(f"{label}: got {_brief(got)}, expected {_brief(expected)}")
        return ok

    def assert_all(self):
        if self.failures:
            raise AssertionError(
                f"{len(self.failures)} check(s) failed:\n  " + "\n  ".join(self.failures)
            )


def _brief(value, limit=90):
    text = repr(value)
    return text if len(text) <= limit else f"{text[:limit]}... ({len(text)} chars)"
