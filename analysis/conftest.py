import sys
from pathlib import Path

import pytest

# Tests run with --import-mode=importlib, which does not put rootdir on sys.path;
# testkit is imported from here.
sys.path.insert(0, str(Path(__file__).resolve().parent))

from testkit import Checker  # noqa: E402


@pytest.fixture
def checker():
    """Collect every check in the test and fail once at the end with all of them."""
    c = Checker()
    yield c
    c.assert_all()
