from pathlib import Path
from shutil import rmtree
from uuid import uuid4

import pytest


@pytest.fixture
def tmp_path():
    """Use an isolated, workspace-local test directory (no OS temp dependency)."""
    root = Path(__file__).resolve().parent / ".test-work"
    root.mkdir(exist_ok=True)
    path = root / uuid4().hex
    path.mkdir()
    try:
        yield path
    finally:
        if path.parent.resolve() == root.resolve() and path.exists():
            rmtree(path)
