"""
Pytest configuration file for integration tests.
"""
import os
import sys
import tempfile
import pytest
from pathlib import Path


@pytest.fixture(scope="session")
def test_data_dir():
    """Return the path to the test data directory."""
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")


@pytest.fixture(scope="function")
def temp_dir():
    """Create a temporary directory for test outputs."""
    with tempfile.TemporaryDirectory() as tmpdir:
        yield tmpdir
