"""
Integration tests for the buscolite CLI.
"""

import os
import subprocess

import pytest


class TestBUSCOliteCLI:
    """Test the buscolite command-line interface."""

    @pytest.fixture(scope="class")
    def test_data_dir(self):
        """Return the path to the test data directory."""
        return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

    @pytest.mark.skip(reason="buscolite is not installed")
    def test_help_command(self):
        """Test that the help command works."""
        result = subprocess.run(
            ["buscolite", "--help"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "BUSCOlite: simplified BUSCO analysis for genome annotation" in result.stdout

    @pytest.mark.skip(reason="buscolite is not installed")
    def test_version_command(self):
        """Test that the version command works."""
        result = subprocess.run(
            ["buscolite", "--version"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 0
        assert "buscolite v" in result.stdout

    @pytest.mark.skip(reason="buscolite is not installed")
    def test_missing_arguments(self):
        """Test that the command fails with missing arguments."""
        result = subprocess.run(
            ["buscolite"],
            capture_output=True,
            text=True,
        )
        assert result.returncode == 1
        assert "the following arguments are required" in result.stderr

    @pytest.mark.skip(reason="buscolite is not installed")
    def test_invalid_mode(self):
        """Test that the command fails with an invalid mode."""
        result = subprocess.run(
            [
                "buscolite",
                "-i",
                "test.fasta",
                "-o",
                "test",
                "-m",
                "invalid",
                "-l",
                "test",
            ],
            capture_output=True,
            text=True,
        )
        assert result.returncode != 0
        assert "invalid choice" in result.stderr
