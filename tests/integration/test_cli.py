"""
Integration tests for the command-line interface.
"""
import os
import pytest
import subprocess
import sys


def run_command(command):
    """Run a command and return the output."""
    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        shell=True,
        universal_newlines=True,
    )
    stdout, stderr = process.communicate()
    return process.returncode, stdout, stderr


class TestCLI:
    """Tests for the command-line interface."""

    def test_help_command(self):
        """Test the help command."""
        returncode, stdout, stderr = run_command("python -m gfftk --help")
        assert returncode == 0
        assert "usage:" in stdout
        assert "gfftk" in stdout

    def test_version_command(self):
        """Test the version command."""
        returncode, stdout, stderr = run_command("python -m gfftk --version")
        assert returncode == 0
        assert "gfftk" in stdout

    def test_sort_command(self, sample_gff_file, temp_dir):
        """Test the sort command."""
        output_file = os.path.join(temp_dir, "sorted.gff3")
        command = f"python -m gfftk sort -i {sample_gff_file} -o {output_file}"
        returncode, stdout, stderr = run_command(command)
        assert returncode == 0
        assert os.path.exists(output_file)

    def test_stats_command(self, sample_gff_file):
        """Test the stats command."""
        command = f"python -m gfftk stats -i {sample_gff_file}"
        returncode, stdout, stderr = run_command(command)
        assert returncode == 0
        assert "Statistics for" in stdout

    def test_convert_command(self, sample_gff_file, temp_dir):
        """Test the convert command."""
        output_file = os.path.join(temp_dir, "converted.gtf")
        command = (
            f"python -m gfftk convert -i {sample_gff_file} -o {output_file} -f gtf"
        )
        returncode, stdout, stderr = run_command(command)
        assert returncode == 0
        assert os.path.exists(output_file)
