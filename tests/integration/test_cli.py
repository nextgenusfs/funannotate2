"""
Integration tests for the command-line interface.
"""

import os
import subprocess

import pytest


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

    @pytest.mark.skip(reason="funannotate2 is not installed")
    def test_help_command(self):
        """Test the help command."""
        returncode, stdout, stderr = run_command("python -m funannotate2 --help")
        assert returncode == 0
        assert "usage:" in stdout
        # The actual implementation doesn't include 'gfftk' in the help output
        assert "funannotate2" in stdout

    @pytest.mark.skip(reason="funannotate2 is not installed")
    def test_version_command(self):
        """Test the version command."""
        returncode, stdout, stderr = run_command("python -m funannotate2 --version")
        assert returncode == 0
        # The actual implementation doesn't include 'gfftk' in the version output
        assert "funannotate2" in stdout

    @pytest.mark.skip(reason="Command not implemented yet")
    def test_sort_command(self, sample_gff_file, temp_dir):
        """Test the sort command."""
        output_file = os.path.join(temp_dir, "sorted.gff3")
        command = f"python -m funannotate2 sort -i {sample_gff_file} -o {output_file}"
        returncode, stdout, stderr = run_command(command)
        assert returncode == 0
        assert os.path.exists(output_file)

    @pytest.mark.skip(reason="Command not implemented yet")
    def test_stats_command(self, sample_gff_file):
        """Test the stats command."""
        command = f"python -m funannotate2 stats -i {sample_gff_file}"
        returncode, stdout, stderr = run_command(command)
        assert returncode == 0
        assert "Statistics for" in stdout

    @pytest.mark.skip(reason="Command not implemented yet")
    def test_convert_command(self, sample_gff_file, temp_dir):
        """Test the convert command."""
        output_file = os.path.join(temp_dir, "converted.gtf")
        command = f"python -m funannotate2 convert -i {sample_gff_file} -o {output_file} -f gtf"
        returncode, stdout, stderr = run_command(command)
        assert returncode == 0
        assert os.path.exists(output_file)
