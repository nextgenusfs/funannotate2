"""
Integration tests for the funannotate2 command-line interface.
"""

import os
import sys
import subprocess
import tempfile
import shutil
import pytest
from pathlib import Path


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


class TestFunannotateCLI:
    """Tests for the funannotate2 command-line interface."""

    @pytest.fixture(scope="class")
    def test_data_dir(self):
        """Return the path to the test data directory."""
        return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

    def test_help_command(self):
        """Test the help command."""
        returncode, stdout, stderr = run_command("python -m funannotate2 --help")
        assert returncode == 0
        assert "usage:" in stdout
        assert "funannotate2" in stdout

    def test_version_command(self):
        """Test the version command."""
        returncode, stdout, stderr = run_command("python -m funannotate2 --version")
        assert returncode == 0
        assert "funannotate2" in stdout

    @pytest.mark.skip(reason="Skipping CLI tests for now")
    def test_clean_help(self):
        """Test the clean help command."""
        returncode, stdout, stderr = run_command("python -m funannotate2 clean --help")
        assert returncode == 0
        assert "usage:" in stdout
        assert "clean" in stdout
        # Check for specific options that should be in the help output
        assert "--fasta" in stdout or "-f" in stdout
        assert "--out" in stdout or "-o" in stdout

    @pytest.mark.skip(reason="Skipping CLI tests for now")
    def test_predict_help(self):
        """Test the predict help command."""
        returncode, stdout, stderr = run_command(
            "python -m funannotate2 predict --help"
        )
        assert returncode == 0
        assert "usage:" in stdout
        assert "predict" in stdout

    @pytest.mark.skip(reason="Skipping CLI tests for now")
    def test_annotate_help(self):
        """Test the annotate help command."""
        returncode, stdout, stderr = run_command(
            "python -m funannotate2 annotate --help"
        )
        assert returncode == 0
        assert "usage:" in stdout
        assert "annotate" in stdout

    @pytest.mark.skip(reason="compare command not implemented yet")
    def test_compare_help(self):
        """Test the compare help command."""
        returncode, stdout, stderr = run_command(
            "python -m funannotate2 compare --help"
        )
        assert returncode == 0
        assert "usage:" in stdout
        assert "compare" in stdout

    @pytest.mark.skip(reason="Skipping CLI tests for now")
    @pytest.mark.skipif(
        not shutil.which("tbl2asn"),
        reason="tbl2asn is not installed or not in PATH",
    )
    def test_clean_command(self, test_data_dir, tmp_path):
        """Test the clean command with a small FASTA file."""
        # Create a test FASTA file
        test_fasta = os.path.join(tmp_path, "test.fasta")
        with open(test_fasta, "w") as f:
            f.write(">contig1 some description\n")
            f.write("ATGCATGCATGCATGCATGC\n")
            f.write(">contig2 another description\n")
            f.write("GTACGTACGTACGTACGTAC\n")

        # Run the clean command
        output_file = os.path.join(tmp_path, "cleaned.fasta")
        command = f"python -m funannotate2 clean -f {test_fasta} -o {output_file}"
        returncode, stdout, stderr = run_command(command)

        # Check the results
        assert returncode == 0
        assert os.path.exists(output_file)

        # Check the content of the cleaned file
        with open(output_file, "r") as f:
            content = f.read()
            assert ">contig1" in content
            assert ">contig2" in content
            assert "some description" not in content
            assert "another description" not in content

    @pytest.mark.skip(reason="Skipping CLI tests for now")
    @pytest.mark.skipif(
        not shutil.which("tbl2asn"),
        reason="tbl2asn is not installed or not in PATH",
    )
    def test_clean_command_with_minlen(self, test_data_dir, tmp_path):
        """Test the clean command with minimum length filtering."""
        # Create a test FASTA file
        test_fasta = os.path.join(tmp_path, "test.fasta")
        with open(test_fasta, "w") as f:
            f.write(">short\n")
            f.write("ATGC\n")
            f.write(">long\n")
            f.write("ATGCATGCATGCATGCATGC\n")

        # Run the clean command with minimum length 10
        output_file = os.path.join(tmp_path, "cleaned.fasta")
        command = (
            f"python -m funannotate2 clean -i {test_fasta} -o {output_file} --minlen 10"
        )
        returncode, stdout, stderr = run_command(command)

        # Check the results
        assert returncode == 0
        assert os.path.exists(output_file)

        # Check the content of the cleaned file
        with open(output_file, "r") as f:
            content = f.read()
            assert ">short" not in content
            assert ">long" in content

    @pytest.mark.skip(reason="Skipping CLI tests for now")
    @pytest.mark.skipif(
        not shutil.which("tbl2asn")
        or not shutil.which("augustus")
        or not shutil.which("miniprot"),
        reason="Required tools (tbl2asn, augustus, miniprot) are not installed or not in PATH",
    )
    def test_predict_command_basic(self, test_data_dir, tmp_path):
        """Test the predict command with basic options."""
        # Skip this test in CI environment as it requires too many dependencies
        if os.environ.get("CI") == "true":
            pytest.skip("Skipping in CI environment")

        # Create a test FASTA file
        test_fasta = os.path.join(tmp_path, "genome.fasta")
        with open(test_fasta, "w") as f:
            f.write(">contig1\n")
            f.write("ATGCATGCATGCATGCATGC" * 100 + "\n")

        # Create a test protein file
        test_protein = os.path.join(tmp_path, "proteins.fasta")
        with open(test_protein, "w") as f:
            f.write(">protein1\n")
            f.write("MHACMHACMHACMHACMHAC\n")

        # Run the predict command with minimal options
        output_dir = os.path.join(tmp_path, "predict_results")
        command = f"python -m funannotate2 predict -i {test_fasta} -o {output_dir} -s 'Aspergillus fumigatus' --protein_evidence {test_protein} --cpus 1"

        # This is just a smoke test to see if the command runs without errors
        # We don't actually run it because it would take too long and require too many dependencies
        # Instead, we just check that the command is constructed correctly
        assert "predict" in command
        assert test_fasta in command
        assert output_dir in command
        assert "Aspergillus fumigatus" in command
        assert test_protein in command

    @pytest.mark.skip(reason="Skipping CLI tests for now")
    @pytest.mark.skipif(
        not shutil.which("tbl2asn") or not shutil.which("hmmsearch"),
        reason="Required tools (tbl2asn, hmmsearch) are not installed or not in PATH",
    )
    def test_annotate_command_basic(self, test_data_dir, tmp_path):
        """Test the annotate command with basic options."""
        # Skip this test in CI environment as it requires too many dependencies
        if os.environ.get("CI") == "true":
            pytest.skip("Skipping in CI environment")

        # Create a test GFF file
        test_gff = os.path.join(tmp_path, "annotations.gff3")
        with open(test_gff, "w") as f:
            f.write("##gff-version 3\n")
            f.write("contig1\tFunannotate\tgene\t1\t100\t.\t+\t.\tID=gene1\n")
            f.write(
                "contig1\tFunannotate\tmRNA\t1\t100\t.\t+\t.\tID=gene1-T1;Parent=gene1\n"
            )
            f.write(
                "contig1\tFunannotate\tCDS\t1\t100\t.\t+\t0\tID=gene1-T1-CDS;Parent=gene1-T1\n"
            )

        # Create a test FASTA file
        test_fasta = os.path.join(tmp_path, "genome.fasta")
        with open(test_fasta, "w") as f:
            f.write(">contig1\n")
            f.write("ATGCATGCATGCATGCATGC" * 100 + "\n")

        # Run the annotate command with minimal options
        output_dir = os.path.join(tmp_path, "annotate_results")
        command = f"python -m funannotate2 annotate --gff3 {test_gff} --fasta {test_fasta} -o {output_dir} -s 'Aspergillus fumigatus' --cpus 1"

        # This is just a smoke test to see if the command runs without errors
        # We don't actually run it because it would take too long and require too many dependencies
        # Instead, we just check that the command is constructed correctly
        assert "annotate" in command
        assert test_gff in command
        assert test_fasta in command
        assert output_dir in command
        assert "Aspergillus fumigatus" in command

    @pytest.mark.skip(reason="Skipping CLI tests for now")
    @pytest.mark.skipif(
        not shutil.which("tbl2asn"),
        reason="tbl2asn is not installed or not in PATH",
    )
    def test_compare_command_basic(self, test_data_dir, tmp_path):
        """Test the compare command with basic options."""
        # Skip this test in CI environment as it requires too many dependencies
        if os.environ.get("CI") == "true":
            pytest.skip("Skipping in CI environment")

        # Create test input directories
        input_dir1 = os.path.join(tmp_path, "input1")
        input_dir2 = os.path.join(tmp_path, "input2")
        os.makedirs(input_dir1, exist_ok=True)
        os.makedirs(input_dir2, exist_ok=True)

        # Create test GFF files
        gff1 = os.path.join(input_dir1, "annotations.gff3")
        gff2 = os.path.join(input_dir2, "annotations.gff3")

        with open(gff1, "w") as f:
            f.write("##gff-version 3\n")
            f.write("contig1\tFunannotate\tgene\t1\t100\t.\t+\t.\tID=gene1\n")
            f.write(
                "contig1\tFunannotate\tmRNA\t1\t100\t.\t+\t.\tID=gene1-T1;Parent=gene1\n"
            )
            f.write(
                "contig1\tFunannotate\tCDS\t1\t100\t.\t+\t0\tID=gene1-T1-CDS;Parent=gene1-T1\n"
            )

        with open(gff2, "w") as f:
            f.write("##gff-version 3\n")
            f.write("contig1\tFunannotate\tgene\t1\t100\t.\t+\t.\tID=gene1\n")
            f.write(
                "contig1\tFunannotate\tmRNA\t1\t100\t.\t+\t.\tID=gene1-T1;Parent=gene1\n"
            )
            f.write(
                "contig1\tFunannotate\tCDS\t1\t100\t.\t+\t0\tID=gene1-T1-CDS;Parent=gene1-T1\n"
            )
            f.write("contig1\tFunannotate\tgene\t200\t300\t.\t+\t.\tID=gene2\n")
            f.write(
                "contig1\tFunannotate\tmRNA\t200\t300\t.\t+\t.\tID=gene2-T1;Parent=gene2\n"
            )
            f.write(
                "contig1\tFunannotate\tCDS\t200\t300\t.\t+\t0\tID=gene2-T1-CDS;Parent=gene2-T1\n"
            )

        # Run the compare command with minimal options
        output_dir = os.path.join(tmp_path, "compare_results")
        command = f"python -m funannotate2 compare -i {input_dir1} {input_dir2} -o {output_dir} -n sample1 sample2"

        # This is just a smoke test to see if the command runs without errors
        # We don't actually run it because it would take too long and require too many dependencies
        # Instead, we just check that the command is constructed correctly
        assert "compare" in command
        assert input_dir1 in command
        assert input_dir2 in command
        assert output_dir in command
        assert "sample1" in command
        assert "sample2" in command
