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


@pytest.fixture(scope="function")
def sample_gff_file(temp_dir):
    """Create a sample GFF file for testing."""
    # Create a sample GFF file
    gff_content = """##gff-version 3
    contig1\tprediction\tgene\t1\t1000\t.\t+\t.\tID=gene1;Name=test_gene
    contig1\tprediction\tmRNA\t1\t1000\t.\t+\t.\tID=mRNA1;Parent=gene1;Name=test_mrna
    contig1\tprediction\texon\t1\t500\t.\t+\t.\tID=exon1;Parent=mRNA1
    contig1\tprediction\texon\t600\t1000\t.\t+\t.\tID=exon2;Parent=mRNA1
    contig1\tprediction\tCDS\t1\t500\t.\t+\t0\tID=cds1;Parent=mRNA1
    contig1\tprediction\tCDS\t600\t1000\t.\t+\t0\tID=cds2;Parent=mRNA1
    """

    # Create a sample FASTA file
    fasta_content = """>contig1
    ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
    ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
    """

    # Write the GFF file
    gff_file = os.path.join(temp_dir, "sample.gff3")
    with open(gff_file, "w") as f:
        f.write(gff_content)

    # Write the FASTA file
    fasta_file = os.path.join(temp_dir, "sample.fasta")
    with open(fasta_file, "w") as f:
        f.write(fasta_content)

    return gff_file
