"""
Configuration file for pytest.
"""

import os
import shutil
import tempfile

import pytest


@pytest.fixture
def test_data_dir():
    """Return the path to the test data directory."""
    return os.path.join(os.path.dirname(os.path.dirname(__file__)), "data")


@pytest.fixture
def temp_dir():
    """Create a temporary directory for test files."""
    temp_dir = tempfile.mkdtemp()
    yield temp_dir
    shutil.rmtree(temp_dir)


@pytest.fixture
def sample_gff_content():
    """Return sample GFF content for testing."""
    return """##gff-version 3
##sequence-region contig1 1 1000
contig1\tprediction\tgene\t1\t1000\t.\t+\t.\tID=gene1;Name=test_gene
contig1\tprediction\tmRNA\t1\t1000\t.\t+\t.\tID=mRNA1;Parent=gene1;Name=test_mrna
contig1\tprediction\texon\t1\t500\t.\t+\t.\tID=exon1;Parent=mRNA1
contig1\tprediction\texon\t600\t1000\t.\t+\t.\tID=exon2;Parent=mRNA1
contig1\tprediction\tCDS\t1\t500\t.\t+\t0\tID=cds1;Parent=mRNA1
contig1\tprediction\tCDS\t600\t1000\t.\t+\t0\tID=cds2;Parent=mRNA1
"""


@pytest.fixture
def sample_fasta_content():
    """Return sample FASTA content for testing."""
    return """>contig1
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC
"""


@pytest.fixture
def sample_gff_file(temp_dir, sample_gff_content):
    """Create a sample GFF file for testing."""
    file_path = os.path.join(temp_dir, "sample.gff3")
    with open(file_path, "w") as f:
        f.write(sample_gff_content)
    return file_path


@pytest.fixture
def sample_fasta_file(temp_dir, sample_fasta_content):
    """Create a sample FASTA file for testing."""
    file_path = os.path.join(temp_dir, "sample.fasta")
    with open(file_path, "w") as f:
        f.write(sample_fasta_content)
    return file_path
