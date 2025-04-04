import pytest
from unittest.mock import patch, MagicMock
from funannotate2.clean import is_duplicated, load_genome, clean
import tempfile
import os

# Test data for load_genome
FASTA_CONTENT = """>contig1
ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
>contig2
ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
>contig3
ATGCGTACGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC
"""


@pytest.fixture
def fasta_file():
    with tempfile.NamedTemporaryFile(delete=False, mode="w") as f:
        f.write(FASTA_CONTENT)
        f.flush()
        yield f.name
    os.remove(f.name)


def test_load_genome(fasta_file):
    genome, n50 = load_genome(fasta_file)
    assert len(genome) == 3
    assert n50 > 0
    assert all(isinstance(item, dict) for item in genome)
    assert all(
        "header" in item and "length" in item and "sequence" in item for item in genome
    )


@patch("funannotate2.clean.mp.Aligner")
@patch("funannotate2.clean.NTF")
def test_is_duplicated(mock_ntf, mock_aligner, fasta_file):
    # Mock the NamedTemporaryFile
    mock_ntf_instance = MagicMock()
    mock_ntf_instance.__enter__.return_value = mock_ntf_instance
    mock_ntf_instance.name = "mock_temp_file.fa"
    mock_ntf.return_value = mock_ntf_instance

    # Mock the Aligner
    mock_aligner_instance = MagicMock()

    # Create a mock hit with the properties we need
    mock_hit = MagicMock()
    mock_hit.mlen = 95  # Match length
    mock_hit.blen = 100  # Block length
    mock_hit.ctg_len = 100  # Contig length

    # Set up the map method to return our mock hit
    mock_aligner_instance.map.return_value = [mock_hit]
    mock_aligner.return_value = mock_aligner_instance

    # Mock os.remove to prevent it from trying to remove our mock file
    with patch("funannotate2.clean.os.remove"):
        data, _ = load_genome(fasta_file)
        with tempfile.TemporaryDirectory() as tmpdir:
            header, result, pident, coverage, _ = is_duplicated(data, 0, tmpdir, 90, 50)

            # Verify the results
            assert header == data[0]["header"]
            assert (
                result is True
            )  # Should be True because pident and coverage exceed thresholds
            assert pident == 95.0  # mlen/blen * 100 = 95/100 * 100 = 95.0
            assert coverage == 100.0  # blen/ctg_len * 100 = 100/100 * 100 = 100.0


@patch("funannotate2.clean.startLogging")
@patch("funannotate2.clean.runThreadJob")
@patch("funannotate2.clean.load_genome")
def test_clean(mock_load_genome, mock_run_thread_job, mock_start_logging, fasta_file):
    # Mock the logger
    mock_logger = MagicMock()
    mock_start_logging.return_value = mock_logger

    # Mock the genome data and N50 value
    mock_genome = [
        {"header": f"contig{i}", "length": 60, "sequence": "A" * 60, "status": None}
        for i in range(1, 4)  # contig1, contig2, contig3
    ]
    mock_n50 = 60
    mock_load_genome.return_value = (mock_genome, mock_n50)

    # Mock the thread job results
    mock_run_thread_job.return_value = [
        MagicMock(result=lambda: (f"contig{i}", False, 0, 0, 60)) for i in range(1, 4)
    ]

    class Args:
        fasta = fasta_file
        minlen = 10
        pident = 90
        cov = 50
        cpus = 1
        rename = None
        out = "output.fa"
        tmpdir = None
        logfile = "logfile.log"
        exhaustive = False

    args = Args()
    clean(args)

    # Check that the expected log messages were generated
    mock_logger.info.assert_any_call(
        "Loaded 3 contigs; 3 are larger than 10; N50 is 60 bp"
    )
    mock_logger.info.assert_any_call("Wrote 3 contigs to output.fa")

    # Check that the output file was created
    assert os.path.exists("output.fa")
    os.remove("output.fa")
