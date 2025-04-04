"""
Comprehensive unit tests for the search module.
"""

import os
import tempfile
import pytest
from unittest.mock import patch, MagicMock, mock_open
import funannotate2.search as search


class TestSearchComprehensive:
    """Comprehensive tests for the search module."""

    @pytest.fixture
    def test_fasta(self, tmp_path):
        """Create a test FASTA file."""
        fasta_file = tmp_path / "test.fasta"
        with open(fasta_file, "w") as f:
            f.write(">seq1\n")
            f.write("ATGCATGCATGCATGCATGC\n")
            f.write(">seq2\n")
            f.write("GTACGTACGTACGTACGTAC\n")
        return str(fasta_file)

    @pytest.fixture
    def test_protein_fasta(self, tmp_path):
        """Create a test protein FASTA file."""
        fasta_file = tmp_path / "test_protein.fasta"
        with open(fasta_file, "w") as f:
            f.write(">seq1\n")
            f.write("MHACMHACMHACMHACMHAC\n")
            f.write(">seq2\n")
            f.write("VRTVRTVRTVRTVRTVRT\n")
        return str(fasta_file)

    @pytest.fixture
    def test_chunks_dir(self, tmp_path):
        """Create a test directory for FASTA chunks."""
        chunks_dir = tmp_path / "chunks"
        chunks_dir.mkdir()

        # Create chunk files
        chunk1 = chunks_dir / "chunk1.fasta"
        with open(chunk1, "w") as f:
            f.write(">seq1\n")
            f.write("MHACMHACMHACMHACMHAC\n")

        chunk2 = chunks_dir / "chunk2.fasta"
        with open(chunk2, "w") as f:
            f.write(">seq2\n")
            f.write("VRTVRTVRTVRTVRTVRT\n")

        return str(chunks_dir)

    def test_digitize_sequences(self, test_protein_fasta):
        """Test the digitize_sequences function."""
        # Call the function
        result = search.digitize_sequences(test_protein_fasta)

        # Check the result
        assert isinstance(result, dict)
        assert len(result) == 2
        assert "seq1" in result
        assert "seq2" in result
        assert isinstance(result["seq1"], int)
        assert isinstance(result["seq2"], int)

    @patch("funannotate2.search.os.path.isfile")
    @patch("funannotate2.search.os.path.getsize")
    @patch("funannotate2.search.subprocess.run")
    def test_pfam_search(
        self, mock_run, mock_getsize, mock_isfile, test_protein_fasta, tmp_path
    ):
        """Test the pfam_search function."""
        # Set up mocks
        mock_isfile.return_value = True
        mock_getsize.return_value = 1000  # Non-zero file size
        mock_run.return_value = MagicMock(returncode=0)

        # Create a mock pfam database
        pfam_db = str(tmp_path / "pfam.hmm")

        # Call the function
        result = search.pfam_search(test_protein_fasta, pfam_db, cpus=1, evalue=1e-5)

        # Check the result
        assert result is not None
        assert isinstance(result, str)
        assert mock_run.call_count == 1

    @patch("funannotate2.search.os.path.isfile")
    @patch("funannotate2.search.os.path.getsize")
    @patch("funannotate2.search.subprocess.run")
    def test_dbcan_search(
        self, mock_run, mock_getsize, mock_isfile, test_protein_fasta, tmp_path
    ):
        """Test the dbcan_search function."""
        # Set up mocks
        mock_isfile.return_value = True
        mock_getsize.return_value = 1000  # Non-zero file size
        mock_run.return_value = MagicMock(returncode=0)

        # Create a mock dbCAN database
        dbcan_db = str(tmp_path / "dbcan.hmm")

        # Call the function
        result = search.dbcan_search(test_protein_fasta, dbcan_db, cpus=1, evalue=1e-5)

        # Check the result
        assert result is not None
        assert isinstance(result, str)
        assert mock_run.call_count == 1

    @patch("funannotate2.search.os.path.isfile")
    @patch("funannotate2.search.os.path.getsize")
    @patch("funannotate2.search.subprocess.run")
    def test_merops_blast(
        self, mock_run, mock_getsize, mock_isfile, test_protein_fasta, tmp_path
    ):
        """Test the merops_blast function."""
        # Set up mocks
        mock_isfile.return_value = True
        mock_getsize.return_value = 1000  # Non-zero file size
        mock_run.return_value = MagicMock(returncode=0)

        # Create a mock MEROPS database
        merops_db = str(tmp_path / "merops.dmnd")

        # Call the function
        result = search.merops_blast(test_protein_fasta, merops_db, cpus=1, evalue=1e-5)

        # Check the result
        assert result is not None
        assert isinstance(result, str)
        assert mock_run.call_count == 1

    @patch("funannotate2.search.os.path.isfile")
    @patch("funannotate2.search.os.path.getsize")
    @patch("funannotate2.search.subprocess.run")
    def test_swissprot_blast(
        self, mock_run, mock_getsize, mock_isfile, test_protein_fasta, tmp_path
    ):
        """Test the swissprot_blast function."""
        # Set up mocks
        mock_isfile.return_value = True
        mock_getsize.return_value = 1000  # Non-zero file size
        mock_run.return_value = MagicMock(returncode=0)

        # Create a mock SwissProt database
        swissprot_db = str(tmp_path / "swissprot.dmnd")

        # Call the function
        result = search.swissprot_blast(
            test_protein_fasta, swissprot_db, cpus=1, evalue=1e-5
        )

        # Check the result
        assert result is not None
        assert isinstance(result, str)
        assert mock_run.call_count == 1

    @patch("funannotate2.search.os.path.isfile")
    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="seq1\tPF00001\t1\t100\t1e-10\n",
    )
    def test_pfam2tsv(self, mock_file, mock_isfile):
        """Test the pfam2tsv function."""
        # Set up mocks
        mock_isfile.return_value = True

        # Call the function
        result = search.pfam2tsv("pfam_results.txt")

        # Check the result
        assert isinstance(result, dict)
        assert "seq1" in result
        assert result["seq1"] == ["PF00001"]

    @patch("funannotate2.search.os.path.isfile")
    @patch(
        "builtins.open", new_callable=mock_open, read_data="seq1\tGT2\t1\t100\t1e-10\n"
    )
    def test_dbcan2tsv(self, mock_file, mock_isfile):
        """Test the dbcan2tsv function."""
        # Set up mocks
        mock_isfile.return_value = True

        # Call the function
        result = search.dbcan2tsv("dbcan_results.txt")

        # Check the result
        assert isinstance(result, dict)
        assert "seq1" in result
        assert result["seq1"] == ["GT2"]

    @patch("funannotate2.search.os.path.isfile")
    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="seq1\tM12.345\t1\t100\t1e-10\t50\n",
    )
    def test_merops2tsv(self, mock_file, mock_isfile):
        """Test the merops2tsv function."""
        # Set up mocks
        mock_isfile.return_value = True

        # Call the function
        result = search.merops2tsv("merops_results.txt")

        # Check the result
        assert isinstance(result, dict)
        assert "seq1" in result
        assert result["seq1"] == ["M12.345"]

    @patch("funannotate2.search.os.path.isfile")
    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="seq1\tP12345\t1\t100\t1e-10\t50\n",
    )
    def test_swissprot2tsv(self, mock_file, mock_isfile):
        """Test the swissprot2tsv function."""
        # Set up mocks
        mock_isfile.return_value = True

        # Call the function
        result = search.swissprot2tsv("swissprot_results.txt")

        # Check the result
        assert isinstance(result, dict)
        assert "seq1" in result
        assert result["seq1"] == ["P12345"]

    @patch("funannotate2.search.os.path.isdir")
    @patch("funannotate2.search.os.path.isfile")
    @patch("funannotate2.search.subprocess.run")
    def test_busco_search(
        self, mock_run, mock_isfile, mock_isdir, test_protein_fasta, tmp_path
    ):
        """Test the busco_search function."""
        # Set up mocks
        mock_isdir.return_value = True
        mock_isfile.return_value = True
        mock_run.return_value = MagicMock(returncode=0)

        # Create a mock BUSCO database
        busco_db = str(tmp_path / "busco_db")
        os.makedirs(busco_db, exist_ok=True)

        # Call the function
        result = search.busco_search(test_protein_fasta, busco_db, "proteins", cpus=1)

        # Check the result
        assert result is not None
        assert isinstance(result, str)
        assert mock_run.call_count == 1

    @patch("funannotate2.search.os.path.isfile")
    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="seq1\tBUSCOxyz\tcomplete\t1\t100\t+\t1e-10\n",
    )
    def test_busco2tsv(self, mock_file, mock_isfile):
        """Test the busco2tsv function."""
        # Set up mocks
        mock_isfile.return_value = True

        # Call the function
        result = search.busco2tsv("busco_results.txt")

        # Check the result
        assert isinstance(result, dict)
        assert "seq1" in result
        assert result["seq1"] == ["BUSCOxyz"]
