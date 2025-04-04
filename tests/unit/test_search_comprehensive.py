"""
Comprehensive unit tests for the search module.
"""

import os
import sys
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
        assert result is not None
        assert isinstance(result, list)
        assert len(result) > 0
        # Check that each item is a DigitalSequence object
        for seq in result:
            assert hasattr(seq, "name")

    @patch("funannotate2.search.env")
    @patch("funannotate2.search.checkfile")
    @patch("funannotate2.search.hmmer_search")
    def test_pfam_search(self, mock_hmmer_search, mock_checkfile, mock_env):
        """Test the pfam_search function."""
        # Set up mocks
        mock_env.get.return_value = "/path/to/db"
        mock_checkfile.return_value = True
        mock_hmmer_search.return_value = [
            {
                "id": "seq1",
                "name": "PF00001",
                "accession": "PF00001.1",
                "description": "Test domain",
                "bitscore": 100.0,
                "evalue": 1e-10,
                "hmm_coverage": 0.8,
            }
        ]

        # Create mock sequences
        sequences = (
            MagicMock()
        )  # This would be a pyhmmer.easel.SequenceFile object in real usage

        # Call the function
        result = search.pfam_search(sequences, cpus=1)

        # Check the result
        assert result is not None
        assert isinstance(result, list)
        assert len(result) > 0
        assert result[0]["id"] == "seq1"
        assert result[0]["accession"] == "PF00001.1"
        assert mock_hmmer_search.call_count == 1
        mock_hmmer_search.assert_called_once_with(
            os.path.join("/path/to/db", "Pfam-A.hmm.h3m"),
            sequences,
            cpus=1,
            bit_cutoffs="gathering",
        )

    @patch("funannotate2.search.env")
    @patch("funannotate2.search.checkfile")
    @patch("funannotate2.search.hmmer_scan")
    def test_dbcan_search(self, mock_hmmer_scan, mock_checkfile, mock_env):
        """Test the dbcan_search function."""
        # Set up mocks
        mock_env.get.return_value = "/path/to/db"
        mock_checkfile.return_value = True
        mock_hmmer_scan.return_value = [
            {
                "id": "seq1",
                "name": "GT2",
                "accession": "GT2.1",
                "description": "Test domain",
                "bitscore": 100.0,
                "evalue": 1e-10,
                "hmm_coverage": 0.8,
            }
        ]

        # Create mock sequences
        sequences = (
            MagicMock()
        )  # This would be a pyhmmer.easel.SequenceFile object in real usage

        # Call the function
        result = search.dbcan_search(sequences, cpus=1, evalue=1e-15)

        # Check the result
        assert result is not None
        assert isinstance(result, list)
        assert len(result) > 0
        assert result[0]["id"] == "seq1"
        assert result[0]["name"] == "GT2"
        assert mock_hmmer_scan.call_count == 1
        mock_hmmer_scan.assert_called_once_with(
            os.path.join("/path/to/db", "dbCAN.hmm.h3m"),
            sequences,
            cpus=1,
            evalue=1e-15,
        )

    @patch("funannotate2.search.env")
    @patch("funannotate2.search.checkfile")
    @patch("funannotate2.search.diamond_blast")
    def test_merops_blast(self, mock_diamond_blast, mock_checkfile, mock_env):
        """Test the merops_blast function."""
        # Set up mocks
        mock_env.get.return_value = "/path/to/db"
        mock_checkfile.return_value = True
        mock_diamond_blast.return_value = [
            {
                "qseqid": "seq1",
                "qtitle": "seq1 test protein",
                "qlen": 100,
                "sseqid": "MER0001",
                "stitle": "MER0001 M12.345",
                "slen": 200,
                "pident": 80.0,
                "length": 90,
                "qstart": 1,
                "qend": 90,
                "sstart": 1,
                "send": 90,
                "evalue": 1e-10,
                "bitscore": 100.0,
            }
        ]

        # Call the function
        result = search.merops_blast("test_protein.fasta", evalue=1e-5, cpus=1)

        # Check the result
        assert result is not None
        assert isinstance(result, list)
        assert len(result) > 0
        assert result[0]["qseqid"] == "seq1"
        assert result[0]["family"] == "M12.345"
        assert mock_diamond_blast.call_count == 1
        mock_diamond_blast.assert_called_once_with(
            os.path.join("/path/to/db", "merops.dmnd"),
            "test_protein.fasta",
            cpus=1,
            evalue=1e-5,
            max_target_seqs=1,
        )

    @patch("funannotate2.search.env")
    @patch("funannotate2.search.checkfile")
    @patch("funannotate2.search.diamond_blast")
    @patch("funannotate2.search.parse_swissprot_headers")
    def test_swissprot_blast(
        self, mock_parse_headers, mock_diamond_blast, mock_checkfile, mock_env
    ):
        """Test the swissprot_blast function."""
        # Set up mocks
        mock_env.get.return_value = "/path/to/db"
        mock_checkfile.return_value = True
        mock_diamond_blast.return_value = [
            {
                "qseqid": "seq1",
                "qtitle": "seq1 test protein",
                "qlen": 100,
                "sseqid": "sp|P12345|TEST_HUMAN",
                "stitle": "sp|P12345|TEST_HUMAN Test protein OS=Human OX=9606 GN=TEST",
                "slen": 200,
                "pident": 80.0,
                "length": 90,
                "qstart": 1,
                "qend": 90,
                "sstart": 1,
                "send": 90,
                "evalue": 1e-10,
                "bitscore": 100.0,
            }
        ]
        mock_parse_headers.return_value = {
            "description": "Test protein",
            "accession": "P12345",
            "sp_name": "TEST_HUMAN",
            "GN": "TEST",
        }

        # Call the function
        result = search.swissprot_blast("test_protein.fasta", evalue=1e-5, cpus=1)

        # Check the result
        assert result is not None
        assert isinstance(result, list)
        assert len(result) > 0
        assert result[0]["query"] == "seq1"
        assert result[0]["accession"] == "P12345"
        assert result[0]["sp_name"] == "TEST_HUMAN"
        assert mock_diamond_blast.call_count == 1
        mock_diamond_blast.assert_called_once_with(
            os.path.join("/path/to/db", "uniprot.dmnd"),
            "test_protein.fasta",
            cpus=1,
            evalue=1e-5,
            max_target_seqs=1,
        )

    @patch("funannotate2.search.json.dump")
    def test_pfam2tsv(self, mock_json_dump):
        """Test the pfam2tsv function."""
        # Create mock results
        results = [
            {
                "id": "seq1",
                "name": "PF00001",
                "accession": "PF00001.1",
                "description": "Test domain",
            }
        ]

        # Call the function
        result = search.pfam2tsv(results, "output.json", "annotations.txt")

        # Check the result
        assert isinstance(result, dict)
        assert "seq1" in result
        assert result["seq1"]["db_xref"] == ["PFAM:PF00001.1"]
        assert mock_json_dump.call_count == 1

    @patch("funannotate2.search.json.dump")
    def test_dbcan2tsv(self, mock_json_dump):
        """Test the dbcan2tsv function."""
        # Create mock results
        results = [
            {
                "id": "seq1",
                "name": "GT2",
                "accession": "GT2.1",
                "description": "Test domain",
            }
        ]

        # Call the function
        result = search.dbcan2tsv(results, "output.json", "annotations.txt")

        # Check the result
        assert isinstance(result, dict)
        assert "seq1" in result
        assert "note" in result["seq1"]
        assert result["seq1"]["note"] == ["CAZy:GT2"]
        assert mock_json_dump.call_count == 1

    @patch("funannotate2.search.json.dump")
    def test_merops2tsv(self, mock_json_dump):
        """Test the merops2tsv function."""
        # Create mock results
        results = [
            {
                "qseqid": "seq1",
                "sseqid": "MER0001",
                "family": "M12.345",
                "pident": 80.0,
                "evalue": 1e-10,
            }
        ]

        # Call the function
        result = search.merops2tsv(results, "output.json", "annotations.txt")

        # Check the result
        assert isinstance(result, dict)
        assert "seq1" in result
        assert "note" in result["seq1"]
        assert result["seq1"]["note"] == ["MEROPS:MER0001 M12.345"]
        assert mock_json_dump.call_count == 1

    @patch("funannotate2.search.json.dump")
    def test_swissprot2tsv(self, mock_json_dump):
        """Test the swissprot2tsv function."""
        # Create mock results
        results = [
            {
                "query": "seq1",
                "accession": "P12345",
                "sp_name": "TEST_HUMAN",
                "description": "Test protein",
                "pident": 80.0,
                "evalue": 1e-10,
                "GN": "TEST",
            }
        ]

        # Call the function
        result = search.swissprot2tsv(results, "output.json", "annotations.txt")

        # Check the result
        assert isinstance(result, dict)
        assert "seq1" in result
        assert "db_xref" in result["seq1"]
        assert result["seq1"]["db_xref"] == ["UniProtKB/Swiss-Prot:P12345"]
        assert "note" in result["seq1"]
        assert "80.0% identical to TEST_HUMAN Test protein" in result["seq1"]["note"][0]
        assert mock_json_dump.call_count == 1

    @patch("funannotate2.search.os.path.isdir")
    @patch("funannotate2.search.os.path.isfile")
    @patch("funannotate2.search.runbusco")
    def test_busco_search(self, mock_runbusco, mock_isfile, mock_isdir):
        """Test the busco_search function."""
        # Set up mocks
        mock_isdir.return_value = True
        mock_isfile.return_value = True
        mock_runbusco.return_value = (
            "/path/to/busco/results/full_table.tsv",
            [],
            {},
            {},
        )

        # Call the function
        result = search.busco_search("test_protein.fasta", "fungi", cpus=1)

        # Check the result
        assert result is not None
        assert isinstance(result, str)
        assert result == "/path/to/busco/results/full_table.tsv"
        assert mock_runbusco.call_count == 1
        mock_runbusco.assert_called_once_with(
            "test_protein.fasta",
            "fungi",
            mode="proteins",
            cpus=1,
            logger=sys.stderr,
            verbosity=0,
        )

    @patch("funannotate2.search.os.path.basename")
    @patch("funannotate2.search.os.listdir")
    @patch("funannotate2.search.os.path.join")
    @patch("funannotate2.search.json.dump")
    @patch(
        "builtins.open",
        new_callable=mock_open,
        read_data="BUSCOxyz\tTest BUSCO\thttp://example.com\n",
    )
    def test_busco2tsv(
        self, mock_file, mock_json_dump, mock_path_join, mock_listdir, mock_basename
    ):
        """Test the busco2tsv function."""
        # Set up mocks
        mock_basename.return_value = "test_db"
        mock_listdir.return_value = ["links_to_ODB.txt"]
        mock_path_join.return_value = "/path/to/db/links_to_ODB.txt"

        # Create mock results
        results = {
            "BUSCOxyz": {
                "name": "BUSCOxyz",
                "hit": "seq1",
                "bitscore": 100.0,
                "evalue": 1e-10,
                "domains": [],
                "length": 100,
                "status": "complete",
            }
        }

        # Call the function
        result = search.busco2tsv(
            results, "/path/to/db", "output.json", "annotations.txt"
        )

        # Check the result
        assert isinstance(result, dict)
        assert "seq1" in result
        assert "note" in result["seq1"]
        assert "BUSCO:BUSCOxyz [test_db] Test BUSCO" in result["seq1"]["note"][0]
        assert mock_json_dump.call_count == 1

        # Verify that the files were opened
        mock_file.assert_any_call("/path/to/db/links_to_ODB.txt", "r")
        mock_file.assert_any_call("output.json", "w")
        mock_file.assert_any_call("annotations.txt", "w")
