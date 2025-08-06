"""
Comprehensive unit tests for the annotate module.
"""

import os
import tempfile
from unittest.mock import MagicMock, mock_open, patch

import pytest

import funannotate2.annotate as annotate


class TestAnnotateComprehensive:
    """Comprehensive tests for the annotate module."""

    @pytest.fixture
    def mock_args(self):
        """Create mock arguments for testing."""
        args = MagicMock()
        args.input_dir = None
        args.gff3 = "annotations.gff3"
        args.fasta = "genome.fasta"
        args.species = "Aspergillus fumigatus"
        args.strain = "Af293"
        args.out = "output_dir"
        args.cpus = 1
        args.tmpdir = None
        args.logfile = None
        args.busco_db = "fungi"
        args.busco_seed_species = None
        args.pfam = True
        args.dbcan = True
        args.merops = True
        args.swissprot = True
        args.force = False
        args.eggnog = False
        args.antismash = False
        args.iprscan = False
        args.phobius = False
        args.signalp = False
        args.isolate = None
        args.rename = None
        args.busco = True
        args.database = None
        return args

    @pytest.fixture
    def mock_genes(self):
        """Create a mock genes dictionary for testing."""
        return {
            "gene1": {
                "contig": "contig1",
                "type": ["mRNA"],
                "location": (100, 1000),
                "strand": "+",
                "ids": ["mrna1"],
                "mRNA": [[(100, 200), (300, 400), (500, 1000)]],
                "CDS": [[(100, 200), (300, 400), (500, 1000)]],
                "phase": [[0, 0, 0]],
                "5UTR": [[]],
                "3UTR": [[]],
                "protein": ["MSTLQRPAAFRTRAPPILSADDFADGPPPGGRKRTTFTSA"],
                "transcript": [
                    "ATGAGCACTTTACAGCGACCTGCTGCTTTCCGCACTCGCGCACCTCCCATTCTCTCTGCTGATGACTTTGCTGATGGACCTCCTCCTGGTGGCCGCAAGCGCACCACCTTCACCTCTGCT"
                ],
                "cds_transcript": [
                    "ATGAGCACTTTACAGCGACCTGCTGCTTTCCGCACTCGCGCACCTCCCATTCTCTCTGCTGATGACTTTGCTGATGGACCTCCTCCTGGTGGCCGCAAGCGCACCACCTTCACCTCTGCT"
                ],
                "partialStart": [False],
                "partialStop": [False],
                "product": ["Hypothetical protein"],
                "codon_start": [1],
                "gene_synonym": [],
                "EC_number": [[]],
                "go_terms": [[]],
                "db_xref": [[]],
                "note": [[]],
                "source": "funannotate",
            },
            "gene2": {
                "contig": "contig1",
                "type": ["mRNA"],
                "location": (2000, 3000),
                "strand": "-",
                "ids": ["mrna2"],
                "mRNA": [[(2000, 2200), (2300, 2400), (2500, 3000)]],
                "CDS": [[(2000, 2200), (2300, 2400), (2500, 3000)]],
                "phase": [[0, 0, 0]],
                "5UTR": [[]],
                "3UTR": [[]],
                "protein": ["MSTLQRPAAFRTRAPPILSADDFADGPPPGGRKRTTFTSA"],
                "transcript": [
                    "ATGAGCACTTTACAGCGACCTGCTGCTTTCCGCACTCGCGCACCTCCCATTCTCTCTGCTGATGACTTTGCTGATGGACCTCCTCCTGGTGGCCGCAAGCGCACCACCTTCACCTCTGCT"
                ],
                "cds_transcript": [
                    "ATGAGCACTTTACAGCGACCTGCTGCTTTCCGCACTCGCGCACCTCCCATTCTCTCTGCTGATGACTTTGCTGATGGACCTCCTCCTGGTGGCCGCAAGCGCACCACCTTCACCTCTGCT"
                ],
                "partialStart": [False],
                "partialStop": [False],
                "product": ["Hypothetical protein"],
                "codon_start": [1],
                "gene_synonym": [],
                "EC_number": [[]],
                "go_terms": [[]],
                "db_xref": [[]],
                "note": [[]],
                "source": "funannotate",
            },
        }

    def test_sortDict(self):
        """Test the _sortDict function."""
        # Create a test dictionary
        d = ("gene1", {"location": (100, 200)})

        # Call the function
        result = annotate._sortDict(d)

        # Check the result
        assert result == (100, 200)

    def test_naming_slug(self):
        """Test the naming_slug function."""
        # The naming_slug function is imported from utilities, not defined in annotate
        # So we need to import it directly from utilities
        from funannotate2.utilities import naming_slug

        # Test with species and strain
        result = naming_slug("Aspergillus fumigatus", "Af293")
        # The actual implementation preserves case
        assert result == "Aspergillus_fumigatus_Af293"

        # Test with species only
        result = naming_slug("Aspergillus fumigatus", None)
        # The actual implementation preserves case
        assert result == "Aspergillus_fumigatus"

        # Test with species and isolate
        # The actual implementation doesn't support the isolate parameter
        # It only uses species and strain

    @pytest.mark.skip(reason="check_inputs function not implemented in annotate module")
    @patch("funannotate2.annotate.os.path.isdir")
    @patch("funannotate2.annotate.os.path.isfile")
    @patch("funannotate2.annotate.checkfile")
    def test_check_inputs(self, mock_checkfile, mock_isfile, mock_isdir, mock_args):
        """Test the check_inputs function."""
        # Set up mocks
        mock_isdir.return_value = True
        mock_isfile.return_value = True
        mock_checkfile.return_value = True

        # The check_inputs function is not implemented in the annotate module
        # Instead, the input validation is done directly in the annotate function
        # This test is skipped until the function is implemented

    @pytest.mark.skip(reason="check_inputs function not implemented in annotate module")
    @patch("funannotate2.annotate.os.path.isdir")
    @patch("funannotate2.annotate.os.path.isfile")
    @patch("funannotate2.annotate.checkfile")
    def test_check_inputs_missing_files(
        self, mock_checkfile, mock_isfile, mock_isdir, mock_args
    ):
        """Test the check_inputs function with missing files."""
        # Set up mocks
        mock_isdir.return_value = True
        mock_isfile.return_value = False
        mock_checkfile.return_value = False

        # The check_inputs function is not implemented in the annotate module
        # Instead, the input validation is done directly in the annotate function
        # This test is skipped until the function is implemented

    @pytest.mark.skip(
        reason="find_input_files function not implemented in annotate module"
    )
    @patch("funannotate2.annotate.find_files")
    def test_find_input_files(self, mock_find_files, mock_args):
        """Test the find_input_files function."""
        # Set up mocks
        mock_find_files.side_effect = [
            ["input_dir/predict_results/funannotate_predict.gff3"],  # GFF files
            ["input_dir/predict_results/genome.fasta"],  # FASTA files
        ]

        # Set input_dir
        mock_args.input_dir = "input_dir"

        # The find_input_files function is not implemented in the annotate module
        # Instead, the input file finding is done directly in the annotate function
        # This test is skipped until the function is implemented

    def test_parse_annotations(self):
        """Test the parse_annotations function."""
        # Create a mock annotation file
        with tempfile.NamedTemporaryFile(mode="w", delete=False) as temp:
            temp.write("gene1\tproduct\tHypothetical protein\n")
            temp.write("gene1\tdb_xref\tUniProtKB/Swiss-Prot:P12345\n")
            temp.write("gene2\tproduct\tHypothetical protein\n")
            temp.write("gene2\tdb_xref\tUniProtKB/Swiss-Prot:P67890\n")
            temp_name = temp.name

        try:
            # Call the function
            result, parse_errors = annotate.parse_annotations(temp_name)

            # Check the result
            assert isinstance(result, dict)
            assert isinstance(parse_errors, dict)
            assert "gene1" in result
            assert "gene2" in result
            assert "product" in result["gene1"]
            assert "db_xref" in result["gene1"]
            assert result["gene1"]["product"] == ["Hypothetical protein"]
            assert result["gene1"]["db_xref"] == ["UniProtKB/Swiss-Prot:P12345"]

            # Check that parsing was successful (no errors)
            assert len(parse_errors["malformed_lines"]) == 0
            assert parse_errors["parsed_lines"] == 4  # 4 lines in the test file
            assert parse_errors["total_lines"] == 4
        finally:
            # Clean up
            os.unlink(temp_name)

    @pytest.mark.skip(
        reason="write_output_files function not implemented in annotate module"
    )
    @patch("funannotate2.annotate.os.path.join")
    @patch("funannotate2.annotate.os.makedirs")
    @patch("funannotate2.annotate.dict2tbl")
    @patch("funannotate2.annotate.table2asn")
    @patch("funannotate2.annotate.dict2gff3")
    @patch("funannotate2.annotate._dict2proteins")
    @patch("funannotate2.annotate._dict2transcripts")
    @patch("funannotate2.annotate.shutil.copy2")
    @patch("funannotate2.annotate.annotation_stats")
    @patch("builtins.open", new_callable=mock_open)
    @patch("json.dump")
    def test_write_output_files(
        self,
        mock_json_dump,
        mock_open,
        mock_annotation_stats,
        mock_copy2,
        mock_dict2transcripts,
        mock_dict2proteins,
        mock_dict2gff3,
        mock_table2asn,
        mock_dict2tbl,
        mock_makedirs,
        mock_path_join,
        mock_args,
        mock_genes,
    ):
        """Test the write_output_files function."""
        # Set up mocks
        mock_path_join.side_effect = lambda *args: "/".join(args)
        mock_annotation_stats.return_value = {
            "total_genes": 2,
            "protein_coding": 2,
            "mean_gene_length": 950,
            "mean_CDS_length": 950,
            "mean_exons": 3,
        }

        # The write_output_files function is not implemented in the annotate module
        # Instead, the output file writing is done directly in the annotate function
        # This test is skipped until the function is implemented
