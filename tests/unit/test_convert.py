"""
Unit tests for the convert module.
"""

import pytest

# Skipping tests for now as the functions don't exist yet
# from gfftk.convert import (
#     gff2bed,
#     gff2gtf,
#     gff2tbl,
# )


class TestConvert:
    """Tests for conversion functions."""

    @pytest.mark.skip(reason="Function not implemented yet")
    def test_gff2bed(self, sample_gff_file, temp_dir):
        """Test converting GFF to BED format."""
        pass

    @pytest.mark.skip(reason="Function not implemented yet")
    def test_gff2gtf(self, sample_gff_file, temp_dir):
        """Test converting GFF to GTF format."""
        pass

    @pytest.mark.skip(reason="Function not implemented yet")
    def test_gff2tbl(self, sample_gff_file, sample_fasta_file, temp_dir):
        """Test converting GFF to TBL format."""
        pass
