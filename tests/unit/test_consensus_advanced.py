"""
Advanced unit tests for the consensus module.
"""
import os
import tempfile
from gfftk.consensus import (
    consensus_gene_models,
    consensus_transcripts,
    consensus_exons,
    get_overlap,
    contained,
)


class TestConsensusAdvanced:
    """Advanced tests for consensus functions."""

    def test_consensus_exons(self):
        """Test generating consensus exons."""
        # Create sample exons
        exons1 = [[1, 100], [200, 300], [400, 500]]
        exons2 = [[1, 100], [210, 310], [400, 500]]

        # Generate consensus exons
        result = consensus_exons([exons1, exons2])

        # Check the result
        assert len(result) == 3
        assert result[0] == [1, 100]
        assert result[2] == [400, 500]
        # The middle exon should be a consensus of [200, 300] and [210, 310]
        assert result[1][0] >= 200 and result[1][0] <= 210
        assert result[1][1] >= 300 and result[1][1] <= 310

    def test_consensus_transcripts(self):
        """Test generating consensus transcripts."""
        # Create sample transcripts
        transcript1 = {
            "exon": [[1, 100], [200, 300], [400, 500]],
            "CDS": [[1, 100], [200, 300], [400, 450]],
            "strand": "+",
            "score": 1.0,
        }

        transcript2 = {
            "exon": [[1, 100], [210, 310], [400, 500]],
            "CDS": [[1, 100], [210, 310], [400, 460]],
            "strand": "+",
            "score": 2.0,
        }

        # Generate consensus transcript
        result = consensus_transcripts([transcript1, transcript2])

        # Check the result
        assert "exon" in result
        assert "CDS" in result
        assert "strand" in result
        assert result["strand"] == "+"
        assert len(result["exon"]) == 3
        assert len(result["CDS"]) == 3

        # The consensus should favor transcript2 due to higher score
        assert result["CDS"][2][1] == 460

    def test_consensus_gene_models_basic(self):
        """Test basic consensus gene model generation."""
        # Create a GFF file for source 1
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".gff3"
        ) as gff1:
            gff1.write("##gff-version 3\n")
            gff1.write(
                "contig1\tsource1\tgene\t1\t1000\t.\t+\t.\tID=gene1;Name=test_gene\n"
            )
            gff1.write(
                "contig1\tsource1\tmRNA\t1\t1000\t.\t+\t.\tID=mRNA1;Parent=gene1;Name=test_mrna\n"
            )
            gff1.write(
                "contig1\tsource1\texon\t1\t500\t.\t+\t.\tID=exon1;Parent=mRNA1\n"
            )
            gff1.write(
                "contig1\tsource1\texon\t600\t1000\t.\t+\t.\tID=exon2;Parent=mRNA1\n"
            )
            gff1.write("contig1\tsource1\tCDS\t1\t500\t.\t+\t0\tID=cds1;Parent=mRNA1\n")
            gff1.write(
                "contig1\tsource1\tCDS\t600\t1000\t.\t+\t0\tID=cds2;Parent=mRNA1\n"
            )
            gff1_name = gff1.name

        # Create a GFF file for source 2
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".gff3"
        ) as gff2:
            gff2.write("##gff-version 3\n")
            gff2.write(
                "contig1\tsource2\tgene\t1\t1000\t.\t+\t.\tID=gene2;Name=test_gene\n"
            )
            gff2.write(
                "contig1\tsource2\tmRNA\t1\t1000\t.\t+\t.\tID=mRNA2;Parent=gene2;Name=test_mrna\n"
            )
            gff2.write(
                "contig1\tsource2\texon\t1\t500\t.\t+\t.\tID=exon3;Parent=mRNA2\n"
            )
            gff2.write(
                "contig1\tsource2\texon\t600\t1000\t.\t+\t.\tID=exon4;Parent=mRNA2\n"
            )
            gff2.write("contig1\tsource2\tCDS\t1\t500\t.\t+\t0\tID=cds3;Parent=mRNA2\n")
            gff2.write(
                "contig1\tsource2\tCDS\t600\t1000\t.\t+\t0\tID=cds4;Parent=mRNA2\n"
            )
            gff2_name = gff2.name

        # Create a temporary FASTA file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".fasta"
        ) as temp_fasta:
            temp_fasta.write(">contig1\n")
            temp_fasta.write("A" * 1000 + "\n")
            temp_fasta_name = temp_fasta.name

        # Create a temporary output file
        with tempfile.NamedTemporaryFile(
            mode="w", delete=False, suffix=".gff3"
        ) as temp_out:
            temp_out_name = temp_out.name

        try:
            # Set up weights
            weights = {"source1": 1, "source2": 2}

            # Generate consensus gene models
            result = consensus_gene_models(
                [gff1_name, gff2_name],
                temp_fasta_name,
                output=temp_out_name,
                weights=weights,
                threshold=3,
            )

            # Check that the output file exists
            assert os.path.exists(temp_out_name)

            # Check the result dictionary
            assert isinstance(result, dict)
            assert len(result) > 0

            # Check the content of the output file
            with open(temp_out_name, "r") as f:
                content = f.read()

            # Basic checks on the content
            assert "##gff-version 3" in content
            assert "gene" in content
            assert "mRNA" in content
            assert "exon" in content
            assert "CDS" in content
        finally:
            # Clean up
            for filename in [gff1_name, gff2_name, temp_fasta_name, temp_out_name]:
                if os.path.exists(filename):
                    os.unlink(filename)

    def test_get_overlap_edge_cases(self):
        """Test get_overlap function with edge cases."""
        # Test with identical ranges
        a = [10, 20]
        b = [10, 20]
        overlap = get_overlap(a, b)
        assert overlap == 10  # Overlap is 10-20 = 10

        # Test with zero-length ranges
        a = [10, 10]
        b = [10, 10]
        overlap = get_overlap(a, b)
        assert overlap == 0  # Zero-length ranges don't overlap

        # Test with negative coordinates
        a = [-10, -5]
        b = [-8, -3]
        overlap = get_overlap(a, b)
        assert overlap == 3  # Overlap is -8 to -5 = 3

        # Test with invalid ranges (start > end)
        a = [20, 10]
        b = [15, 25]
        overlap = get_overlap(a, b)
        assert overlap == 0  # Invalid ranges don't overlap

    def test_contained_edge_cases(self):
        """Test contained function with edge cases."""
        # Test with identical ranges
        a = [10, 20]
        b = [10, 20]
        result = contained(a, b)
        assert result is True  # Identical ranges are contained

        # Test with zero-length ranges
        a = [10, 10]
        b = [10, 10]
        result = contained(a, b)
        assert result is True  # Identical zero-length ranges are contained

        # Test with negative coordinates
        a = [-10, -5]
        b = [-12, -3]
        result = contained(a, b)
        assert result is True  # [-10, -5] is contained in [-12, -3]

        # Test with invalid ranges (start > end)
        a = [20, 10]
        b = [5, 25]
        result = contained(a, b)
        assert result is False  # Invalid ranges are not contained
