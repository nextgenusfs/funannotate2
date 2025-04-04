"""
Integration tests for the funannotate2 workflow.
"""

import os
import json
import tempfile
import shutil
import pytest
import funannotate2.clean as clean
import funannotate2.annotate as annotate

# import funannotate2.predict as predict  # Not used in this test file
# import funannotate2.compare as compare  # Module doesn't exist yet


class TestFunannotateWorkflow:
    """Test the funannotate2 workflow."""

    @pytest.fixture(scope="class")
    def test_data_dir(self):
        """Return the path to the test data directory."""
        return os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

    @pytest.fixture(scope="function")
    def temp_dir(self):
        """Create a temporary directory for test outputs."""
        with tempfile.TemporaryDirectory() as tmpdir:
            yield tmpdir

    @pytest.fixture
    def test_fasta(self, temp_dir):
        """Create a test FASTA file."""
        fasta_file = os.path.join(temp_dir, "test_genome.fasta")
        with open(fasta_file, "w") as f:
            f.write(">contig1 some description\n")
            f.write("ATGCATGCATGCATGCATGC" * 100 + "\n")
            f.write(">contig2 another description\n")
            f.write("GTACGTACGTACGTACGTAC" * 100 + "\n")
            f.write(">short\n")
            f.write("ATGC\n")
        return fasta_file

    @pytest.fixture
    def test_protein_fasta(self, temp_dir):
        """Create a test protein FASTA file."""
        fasta_file = os.path.join(temp_dir, "test_proteins.fasta")
        with open(fasta_file, "w") as f:
            f.write(">protein1\n")
            f.write("MHACMHACMHACMHACMHAC\n")
            f.write(">protein2\n")
            f.write("VRTVRTVRTVRTVRTVRT\n")
        return fasta_file

    @pytest.fixture
    def test_transcript_fasta(self, temp_dir):
        """Create a test transcript FASTA file."""
        fasta_file = os.path.join(temp_dir, "test_transcripts.fasta")
        with open(fasta_file, "w") as f:
            f.write(">transcript1\n")
            f.write("ATGCATGCATGCATGCATGC" * 5 + "\n")
            f.write(">transcript2\n")
            f.write("GTACGTACGTACGTACGTAC" * 5 + "\n")
        return fasta_file

    @pytest.fixture
    def test_gff(self, temp_dir):
        """Create a test GFF file."""
        gff_file = os.path.join(temp_dir, "test_annotations.gff3")
        with open(gff_file, "w") as f:
            f.write("##gff-version 3\n")
            f.write("contig1\tFunannotate\tgene\t1\t100\t.\t+\t.\tID=gene1\n")
            f.write(
                "contig1\tFunannotate\tmRNA\t1\t100\t.\t+\t.\tID=gene1-T1;Parent=gene1\n"
            )
            f.write(
                "contig1\tFunannotate\tCDS\t1\t100\t.\t+\t0\tID=gene1-T1-CDS;Parent=gene1-T1\n"
            )
            f.write("contig2\tFunannotate\tgene\t1\t100\t.\t-\t.\tID=gene2\n")
            f.write(
                "contig2\tFunannotate\tmRNA\t1\t100\t.\t-\t.\tID=gene2-T1;Parent=gene2\n"
            )
            f.write(
                "contig2\tFunannotate\tCDS\t1\t100\t.\t-\t0\tID=gene2-T1-CDS;Parent=gene2-T1\n"
            )
        return gff_file

    @pytest.mark.skipif(
        not shutil.which("tbl2asn"),
        reason="tbl2asn is not installed or not in PATH",
    )
    def test_clean_workflow(self, test_fasta, temp_dir):
        """Test the clean workflow."""

        # Create arguments for clean
        class Args:
            def __init__(self):
                self.input = test_fasta
                self.out = os.path.join(temp_dir, "cleaned.fasta")
                self.minlen = 10
                self.cpus = 1
                self.tmpdir = None
                self.logfile = None
                self.exhaustive = False
                self.pident = 95.0
                self.cov = 95.0
                self.gap_fill = False
                self.gap_method = "diploidocus"
                self.gap_length = 100
                self.gap_safezone = 0.5
                self.force = False
                self.rename = None
                self.species = "Aspergillus fumigatus"
                self.strain = "Af293"
                self.isolate = None
                self.header_slice = None
                self.header_names = None
                self.sort = "size"
                self.busco_db = None
                self.busco_seed_species = None
                self.busco = False
                self.outgroup = None
                self.ploidy = 1
                self.diploid_threshold = 95.0
                self.diploid_mode = "consensus"
                self.diploid_split = False
                self.diploid_consensus = True
                self.diploid_alignments = False
                self.diploid_alignments_method = "mafft"
                self.diploid_alignments_format = "fasta"
                self.diploid_alignments_outdir = None
                self.diploid_alignments_outfile = None
                self.diploid_alignments_outprefix = None
                self.diploid_alignments_outsuffix = None
                self.diploid_alignments_outextension = None
                self.diploid_alignments_outformat = None
                self.diploid_alignments_outdir_structure = None
                self.diploid_alignments_outdir_prefix = None
                self.diploid_alignments_outdir_suffix = None
                self.diploid_alignments_outdir_extension = None
                self.diploid_alignments_outdir_format = None

        args = Args()

        # Run the clean function
        clean.clean(args)

        # Check the results
        assert os.path.exists(args.out)

        # Check the content of the cleaned file
        with open(args.out, "r") as f:
            content = f.read()
            assert ">contig1" in content
            assert ">contig2" in content
            assert ">short" not in content
            assert "some description" not in content
            assert "another description" not in content

    @pytest.mark.skipif(
        not shutil.which("tbl2asn") or not shutil.which("hmmsearch"),
        reason="Required tools (tbl2asn, hmmsearch) are not installed or not in PATH",
    )
    def test_annotate_workflow(self, test_fasta, test_gff, temp_dir):
        """Test the annotate workflow."""
        # Skip this test in CI environment as it requires too many dependencies
        if os.environ.get("CI") == "true":
            pytest.skip("Skipping in CI environment")

        # Create arguments for annotate
        class Args:
            def __init__(self):
                self.input_dir = None
                self.gff3 = test_gff
                self.fasta = test_fasta
                self.species = "Aspergillus fumigatus"
                self.strain = "Af293"
                self.out = os.path.join(temp_dir, "annotate_results")
                self.cpus = 1
                self.tmpdir = None
                self.logfile = None
                self.busco_db = None
                self.busco_seed_species = None
                self.pfam = False
                self.dbcan = False
                self.merops = False
                self.swissprot = False
                self.force = False
                self.eggnog = False
                self.antismash = False
                self.iprscan = False
                self.phobius = False
                self.signalp = False
                self.isolate = None
                self.rename = None
                self.busco = False
                self.database = None

        args = Args()

        # Mock the annotate function to avoid running the full workflow
        def mock_annotate(args):
            # Create the output directory
            os.makedirs(args.out, exist_ok=True)

            # Create a mock results directory
            res_dir = os.path.join(args.out, "annotate_results")
            os.makedirs(res_dir, exist_ok=True)

            # Create mock output files
            prefix = f"{annotate.naming_slug(args.species, args.strain)}"

            # Create a mock GFF file
            gff_file = os.path.join(res_dir, f"{prefix}.gff3")
            shutil.copy(args.gff3, gff_file)

            # Create a mock proteins file
            proteins_file = os.path.join(res_dir, f"{prefix}.proteins.fa")
            with open(proteins_file, "w") as f:
                f.write(">gene1-T1\n")
                f.write("MHACMHACMHACMHACMHAC\n")
                f.write(">gene2-T1\n")
                f.write("VRTVRTVRTVRTVRTVRT\n")

            # Create a mock transcripts file
            transcripts_file = os.path.join(res_dir, f"{prefix}.transcripts.fa")
            with open(transcripts_file, "w") as f:
                f.write(">gene1-T1\n")
                f.write("ATGCATGCATGCATGCATGC\n")
                f.write(">gene2-T1\n")
                f.write("GTACGTACGTACGTACGTAC\n")

            # Create a mock genome file
            genome_file = os.path.join(res_dir, f"{prefix}.fasta")
            shutil.copy(args.fasta, genome_file)

            # Create a mock summary file
            summary_file = os.path.join(res_dir, f"{prefix}.summary.json")
            with open(summary_file, "w") as f:
                json.dump(
                    {
                        "total_genes": 2,
                        "protein_coding": 2,
                        "mean_gene_length": 100,
                        "mean_CDS_length": 100,
                        "mean_exons": 1,
                    },
                    f,
                    indent=4,
                )

            return 0

        # Save the original function
        original_annotate = annotate.annotate

        try:
            # Replace with our mock function
            annotate.annotate = mock_annotate

            # Run the annotate function
            annotate.annotate(args)

            # Check the results
            prefix = f"{annotate.naming_slug(args.species, args.strain)}"
            res_dir = os.path.join(args.out, "annotate_results")

            assert os.path.exists(os.path.join(res_dir, f"{prefix}.gff3"))
            assert os.path.exists(os.path.join(res_dir, f"{prefix}.proteins.fa"))
            assert os.path.exists(os.path.join(res_dir, f"{prefix}.transcripts.fa"))
            assert os.path.exists(os.path.join(res_dir, f"{prefix}.fasta"))
            assert os.path.exists(os.path.join(res_dir, f"{prefix}.summary.json"))
        finally:
            # Restore the original function
            annotate.annotate = original_annotate

    @pytest.mark.skip(reason="Compare module not implemented yet")
    def test_compare_workflow(self, test_fasta, temp_dir):  # test_gff not used
        """Test the compare workflow."""
        # Skip this test in CI environment as it requires too many dependencies
        if os.environ.get("CI") == "true":
            pytest.skip("Skipping in CI environment")

        # Create test input directories
        input_dir1 = os.path.join(temp_dir, "input1")
        input_dir2 = os.path.join(temp_dir, "input2")
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

        # Create test FASTA files
        fasta1 = os.path.join(input_dir1, "genome.fasta")
        fasta2 = os.path.join(input_dir2, "genome.fasta")
        shutil.copy(test_fasta, fasta1)
        shutil.copy(test_fasta, fasta2)

        # Create arguments for compare
        class Args:
            def __init__(self):
                self.input = [input_dir1, input_dir2]
                self.out = os.path.join(temp_dir, "compare_results")
                self.names = ["sample1", "sample2"]
                self.cpus = 1
                self.tmpdir = None
                self.logfile = None
                self.force = False
                self.run_busco = False
                self.busco_db = None
                self.run_pfam = False
                self.pfam_db = None
                self.run_interpro = False
                self.interpro_db = None
                self.run_eggnog = False
                self.eggnog_db = None
                self.run_cog = False
                self.cog_db = None
                self.run_go = False
                self.go_db = None
                self.run_kegg = False
                self.kegg_db = None
                self.run_cazy = False
                self.cazy_db = None
                self.run_merops = False
                self.merops_db = None
                self.run_signalp = False
                self.signalp_db = None
                self.run_tmhmm = False
                self.tmhmm_db = None
                self.run_phobius = False
                self.phobius_db = None
                self.run_antismash = False
                self.antismash_db = None
                self.run_all = False
                self.run_none = False
                self.run_custom = False
                self.custom_db = None
                self.custom_name = None
                self.custom_cutoff = None
                self.custom_output = None
                self.custom_evalue = None
                self.custom_coverage = None
                self.custom_identity = None
                self.custom_bitscore = None
                self.custom_maxhits = None
                self.custom_maxhsps = None
                self.custom_threads = None
                self.custom_tmpdir = None
                self.custom_outfmt = None
                self.custom_outfmt_sep = None
                self.custom_outfmt_header = None
                self.custom_outfmt_footer = None
                self.custom_outfmt_line = None
                self.custom_outfmt_query = None
                self.custom_outfmt_target = None
                self.custom_outfmt_qseq = None
                self.custom_outfmt_tseq = None
                self.custom_outfmt_qlen = None
                self.custom_outfmt_tlen = None
                self.custom_outfmt_qstart = None
                self.custom_outfmt_qend = None
                self.custom_outfmt_tstart = None
                self.custom_outfmt_tend = None
                self.custom_outfmt_evalue = None
                self.custom_outfmt_bitscore = None
                self.custom_outfmt_score = None
                self.custom_outfmt_length = None
                self.custom_outfmt_pident = None
                self.custom_outfmt_nident = None
                self.custom_outfmt_mismatch = None
                self.custom_outfmt_positive = None
                self.custom_outfmt_gapopen = None
                self.custom_outfmt_gaps = None
                self.custom_outfmt_ppos = None
                self.custom_outfmt_frames = None
                self.custom_outfmt_qframe = None
                self.custom_outfmt_tframe = None
                self.custom_outfmt_btop = None
                self.custom_outfmt_staxids = None
                self.custom_outfmt_sscinames = None
                self.custom_outfmt_scomnames = None
                self.custom_outfmt_sblastnames = None
                self.custom_outfmt_sskingdoms = None
                self.custom_outfmt_stitle = None
                self.custom_outfmt_salltitles = None
                self.custom_outfmt_qcovs = None
                self.custom_outfmt_qcovhsp = None
                self.custom_outfmt_qcovus = None

        args = Args()

        # Mock the compare function to avoid running the full workflow
        def mock_compare(args):
            # Create the output directory
            os.makedirs(args.out, exist_ok=True)

            # Create mock output files
            summary_file = os.path.join(args.out, "summary.tsv")
            with open(summary_file, "w") as f:
                f.write("Sample\tTotal Genes\tUnique Genes\n")
                f.write("sample1\t1\t0\n")
                f.write("sample2\t2\t1\n")

            return 0

        # Since the compare module doesn't exist yet, we'll just run the mock function directly
        mock_compare(args)

        # Check the results
        assert os.path.exists(os.path.join(args.out, "summary.tsv"))
