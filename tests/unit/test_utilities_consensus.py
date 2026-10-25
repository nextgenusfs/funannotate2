"""
Unit tests for the consensus re-use helpers in the utilities module.
"""

import json

from funannotate2.utilities import (
    consensus_needs_rebuild,
    remove_files,
    utr_refinement_needs_rebuild,
    write_consensus_weights,
    write_utr_refinement_record,
)

WEIGHTS = ["augustus-hiq:4", "augustus:2", "helixer:2", "glimmerhmm:1", "snap:1"]


def make_consensus(tmp_path, content="##gff-version 3\n"):
    consensus = tmp_path / "consensus.predictions.gff3"
    consensus.write_text(content)
    return str(consensus)


class TestConsensusNeedsRebuild:
    """Tests for consensus_needs_rebuild."""

    def test_no_consensus(self, tmp_path):
        """A missing consensus file needs to be built."""
        rebuild, reason = consensus_needs_rebuild(
            str(tmp_path / "consensus.predictions.gff3"),
            str(tmp_path / "consensus.predictions.weights.json"),
            WEIGHTS,
        )
        assert rebuild
        assert "no existing consensus" in reason

    def test_no_record(self, tmp_path):
        """A consensus without a weights record is rebuilt."""
        consensus = make_consensus(tmp_path)
        rebuild, reason = consensus_needs_rebuild(
            consensus, str(tmp_path / "consensus.predictions.weights.json"), WEIGHTS
        )
        assert rebuild
        assert "no record" in reason

    def test_same_weights(self, tmp_path):
        """A consensus built with the same weights is re-used."""
        consensus = make_consensus(tmp_path)
        record = str(tmp_path / "consensus.predictions.weights.json")
        write_consensus_weights(record, WEIGHTS)
        rebuild, _ = consensus_needs_rebuild(consensus, record, WEIGHTS)
        assert not rebuild

    def test_same_weights_different_order(self, tmp_path):
        """Weight order does not matter."""
        consensus = make_consensus(tmp_path)
        record = str(tmp_path / "consensus.predictions.weights.json")
        write_consensus_weights(record, WEIGHTS)
        rebuild, _ = consensus_needs_rebuild(consensus, record, list(reversed(WEIGHTS)))
        assert not rebuild

    def test_changed_weight(self, tmp_path):
        """A changed weight triggers a rebuild and is named in the reason."""
        consensus = make_consensus(tmp_path)
        record = str(tmp_path / "consensus.predictions.weights.json")
        write_consensus_weights(record, WEIGHTS)
        new_weights = [w if w != "helixer:2" else "helixer:4" for w in WEIGHTS]
        rebuild, reason = consensus_needs_rebuild(consensus, record, new_weights)
        assert rebuild
        assert "helixer: 2 -> 4" in reason

    def test_added_and_removed_tools(self, tmp_path):
        """Adding or removing a tool triggers a rebuild."""
        consensus = make_consensus(tmp_path)
        record = str(tmp_path / "consensus.predictions.weights.json")
        write_consensus_weights(record, WEIGHTS)
        new_weights = [w for w in WEIGHTS if not w.startswith("snap")] + ["miniprot-gene:1"]
        rebuild, reason = consensus_needs_rebuild(consensus, record, new_weights)
        assert rebuild
        assert "snap: 1 -> none" in reason
        assert "miniprot-gene: none -> 1" in reason

    def test_malformed_record(self, tmp_path):
        """An unreadable record is treated as missing rather than raising."""
        consensus = make_consensus(tmp_path)
        record = tmp_path / "consensus.predictions.weights.json"
        for content in ["{not json", "[1, 2]", '{"augustus": "high"}']:
            record.write_text(content)
            rebuild, reason = consensus_needs_rebuild(consensus, str(record), WEIGHTS)
            assert rebuild
            assert "no record" in reason

    def test_empty_consensus(self, tmp_path):
        """An empty consensus file is rebuilt even with a matching record."""
        consensus = make_consensus(tmp_path, content="")
        record = str(tmp_path / "consensus.predictions.weights.json")
        write_consensus_weights(record, WEIGHTS)
        rebuild, _ = consensus_needs_rebuild(consensus, record, WEIGHTS)
        assert rebuild

    def test_changed_max_intron(self, tmp_path):
        """A change in max_intron triggers a rebuild."""
        consensus = make_consensus(tmp_path)
        record = str(tmp_path / "consensus.predictions.weights.json")
        write_consensus_weights(record, WEIGHTS, max_intron=3000)
        rebuild, _ = consensus_needs_rebuild(consensus, record, WEIGHTS, max_intron=3000)
        assert not rebuild
        rebuild, reason = consensus_needs_rebuild(consensus, record, WEIGHTS, max_intron=5000)
        assert rebuild
        assert "max_intron changed (3000 -> 5000)" in reason

    def test_abinitio_preds_added_removed_modified(self, tmp_path):
        """Adding, removing, or modifying a prediction GFF3 triggers a rebuild."""
        consensus = make_consensus(tmp_path)
        record = str(tmp_path / "consensus.predictions.weights.json")
        pred_a = tmp_path / "predictions.augustus.gff3"
        pred_b = tmp_path / "predictions.snap.gff3"
        pred_a.write_text("##gff-version 3\nchr1\taugustus\tgene\t1\t100\t.\t+\t.\tID=g1\n")
        pred_b.write_text("##gff-version 3\nchr1\tsnap\tgene\t1\t100\t.\t+\t.\tID=g2\n")

        write_consensus_weights(record, WEIGHTS, abinitio_preds=[str(pred_a), str(pred_b)])
        rebuild, _ = consensus_needs_rebuild(
            consensus, record, WEIGHTS, abinitio_preds=[str(pred_a), str(pred_b)]
        )
        assert not rebuild

        # Modify pred_a content (e.g. external prediction re-run with updated GFF3)
        pred_a.write_text("##gff-version 3\nchr1\taugustus\tgene\t1\t200\t.\t+\t.\tID=g1\n")
        rebuild, reason = consensus_needs_rebuild(
            consensus, record, WEIGHTS, abinitio_preds=[str(pred_a), str(pred_b)]
        )
        assert rebuild
        assert "modified predictions.augustus.gff3" in reason

        # Remove pred_b and add pred_c
        pred_c = tmp_path / "predictions.helixer.gff3"
        pred_c.write_text("##gff-version 3\n")
        rebuild, reason = consensus_needs_rebuild(
            consensus, record, WEIGHTS, abinitio_preds=[str(pred_a), str(pred_c)]
        )
        assert rebuild
        assert "removed predictions.snap.gff3" in reason
        assert "added predictions.helixer.gff3" in reason

    def test_protein_and_transcript_alignments_changed(self, tmp_path):
        """Changes to protein or transcript alignment files trigger a rebuild."""
        consensus = make_consensus(tmp_path)
        record = str(tmp_path / "consensus.predictions.weights.json")
        prot = tmp_path / "protein-alignments.gff3"
        tran = tmp_path / "transcript-alignments.gff3"
        prot.write_text("##gff-version 3\nchr1\tminiprot\tmRNA\t1\t100\t.\t+\t.\tID=mp1\n")
        tran.write_text("##gff-version 3\nchr1\tgapmm2\tcDNA_match\t1\t100\t.\t+\t.\tID=tr1\n")

        write_consensus_weights(
            record,
            WEIGHTS,
            protein_alignments=[str(prot)],
            transcript_alignments=[str(tran)],
        )
        rebuild, _ = consensus_needs_rebuild(
            consensus,
            record,
            WEIGHTS,
            protein_alignments=[str(prot)],
            transcript_alignments=[str(tran)],
        )
        assert not rebuild

        tran.write_text("##gff-version 3\nchr1\tgapmm2\tcDNA_match\t1\t250\t.\t+\t.\tID=tr1\n")
        rebuild, reason = consensus_needs_rebuild(
            consensus,
            record,
            WEIGHTS,
            protein_alignments=[str(prot)],
            transcript_alignments=[str(tran)],
        )
        assert rebuild
        assert "transcript alignments changed (modified transcript-alignments.gff3)" in reason


class TestWriteConsensusWeights:
    """Tests for write_consensus_weights."""

    def test_record_contents(self, tmp_path):
        """The record stores weights as a tool to integer mapping."""
        record = tmp_path / "consensus.predictions.weights.json"
        write_consensus_weights(str(record), ["augustus:2", "snap:1"])
        assert json.loads(record.read_text()) == {"augustus": 2, "snap": 1}

    def test_record_contents_with_extended_inputs(self, tmp_path):
        """When additional consensus inputs are provided, a structured manifest is written."""
        record = tmp_path / "consensus.predictions.weights.json"
        pred = tmp_path / "predictions.snap.gff3"
        pred.write_text("##gff-version 3\n")
        write_consensus_weights(
            str(record),
            ["augustus:2", "snap:1"],
            max_intron=3000,
            abinitio_preds=[str(pred)],
            protein_alignments=[],
            transcript_alignments=[],
        )
        data = json.loads(record.read_text())
        assert data["weights"] == {"augustus": 2, "snap": 1}
        assert data["max_intron"] == 3000
        assert "predictions.snap.gff3" in data["abinitio_preds"]
        assert data["protein_alignments"] == {}
        assert data["transcript_alignments"] == {}


class TestUtrRefinementNeedsRebuild:
    """Tests for utr_refinement_needs_rebuild and write_utr_refinement_record."""

    def test_utr_refinement_lifecycle(self, tmp_path):
        """UTR refinement is re-used when consensus, BAM file, and library match, and rebuilt otherwise."""
        consensus = make_consensus(tmp_path)
        refined = tmp_path / "consensus.predictions.utr_refined.gff3"
        record = str(tmp_path / "consensus.predictions.utr_refined.json")
        bam = tmp_path / "rnaseq.bam"
        bam.write_bytes(b"BAM\x01mockcontent")

        # Missing refined file -> needs rebuild
        rebuild, reason = utr_refinement_needs_rebuild(
            str(refined), record, consensus, str(bam), "RF"
        )
        assert rebuild
        assert "no existing UTR-refined consensus" in reason

        # Write refined file and record -> re-used
        refined.write_text("##gff-version 3\n")
        write_utr_refinement_record(record, consensus, str(bam), "RF")
        rebuild, _ = utr_refinement_needs_rebuild(
            str(refined), record, consensus, str(bam), "RF"
        )
        assert not rebuild

        # Changing bam_library triggers rebuild
        rebuild, reason = utr_refinement_needs_rebuild(
            str(refined), record, consensus, str(bam), "FR"
        )
        assert rebuild
        assert "bam_library changed (RF -> FR)" in reason

        # Changing base consensus triggers rebuild
        consensus = make_consensus(tmp_path, content="##gff-version 3\n# updated\n")
        rebuild, reason = utr_refinement_needs_rebuild(
            str(refined), record, consensus, str(bam), "RF"
        )
        assert rebuild
        assert "base consensus predictions changed" in reason


class TestRemoveFiles:
    """Tests for remove_files."""

    def test_removes_existing_files(self, tmp_path):
        """All existing files are removed and returned."""
        paths = []
        for name in ["a.gff3", "b.gff3", "c.json"]:
            f = tmp_path / name
            f.write_text("x")
            paths.append(str(f))
        assert remove_files(paths) == paths
        assert not any((tmp_path / name).exists() for name in ["a.gff3", "b.gff3", "c.json"])

    def test_consensus_missing_refined_present(self, tmp_path):
        """A stale refined file is removed even when the consensus was deleted by hand."""
        consensus = tmp_path / "consensus.predictions.gff3"
        refined = tmp_path / "consensus.predictions.utr_refined.gff3"
        refined.write_text("##gff-version 3\n")
        removed = remove_files([str(consensus), str(refined)])
        assert removed == [str(refined)]
        assert not refined.exists()

    def test_nothing_to_remove(self, tmp_path):
        """Missing files are skipped without error."""
        assert remove_files([str(tmp_path / "missing.gff3")]) == []
