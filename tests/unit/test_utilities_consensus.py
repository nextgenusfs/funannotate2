"""
Unit tests for the consensus re-use helpers in the utilities module.
"""

import json

from funannotate2.utilities import (
    consensus_needs_rebuild,
    remove_files,
    write_consensus_weights,
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


class TestWriteConsensusWeights:
    """Tests for write_consensus_weights."""

    def test_record_contents(self, tmp_path):
        """The record stores weights as a tool to integer mapping."""
        record = tmp_path / "consensus.predictions.weights.json"
        write_consensus_weights(str(record), ["augustus:2", "snap:1"])
        assert json.loads(record.read_text()) == {"augustus": 2, "snap": 1}


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
