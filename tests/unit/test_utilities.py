"""
Unit tests for the utilities module.
"""

import json
import os
from unittest.mock import MagicMock, patch

from funannotate2.utilities import (
    create_tmpdir,
    merge_coordinates,
    naming_slug,
    readBlocks,
    runSubprocess,
)


class TestMergeCoordinates:
    """Tests for the merge_coordinates function."""

    def test_empty_list(self):
        """Test with an empty list."""
        assert merge_coordinates([]) == []

    def test_single_interval(self):
        """Test with a single interval."""
        assert merge_coordinates([[1, 3]]) == [[1, 3]]

    def test_non_overlapping_intervals(self):
        """Test with non-overlapping intervals."""
        assert merge_coordinates([[1, 3], [5, 7], [9, 11]]) == [[1, 3], [5, 7], [9, 11]]

    def test_overlapping_intervals(self):
        """Test with overlapping intervals."""
        assert merge_coordinates([[1, 3], [2, 6], [8, 10], [15, 18]]) == [
            [1, 6],
            [8, 10],
            [15, 18],
        ]

    def test_completely_overlapping_intervals(self):
        """Test with completely overlapping intervals."""
        assert merge_coordinates([[1, 10], [2, 5], [3, 7]]) == [[1, 10]]

    def test_adjacent_intervals(self):
        """Test with adjacent intervals."""
        assert merge_coordinates([[1, 3], [3, 5], [7, 9]]) == [[1, 5], [7, 9]]

    def test_unsorted_intervals(self):
        """Test with unsorted intervals."""
        assert merge_coordinates([[8, 10], [1, 3], [2, 6], [15, 18]]) == [
            [1, 6],
            [8, 10],
            [15, 18],
        ]


class TestNamingSlug:
    """Tests for the naming_slug function."""

    def test_basic_slug(self):
        """Test basic slug generation."""
        assert naming_slug("aspergillus fumigatus", "Af293") == "Aspergillus_fumigatus_Af293"

    def test_lowercase_slug(self):
        """Test lowercase slug generation."""
        assert (
            naming_slug("Aspergillus fumigatus", "Af293", lowercase=True)
            == "aspergillus_fumigatus_Af293"
        )

    def test_no_strain(self):
        """Test slug generation without strain."""
        assert naming_slug("Aspergillus fumigatus", None) == "Aspergillus_fumigatus"

    def test_strain_with_spaces(self):
        """Test slug generation with spaces in strain."""
        assert naming_slug("Aspergillus fumigatus", "Af 293") == "Aspergillus_fumigatus_Af293"


class TestCreateTmpdir:
    """Tests for the create_tmpdir function."""

    def test_create_tmpdir_with_outdir(self, temp_dir):
        """Test creating a temporary directory with an output directory."""
        tmpdir = create_tmpdir(temp_dir, base="test")
        assert os.path.exists(tmpdir)
        assert os.path.isdir(tmpdir)
        assert tmpdir.startswith(os.path.abspath(temp_dir))

    def test_create_tmpdir_with_tmp(self):
        """Test creating a temporary directory in /tmp."""
        tmpdir = create_tmpdir("/tmp", base="test")
        assert os.path.exists(tmpdir)
        assert os.path.isdir(tmpdir)
        assert tmpdir.startswith("/tmp")

    def test_create_tmpdir_without_outdir(self):
        """Test creating a temporary directory without an output directory."""
        tmpdir = create_tmpdir(None, base="test")
        assert os.path.exists(tmpdir)
        assert os.path.isdir(tmpdir)
        # Clean up
        os.rmdir(tmpdir)


class TestReadBlocks:
    """Tests for the readBlocks function."""

    def test_read_blocks(self):
        """Test reading blocks from a source."""
        source = [
            "# Block 1",
            "Line 1",
            "Line 2",
            "# Block 2",
            "Line 3",
            "Line 4",
            "# Block 3",
            "Line 5",
        ]
        blocks = list(readBlocks(source, "#"))
        assert len(blocks) == 3
        assert blocks[0] == ["# Block 1", "Line 1", "Line 2"]
        assert blocks[1] == ["# Block 2", "Line 3", "Line 4"]
        assert blocks[2] == ["# Block 3", "Line 5"]

    def test_read_blocks_empty_source(self):
        """Test reading blocks from an empty source."""
        blocks = list(readBlocks([], "#"))
        assert len(blocks) == 1
        assert blocks[0] == []

    def test_read_blocks_no_pattern(self):
        """Test reading blocks with no pattern match."""
        source = ["Line 1", "Line 2", "Line 3"]
        blocks = list(readBlocks(source, "#"))
        assert len(blocks) == 1
        assert blocks[0] == ["Line 1", "Line 2", "Line 3"]


class TestRunSubprocess:
    """Regression tests for subprocess logging behavior."""

    @patch("funannotate2.utilities.subprocess.Popen")
    @patch("funannotate2.memory.MemoryMonitor")
    def test_monitoring_keeps_memory_report_out_of_callable_logs(
        self, mock_monitor_cls, mock_popen, tmp_path, monkeypatch
    ):
        """Memory monitoring should not emit formatted reports via callable logs."""
        process = MagicMock()
        process.pid = 12345
        process.returncode = 0
        process.stdout = ""
        process.stderr = ""
        mock_popen.return_value = process

        mock_monitor = mock_monitor_cls.return_value
        mock_monitor.monitor_process.return_value = {
            "process_name": "snap-scaffold_1.fasta",
            "duration_seconds": 1.5,
            "peak_rss_mb": 10.0,
            "peak_vms_mb": 20.0,
            "avg_rss_mb": 9.0,
            "avg_vms_mb": 19.0,
            "sample_count": 2,
            "samples": [],
        }

        messages = []
        monkeypatch.setenv("FUNANNOTATE2_OUTPUT_DIR", str(tmp_path))

        runSubprocess(
            ["snap", "input.fasta"],
            messages.append,
            cwd=str(tmp_path),
            monitor_memory=True,
            process_name="snap-scaffold_1.fasta",
        )

        assert messages == []

        memory_log = tmp_path / "logfiles" / "predict-abinitio-memory-monitoring.jsonl"
        assert memory_log.exists()

        records = [json.loads(line) for line in memory_log.read_text().splitlines()]
        assert len(records) == 1
        assert records[0]["process_name"] == "snap-scaffold_1.fasta"
        assert records[0]["tool_name"] == "snap"
        assert records[0]["memory_stats"]["peak_rss_mb"] == 10.0

    @patch("funannotate2.utilities.subprocess.Popen")
    @patch("funannotate2.memory.MemoryMonitor")
    def test_monitoring_does_not_debug_log_formatted_memory_report(
        self, mock_monitor_cls, mock_popen, tmp_path, monkeypatch
    ):
        """Detailed memory reports should stay out of the normal logfile path."""
        process = MagicMock()
        process.pid = 12345
        process.returncode = 0
        process.stdout = ""
        process.stderr = ""
        mock_popen.return_value = process

        mock_monitor = mock_monitor_cls.return_value
        mock_monitor.monitor_process.return_value = {
            "process_name": "snap-scaffold_1.fasta",
            "duration_seconds": 1.5,
            "peak_rss_mb": 10.0,
            "peak_vms_mb": 20.0,
            "avg_rss_mb": 9.0,
            "avg_vms_mb": 19.0,
            "sample_count": 2,
            "samples": [],
        }

        logger = MagicMock()
        monkeypatch.setenv("FUNANNOTATE2_OUTPUT_DIR", str(tmp_path))

        runSubprocess(
            ["snap", "input.fasta"],
            logger,
            cwd=str(tmp_path),
            monitor_memory=True,
            process_name="snap-scaffold_1.fasta",
        )

        debug_messages = [call.args[0] for call in logger.debug.call_args_list]
        assert "snap input.fasta" in debug_messages
        assert not any("Memory usage for" in message for message in debug_messages)
