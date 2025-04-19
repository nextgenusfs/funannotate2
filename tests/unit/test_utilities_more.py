"""
Additional unit tests for the utilities module.
"""

import os
import tempfile

import pytest

from funannotate2.utilities import checkfile, process_handle


class TestCheckfile:
    """Tests for the checkfile function."""

    def test_file_exists_and_not_empty(self):
        """Test with a file that exists and is not empty."""
        with tempfile.NamedTemporaryFile(delete=False) as temp:
            temp.write(b"test content")
            temp_name = temp.name

        try:
            assert checkfile(temp_name) is True
        finally:
            os.unlink(temp_name)

    def test_file_exists_but_empty(self):
        """Test with a file that exists but is empty."""
        with tempfile.NamedTemporaryFile(delete=False) as temp:
            temp_name = temp.name

        try:
            assert checkfile(temp_name) is False
        finally:
            os.unlink(temp_name)

    def test_file_does_not_exist(self):
        """Test with a file that does not exist."""
        assert checkfile("nonexistent_file.txt") is False

    def test_symlink(self):
        """Test with a symbolic link."""
        with tempfile.NamedTemporaryFile(delete=False) as temp:
            temp.write(b"test content")
            temp_name = temp.name

        symlink_name = temp_name + ".link"
        try:
            os.symlink(temp_name, symlink_name)
            assert checkfile(symlink_name) is True
        finally:
            os.unlink(temp_name)
            if os.path.exists(symlink_name):
                os.unlink(symlink_name)


class TestProcessHandle:
    """Tests for the process_handle context manager."""

    def test_pipe(self):
        """Test with handle=True (pipe)."""
        import subprocess

        with process_handle(True) as handle:
            assert handle == subprocess.PIPE

    def test_devnull(self):
        """Test with handle=False (devnull)."""
        import subprocess

        with process_handle(False) as handle:
            assert handle == subprocess.DEVNULL

    def test_stdout(self):
        """Test with handle='STDOUT'."""
        import subprocess

        with process_handle("STDOUT") as handle:
            assert handle == subprocess.STDOUT

    def test_none(self):
        """Test with handle=None."""
        with process_handle(None) as handle:
            assert handle is None

    def test_file_write(self):
        """Test with handle as a file path (write mode)."""
        with tempfile.NamedTemporaryFile(delete=False) as temp:
            temp_name = temp.name

        try:
            with process_handle(temp_name, mode="w") as handle:
                assert hasattr(handle, "write")
                handle.write("test content")

            with open(temp_name, "r") as f:
                content = f.read()
            assert content == "test content"
        finally:
            os.unlink(temp_name)

    def test_file_append(self):
        """Test with handle as a file path (append mode)."""
        with tempfile.NamedTemporaryFile(delete=False) as temp:
            temp.write(b"initial content\n")
            temp_name = temp.name

        try:
            with process_handle(temp_name, mode="a") as handle:
                assert hasattr(handle, "write")
                handle.write("appended content")

            with open(temp_name, "r") as f:
                content = f.read()
            assert content == "initial content\nappended content"
        finally:
            os.unlink(temp_name)

    def test_invalid_mode(self):
        """Test with an invalid mode."""
        with pytest.raises(ValueError):
            with process_handle("test.txt", mode="invalid"):
                pass
