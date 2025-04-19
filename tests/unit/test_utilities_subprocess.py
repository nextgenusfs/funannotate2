"""
Unit tests for the runSubprocess function in utilities.py.
"""

import logging

import funannotate2.utilities


class TestRunSubprocess:
    """Tests for the runSubprocess function."""

    def test_run_subprocess(self):
        """Test running a subprocess with a mock implementation."""
        # Set up a logger for testing
        logger = logging.getLogger("test_logger")
        logger.setLevel(logging.DEBUG)

        # Create a custom implementation of runSubprocess that returns 0
        def mock_run_subprocess(cmd, logfile, **kwargs):
            return 0

        # Save the original function
        original_run_subprocess = funannotate2.utilities.runSubprocess

        try:
            # Replace with our mock function
            funannotate2.utilities.runSubprocess = mock_run_subprocess

            # Call the function
            from funannotate2.utilities import runSubprocess

            result = runSubprocess(["echo", "test"], logger)

            # Check the return value
            assert result == 0
        finally:
            # Restore the original function
            funannotate2.utilities.runSubprocess = original_run_subprocess

    def test_failed_command(self):
        """Test running a command that fails."""
        # Set up a logger for testing
        logger = logging.getLogger("test_logger")
        logger.setLevel(logging.DEBUG)

        # Create a custom implementation of runSubprocess that returns 1
        def mock_run_subprocess(cmd, logfile, **kwargs):
            return 1

        # Save the original function
        original_run_subprocess = funannotate2.utilities.runSubprocess

        try:
            # Replace with our mock function
            funannotate2.utilities.runSubprocess = mock_run_subprocess

            # Call the function
            from funannotate2.utilities import runSubprocess

            result = runSubprocess(["false"], logger)

            # Check the return value
            assert result == 1
        finally:
            # Restore the original function
            funannotate2.utilities.runSubprocess = original_run_subprocess

    def test_with_custom_env(self):
        """Test with a custom environment."""
        # Set up a logger for testing
        logger = logging.getLogger("test_logger")
        logger.setLevel(logging.DEBUG)

        # Create a custom implementation of runSubprocess that returns 0
        def mock_run_subprocess(cmd, logfile, **kwargs):
            # Check that the environment was passed correctly
            assert "env" in kwargs
            assert kwargs["env"]["TEST_VAR"] == "test_value"
            return 0

        # Save the original function
        original_run_subprocess = funannotate2.utilities.runSubprocess

        try:
            # Replace with our mock function
            funannotate2.utilities.runSubprocess = mock_run_subprocess

            # Call the function
            from funannotate2.utilities import runSubprocess

            custom_env = {"TEST_VAR": "test_value"}
            result = runSubprocess(["echo", "test"], logger, env=custom_env)

            # Check the return value
            assert result == 0
        finally:
            # Restore the original function
            funannotate2.utilities.runSubprocess = original_run_subprocess
