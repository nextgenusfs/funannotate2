# BUSCOlite Integration Tests

This directory contains integration tests for the BUSCOlite package. These tests verify that the package works correctly as a whole, testing the entire workflow from start to finish.

## Test Structure

The integration tests are organized as follows:

- `test_buscolite_cli.py`: Tests for the command-line interface
- `test_buscolite_workflow.py`: Tests for the workflow using the Python API
- `test_buscolite_api.py`: Tests for the Python API functions
- `test_buscolite_end_to_end.py`: End-to-end tests using subprocess

## Test Data

The test data is stored in the `data` directory:

- `test_genome.fasta`: A small test genome
- `test_proteins.fasta`: A small test proteome
- `mock_lineage`: A mock BUSCO lineage directory

## Running the Tests

To run the integration tests, use the following command:

```bash
# Run all integration tests
pytest tests/integration

# Run a specific test file
pytest tests/integration/test_buscolite_cli.py

# Run a specific test
pytest tests/integration/test_buscolite_cli.py::TestBUSCOliteCLI::test_help_command
```

## Requirements

The integration tests require the following dependencies:

- pytest
- buscolite (installed in development mode)
- augustus (for genome mode tests)
- miniprot (for genome mode tests)

## Notes

- Some tests are skipped if the required dependencies are not installed
- The tests create temporary directories for output files, which are automatically cleaned up after the tests finish
