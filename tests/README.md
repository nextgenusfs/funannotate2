# Testing funannotate2

This directory contains tests for the funannotate2 package.

## Test Structure

- `unit/`: Unit tests for individual functions and classes
- `functional/`: Functional tests for larger components and workflows
- `data/`: Test data files used by the tests

## Running Tests

To run the tests, use pytest:

```bash
# Run all tests
pytest

# Run specific test file
pytest tests/unit/test_utilities.py

# Run specific test
pytest tests/unit/test_utilities.py::TestMergeCoordinates::test_overlapping_intervals

# Run with verbose output
pytest -v

# Run with coverage report
pytest --cov=funannotate2
```

## Writing Tests

When writing tests:

1. Place unit tests in the `unit/` directory
2. Name test files with the prefix `test_`
3. Name test classes with the prefix `Test`
4. Name test methods with the prefix `test_`
5. Use fixtures from `conftest.py` when possible
6. Add test data to the `data/` directory when needed

## Test Dependencies

The tests require the following packages:

- pytest
- pytest-cov (optional, for coverage reports)

You can install these with:

```bash
pip install pytest pytest-cov
```
