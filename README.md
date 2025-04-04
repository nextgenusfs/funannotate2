[![Latest Github release](https://img.shields.io/github/release/nextgenusfs/funannotate2.svg)](https://github.com/nextgenusfs/funannotate2/releases/latest)
![Conda](https://img.shields.io/conda/dn/bioconda/funannotate2)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Tests](https://github.com/nextgenusfs/funannotate2/actions/workflows/tests.yml/badge.svg)](https://github.com/nextgenusfs/funannotate2/actions/workflows/tests.yml)

# funannotate2: eukaryotic genome annotation pipeline

Funannotate2 is a comprehensive eukaryotic genome annotation pipeline that provides a complete workflow for annotating fungal, plant, and other eukaryotic genomes. It integrates various tools and databases to produce high-quality gene predictions and functional annotations.

#### This is work in progress. Do not expect it to work until a release is tagged.

Setting up a conda environment with these packages is necessary, then it can be installed with pip.

```shell
mamba create -n funannotate2 "python<=3.10" biopython "evidencemodeler>=2" minimap2 miniprot snap "augustus==3.5.0" glimmerhmm diamond trnascan-se table2asn
conda activate funannotate2
python -m pip install git+https://github.com/nextgenusfs/funannotate2.git
```

Additional tools like genemarkHMM must be installed manually due to licensing.

To install release versions use the pip package manager, like so:

```shell
pip install funannotate2
```

## Development

### Testing

Funannotate2 includes both unit tests and integration tests to ensure the code works correctly.

#### Running Tests

To run the tests, you need to install pytest and the package in development mode:

```bash
# Install pytest and coverage tools
pip install pytest pytest-cov

# Install funannotate2 in development mode
pip install -e .

# Run all tests
pytest

# Run with coverage report
pytest --cov=funannotate2

# Generate HTML coverage report
python scripts/run_coverage.py
```

For more information about testing, see the [TESTING.md](TESTING.md) file.

### Development Dependencies

To work on funannotate2 development, you'll need to install the development dependencies:

```shell
pip install pytest pytest-cov
```

### Documentation

Funannotate2 includes comprehensive documentation that covers installation, usage, API reference, and more. To build the documentation:

```bash
# Install Sphinx and the theme
pip install sphinx sphinx_rtd_theme

# Build the documentation
cd docs
make html
```

The built documentation will be in the `docs/_build/html` directory.

For more information about the documentation, see the [docs/README.md](docs/README.md) file.

### Running Tests

After installing the development dependencies, you can run the tests with:

```shell
python -m pytest
```

To run tests with coverage reporting:

```shell
python -m pytest --cov=funannotate2 --cov-report=term-missing
```

Or use the provided script to generate an HTML coverage report:

```shell
python scripts/run_coverage.py
```

To install the most up to date code from this repo, you can run:
```
python -m pip install git+https://github.com/nextgenusfs/funannotate2.git --upgrade --force --no-deps
```
