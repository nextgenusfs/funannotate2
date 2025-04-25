[![Latest Github release](https://img.shields.io/github/release/nextgenusfs/funannotate2.svg)](https://github.com/nextgenusfs/funannotate2/releases/latest)
![Conda](https://img.shields.io/conda/dn/bioconda/funannotate2)
[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)
[![Tests](https://github.com/nextgenusfs/funannotate2/actions/workflows/tests.yml/badge.svg)](https://github.com/nextgenusfs/funannotate2/actions/workflows/tests.yml)

# funannotate2: eukaryotic genome annotation pipeline

Funannotate2 is a comprehensive eukaryotic genome annotation pipeline that provides a complete workflow for annotating
eukaryotic genomes. It integrates various tools and databases to produce high-quality gene predictions and functional annotations.


### Quick start: Installation

#### Linux systems

Unit this gets pushed to bioconda, can try this:
```shell
mamba create -n funannotate2 gfftk gapmm2 minimap2 miniprot snap "augustus==3.5.0" glimmerhmm diamond trnascan-se table2asn gb-io buscolite
conda activate funannotate2
python -m pip install git+https://github.com/nextgenusfs/funannotate2.git
```

#### Apple Silicon (M series)
Installation on apple silicon (M series) is a little bit more involved due to some dependency issues and non-native builds of some software.  I've not been able to find or build a version of `augustus` that will run, so instead I've been running `augustus` and `genemark` locally with Docker.  I've setup two repos with instructions on how to get this working (Need Docker Desktop installed) and then will need to put the bash wrapper files in your PATH to mimic the CLI interface.

https://github.com/nextgenusfs/dockerized-augustus

https://github.com/nextgenusfs/dockerized-genemark

Once that is working, you can then install most of the remaining dependencies with conda, although we need to leave out both `buscolite` and `funannotate2` because they have `augustus` as a dependency, instead we will install those python packages with pip.  The conda mkl<2022 is to avoid an annoying warning on apple silicon with the intel mkl package.

```shell
# first install most of the dependencies
mamba create -n funannotate2 --platform osx-64 "python>=3.7,<3.13" gfftk gapmm2 minimap2 miniprot snap glimmerhmm diamond trnascan-se gb-io pyhmmer pyfastx requests json-repair pytantan "mkl<2022"

# we can then add the required FUNANNOTATE2_DB env variable to the conda environment, note need to reactivate to use it
conda activate funannotate2
conda env config vars set FUNANNOTATE2_DB=/path/to/funannotate2-db
conda env config vars set AUGUSTUS_CONFIG_PATH=/path/to/augustus-3.5.0/config
conda deactivate

# now reactivate environment, and install the remaining python dependencies with pip
conda activate funannotate2
python -m pip install buscolite git+https://github.com/nextgenusfs/funannotate2.git

# now we can install the databases
funannotate2 install -d all
```

#### Other/Manual Installation

Additional tools like genemarkHMM must be installed manually due to licensing.

`funannotate2` is a python package, to install release versions use the pip package manager, like so:

```shell
pip install funannotate2
```
Or to install the bleeding edge version from github repo:

```shell
python -m pip install git+https://github.com/nextgenusfs/funannotate2.git
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

### Citation

Funannotate2 includes a CITATION.cff file that provides citation information for the software. The version and release date in this file are automatically updated when a new release is created.

To cite funannotate2 in your work, you can use the citation information from the CITATION.cff file or generate a citation in your preferred format using tools like [citeas.org](https://citeas.org/).

