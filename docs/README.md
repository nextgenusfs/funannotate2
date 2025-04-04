# Funannotate2 Documentation

This directory contains the documentation for Funannotate2.

## Building the Documentation

To build the documentation, you need to have Sphinx installed:

```bash
pip install sphinx sphinx_rtd_theme
```

Then, you can build the documentation:

```bash
cd docs
make html
```

The built documentation will be in the `_build/html` directory.

## Documentation Structure

- `index.rst`: Main index file
- `installation.rst`: Installation instructions
- `usage.rst`: Usage guide
- `modules.rst`: Module reference
- `API.rst`: API reference
- `tutorial.rst`: Tutorial
- `faq.rst`: Frequently asked questions
- `changelog.rst`: Changelog
- `API/`: Directory containing API documentation for each module
  - `clean.rst`: Documentation for the clean module
  - `predict.rst`: Documentation for the predict module
  - `annotate.rst`: Documentation for the annotate module
  - `compare.rst`: Documentation for the compare module
  - `search.rst`: Documentation for the search module
  - `utilities.rst`: Documentation for the utilities module
  - `log.rst`: Documentation for the log module
