# Git Hooks for CITATION.cff Management

This directory contains Git hooks to help manage the CITATION.cff file in the repository.

## Pre-commit Hook

The pre-commit hook automatically updates the CITATION.cff file when changes to version information in pyproject.toml are committed.

### Installation

To install the pre-commit hook, run:

```bash
git config core.hooksPath .githooks
chmod +x .githooks/pre-commit
```

This will configure Git to use the hooks in this directory and make the pre-commit hook executable.

## Manual Update

If you need to manually update the CITATION.cff file, you can use the GitHub Actions workflow by creating a new release/tag, or modify the file directly.
