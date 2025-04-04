#!/usr/bin/env python3
"""
Script to run tests with coverage and generate HTML reports.
"""
import os
import subprocess
import sys
import webbrowser
from pathlib import Path


def main():
    """Run tests with coverage and generate HTML reports."""
    # Create the htmlcov directory if it doesn't exist
    htmlcov_dir = Path("htmlcov")
    htmlcov_dir.mkdir(exist_ok=True)

    # Run pytest with coverage
    cmd = [
        "python",
        "-m",
        "pytest",
        "--cov=funannotate2",
        "--cov-report=term-missing",
        "--cov-report=html:htmlcov",
    ]

    # Add any additional arguments passed to this script
    cmd.extend(sys.argv[1:])

    # Run the command
    result = subprocess.run(cmd, check=False)

    # Open the HTML report in the default browser
    if result.returncode == 0:
        index_path = htmlcov_dir / "index.html"
        if index_path.exists():
            print(f"\nOpening coverage report in browser: {index_path}")
            webbrowser.open(f"file://{index_path.absolute()}")
        else:
            print(f"Error: Could not find coverage report at {index_path}")

    # Return the pytest exit code
    return result.returncode


if __name__ == "__main__":
    sys.exit(main())
