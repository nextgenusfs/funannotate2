#!/usr/bin/env python3

import os
import shutil
import tempfile
import unittest
from unittest import mock

from funannotate2.name_cleaner import (
    NameCleaner,
    write_new_valid_annotations,
    write_problematic_annotations,
)


class TestAnnotateIntegration(unittest.TestCase):
    """Test the integration of name cleaning in the annotate module."""

    def setUp(self):
        """Set up test environment."""
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()

        # Create a mock curated gene products file
        self.curated_file = os.path.join(self.test_dir, "ncbi_cleaned_gene_products.txt")
        with open(self.curated_file, "w") as f:
            f.write("# Curated gene names and products\n")
            f.write("ACT1\tActin\n")
            f.write("CDC42\tCell division control protein 42\n")

        # Mock the environment variable
        self.patcher = mock.patch.dict(os.environ, {"FUNANNOTATE2_DB": self.test_dir})
        self.patcher.start()

        # Create sample merged annotations
        self.merged = {
            "gene1": {"name": ["ACT1", "actin"], "product": ["Actin-like protein"]},
            "gene2": {"name": ["YPT7"], "product": ["GTPase YPT7"]},
            "gene3": {"name": ["orf19.123"], "product": ["Hypothetical protein"]},
            "gene4": {
                "name": ["CDC42"],
                "product": ["Required for cell division and establishment of cell polarity"],
            },
        }

    def tearDown(self):
        """Clean up test environment."""
        # Stop the patcher
        self.patcher.stop()

        # Remove the temporary directory
        shutil.rmtree(self.test_dir)

    def test_write_problematic_annotations(self):
        """Test writing problematic annotations to file."""
        output_file = os.path.join(self.test_dir, "problematic.txt")
        count = write_problematic_annotations(self.merged, output_file)

        # Check that the correct number of problematic annotations were found
        self.assertEqual(count, 1)

        # Check that the output file was created
        self.assertTrue(os.path.exists(output_file))

        # Check the content of the output file
        with open(output_file, "r") as f:
            lines = f.readlines()

        self.assertEqual(len(lines), 2)  # Header + 1 problematic annotation
        self.assertIn("gene4", lines[1])
        self.assertIn("CDC42", lines[1])
        self.assertIn("Required for cell division", lines[1])

    def test_write_new_valid_annotations(self):
        """Test writing new valid annotations to file."""
        output_file = os.path.join(self.test_dir, "new_valid.txt")
        count = write_new_valid_annotations(self.merged, output_file)

        # Check that the correct number of new valid annotations were found
        self.assertEqual(count, 1)

        # Check that the output file was created
        self.assertTrue(os.path.exists(output_file))

        # Check the content of the output file
        with open(output_file, "r") as f:
            lines = f.readlines()

        self.assertEqual(len(lines), 2)  # Header + 1 new valid annotation
        self.assertIn("YPT7", lines[1])
        self.assertIn("GTPase YPT7", lines[1])

    def test_name_cleaner_integration(self):
        """Test the integration of NameCleaner with the annotate module."""
        # Create a NameCleaner instance
        cleaner = NameCleaner()

        # Clean the merged annotations
        cleaned_merged = {}
        for gene_id, annot in self.merged.items():
            cleaned_annot = cleaner.process_annotation(annot)
            cleaned_merged[gene_id] = cleaned_annot

        # Check that curated names and products are preserved
        self.assertEqual(cleaned_merged["gene1"]["name"], ["ACT1"])
        self.assertEqual(cleaned_merged["gene1"]["product"], ["Actin"])

        # Check that valid non-curated names and products are preserved
        self.assertEqual(cleaned_merged["gene2"]["name"], ["YPT7"])
        self.assertEqual(cleaned_merged["gene2"]["product"], ["GTPase YPT7"])

        # Check that invalid names are removed
        self.assertNotIn("name", cleaned_merged["gene3"])

        # Check that problematic products are replaced with curated ones
        self.assertEqual(cleaned_merged["gene4"]["name"], ["CDC42"])
        self.assertEqual(cleaned_merged["gene4"]["product"], ["Cell division control protein 42"])


if __name__ == "__main__":
    unittest.main()
