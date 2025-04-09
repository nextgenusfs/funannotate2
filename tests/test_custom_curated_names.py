#!/usr/bin/env python3

import os
import tempfile
import unittest
from unittest import mock
import shutil

from funannotate2.name_cleaner import NameCleaner


class TestCustomCuratedNames(unittest.TestCase):
    """Test the custom curated names functionality."""

    def setUp(self):
        """Set up test environment."""
        # Create a temporary directory
        self.test_dir = tempfile.mkdtemp()

        # Create a mock curated gene products file
        self.curated_file = os.path.join(
            self.test_dir, "ncbi_cleaned_gene_products.txt"
        )
        with open(self.curated_file, "w") as f:
            f.write("# Curated gene names and products\n")
            f.write("ACT1\tActin\n")
            f.write("CDC42\tCell division control protein 42\n")

        # Create a mock custom annotations file
        self.custom_file = os.path.join(self.test_dir, "custom_annotations.txt")
        with open(self.custom_file, "w") as f:
            f.write("# Custom annotations for specific genes/transcripts\n")
            f.write("gene_ypt7\tname\tYPT7\n")
            f.write("gene_ypt7\tproduct\tCustom GTPase YPT7\n")
            f.write("gene_act1\tname\tACT1\n")
            f.write(
                "gene_act1\tproduct\tCustom Actin\n"
            )  # This should override the standard database

        # Mock the environment variable
        self.patcher = mock.patch.dict(os.environ, {"FUNANNOTATE2_DB": self.test_dir})
        self.patcher.start()

    def tearDown(self):
        """Clean up test environment."""
        # Stop the patcher
        self.patcher.stop()

        # Remove the temporary directory
        shutil.rmtree(self.test_dir)

    def test_standard_database(self):
        """Test loading the standard database."""
        cleaner = NameCleaner()

        # Check that standard database was loaded
        self.assertIn("ACT1", cleaner.curated_names)
        self.assertEqual(cleaner.curated_names["ACT1"], "Actin")
        self.assertIn("CDC42", cleaner.curated_names)
        self.assertEqual(
            cleaner.curated_names["CDC42"], "Cell division control protein 42"
        )

        # Check that custom names are not present
        self.assertNotIn("YPT7", cleaner.curated_names)

    def test_custom_database(self):
        """Test loading a custom database."""
        cleaner = NameCleaner(custom_file=self.custom_file)

        # Check that standard database was loaded
        self.assertIn("ACT1", cleaner.curated_names)
        self.assertIn("CDC42", cleaner.curated_names)

        # Check that custom annotations were loaded
        self.assertIn("gene_ypt7", cleaner.custom_annotations)
        self.assertIn("gene_act1", cleaner.custom_annotations)

        # Check specific custom annotations
        self.assertEqual(cleaner.custom_annotations["gene_ypt7"]["name"], ["YPT7"])
        self.assertEqual(
            cleaner.custom_annotations["gene_ypt7"]["product"], ["Custom GTPase YPT7"]
        )
        self.assertEqual(cleaner.custom_annotations["gene_act1"]["name"], ["ACT1"])
        self.assertEqual(
            cleaner.custom_annotations["gene_act1"]["product"], ["Custom Actin"]
        )

    def test_process_annotation_with_custom_database(self):
        """Test processing annotations with a custom database."""
        cleaner = NameCleaner(custom_file=self.custom_file)

        # Test with a gene that has custom annotations
        annotation = {"name": ["Original Name"], "product": ["Original Product"]}
        processed = cleaner.process_annotation(annotation, gene_id="gene_ypt7")

        # Check that the custom annotations were applied
        self.assertEqual(processed["name"], ["YPT7"])
        self.assertEqual(processed["product"], ["Custom GTPase YPT7"])

        # Test with another gene that has custom annotations
        annotation = {"name": ["Original Name"], "product": ["Original Product"]}
        processed = cleaner.process_annotation(annotation, gene_id="gene_act1")

        # Check that the custom annotations were applied
        self.assertEqual(processed["name"], ["ACT1"])
        self.assertEqual(processed["product"], ["Custom Actin"])

        # Test with a gene that doesn't have custom annotations
        annotation = {"name": ["ACT1", "actin"], "product": ["Actin-like protein"]}
        processed = cleaner.process_annotation(annotation, gene_id="gene_unknown")

        # Check that standard processing was applied
        self.assertEqual(processed["name"], ["ACT1"])
        # The product should be cleaned but not replaced with a custom value
        self.assertNotEqual(processed["product"], ["Custom Actin"])


if __name__ == "__main__":
    unittest.main()
