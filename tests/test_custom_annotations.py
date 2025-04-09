#!/usr/bin/env python3

import os
import tempfile
import unittest
from unittest import mock
import shutil

from funannotate2.name_cleaner import NameCleaner


class TestCustomAnnotations(unittest.TestCase):
    """Test the custom annotations functionality."""

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
            f.write("gene123\tname\tACT1\n")
            f.write("gene123\tproduct\tCustom Actin\n")
            f.write("gene123\tnote\tManually curated annotation\n")
            f.write("gene456\tname\tCDC42\n")
            f.write("gene456\tproduct\tCustom Cell division control protein 42\n")
            f.write("gene789\tgo_term\tGO:0005524\n")
            f.write("gene789\tec_number\t3.6.4.13\n")

        # Mock the environment variable
        self.patcher = mock.patch.dict(os.environ, {"FUNANNOTATE2_DB": self.test_dir})
        self.patcher.start()

    def tearDown(self):
        """Clean up test environment."""
        # Stop the patcher
        self.patcher.stop()

        # Remove the temporary directory
        shutil.rmtree(self.test_dir)

    def test_load_custom_annotations(self):
        """Test loading custom annotations."""
        cleaner = NameCleaner(custom_file=self.custom_file)

        # Check that custom annotations were loaded
        self.assertIn("gene123", cleaner.custom_annotations)
        self.assertIn("gene456", cleaner.custom_annotations)
        self.assertIn("gene789", cleaner.custom_annotations)

        # Check specific annotations
        self.assertEqual(cleaner.custom_annotations["gene123"]["name"], ["ACT1"])
        self.assertEqual(
            cleaner.custom_annotations["gene123"]["product"], ["Custom Actin"]
        )
        self.assertEqual(
            cleaner.custom_annotations["gene123"]["note"],
            ["Manually curated annotation"],
        )

        self.assertEqual(cleaner.custom_annotations["gene456"]["name"], ["CDC42"])
        self.assertEqual(
            cleaner.custom_annotations["gene456"]["product"],
            ["Custom Cell division control protein 42"],
        )

        self.assertEqual(
            cleaner.custom_annotations["gene789"]["go_term"], ["GO:0005524"]
        )
        self.assertEqual(
            cleaner.custom_annotations["gene789"]["ec_number"], ["3.6.4.13"]
        )

    def test_process_annotation_with_custom_annotations(self):
        """Test processing annotations with custom annotations."""
        cleaner = NameCleaner(custom_file=self.custom_file)

        # Test with a gene that has custom annotations
        annotation = {"name": ["YPT7"], "product": ["GTPase YPT7"]}
        processed = cleaner.process_annotation(annotation, gene_id="gene123")

        # Check that custom annotations were applied
        self.assertEqual(processed["name"], ["ACT1"])
        self.assertEqual(processed["product"], ["Custom Actin"])
        self.assertEqual(processed["note"], ["Manually curated annotation"])

        # Test with a gene that doesn't have custom annotations
        annotation = {"name": ["YPT7"], "product": ["GTPase YPT7"]}
        processed = cleaner.process_annotation(annotation, gene_id="gene999")

        # Check that standard processing was applied
        self.assertEqual(processed["name"], ["YPT7"])
        self.assertEqual(processed["product"], ["GTPase YPT7"])

        # Test with a gene that has GO terms and EC numbers
        annotation = {
            "name": ["YPT7"],
            "product": ["GTPase YPT7"],
            "go_term": ["GO:0003700", "GO:0006355"],
            "ec_number": ["1.2.3.4"],
        }
        processed = cleaner.process_annotation(annotation, gene_id="gene789")

        # Check that custom annotations were added to existing ones
        self.assertEqual(processed["name"], ["YPT7"])
        self.assertEqual(processed["product"], ["GTPase YPT7"])
        self.assertEqual(
            set(processed["go_term"]), {"GO:0003700", "GO:0006355", "GO:0005524"}
        )
        self.assertEqual(set(processed["ec_number"]), {"1.2.3.4", "3.6.4.13"})


if __name__ == "__main__":
    unittest.main()
