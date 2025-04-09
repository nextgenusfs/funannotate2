#!/usr/bin/env python3

import os
import tempfile
import unittest
from unittest import mock
import shutil

from funannotate2.name_cleaner import (
    NameCleaner,
    clean_annotations,
    write_problematic_annotations,
    write_new_valid_annotations,
    number_present,
    morethanXnumbers,
    capfirst,
)


class TestHelperFunctions(unittest.TestCase):
    """Test helper functions in the name_cleaner module."""

    def test_number_present(self):
        """Test the number_present function."""
        self.assertTrue(number_present("abc123"))
        self.assertTrue(number_present("123"))
        self.assertFalse(number_present("abc"))
        self.assertFalse(number_present(""))

    def test_morethanXnumbers(self):
        """Test the morethanXnumbers function."""
        self.assertTrue(morethanXnumbers("abc123def456", 2))
        self.assertTrue(morethanXnumbers("123456", 3))
        self.assertFalse(morethanXnumbers("abc123", 2))
        self.assertFalse(morethanXnumbers("", 1))

    def test_capfirst(self):
        """Test the capfirst function."""
        self.assertEqual(capfirst("abc"), "Abc")
        self.assertEqual(capfirst("ABC"), "ABC")
        self.assertEqual(capfirst("a"), "A")
        self.assertEqual(capfirst(""), "")


class TestNameCleaner(unittest.TestCase):
    """Test the NameCleaner class."""

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
            f.write("HSP90\tHeat shock protein 90\n")
            f.write("TUB1\tTubulin alpha-1 chain\n")

        # Mock the environment variable
        self.patcher = mock.patch.dict(os.environ, {"FUNANNOTATE2_DB": self.test_dir})
        self.patcher.start()

        # Create a NameCleaner instance
        self.cleaner = NameCleaner()

    def tearDown(self):
        """Clean up test environment."""
        # Stop the patcher
        self.patcher.stop()

        # Remove the temporary directory
        shutil.rmtree(self.test_dir)

    def test_load_curated_names(self):
        """Test loading curated names from file."""
        # Check that curated names were loaded correctly
        self.assertIn("ACT1", self.cleaner.curated_names)
        self.assertIn("act1", self.cleaner.curated_names)
        self.assertEqual(self.cleaner.curated_names["ACT1"], "Actin")
        self.assertEqual(self.cleaner.curated_names["act1"], "Actin")

    def test_clean_name(self):
        """Test cleaning gene names."""
        # Valid names
        self.assertEqual(self.cleaner.clean_name("ACT1"), "ACT1")
        self.assertEqual(self.cleaner.clean_name("CDC42"), "CDC42")
        self.assertEqual(self.cleaner.clean_name("YPT7"), "YPT7")

        # Invalid names
        self.assertIsNone(self.cleaner.clean_name("orf19.123"))
        self.assertIsNone(self.cleaner.clean_name("123abc"))
        self.assertIsNone(self.cleaner.clean_name("A1B2C3D4"))
        self.assertIsNone(self.cleaner.clean_name(""))
        self.assertIsNone(self.cleaner.clean_name(None))

    def test_clean_product(self):
        """Test cleaning product descriptions."""
        # Basic cleaning
        self.assertEqual(
            self.cleaner.clean_product("potential actin binding protein"),
            "putative actin binding protein",
        )
        self.assertEqual(
            self.cleaner.clean_product("conserved hypothetical protein"),
            "hypothetical protein",
        )

        # Remove problematic phrases
        self.assertEqual(
            self.cleaner.clean_product("protein similar to CDC42"), "protein CDC42"
        )
        self.assertEqual(
            self.cleaner.clean_product("domain-containing protein"), "domain protein"
        )

        # Default for empty
        self.assertEqual(self.cleaner.clean_product(""), "hypothetical protein")

        # With gene name
        self.assertEqual(self.cleaner.clean_product("ACT1 protein", "ACT1"), "act1p")

        # Problematic descriptions
        self.assertEqual(
            self.cleaner.clean_product(
                "Required for actin cytoskeleton organization", "ACT1"
            ),
            "Act1p",
        )
        self.assertEqual(
            self.cleaner.clean_product("Involved in cell division", "CDC42"), "Cdc42p"
        )

    def test_get_curated_product(self):
        """Test getting curated product descriptions."""
        self.assertEqual(self.cleaner.get_curated_product("ACT1"), "Actin")
        self.assertEqual(self.cleaner.get_curated_product("act1"), "Actin")
        self.assertEqual(
            self.cleaner.get_curated_product("CDC42"),
            "Cell division control protein 42",
        )
        self.assertIsNone(self.cleaner.get_curated_product("YPT7"))

    def test_process_annotation(self):
        """Test processing annotations."""
        # Test with curated name
        annotation = {"name": ["ACT1", "actin"], "product": ["Actin-like protein"]}
        processed = self.cleaner.process_annotation(annotation)
        self.assertEqual(processed["name"], ["ACT1"])
        self.assertEqual(processed["product"], ["Actin"])

        # Test with valid but non-curated name
        annotation = {"name": ["YPT7"], "product": ["GTPase YPT7"]}
        processed = self.cleaner.process_annotation(annotation)
        self.assertEqual(processed["name"], ["YPT7"])
        self.assertEqual(processed["product"], ["GTPase YPT7"])

        # Test with invalid name
        annotation = {"name": ["orf19.123"], "product": ["Hypothetical protein"]}
        processed = self.cleaner.process_annotation(annotation)
        self.assertNotIn("name", processed)
        self.assertEqual(processed["product"], ["Hypothetical protein"])

        # Test with no name
        annotation = {"product": ["Hypothetical protein"]}
        processed = self.cleaner.process_annotation(annotation)
        self.assertEqual(processed["product"], ["Hypothetical protein"])


class TestAnnotationFunctions(unittest.TestCase):
    """Test functions for working with collections of annotations."""

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

        # Mock the environment variable
        self.patcher = mock.patch.dict(os.environ, {"FUNANNOTATE2_DB": self.test_dir})
        self.patcher.start()

        # Create sample annotations
        self.annotations = {
            "gene1": {"name": ["ACT1", "actin"], "product": ["Actin-like protein"]},
            "gene2": {"name": ["YPT7"], "product": ["GTPase YPT7"]},
            "gene3": {"name": ["orf19.123"], "product": ["Hypothetical protein"]},
            "gene4": {
                "name": ["CDC42"],
                "product": [
                    "Required for cell division and establishment of cell polarity"
                ],
            },
            "gene5": {"product": ["Hypothetical protein"]},
        }

    def tearDown(self):
        """Clean up test environment."""
        # Stop the patcher
        self.patcher.stop()

        # Remove the temporary directory
        shutil.rmtree(self.test_dir)

    def test_clean_annotations(self):
        """Test cleaning a collection of annotations."""
        cleaned = clean_annotations(self.annotations)

        # Check that curated names and products are preserved
        self.assertEqual(cleaned["gene1"]["name"], ["ACT1"])
        self.assertEqual(cleaned["gene1"]["product"], ["Actin"])

        # Check that valid non-curated names and products are preserved
        self.assertEqual(cleaned["gene2"]["name"], ["YPT7"])
        self.assertEqual(cleaned["gene2"]["product"], ["GTPase YPT7"])

        # Check that invalid names are removed
        self.assertNotIn("name", cleaned["gene3"])

        # Check that problematic products are replaced
        self.assertEqual(cleaned["gene4"]["name"], ["CDC42"])
        self.assertEqual(
            cleaned["gene4"]["product"], ["Cell division control protein 42"]
        )

    def test_write_problematic_annotations(self):
        """Test writing problematic annotations to file."""
        output_file = os.path.join(self.test_dir, "problematic.txt")
        count = write_problematic_annotations(self.annotations, output_file)

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
        count = write_new_valid_annotations(self.annotations, output_file)

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


if __name__ == "__main__":
    unittest.main()
