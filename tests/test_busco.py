#!/usr/bin/env python3

import os
import unittest
import tempfile
from buscolite.busco import load_config, load_cutoffs, check_lineage


class TestBuscoFunctions(unittest.TestCase):
    """Test basic functions in the busco module."""

    def setUp(self):
        """Set up test fixtures."""
        # Create a temporary directory to simulate a BUSCO lineage
        self.temp_dir = tempfile.TemporaryDirectory()
        self.lineage_dir = self.temp_dir.name

        # Create required directories
        os.makedirs(os.path.join(self.lineage_dir, "hmms"), exist_ok=True)
        os.makedirs(os.path.join(self.lineage_dir, "prfl"), exist_ok=True)

        # Create required files
        with open(os.path.join(self.lineage_dir, "dataset.cfg"), "w") as f:
            f.write("name=test_lineage\n")
            f.write("species=test_species\n")
            f.write("domain=test_domain\n")

        with open(os.path.join(self.lineage_dir, "scores_cutoff"), "w") as f:
            f.write("busco1\t100.0\n")
            f.write("busco2\t200.0\n")

        with open(os.path.join(self.lineage_dir, "lengths_cutoff"), "w") as f:
            f.write("busco1\t0\t1.5\t300\n")
            f.write("busco2\t0\t0.0\t400\n")

        # Create empty files for the remaining required files
        open(os.path.join(self.lineage_dir, "ancestral"), "w").close()
        open(os.path.join(self.lineage_dir, "ancestral_variants"), "w").close()

    def tearDown(self):
        """Tear down test fixtures."""
        self.temp_dir.cleanup()

    def test_load_config(self):
        """Test loading configuration from dataset.cfg."""
        config = load_config(self.lineage_dir)
        self.assertIsInstance(config, dict)
        self.assertEqual(config["name"], "test_lineage")
        self.assertEqual(config["species"], "test_species")
        self.assertEqual(config["domain"], "test_domain")

    def test_load_cutoffs(self):
        """Test loading cutoffs from scores_cutoff and lengths_cutoff."""
        cutoffs = load_cutoffs(self.lineage_dir)
        self.assertIsInstance(cutoffs, dict)
        self.assertIn("busco1", cutoffs)
        self.assertIn("busco2", cutoffs)

        # Check busco1 values
        self.assertEqual(cutoffs["busco1"]["score"], 100.0)
        self.assertEqual(cutoffs["busco1"]["sigma"], 1.5)
        self.assertEqual(cutoffs["busco1"]["length"], 300)

        # Check busco2 values
        self.assertEqual(cutoffs["busco2"]["score"], 200.0)
        self.assertEqual(
            cutoffs["busco2"]["sigma"], 1
        )  # Should be 1 because sigma was 0.0
        self.assertEqual(cutoffs["busco2"]["length"], 400)

    def test_check_lineage_valid(self):
        """Test checking a valid lineage directory."""
        valid, message = check_lineage(self.lineage_dir)
        self.assertTrue(valid)
        self.assertEqual(message, "")

    def test_check_lineage_invalid_dir(self):
        """Test checking an invalid directory."""
        invalid_dir = os.path.join(self.lineage_dir, "nonexistent")
        valid, message = check_lineage(invalid_dir)
        self.assertFalse(valid)
        self.assertIn("is not a directory", message)

    def test_check_lineage_missing_dir(self):
        """Test checking a lineage with a missing directory."""
        # Remove the hmms directory
        os.rmdir(os.path.join(self.lineage_dir, "hmms"))

        valid, message = check_lineage(self.lineage_dir)
        self.assertFalse(valid)
        self.assertIn("hmms directory was not found", message)

    def test_check_lineage_missing_file(self):
        """Test checking a lineage with a missing file."""
        # Remove the scores_cutoff file
        os.remove(os.path.join(self.lineage_dir, "scores_cutoff"))

        valid, message = check_lineage(self.lineage_dir)
        self.assertFalse(valid)
        self.assertIn("scores_cutoff file is missing", message)


if __name__ == "__main__":
    unittest.main()
