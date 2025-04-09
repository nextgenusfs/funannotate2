#!/usr/bin/env python3

import os
import tempfile
import unittest
from unittest import mock
import shutil
import json

from funannotate2.name_cleaner import NameCleaner


class TestAnnotateNameCleaning(unittest.TestCase):
    """Test the integration of name cleaning in the annotate module."""

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

        # Create sample merged annotations
        self.merged = {
            "gene1": {
                "name": ["ACT1", "actin"],
                "product": ["Actin-like protein"],
                "note": ["Cytoskeletal protein"],
                "db_xref": ["UniProtKB:P60010"],
                "ec_number": ["3.6.4.-"],
                "go_terms": ["GO:0005524", "GO:0005856"],
            },
            "gene2": {
                "name": ["YPT7"],
                "product": ["GTPase YPT7"],
                "note": ["Involved in vesicle transport"],
                "db_xref": ["UniProtKB:P32939"],
                "ec_number": ["3.6.5.-"],
                "go_terms": ["GO:0005525", "GO:0006886"],
            },
            "gene3": {
                "name": ["orf19.123"],
                "product": ["Hypothetical protein"],
                "note": ["Predicted protein"],
                "db_xref": ["UniProtKB:Q5A7N1"],
            },
            "gene4": {
                "name": ["CDC42"],
                "product": [
                    "Required for cell division and establishment of cell polarity"
                ],
                "note": ["GTP-binding protein"],
                "db_xref": ["UniProtKB:P19073"],
                "ec_number": ["3.6.5.-"],
                "go_terms": ["GO:0005525", "GO:0007049"],
            },
            "gene5": {
                "product": ["Hypothetical protein"],
                "note": ["Predicted protein"],
            },
        }

        # Create sample Genes dictionary
        self.genes = {
            "g1": {
                "ids": ["gene1"],
                "product": ["Unknown"],
                "db_xref": [[]],
                "note": [[]],
                "ec_number": [[]],
                "go_terms": [[]],
            },
            "g2": {
                "ids": ["gene2"],
                "product": ["Unknown"],
                "db_xref": [[]],
                "note": [[]],
                "ec_number": [[]],
                "go_terms": [[]],
            },
            "g3": {
                "ids": ["gene3"],
                "product": ["Unknown"],
                "db_xref": [[]],
                "note": [[]],
                "ec_number": [[]],
                "go_terms": [[]],
            },
            "g4": {
                "ids": ["gene4"],
                "product": ["Unknown"],
                "db_xref": [[]],
                "note": [[]],
                "ec_number": [[]],
                "go_terms": [[]],
            },
            "g5": {
                "ids": ["gene5"],
                "product": ["Unknown"],
                "db_xref": [[]],
                "note": [[]],
                "ec_number": [[]],
                "go_terms": [[]],
            },
        }

        # Create a NameCleaner instance
        self.cleaner = NameCleaner()

    def tearDown(self):
        """Clean up test environment."""
        # Stop the patcher
        self.patcher.stop()

        # Remove the temporary directory
        shutil.rmtree(self.test_dir)

    def test_clean_merged_annotations(self):
        """Test cleaning merged annotations."""
        # Clean the merged annotations
        cleaned_merged = {}
        for gene_id, annot in self.merged.items():
            cleaned_annot = self.cleaner.process_annotation(annot)
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
        self.assertEqual(
            cleaned_merged["gene4"]["product"], ["Cell division control protein 42"]
        )

        # Check that other fields are preserved
        self.assertEqual(cleaned_merged["gene1"]["note"], ["Cytoskeletal protein"])
        self.assertEqual(cleaned_merged["gene1"]["db_xref"], ["UniProtKB:P60010"])
        self.assertEqual(cleaned_merged["gene1"]["ec_number"], ["3.6.4.-"])
        self.assertEqual(
            cleaned_merged["gene1"]["go_terms"], ["GO:0005524", "GO:0005856"]
        )

    def test_apply_annotations_to_genes(self):
        """Test applying cleaned annotations to genes."""
        # Clean the merged annotations
        cleaned_merged = {}
        for gene_id, annot in self.merged.items():
            cleaned_annot = self.cleaner.process_annotation(annot)
            cleaned_merged[gene_id] = cleaned_annot

        # Apply annotations to genes
        annotation = {}
        for k, v in self.genes.items():
            n = v.copy()
            n["gene_synonym"] = []  # Initialize gene_synonym

            for i, x in enumerate(n["ids"]):
                if x in cleaned_merged:  # then functional annotation to add
                    fa = cleaned_merged.get(x)
                    if "product" in fa:
                        n["product"][i] = fa["product"][0]
                    if "db_xref" in fa:
                        n["db_xref"][i] = fa["db_xref"]
                    if "name" in fa:
                        n["name"] = fa["name"][0]
                        if len(fa["name"]) > 1:
                            n["gene_synonym"] += fa["name"][1:]
                    if "note" in fa:
                        n["note"][i] = fa["note"]
                    if "ec_number" in fa:
                        n["ec_number"][i] = fa["ec_number"]
                    if "go_terms" in fa:
                        n["go_terms"][i] = fa["go_terms"]
            annotation[k] = n

        # Check that annotations were applied correctly
        self.assertEqual(annotation["g1"]["product"][0], "Actin")
        self.assertEqual(annotation["g1"]["name"], "ACT1")
        self.assertEqual(annotation["g1"]["gene_synonym"], [])

        self.assertEqual(annotation["g2"]["product"][0], "GTPase YPT7")
        self.assertEqual(annotation["g2"]["name"], "YPT7")

        self.assertEqual(annotation["g3"]["product"][0], "Hypothetical protein")
        self.assertNotIn("name", annotation["g3"])

        self.assertEqual(
            annotation["g4"]["product"][0], "Cell division control protein 42"
        )
        self.assertEqual(annotation["g4"]["name"], "CDC42")

        self.assertEqual(annotation["g5"]["product"][0], "Hypothetical protein")
        self.assertNotIn("name", annotation["g5"])


if __name__ == "__main__":
    unittest.main()
