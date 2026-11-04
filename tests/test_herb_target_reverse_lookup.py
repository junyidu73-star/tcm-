import unittest

from herb_target_reverse_lookup import _parse_target_relation


class HerbTargetReverseLookupTests(unittest.TestCase):
    def test_official_herb_detail_field_names_preserve_fdr(self):
        parsed = _parse_target_relation({
            "Herb id": {"title": "HERB000001"},
            "Herb name": "Example herb",
            "P-value": "0.002",
            "FDR_BH": "0.018",
        }, "HERB herb_target", "MMP9")

        self.assertEqual(parsed["entity_id"], "HERB000001")
        self.assertEqual(parsed["entity_type"], "Herb")
        self.assertEqual(parsed["p_value"], "0.002")
        self.assertEqual(parsed["fdr_bh"], "0.018")

    def test_paper_target_rows_keep_pubmed_provenance(self):
        parsed = _parse_target_relation({
            "Ingredient id": {"title": "HBIN000001"},
            "Ingredient name": "Example compound",
            "PubMed id": "12345678",
            "Paper id": "HBREF000001",
        }, "HERB paper_ingredient_target", "MMP9")

        self.assertEqual(parsed["entity_type"], "Ingredient")
        self.assertEqual(parsed["references"][0]["pmid"], "12345678")
        self.assertEqual(parsed["references"][0]["reference_id"], "HBREF000001")


if __name__ == "__main__":
    unittest.main()
