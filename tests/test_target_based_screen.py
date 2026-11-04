import json
import unittest
from pathlib import Path

from target_based_screen import screen_herbs_by_shared_targets, screen_from_reverse_target_links


class TargetBasedScreenTests(unittest.TestCase):
    def test_ingredient_target_overlap_is_reported_as_its_own_candidate_type(self):
        report = {
            "multi_cohort_cad_proxy": {"meta_analysis": {"significant_n": 0, "tested_genes": 1, "results": [
                {"gene": "MMP9", "study_n": 2, "direction_consistency": 2, "pvalue": 0.01, "FDR": 0.5}
            ]}},
            "multi_source_target_evidence": {},
            "external_tcm_database": {"herb_candidates": [], "formula_candidates": []},
        }
        reverse = {"state": "complete", "requested_genes": 1, "completed_genes": 1, "errors": [], "by_gene": {
            "MMP9": {"complete": True, "herb": [], "formula": [], "ingredient": [
                {"entity_id": "HBIN000001", "entity_type": "Ingredient", "name": "example compound",
                 "gene": "MMP9", "relation_source": "HERB ingredient_target"}
            ]}
        }}
        result = screen_from_reverse_target_links(report, reverse)
        self.assertEqual(result["ingredient_candidate_count"], 0)
        self.assertEqual(result["state"], "no_reverse_target_overlap")

    def test_proxy_and_nominal_meta_genes_do_not_enter_herb_target_pool(self):
        report = {
            "multi_cohort_cad_proxy": {"meta_analysis": {"tested_genes": 1, "significant_n": 0, "results": [
                {"gene": "MMP9", "study_n": 3, "direction_consistency": 3,
                 "pvalue": 0.001, "FDR": 0.8, "log2FC": 0.2}
            ]}}
        }
        from target_based_screen import reverse_lookup_genes
        self.assertEqual(reverse_lookup_genes(report), [])

    def test_candidate_requires_shared_target_and_expression_is_supporting_only(self):
        report = {
            "multi_source_target_evidence": {"GWAS_Catalog": {"association_sets": [{"scope": "exact", "associations": [
                {"mapped_genes": ["MMP9"], "p_value": 1e-8, "accession_id": "GCST1"}
            ]}]}},
            "multi_cohort_cad_proxy": {
                "meta_analysis": {
                    "significant_n": 0,
                    "tested_genes": 3,
                    "results": [
                        {"gene": "MMP9", "study_n": 3, "direction_consistency": 3, "pvalue": 0.01, "FDR": 0.7, "log2FC": 0.2},
                        {"gene": "MYC", "study_n": 2, "direction_consistency": 2, "pvalue": 0.02, "FDR": 0.8, "log2FC": -0.1},
                        {"gene": "FADD", "study_n": 1, "pvalue": 0.001, "FDR": 0.9, "log2FC": 0.3},
                    ],
                }
            },
            "external_tcm_database": {
                "herb_candidates": [
                    {"entity_id": "HERB0001", "Herb cn name": "候选甲"},
                    {"entity_id": "HERB0002", "Herb cn name": "无交集药材"},
                ],
                "formula_candidates": [],
            },
        }
        signatures = [
            {
                "candidate_id": "HERB0001",
                "source_accession": "GSE_TEST",
                "genes": [
                    {"gene": "MMP9", "log2FC": -1.2, "fdr": 0.03},
                    {"gene": "NOT_SHARED", "log2FC": 2, "fdr": 0.01},
                ],
            }
        ]
        cache = Path(__file__).resolve().parent / "fixtures" / "target_based_screen"
        result = screen_herbs_by_shared_targets(report, cache, signatures, top_n=2)

        self.assertEqual(result["state"], "multi_source_proxy_target_screen")
        self.assertEqual(len(result["candidates"]), 1)
        candidate = result["candidates"][0]
        self.assertEqual(candidate["shared_target_count"], 1)
        self.assertEqual(candidate["shared_targets"][0]["gene"], "MMP9")
        measured = candidate["shared_targets"][0]["measured_perturbation_at_target"]
        self.assertEqual(len(measured), 1)
        self.assertEqual(measured[0]["log2FC"], -1.2)


if __name__ == "__main__":
    unittest.main()
