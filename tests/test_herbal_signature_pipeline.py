import unittest

from herbal_signature_pipeline import score_reversal


class HerbalSignaturePipelineTests(unittest.TestCase):
    def test_inverse_signature_gets_positive_exploratory_score(self):
        disease = [{"gene": f"G{i}", "log2FC": float(i)} for i in range(1, 6)]
        treatment = [{"gene": f"G{i}", "log2FC": -float(i)} for i in range(1, 6)]
        result = score_reversal(disease, [{
            "candidate_id": "HERB1", "name": "test herb", "source_id": "EXP1",
            "genes": treatment, "context": {"species": "Homo sapiens", "tissue": "PBMC"},
        }], min_shared_genes=5, disease_context={"tissue": "PBMC"})
        self.assertEqual(result["state"], "scored_exploratory")
        self.assertGreater(result["scores"][0]["connectivity_score"], 0.99)
        self.assertEqual(result["scores"][0]["context_match"], "same_context")
        self.assertIsNone(result["p_value"])
        self.assertIsNone(result["fdr"])

    def test_nonhuman_and_insufficient_signatures_are_rejected(self):
        disease = [{"gene": "A", "log2FC": 1.0}]
        result = score_reversal(disease, [
            {"source_id": "mouse", "genes": [{"gene": "A", "log2FC": -1}],
             "context": {"species": "Mus musculus"}},
            {"source_id": "few", "genes": [{"gene": "A", "log2FC": -1}],
             "context": {"species": "Homo sapiens"}},
        ], min_shared_genes=2)
        self.assertEqual(result["state"], "no_eligible_measured_signatures")
        self.assertEqual(len(result["rejected"]), 2)

    def test_weak_disease_signature_is_diagnostic_only(self):
        disease = [{"gene": f"G{i}", "log2FC": float(i)} for i in range(1, 6)]
        treatment = [{"gene": f"G{i}", "log2FC": -float(i)} for i in range(1, 6)]
        result = score_reversal(disease, [{
            "candidate_id": "HERB1", "name": "test herb", "source_id": "EXP1",
            "genes": treatment, "context": {"species": "Homo sapiens", "tissue": "PBMC"},
        }], min_shared_genes=5, disease_context={"tissue": "PBMC"},
            disease_significant_n=0, min_disease_significant_genes=2)
        self.assertEqual(result["state"], "insufficient_disease_signature")
        self.assertEqual(result["scores"], [])
        self.assertEqual(len(result["diagnostic_correlations"]), 1)

    def test_low_confidence_mode_emits_ranked_hypothesis_without_calling_it_recommendation(self):
        disease = [{"gene": f"G{i}", "log2FC": float(i)} for i in range(1, 6)]
        treatment = [{"gene": f"G{i}", "log2FC": -float(i)} for i in range(1, 6)]
        result = score_reversal(disease, [{
            "candidate_id": "HERB1", "name": "test herb", "source_id": "EXP1",
            "genes": treatment, "context": {"species": "Homo sapiens", "tissue": "whole blood"},
        }], min_shared_genes=5, disease_context={"tissue": "peripheral blood"},
            disease_significant_n=0, min_disease_significant_genes=2,
            allow_low_confidence=True)
        self.assertEqual(result["state"], "exploratory_low_confidence")
        self.assertEqual(len(result["scores"]), 1)
        self.assertEqual(result["scores"][0]["evidence_label"], "exploratory_proxy_signature_match_low_confidence")
        self.assertIsNone(result["fdr"])


if __name__ == "__main__":
    unittest.main()
