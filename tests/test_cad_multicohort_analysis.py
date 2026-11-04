import unittest

import pandas as pd

from cad_multicohort_analysis import _random_effects_meta


class CadMultiCohortAnalysisTests(unittest.TestCase):
    def test_meta_analysis_requires_cross_study_direction_support(self):
        studies = [
            pd.DataFrame({"gene": ["CONSISTENT", "CONFLICT"], "log2FC": [0.8, 1.0], "SE": [0.1, 0.1]}),
            pd.DataFrame({"gene": ["CONSISTENT", "CONFLICT"], "log2FC": [0.9, -1.0], "SE": [0.1, 0.1]}),
            pd.DataFrame({"gene": ["CONSISTENT", "CONFLICT"], "log2FC": [0.7, 1.0], "SE": [0.1, 0.1]}),
        ]
        result = _random_effects_meta(studies).set_index("gene")
        self.assertEqual(result.loc["CONSISTENT", "study_n"], 3)
        self.assertEqual(result.loc["CONSISTENT", "direction_consistency"], 3)
        self.assertGreater(result.loc["CONSISTENT", "log2FC"], 0)
        self.assertEqual(result.loc["CONFLICT", "direction_consistency"], 2)
        self.assertGreater(result.loc["CONFLICT", "I2"], 0)


if __name__ == "__main__":
    unittest.main()
