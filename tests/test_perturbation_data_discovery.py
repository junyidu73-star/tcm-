import json
import unittest
from urllib.parse import parse_qs, urlparse

from perturbation_data_discovery import discover_candidate_perturbation_studies


class PerturbationDataDiscoveryTests(unittest.TestCase):
    def test_geo_and_biostudies_hits_are_deduplicated_but_sources_retained(self):
        calls = []

        def get(url):
            calls.append(url)
            parsed = urlparse(url)
            query = parse_qs(parsed.query)
            if "esearch.fcgi" in parsed.path:
                return json.dumps({"esearchresult": {"count": "1", "idlist": ["101"]}}).encode()
            if "esummary.fcgi" in parsed.path:
                return json.dumps({"result": {"uids": ["101"], "101": {
                    "accession": "GSE99999", "title": "Herb treatment expression", "summary": "study",
                    "gpl": ["GPL1"], "taxon": "Homo sapiens",
                    "gdstype": "Expression profiling by high throughput sequencing",
                    "samples": [{"accession": "GSM1", "title": "herb extract dose rep1"},
                                {"accession": "GSM2", "title": "herb extract dose rep2"},
                                {"accession": "GSM3", "title": "Vehicle DMSO rep1"},
                                {"accession": "GSM4", "title": "Vehicle DMSO rep2"}]}}}).encode()
            if "biostudies" in parsed.path:
                return json.dumps({"hits": [{"accession": "E-GEOD-99999", "title": "Herb treatment expression",
                                               "content": "study"}]}).encode()
            raise AssertionError(url)

        result = discover_candidate_perturbation_studies(
            get, lambda *_: None,
            [{"entity_id": "H1", "Herb en name": "herb"}],
            [{"entity_id": "F1", "Formula en name": "formula"}], pause=0)
        self.assertEqual(result["candidate_count"], 2)
        self.assertEqual(result["study_record_count"], 2)
        herb = next(record for record in result["records"] if record["candidate_id"] == "H1")
        self.assertEqual(herb["sources"], ["GEO/E-utilities", "BioStudies/ArrayExpress"])
        self.assertEqual(herb["metadata_only"], True)
        self.assertEqual(herb["candidate_sample_count"], 2)
        self.assertEqual(herb["control_like_sample_count"], 2)
        self.assertEqual(herb["sample_design_triage"], "candidate_and_control_labels_found_review_required")
        self.assertEqual(len(calls), 6)

    def test_missing_entity_alias_is_not_searched_as_database_id(self):
        result = discover_candidate_perturbation_studies(lambda _url: b"", lambda *_: None,
                                                          [{"entity_id": "H1", "Herb cn name": "NA"}], pause=0)
        self.assertEqual(result["query_count"], 0)
        self.assertEqual(result["study_record_count"], 0)


if __name__ == "__main__":
    unittest.main()
