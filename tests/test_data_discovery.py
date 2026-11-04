import json
import gzip
import urllib.parse
import unittest
from unittest.mock import patch

import research_engine as research_module
from data_discovery import biostudies_inventory, pride_inventory, search_biostudies, search_pride
from research_engine import ResearchEngine, differential_expression, differential_expression_counts


class DataDiscoveryTests(unittest.TestCase):
    def test_geo_supplementary_counts_are_aligned_to_soft_sample_ids(self):
        accession="GSE123456"
        root="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/"
        soft=("^SERIES = GSE123456\n!Series_overall_design = Raw integer counts\n"
              "^SAMPLE = GSM1\n!Sample_title = Control one\n!Sample_geo_accession = GSM1\n!Sample_platform_id = GPL1\n!Sample_characteristics_ch1 = condition: control\n"
              "^SAMPLE = GSM2\n!Sample_title = Case one\n!Sample_geo_accession = GSM2\n!Sample_platform_id = GPL1\n!Sample_characteristics_ch1 = condition: case\n")
        table="gene_id\tGSM1\tGSM2\n"+"".join(f"GENE{i}\t10\t20\n" for i in range(101))
        responses={
            root+"matrix/":b"<html></html>",
            root+"suppl/":b'<a href="GSE123456_counts.tsv.gz">counts</a>',
            root+"soft/":b'<a href="GSE123456_family.soft.gz">SOFT</a>',
            root+"soft/GSE123456_family.soft.gz":gzip.compress(soft.encode()),
            root+"suppl/GSE123456_counts.tsv.gz":gzip.compress(table.encode()),
        }
        def network(url,*_args,**_kwargs): return responses[url],None
        engine=ResearchEngine({}, {}, None, network, None)
        metadata,data=engine.matrix(accession)
        self.assertEqual(data[0], ["gene_id","GSM1","GSM2"])
        self.assertEqual(len(data),102)
        self.assertEqual(metadata["!Sample_geo_accession"],["GSM1","GSM2"])
        self.assertEqual(engine._last_matrix_input["filename"],"GSE123456_counts.tsv.gz")

    def test_paired_expression_analysis_uses_within_pair_differences(self):
        result = differential_expression(
            [[2.0, 1.0, 4.0, 3.0, 6.0, 5.0], [2.0, 2.0, 3.0, 3.0, 4.0, 4.0]],
            ["GENE1", "GENE2"], [0, 2, 4], [1, 3, 5],
            [([0], [1]), ([2], [3]), ([4], [5])],
        )
        self.assertTrue(result["paired"])
        self.assertEqual(result["case_n"], 3)
        self.assertIn("Paired t-test", result["method"])
        self.assertEqual(result["results"][0]["gene"], "GENE1")
        self.assertAlmostEqual(result["results"][0]["log2FC"], 1.0)

    def test_raw_count_route_rejects_normalized_decimal_values(self):
        with self.assertRaisesRegex(ValueError, "non-negative integer|非负整数"):
            differential_expression_counts([[1.0, 1.5, 2.0, 3.0, 4.0, 5.0]], ["GENE1"],
                                            ["s1","s2","s3","s4","s5","s6"],
                                            [0,1,2],[3,4,5])

    def test_biostudies_cursor_pagination_and_deduplication(self):
        calls = []
        def get(url):
            q = urllib.parse.parse_qs(urllib.parse.urlsplit(url).query)
            calls.append(q)
            if q.get("cursor"):
                return json.dumps({"hits":[{"accession":"E-GEOD-2","title":"Disease B"}],"nextCursor":None,"totalHits":2,"isTotalHitsExact":True}).encode()
            return json.dumps({"hits":[{"accession":"E-GEOD-1","title":"Disease A"}],"nextCursor":"next","totalHits":2,"isTotalHitsExact":True}).encode()
        with patch("data_discovery.time.sleep", lambda *_: None):
            result = search_biostudies(get, lambda *_: None, ["stable angina"])
        self.assertTrue(result["complete"])
        self.assertEqual(result["retrieved"], 2)
        self.assertEqual(calls[0]["query"], ['title:"stable angina"'])

    def test_pride_search_reads_full_pages(self):
        def get(url):
            q = urllib.parse.parse_qs(urllib.parse.urlsplit(url).query)
            page = int(q["page"][0])
            rows = [{"accession":f"PXD{i:06d}","title":f"Study {i}"} for i in range(page*100,(page+1)*100)] if page == 0 else [{"accession":"PXD000100","title":"Study 100"}]
            return json.dumps(rows).encode()
        with patch("data_discovery.time.sleep", lambda *_: None):
            result = search_pride(get, lambda *_: None, ["angina"])
        self.assertTrue(result["complete"])
        self.assertEqual(result["retrieved"], 101)
        self.assertEqual(result["queries"][0]["pages"], 2)

    def test_geo_search_marks_missing_summary_as_incomplete(self):
        def network(url, data=None, headers=None, timeout=None):
            q = urllib.parse.parse_qs(urllib.parse.urlsplit(url).query)
            if "esearch.fcgi" in url:
                obj = {"esearchresult":{"count":"1","idlist":["1"]}}
            else:
                obj = {"result":{"uids":[]}}
            return json.dumps(obj).encode(), "application/json"
        engine = ResearchEngine({}, {}, None, network, None)
        with patch("research_engine.time.sleep", lambda *_: None):
            result = engine.search_geo("disease")
        self.assertFalse(result["complete"])
        self.assertEqual(result["retrieved_ids"], 1)
        self.assertTrue(result["errors"])

    def test_disease_relevance_uses_synonyms_and_angina_boundaries(self):
        direct = ResearchEngine.dataset_relevance({"title":"Type 2 diabetes cohort"}, "diabetes mellitus", ["type 2 diabetes"])
        self.assertEqual(direct["tier"], "direct_phenotype")
        unstable = ResearchEngine.dataset_relevance({"title":"Unstable angina patients"}, "stable angina")
        self.assertEqual(unstable["tier"], "related_proxy")
        unrelated = ResearchEngine.dataset_relevance({"title":"Stable cell line study"}, "stable angina")
        self.assertEqual(unrelated["tier"], "unrelated_or_unverified")

    def test_file_inventories_only_expand_direct_disease_records(self):
        engine = ResearchEngine({}, {}, None, lambda *_args, **_kwargs: json.dumps({"data":{"search":{"hits":[]}}}).encode(), None)
        ae = {"complete":True,"records":[
            {"accession":"E-GEOD-1","title":"Stable angina transcriptome","canonical_study":"GSE1"},
            {"accession":"E-GEOD-2","title":"Unrelated cell study","canonical_study":"GSE2"}]}
        pride = {"complete":True,"records":[
            {"accession":"PXD000001","title":"Stable angina proteomics"},
            {"accession":"PXD000002","title":"Unrelated cell study"}]}
        called=[]
        with patch.object(research_module,"search_biostudies",return_value=ae), \
             patch.object(research_module,"search_pride",return_value=pride), \
             patch.object(research_module,"biostudies_inventory",side_effect=lambda _get, acc: called.append(acc) or {"status":"inventory_only","files":[]}), \
             patch.object(research_module,"pride_inventory",side_effect=lambda _get, acc: called.append(acc) or {"status":"inventory_only","files":[]}):
            result=engine.additional_databases("stable angina",False,["stable angina"])
        self.assertEqual(set(called),{"E-GEOD-1","PXD000001"})
        self.assertEqual(result["PRIDE"]["records"][1]["file_inventory"]["status"],"not_enumerated_relevance")
        self.assertEqual(result["ArrayExpress"]["records"][1]["file_inventory"]["status"],"not_enumerated_relevance")

    def test_biostudies_file_tree_recurses_and_deduplicates(self):
        def get(_):
            return json.dumps({"section":{"files":[[{"path":"matrix.tsv"}]],"subsections":[{"files":[{"path":"matrix.tsv"},{"path":"sample.txt"}]}]}}).encode()
        result = biostudies_inventory(get, "E-MTAB-1")
        self.assertEqual([f["name"] for f in result["files"]], ["matrix.tsv", "sample.txt"])

    def test_pride_inventory_retains_public_download_and_checksum_metadata(self):
        def get(_):
            return json.dumps([{"fileName":"sample.raw","fileSizeBytes":123,"checksum":"abc","publicFileLocations":[{"name":"FTP","value":"ftp://ftp.pride.ebi.ac.uk/pride/data/archive/sample.raw"}]}]).encode()
        result = pride_inventory(get, "PXD000001")
        self.assertEqual(result["files"][0]["checksum"], "abc")
        self.assertEqual(result["files"][0]["size"], 123)
        self.assertEqual(result["files"][0]["locations"][0]["value"], "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/sample.raw")
        self.assertEqual(result["files"][0]["url"], "https://ftp.pride.ebi.ac.uk/pride/data/archive/sample.raw")

    def test_europe_pmc_partial_failure_preserves_prior_pages(self):
        class Engine(ResearchEngine):
            def __init__(self): super().__init__({}, {}, None, None, None); self.n = 0
            def get(self, _):
                self.n += 1
                if self.n == 1:
                    return json.dumps({"hitCount":2,"nextCursorMark":"next","resultList":{"result":[{"pmid":"1","title":"A"}]}}).encode()
                raise TimeoutError("page failed")
        engine = Engine()
        with patch("research_engine.time.sleep", lambda *_: None):
            result = engine.search_europe_pmc("TITLE_ABS:disease")
        self.assertFalse(result["complete"])
        self.assertEqual(result["retrieved"], 1)
        self.assertEqual(result["results"][0]["pmid"], "1")
        self.assertTrue(result["errors"])


if __name__ == "__main__":
    unittest.main()
