import json
from types import SimpleNamespace

from gene_drug_discovery import (attach_gene_drug_discovery, attach_gene_to_herb_candidates,
                                 collect_gene_evidence)


def empty_herb_lookup(genes, cache_dir, **kwargs):
    return {"state": "complete", "requested_genes": len(genes), "completed_genes": len(genes), "errors": [],
            "by_gene": {gene: {"complete": True, "herb": [], "formula": [], "ingredient": []} for gene in genes}}


def test_collects_disease_genes_and_keeps_proxy_scope_separate():
    report = {
        "target_evidence": [{"gene": "IL6", "source": "Open Targets", "mapping_scope": "exact", "score": 0.8}],
        "omics": {"dataset": {"accession": "GSE1"}, "design": {"phenotype_scope": "direct"},
                  "analysis": {"results": [{"gene": "TNF", "FDR": 0.01, "pvalue": 0.001, "log2FC": 1.2},
                                             {"gene": "JUNK", "FDR": 0.8, "pvalue": 0.7}] }},
        "multi_cohort_cad_proxy": {"meta_analysis": {"results": [
            {"gene": "CRP", "study_n": 3, "direction_consistency": 3, "pvalue": 0.02}
        ]}},
    }
    genes = collect_gene_evidence(report)
    assert set(genes) == {"IL6", "TNF", "CRP"}
    assert genes["IL6"]["evidence"][0]["scope"] == "exact"
    assert genes["CRP"]["evidence"][0]["scope"] == "broad_CAD_proxy"


def test_collects_archived_gene_prioritization_key_without_recomputing():
    report = {
        "gene_prioritization": {
            "eligible_genes": [{"gene": "LPA", "tier": "medium", "reason": "exact GWAS",
                                "sources": [{"source": "GWAS Catalog", "scope": "exact", "pvalue": 1e-9}]}],
            "proxy_genes": [{"gene": "CRP", "evidence": [{"source": "GWAS Catalog", "scope": "related_proxy"}]}],
        }
    }
    genes = collect_gene_evidence(report)
    assert set(genes) == {"LPA", "CRP"}
    assert genes["LPA"]["tier"] == "medium"
    assert genes["LPA"]["evidence"][0]["scope"] == "exact"
    assert genes["CRP"]["tier"] == "proxy"


def test_legacy_tiers_are_recomputed_when_raw_data_exists():
    report = {
        "gene_prioritization": {"eligible_genes": [
            {"gene": "NOMINAL", "tier": "medium", "sources": [{"source": "legacy"}]}
        ], "proxy_genes": []},
        "omics": {"dataset": {"accession": "GSE1"}, "design": {"phenotype_scope": "direct"},
                  "analysis": {"results": [{"gene": "NOMINAL", "FDR": 0.3, "pvalue": 0.02, "log2FC": 2.0}]}},
    }
    assert collect_gene_evidence(report) == {}


def test_drug_lookup_uses_archived_eligible_genes(tmp_path):
    report = {"gene_prioritization": {
        "eligible_genes": [{"gene": "LPA", "tier": "medium", "reason": "exact GWAS",
                            "sources": [{"source": "GWAS Catalog", "scope": "exact", "pvalue": 1e-9}]}],
        "proxy_genes": [],
    }}

    def network(url, payload, headers, timeout):
        request = json.loads(payload.decode("utf-8"))
        assert request["variables"]["genes"] == ["LPA"]
        response = {"data": {"genes": {"nodes": [{"name": "LPA", "longName": "Lipoprotein(a)",
                    "interactions": [{"drug": {"name": "Example drug", "conceptId": "C1"},
                                     "sources": [{"sourceDbName": "ExampleDB"}], "publications": [],
                                     "interactionTypes": [], "interactionAttributes": []}]}]}}}
        return json.dumps(response).encode("utf-8"), "application/json"

    result = attach_gene_drug_discovery(report, network, tmp_path, herb_retrieve=empty_herb_lookup)
    assert result["state"] == "completed"
    assert result["candidate_gene_count"] == 1
    assert result["genes_queried"] == 1
    assert result["candidates"][0]["gene"] == "LPA"
    assert result["candidates"][0]["drug"] == "Example drug"


def test_gene_first_lookup_merges_dgidb_with_chembl(tmp_path):
    report = {
        "target_evidence": [{"gene": "IL6", "source": "Open Targets", "mapping_scope": "exact", "score": 0.7}],
        "drug_candidates": [{"gene": "IL6", "molecule": "Example compound", "molecule_chembl_id": "CHEMBL1",
                             "mechanism": "inhibitor", "url": "https://www.ebi.ac.uk/chembl/explore/compound/CHEMBL1"}],
    }

    def network(url, payload, headers, timeout):
        assert url == "https://dgidb.org/api/graphql"
        request = json.loads(payload.decode("utf-8"))
        assert request["variables"]["genes"] == ["IL6"]
        response = {"data": {"genes": {"nodes": [{"name": "IL6", "longName": "Interleukin 6", "interactions": [
            {"drug": {"name": "Example compound", "conceptId": "CHEMBL1"}, "interactionScore": 0.4,
             "interactionTypes": [{"type": "inhibitor", "directionality": "negative"}],
             "interactionAttributes": [], "publications": [{"pmid": "12345"}],
             "sources": [{"sourceDbName": "ExampleDB"}]}
        ]}]}}}
        return json.dumps(response).encode("utf-8"), "application/json"

    result = attach_gene_drug_discovery(report, network, tmp_path, herb_retrieve=empty_herb_lookup)
    assert result["state"] == "completed"
    assert result["candidate_gene_count"] == 1
    assert len(result["candidates"]) == 1
    candidate = result["candidates"][0]
    assert candidate["gene"] == "IL6"
    assert set(candidate["evidence_databases"]) == {"ChEMBL", "DGIdb"}
    assert candidate["pmids"] == ["12345"]
    assert "not disease-specific efficacy" in candidate["evidence_scope"]


def test_dgidb_batches_all_eligible_genes_instead_of_dropping_after_40(tmp_path):
    genes = [f"G{i:03d}" for i in range(45)]
    batches = []

    def network(url, payload, headers, timeout):
        batch = json.loads(payload.decode("utf-8"))["variables"]["genes"]
        batches.append(batch)
        response = {"data": {"genes": {"nodes": [
            {"name": gene, "longName": gene, "interactions": []} for gene in batch
        ]}}}
        return json.dumps(response).encode("utf-8"), "application/json"

    from gene_drug_discovery import retrieve_dgidb_interactions
    result = retrieve_dgidb_interactions(genes, network, tmp_path)
    assert result["state"] == "completed"
    assert result["queried_genes"] == genes
    assert result["mapped_genes"] == genes
    assert [len(batch) for batch in batches] == [40, 5]
    assert result["batch_count"] == 2


def test_dgidb_failure_is_reported_without_losing_existing_chembl(tmp_path):
    report = {
        "target_evidence": [{"gene": "TNF", "source": "GEO", "mapping_scope": "direct_phenotype",
                             "FDR": 0.01, "log2FC": 0.8}],
        "drug_candidates": [{"gene": "TNF", "molecule": "Known molecule"}],
    }

    def network(*args, **kwargs):
        raise TimeoutError("connection timed out")

    result = attach_gene_drug_discovery(report, network, tmp_path, herb_retrieve=empty_herb_lookup)
    assert result["state"] == "partial"
    assert result["source_status"]["DGIdb"]["state"] == "failed"
    assert result["source_status"]["DGIdb"]["error"] == "上游数据库请求超时；已保留其他数据源结果。"
    assert result["candidates"][0]["drug"] == "Known molecule"


def test_gene_first_herb_lookup_uses_shared_disease_targets(tmp_path):
    report = {"target_evidence": [{"gene": "IL6", "source": "Open Targets", "score": 0.7}]}
    calls = []

    def retrieve(genes, cache_dir, **kwargs):
        calls.append((genes, cache_dir, kwargs))
        return {"state": "complete", "requested_genes": 1, "completed_genes": 1, "errors": [],
                "by_gene": {"IL6": {"complete": True, "herb": [{
                    "entity_id": "HBH0001", "entity_type": "Herb", "name": "示例药材",
                    "relation_source": "HERB herb_target", "fdr_bh": 0.01}],
                    "formula": [], "ingredient": []}}}

    result = attach_gene_to_herb_candidates(report, tmp_path, retrieve)
    assert calls[0][0] == ["IL6"]
    assert result["candidate_count"] == 1
    assert result["candidates"][0]["name"] == "示例药材"
    assert result["candidates"][0]["shared_targets"][0]["gene"] == "IL6"
