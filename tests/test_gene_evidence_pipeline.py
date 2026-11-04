from gene_evidence_pipeline import prioritize_gene_evidence


def test_gene_filters_keep_exploratory_and_proxy_signals_out_of_eligible_pool():
    report = {
        "target_evidence": [
            {"gene": "GWAS_PASS", "source_database": "GWAS Catalog", "mapping_scope": "exact", "p_value": 4.9e-8},
            {"gene": "GWAS_BOUNDARY", "source_database": "GWAS Catalog", "mapping_scope": "exact", "p_value": 5e-8},
            {"gene": "OT_PASS", "source_database": "Open Targets", "mapping_scope": "exact", "score": 0.1},
            {"gene": "OT_LOW", "source_database": "Open Targets", "mapping_scope": "exact", "score": 0.099},
            {"gene": "GENERIC", "source_database": "unclassified adapter", "mapping_scope": "exact", "score": 1.0},
        ],
        "omics_cohorts": [
            {"dataset": {"accession": "GSE_DIRECT"},
             "design": {"phenotype_scope": "direct"},
             "analysis": {"results": [
                 {"gene": "DE_PASS", "FDR": 0.049, "pvalue": 0.001, "log2FC": 0.5},
                 {"gene": "NOMINAL_ONLY", "FDR": 0.051, "pvalue": 0.009, "log2FC": 1.2},
                 {"gene": "SMALL_EFFECT", "FDR": 0.01, "pvalue": 0.001, "log2FC": 0.49},
                 {"gene": "WEAK_P", "FDR": 0.8, "pvalue": 0.02, "log2FC": 1.0},
             ]}},
            {"dataset": {"accession": "GSE_PROXY"},
             "design": {"phenotype_scope": "proxy"},
             "analysis": {"results": [
                 {"gene": "PROXY_HIT", "FDR": 0.001, "pvalue": 1e-5, "log2FC": 2.0},
             ]}},
        ],
    }
    result = prioritize_gene_evidence(report)
    eligible = {row["gene"] for row in result["eligible_genes"]}
    exploratory = {row["gene"] for row in result["exploratory_genes"]}
    proxy = {row["gene"] for row in result["proxy_genes"]}
    assert eligible == {"GWAS_PASS", "OT_PASS", "DE_PASS"}
    assert exploratory == {"NOMINAL_ONLY", "SMALL_EFFECT"}
    assert "PROXY_HIT" in proxy
    assert not eligible & (exploratory | proxy)


def test_project_thresholds_are_applied_and_saved_in_prioritization():
    report = {"target_evidence": [
        {"gene": "OT_TARGET", "source_database": "Open Targets", "mapping_scope": "exact", "score": 0.15},
        {"gene": "GWAS_TARGET", "source_database": "GWAS Catalog", "mapping_scope": "exact", "p_value": 2e-8},
    ]}
    filters = {"opentargets_score_min": 0.2, "gwas_p_threshold": 1e-8}
    result = prioritize_gene_evidence(report, filters)
    assert result["eligible_genes"] == []
    assert result["filters"]["opentargets_score_min"] == 0.2
    assert result["filters"]["gwas_p_threshold"] == 1e-8
