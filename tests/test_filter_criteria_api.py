import pytest

from server import validate_gene_filter_criteria


def test_filter_criteria_defaults_match_platform_rules():
    values = validate_gene_filter_criteria({})
    assert values == {
        "fdr_threshold": 0.05,
        "abs_log2fc_min": 0.5,
        "exploratory_p_threshold": 0.01,
        "gwas_p_threshold": 5e-8,
        "opentargets_score_min": 0.1,
        "meta_direction_concordance_min": 0.75,
        "meta_i2_max": 0.75,
    }


def test_gwas_threshold_cannot_be_weakened_past_genome_wide_significance():
    with pytest.raises(ValueError):
        validate_gene_filter_criteria({"gwas_p_threshold": 1e-7})


def test_custom_project_filters_are_kept():
    values = validate_gene_filter_criteria({"fdr_threshold": 0.03, "abs_log2fc_min": 0.8,
                                            "meta_i2_max": 0.5})
    assert values["fdr_threshold"] == 0.03
    assert values["abs_log2fc_min"] == 0.8
    assert values["meta_i2_max"] == 0.5
