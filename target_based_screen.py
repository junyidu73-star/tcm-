"""Target-first herb screening from disease candidate genes and HERB target evidence."""
from __future__ import annotations

import json
from pathlib import Path

DIRECT_SCOPES = {"exact", "direct", "direct_phenotype", "exact_clinical_concept", "direct_disease_candidate"}


def _is_direct_scope(scope):
    return str(scope or "").casefold() in DIRECT_SCOPES


def _paper_supports(reference):
    if not isinstance(reference, list) or not reference:
        return []
    headers = [str(x) for x in reference[0]] if isinstance(reference[0], list) else []
    supports = []
    for row in reference[1:]:
        if not isinstance(row, list) or len(row) != len(headers):
            continue
        item = dict(zip(headers, row))
        supports.append({
            "pmid": str((item.get("PubMed ID") or {}).get("title", "")) if isinstance(item.get("PubMed ID"), dict) else str(item.get("PubMed ID") or ""),
            "title": str(item.get("Reference title") or ""),
            "relationship": str(item.get("Relationship") or ""),
            "grade": str(item.get("Grade") or ""),
            "supporting_sentence": str(item.get("Supporting sentences") or ""),
        })
    return supports


def _herb_target_map(detail):
    targets = {}
    for row in (detail.get("herb_target") or [])[1:]:
        if not isinstance(row, list) or len(row) < 6 or not isinstance(row[1], str):
            continue
        symbol = row[1].strip().upper()
        if not symbol:
            continue
        item = targets.setdefault(symbol, {"gene": symbol, "evidence": []})
        item["evidence"].append({"source": "HERB herb_target", "target_id": str((row[0] or {}).get("title", "")) if isinstance(row[0], dict) else str(row[0]),
                                 "p_value": row[5], "fdr_bh": row[6] if len(row) > 6 else None})
    for row in (detail.get("drug_paper_target") or [])[1:]:
        if not isinstance(row, list) or len(row) < 4 or not isinstance(row[1], str):
            continue
        symbol = row[1].strip().upper()
        if not symbol:
            continue
        item = targets.setdefault(symbol, {"gene": symbol, "evidence": []})
        item["evidence"].append({"source": "HERB drug_paper_target", "target_id": str((row[0] or {}).get("title", "")) if isinstance(row[0], dict) else str(row[0]),
                                 "supports": _paper_supports(row[3])})
    return targets


def _candidate_disease_targets(report, top_n=None):
    """Return only genes admitted by the shared, explicit evidence tiers.

    Earlier versions independently re-filtered raw rows here using nominal
    P<0.05 and a proxy-meta P<0.1 rule. That bypassed the main gene pipeline
    and let exploratory/proxy signals enter herb ranking. Reuse its eligible
    tiers so every downstream path follows the same criteria.
    """
    from gene_evidence_pipeline import prioritize_gene_evidence

    prioritization = (report.get("gene_evidence_prioritization")
                      or report.get("gene_prioritization"))
    raw_evidence_available = bool(report.get("multi_source_target_evidence") or report.get("omics_cohorts")
                                  or report.get("omics") or report.get("target_evidence"))
    if (not isinstance(prioritization, dict)
            or not isinstance(prioritization.get("eligible_genes"), list)
            or (raw_evidence_available and not isinstance(prioritization.get("filters"), dict))):
        prioritization = prioritize_gene_evidence(report)
    merged = {}
    for item in prioritization.get("eligible_genes", []):
        if not isinstance(item, dict) or item.get("tier") not in {"high", "medium"}:
            continue
        symbol = str(item.get("gene") or "").strip().upper()
        if not symbol or symbol in {"NA", "N/A", "NONE", "NULL"}:
            continue
        sources = []
        for evidence in item.get("sources", []) or []:
            source = str(evidence.get("source") or "qualified disease evidence")
            if "open target" in source.casefold():
                database = "Open Targets"
            elif "gwas" in source.casefold():
                database = "GWAS Catalog"
            else:
                database = source
            sources.append({**evidence, "source_database": database,
                            "scope": evidence.get("scope") or "direct_disease_candidate",
                            "tier": item.get("tier"), "reason": item.get("reason")})
        merged[symbol] = {"gene": symbol, "tier": item["tier"], "reason": item.get("reason"),
                          "sources": sources, "meta_analysis": item.get("meta_analysis")}

    meta = prioritization.get("random_effects_meta") or {}
    ranked_n = meta.get("tested_gene_count", 0)
    if not ranked_n:
        proxy_meta = (report.get("multi_cohort_cad_proxy") or {}).get("meta_analysis") or {}
        ranked_n = proxy_meta.get("tested_genes", len(proxy_meta.get("results") or []))
    significant_n = meta.get("high_priority_count", 0)
    if top_n is not None:
        ordered = sorted(merged.values(), key=lambda x: (
            x.get("tier") != "high",
            -len({s.get("source_database") for s in x["sources"]}),
            min((s.get("pvalue") for s in x["sources"] if isinstance(s.get("pvalue"), (int, float))), default=1),
            x["gene"]))[:top_n]
        merged = {x["gene"]: x for x in ordered}
    return merged, ranked_n, significant_n, ranked_n


def screen_herbs_by_shared_targets(report, herb_cache_dir: str | Path, signatures=(), top_n=None):
    """Rank herbs by explicit disease-candidate-target overlap, not global expression correlation."""
    disease_targets, ranked_universe_n, disease_fdr_significant_n, tested_genes = _candidate_disease_targets(report, top_n)
    if not disease_targets:
        return {"state": "no_disease_candidate_targets", "method": "explicit disease-target ∩ herb-target gene-symbol overlap",
                "candidates": [], "disease_targets": [], "limitations": ["No disease target candidates could be constructed from the analyzed cohorts."]}
    signature_index = {}
    for signature in signatures:
        candidate_id = str(signature.get("candidate_id"))
        for item in signature.get("genes", []):
            gene = str(item.get("gene", "")).upper()
            if gene:
                signature_index.setdefault((candidate_id, gene), []).append({
                    "source": signature.get("source_id"), "accession": signature.get("source_accession"),
                    "condition": (signature.get("context") or {}).get("treatment"),
                    "log2FC": item.get("log2FC"), "fdr": item.get("fdr"),
                })
    herb_cache_dir = Path(herb_cache_dir)
    herb_rows = (report.get("external_tcm_database") or {}).get("herb_candidates", [])
    disease_gene_count = len(disease_targets)
    disease_symbols = set(disease_targets)
    universe = {str(x.get("gene", "")).upper() for x in report.get("multi_cohort_cad_proxy", {}).get("meta_analysis", {}).get("results", []) if x.get("gene")}
    universe.update(disease_symbols)
    candidates = []
    missing_detail = []
    for herb in herb_rows:
        entity_id = str(herb.get("entity_id") or "")
        detail_path = herb_cache_dir / f"{entity_id}.json"
        if not detail_path.is_file():
            missing_detail.append(entity_id)
            continue
        try:
            detail = json.loads(detail_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            missing_detail.append(entity_id)
            continue
        target_map = _herb_target_map(detail)
        overlap = sorted(disease_symbols & set(target_map))
        if not overlap:
            continue
        target_details = []
        for symbol in overlap:
            treatment_signatures = signature_index.get((entity_id, symbol), [])
            source_records = disease_targets[symbol].get("sources", [])
            target_details.append({"gene": symbol, "disease_evidence": disease_targets[symbol],
                                   "herb_target_evidence": target_map[symbol]["evidence"],
                                   "disease_source_count": len({x.get("source_database") for x in source_records}),
                                   "disease_source_databases": sorted({x.get("source_database") for x in source_records if x.get("source_database")}),
                                   "exact_disease_source_count": len({x.get("source_database") for x in source_records if _is_direct_scope(x.get("scope"))}),
                                   "measured_perturbation_at_target": treatment_signatures})
        herb_target_genes = set(target_map) & universe
        herb_candidate = {"candidate_id": entity_id, "name": herb.get("Herb cn name") or herb.get("Herb pinyin name") or entity_id,
                          "aliases": [x for x in (herb.get("Herb pinyin name"), herb.get("Herb en name"), herb.get("Herb latin name")) if x and x != "NA"],
                          "disease_association_prior": {"p_value": herb.get("p_value"), "fdr_bh": herb.get("fdr_bh"), "method": herb.get("association_method")},
                          "shared_target_count": len(overlap), "shared_targets": target_details,
                          "multi_database_shared_target_count": sum(1 for x in target_details if x["disease_source_count"] >= 2),
                          "exact_disease_shared_target_count": sum(1 for x in target_details if x["exact_disease_source_count"] >= 1),
                          "supporting_database_count": len({db for x in target_details for db in x["disease_source_databases"]}),
                          "target_set_size_in_disease_universe": len(herb_target_genes), "overlap_p_value": None,
                          "overlap_fdr_bh": None,
                          "evidence_label": "qualified_disease_gene_target_overlap; efficacy not established"}
        candidates.append(herb_candidate)
    candidates.sort(key=lambda x: (-x["exact_disease_shared_target_count"], -x["multi_database_shared_target_count"],
                                   -x["shared_target_count"], x["name"].casefold()))
    source_record = report.get("multi_source_target_evidence") or {}
    used_databases = sorted({source.get("source_database") for gene in disease_targets.values() for source in gene.get("sources", []) if source.get("source_database")})
    exact_supported_gene_count = sum(
        1 for gene in disease_targets.values()
        if any(_is_direct_scope(source.get("scope")) for source in gene.get("sources", []))
    )
    return {"state": "multi_source_proxy_target_screen" if candidates else "no_shared_targets",
            "method": "single workflow: merge retrieved disease-target evidence from independent databases and analyzed omics, intersect with herb target evidence, then rank by exact-phenotype support, multi-database support, and unique target overlap",
            "selection_rule": {"disease_sources_used": used_databases,
                               "disease_target_set": f"all genes supported by source evidence rules; no arbitrary top-N cutoff; {disease_gene_count} distinct genes",
                               "OpenTargets_complete": (source_record.get("OpenTargets") or {}).get("complete"),
                               "GWAS_Catalog_complete": (source_record.get("GWAS_Catalog") or {}).get("complete"),
                               "herb_target_sources": ["HERB herb_target", "HERB drug_paper_target"],
                               "primary_rank": "exact-disease shared targets, then targets supported by multiple disease sources, then total shared targets",
                               "expression_signature_role": "only displayed at shared target genes as supporting perturbation evidence; not used for primary candidate selection"},
            "disease_candidate_target_count": disease_gene_count, "exact_disease_target_count": exact_supported_gene_count,
            "ranked_disease_gene_universe": ranked_universe_n,
            "disease_fdr_significant_gene_count": disease_fdr_significant_n, "tested_gene_universe": tested_genes,
            "disease_targets": list(disease_targets.values()), "candidate_count": len(candidates), "candidates": candidates,
            "details_missing": missing_detail,
            "unscored_formula_count": len((report.get("external_tcm_database") or {}).get("formula_candidates", [])),
                            "limitations": ["Only genes admitted by the high/medium evidence tiers enter this target overlap; proxy and exploratory genes stay separate.",
                            "GWAS Catalog provides curated top associations and does not expose the complete genome-wide summary-statistics archive through its REST API.",
                            "HERB target links can include database-inferred or model-specific literature evidence; inspect each evidence source and experimental context.",
                            "A shared target is a mechanistic hypothesis, not proof of binding, causal mediation, efficacy, dose, or safety.",
                            "Formula target profiles are not yet aggregated from constituent herbs; this result ranks single herbs only."]}


def reverse_lookup_genes(report):
    """Return only high/medium genes from the shared disease evidence tiers."""
    merged, _, _, _ = _candidate_disease_targets(report)
    return sorted(merged)


def _candidate_gene_selection_audit(report):
    from gene_evidence_pipeline import prioritize_gene_evidence

    merged, _, _, _ = _candidate_disease_targets(report)
    prioritization = (report.get("gene_evidence_prioritization")
                      or report.get("gene_prioritization")
                      or prioritize_gene_evidence(report))
    selected = set(reverse_lookup_genes(report))
    source_counts = {}
    for db in ("Open Targets", "GWAS Catalog", "GEO", "GEO multi-cohort Meta"):
        source_counts[db] = len({gene for gene, row in merged.items() if gene in selected
                                 and any(x.get("source_database") == db for x in row.get("sources", []))})
    return {"merged_unique_gene_count": len(merged), "qualified_unique_gene_count": len(selected),
            "qualified_by_source_gene_counts": source_counts,
            "rules": {**(prioritization.get("thresholds") or {}),
                      "proxy_or_nominal": "retained in separate exploratory/proxy lists and excluded from herb/drug target lookup"}}


def screen_from_reverse_target_links(report, reverse_lookup, signatures=(), catalog_totals=None):
    """Screen the broad HERB candidate pool by reversing target-to-herb links.

    HERB disease-associated herb rows are retained as a prior field only; they
    do not determine candidate eligibility.
    """
    all_disease_targets, ranked_n, significant_n, tested_n = _candidate_disease_targets(report)
    selection_audit = _candidate_gene_selection_audit(report)
    genes = reverse_lookup.get("by_gene", {})
    disease_targets = {}
    for gene in reverse_lookup_genes(report):
        if gene in genes and genes[gene].get("complete"):
            disease_targets[gene] = all_disease_targets.get(gene, {
                "gene": gene, "sources": [{"source_database": "GEO multi-cohort Meta", "source_id": "GSE12288;GSE20680;GSE20681",
                                             "scope": "broad_CAD_proxy", "evidence_type": "random_effects_meta_nominal_p_lt_0.1"}]
            })

    disease_priors = {str(x.get("entity_id")): x for x in
                      ((report.get("external_tcm_database") or {}).get("herb_candidates") or [])}
    signatures_by_candidate = {}
    for signature in signatures:
        candidate_id = str(signature.get("candidate_id"))
        for row in signature.get("genes", []):
            signature_gene = str(row.get("gene", "")).upper()
            signatures_by_candidate.setdefault((candidate_id, signature_gene), []).append({
                "source": signature.get("source_id"), "accession": signature.get("source_accession"),
                "condition": (signature.get("context") or {}).get("treatment"),
                "log2FC": row.get("log2FC"), "fdr": row.get("fdr"),
            })

    grouped = {"Herb": {}, "Formula": {}, "Ingredient": {}}
    for gene, target_data in genes.items():
        if gene not in disease_targets:
            continue
        for kind, rows in (("Herb", target_data.get("herb", [])),
                           ("Formula", target_data.get("formula", [])),
                           ("Ingredient", target_data.get("ingredient", []))):
            for edge in rows:
                entity_id = edge.get("entity_id")
                if not entity_id:
                    continue
                candidate = grouped[kind].setdefault(entity_id, {"candidate_id": entity_id, "entity_type": kind, "name": edge.get("name") or entity_id,
                    "aliases": edge.get("aliases", []), "shared_targets": {}})
                gene_item = candidate["shared_targets"].setdefault(gene, {
                    "gene": gene, "disease_evidence": disease_targets[gene], "herb_target_evidence": [],
                    "disease_source_databases": sorted({x.get("source_database") for x in disease_targets[gene].get("sources", []) if x.get("source_database")}),
                    "exact_disease_source_count": sum(1 for x in disease_targets[gene].get("sources", []) if _is_direct_scope(x.get("scope"))),
                    "disease_source_count": len({x.get("source_database") for x in disease_targets[gene].get("sources", []) if x.get("source_database")}),
                    "measured_perturbation_at_target": signatures_by_candidate.get((entity_id, gene), []),
                })
                gene_item["herb_target_evidence"].append(edge)

    ranked = {}
    def disease_pvalue(target):
        value = target.get("disease_evidence", {}).get("pvalue")
        try:
            return float(value)
        except (TypeError, ValueError):
            return float("inf")
    # Precompute disease evidence weights once; do not rescan the entire gene
    # universe for every herb/ingredient/formula candidate.
    disease_gene_weights = {}
    for gene, evidence in disease_targets.items():
        weights = []
        for source in evidence.get("sources", []):
            pvalue = source.get("pvalue")
            try:
                if pvalue is not None:
                    weights.append(min(8.0, max(0.0, -__import__("math").log10(max(float(pvalue), 1e-30)))))
            except (TypeError, ValueError):
                pass
            score = source.get("score")
            if isinstance(score, (int, float)):
                weights.append(1.0 + 5.0 * max(0.0, min(1.0, float(score))))
            if _is_direct_scope(source.get("scope")):
                weights.append(8.0)
        disease_gene_weights[gene] = max(weights, default=0.5) or 0.5
    total_target_weight = sum(disease_gene_weights.values())

    for kind, candidates in grouped.items():
        out = []
        for candidate in candidates.values():
            targets = list(candidate["shared_targets"].values())
            support_weight = 0.0
            for gene, target_info in candidate["shared_targets"].items():
                weight = disease_gene_weights.get(gene, 0.5)
                edge_weight = 0.0
                for edge in target_info["herb_target_evidence"]:
                    if edge.get("relation_source") in {"HERB herb_target", "HERB formula_target", "HERB ingredient_target"}:
                        fdr = edge.get("fdr_bh")
                        edge_weight = max(edge_weight, max(0.25, 1.0 - float(fdr)) if isinstance(fdr, (int, float)) else 0.5)
                    elif edge.get("references"):
                        grade = {"Grade A": 1.0, "Grade B": 0.75, "Grade C": 0.5, "Grade D": 0.25}
                        edge_weight = max(edge_weight, max((grade.get(ref.get("grade"), 0.3) for ref in edge["references"]), default=0.3))
                    else:
                        edge_weight = max(edge_weight, 0.3)
                support_weight += weight * edge_weight
            coverage = support_weight / total_target_weight if total_target_weight else 0.0
            candidate.update({
                "shared_targets": sorted(targets, key=lambda x: (disease_pvalue(x), x["gene"])),
                "shared_target_count": len(targets),
                "multi_database_shared_target_count": sum(1 for x in targets if x["disease_source_count"] >= 2),
                "exact_disease_shared_target_count": sum(1 for x in targets if x["exact_disease_source_count"] >= 1),
                "supporting_database_count": len({db for x in targets for db in x["disease_source_databases"]}),
                "weighted_target_coverage": coverage,
                "overlap_p_value": None, "overlap_fdr_bh": None,
                "direct_disease_shared_target_count": sum(1 for x in targets if x["exact_disease_source_count"] >= 1),
                "proxy_shared_target_count": sum(1 for x in targets if x["exact_disease_source_count"] == 0),
                "evidence_label": "qualified_disease_gene_target_overlap; efficacy not established",
            })
            if kind == "Herb":
                prior = disease_priors.get(candidate["candidate_id"], {})
                candidate["disease_association_prior"] = {"p_value": prior.get("p_value"), "fdr_bh": prior.get("fdr_bh"),
                                                           "method": prior.get("association_method")}
            out.append(candidate)
        out.sort(key=lambda x: (-x["exact_disease_shared_target_count"], -x["multi_database_shared_target_count"],
                                -x["weighted_target_coverage"], -x["shared_target_count"], x["name"].casefold()))
        ranked[kind] = out

    # Disease evidence is already stored once in disease_targets. Candidate
    # rows keep the gene key and herb/ingredient edge evidence, avoiding a
    # potentially enormous repeated copy of the same GWAS/Open Targets records.
    for candidate_rows in ranked.values():
        for candidate in candidate_rows:
            for target in candidate.get("shared_targets", []):
                target.pop("disease_evidence", None)

    lookup_incomplete = bool(reverse_lookup.get("errors")) or reverse_lookup.get("completed_genes", 0) < reverse_lookup.get("requested_genes", 0)
    return {
        "state": ("partial_reverse_target_screen" if lookup_incomplete and any(ranked[k] for k in ranked)
                  else "reverse_target_source_unavailable" if lookup_incomplete
                  else "reverse_target_screen_exploratory" if any(ranked[k] for k in ranked) else "no_reverse_target_overlap"),
        "disease_ontology": (report.get("multi_source_target_evidence") or {}).get("ontology_crosswalk"),
        "method": "source-qualified disease genes from Open Targets, GWAS Catalog, GEO and multi-cohort analyses are merged and reverse-matched to HERB herb/formula target relations; direct and proxy evidence remain labeled separately",
        "selection_rule": {
            "disease_target_source": selection_audit["rules"],
            "herb_candidate_source": "reverse lookup of HERB target records across herbs, formulas and ingredients; direct ingredient-target edges are retained as a separate entity class",
            "rank": "direct phenotype support first, then multi-source disease support, source-weighted target coverage and unique shared-target count",
            "candidate_gene_selection": selection_audit,
            "p_value_test": "not computed; no invalid hypergeometric p-values are reported",
            "expression_signature_role": "supporting evidence shown on shared targets only",
        },
        "candidate_universe": catalog_totals or {},
        "reverse_lookup": {k: reverse_lookup.get(k) for k in ("state", "requested_genes", "completed_genes", "errors")},
        "disease_candidate_target_count": len(disease_targets),
        "genes_qualified_for_reverse_lookup": selection_audit["qualified_unique_gene_count"],
        "genes_successfully_reverse_queried": reverse_lookup.get("completed_genes", 0),
        "genes_not_successfully_reverse_queried": max(0, reverse_lookup.get("requested_genes", 0) - reverse_lookup.get("completed_genes", 0)),
        "direct_evidence_candidate_target_count": sum(1 for x in disease_targets.values() if any(_is_direct_scope(s.get("scope")) for s in x.get("sources", []))),
        "proxy_evidence_candidate_target_count": sum(1 for x in disease_targets.values() if not any(_is_direct_scope(s.get("scope")) for s in x.get("sources", []))),
        "exact_disease_target_count": sum(1 for x in disease_targets.values() if any(_is_direct_scope(s.get("scope")) for s in x.get("sources", []))),
        "ranked_disease_gene_universe": ranked_n, "disease_fdr_significant_gene_count": significant_n,
        "tested_gene_universe": tested_n, "disease_targets": list(disease_targets.values()),
        "candidate_count": len(ranked["Herb"]), "candidates": ranked["Herb"],
        "formula_candidate_count": len(ranked["Formula"]), "formula_candidates": ranked["Formula"],
        "ingredient_candidate_count": len(ranked["Ingredient"]), "ingredient_candidates": ranked["Ingredient"],
        "unscored_formula_count": max(0, int((catalog_totals or {}).get("Formula", {}).get("total", 0)) - len(ranked["Formula"])),
        "limitations": [
            "No exact stable-angina disease target set or FDR-significant CAD meta-analysis gene was available; this is broad-CAD proxy prioritization.",
            *( [f"HERB target API did not complete {max(0, reverse_lookup.get('requested_genes', 0) - reverse_lookup.get('completed_genes', 0))} of {reverse_lookup.get('requested_genes', 0)} eligible gene queries; candidates only reflect successfully retrieved cached/API relations."] if lookup_incomplete else [] ),
            "The disease target threshold is exploratory and must be validated on independent cohorts; no causal claim is made.",
            "HERB target evidence retains its source type and source context; shared-target overlap does not establish efficacy, dose, safety, or direct binding.",
            "HERB disease association rows do not restrict eligibility and are displayed only as a separate prior.",
        ],
    }
