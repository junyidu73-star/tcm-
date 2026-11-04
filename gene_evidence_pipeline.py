from __future__ import annotations

import math
import re
from collections import defaultdict
from pathlib import Path
from typing import Any


DEFAULT_GENE_FILTERS = {
    "fdr_threshold": 0.05,
    "abs_log2fc_min": 0.5,
    "exploratory_p_threshold": 0.01,
    "gwas_p_threshold": 5e-8,
    "opentargets_score_min": 0.1,
    "meta_direction_concordance_min": 0.75,
    "meta_i2_max": 0.75,
}


def normalize_gene_filter_criteria(value=None):
    """Normalize pre-registered thresholds and fill omitted values from defaults."""
    value = value if isinstance(value, dict) else {}
    output = dict(DEFAULT_GENE_FILTERS)
    for key, default in DEFAULT_GENE_FILTERS.items():
        try:
            output[key] = float(value.get(key, default))
        except (TypeError, ValueError):
            output[key] = default
    return output


def _accession(record: dict) -> str:
    return str(record.get("accession") or record.get("project_id") or record.get("id") or "").strip()


def _canonical_id(value: str) -> str:
    value = str(value or "").strip().upper()
    match = re.search(r"(?:E-GEOD-|GSE)(\d+)", value)
    return "GSE" + match.group(1) if match else value


def _records(report: dict) -> list[dict]:
    out = []
    datasets = report.get("datasets") or {}
    if isinstance(datasets, dict):
        for row in datasets.get("records", []) or []:
            if isinstance(row, dict): out.append(("GEO", row))
    for source, content in (report.get("databases") or {}).items():
        if not isinstance(content, dict): continue
        for row in content.get("records", []) or []:
            if isinstance(row, dict): out.append((str(source), row))
    return out


def build_dataset_manifest(report: dict) -> dict:
    """Keep every discovered study visible, while deduplicating canonical accessions."""
    seen, rows = {}, []
    analysis_by_id = {}
    analyses = list(report.get("omics_cohorts") or []) + list(report.get("database_omics_analyses") or [])
    if isinstance(report.get("omics"), dict): analyses.append(report["omics"])
    for analysis in analyses:
        accession = _accession(analysis.get("dataset") or analysis)
        if accession: analysis_by_id[_canonical_id(accession)] = analysis
    for source, record in _records(report):
        accession = _accession(record)
        if not accession: continue
        canonical = _canonical_id(record.get("canonical_study") or accession)
        relevance = record.get("relevance") or {}
        title = str(record.get("title") or record.get("name") or "")
        summary = str(record.get("summary") or record.get("description") or "")
        kind_text = f"{title} {summary} {record.get('data_type','')} {record.get('organism','')}".casefold()
        if source.upper().startswith("PRIDE") or re.search(r"proteom|mass spect|phosphoprote", kind_text): modality = "proteomics"
        elif source.casefold() == "expression atlas": modality = "expression_atlas_experiment_metadata"
        elif re.search(r"single[ -]?cell|scrna|sc-rna", kind_text): modality = "single_cell_rna"
        elif re.search(r"mirna|micro-?rna|circrna|lncrna|non.?coding rna", kind_text): modality = "noncoding_rna"
        elif re.search(r"rna.?seq|transcript|expression profil|microarray|gene expression", kind_text): modality = "bulk_transcriptomics"
        else: modality = "unknown_or_other"
        if re.search(r"\b(mouse|mice|murine|rat|rats|zebrafish|drosophila)\b", kind_text): evidence_class = "animal_or_model"
        elif re.search(r"in vitro|cell line|cell culture|cultured cells|organoid|treated cells|stimulated cells", kind_text): evidence_class = "cell_or_perturbation"
        elif relevance.get("tier") == "direct_phenotype": evidence_class = "direct_disease_candidate"
        elif relevance.get("tier") == "related_proxy": evidence_class = "related_disease_proxy"
        else: evidence_class = "unresolved_metadata"
        analysis = analysis_by_id.get(canonical)
        if analysis:
            state = "analyzed_independently"
            design = analysis.get("design") or {}
            if design.get("eligible") is False: state = "analysis_rejected_by_design_gate"
        elif modality == "proteomics": state = "metadata_only_proteomics_quant_parser_unavailable"
        elif modality in {"single_cell_rna", "noncoding_rna"}: state = "retained_not_supported_by_current_bulk_mrna_analyzer"
        elif source.upper().startswith("PRIDE"): state = "metadata_only_not_expression_cohort"
        else: state = "discovered_not_yet_analyzed"
        row = {
            "source": source, "accession": accession, "canonical_accession": canonical,
            "title": title, "url": record.get("url"), "relevance_tier": relevance.get("tier", "unclassified"),
            "relevance_basis": relevance.get("basis") or relevance.get("reason"),
            "evidence_class": evidence_class, "modality": modality,
            "species": record.get("organism") or record.get("species"),
            "sample_size": record.get("n_samples") or record.get("sample_count"),
            "platform": record.get("platform") or record.get("technology"),
            "data_type": record.get("data_type"), "analysis_status": state,
            "analysis_scope": "direct" if evidence_class == "direct_disease_candidate" else "separate_or_unresolved",
            "duplicate_of": seen.get(canonical),
            "file_inventory_state": (record.get("file_inventory") or {}).get("status") if isinstance(record.get("file_inventory"), dict) else None,
        }
        if canonical not in seen: seen[canonical] = accession
        rows.append(row)
    summary = defaultdict(int)
    for row in rows:
        summary[row["analysis_status"]] += 1
    return {
        "policy": "search all returned records; deduplicate canonical study IDs; classify evidence before analysis; analyze supported direct human disease cohorts separately; keep proxy, animal, cell/perturbation, and modality-specific evidence separate",
        "record_count": len(rows), "unique_study_count": len({r["canonical_accession"] for r in rows}),
        "status_counts": dict(summary), "records": rows,
        "limitations": ["A metadata record is not proof that a matrix was downloaded or analyzed.",
                        "PRIDE quantification and single-cell/noncoding RNA require modality-specific analyzers; records remain in the manifest with explicit status.",
                        "Public metadata rarely exposes participant-level identifiers; canonical accessions are deduplicated, but cross-study participant overlap may remain unknown."],
    }


def _bh(pvalues: list[float]) -> list[float]:
    n=len(pvalues); order=sorted(range(n), key=lambda i:pvalues[i]); out=[1.0]*n; running=1.0
    for rank in range(n,0,-1):
        i=order[rank-1]; running=min(running,pvalues[i]*n/rank); out[i]=min(1.0,running)
    return out


def random_effects_meta(report: dict, criteria=None) -> dict:
    """DerSimonian-Laird meta-analysis from direct, unique cohort effects and SEs."""
    filters=normalize_gene_filter_criteria(criteria or report.get("gene_filter_criteria"))
    cohorts=[]; seen=set()
    all_analyses=list(report.get("omics_cohorts") or [])
    if isinstance(report.get("omics"),dict): all_analyses.append(report["omics"])
    manifest={r["canonical_accession"]:r for r in (report.get("dataset_manifest") or {}).get("records",[])}
    for analysis in all_analyses:
        ds=analysis.get("dataset") or {}; accession=_accession(ds)
        canonical=_canonical_id(accession)
        if not canonical or canonical in seen: continue
        seen.add(canonical)
        entry=manifest.get(canonical,{})
        design=analysis.get("design") or {}
        if (entry.get("evidence_class") != "direct_disease_candidate" and design.get("phenotype_scope") != "direct") or design.get("eligible") is False: continue
        cohorts.append((canonical, analysis))
    gene_rows=defaultdict(list)
    for accession, analysis in cohorts:
        for row in ((analysis.get("analysis") or {}).get("results") or []):
            if not isinstance(row,dict): continue
            gene=str(row.get("gene") or "").strip().upper()
            effect=row.get("log2FC"); se=row.get("SE") or row.get("se")
            try: effect=float(effect); se=float(se)
            except (TypeError,ValueError): continue
            if gene and math.isfinite(effect) and math.isfinite(se) and se>0:
                gene_rows[gene].append({"study_id":accession,"effect":effect,"SE":se,
                                        "pvalue":row.get("pvalue",row.get("p")),"FDR":row.get("FDR")})
    pooled=[]
    for gene, rows in gene_rows.items():
        if len(rows)<2: continue
        weights=[1/(r["SE"]**2) for r in rows]
        fixed=sum(w*r["effect"] for w,r in zip(weights,rows))/sum(weights)
        q=sum(w*(r["effect"]-fixed)**2 for w,r in zip(weights,rows)); df=len(rows)-1
        c=sum(weights)-sum(w*w for w in weights)/sum(weights)
        tau2=max(0.0,(q-df)/c) if c>0 else 0.0
        random_w=[1/(r["SE"]**2+tau2) for r in rows]
        effect=sum(w*r["effect"] for w,r in zip(random_w,rows))/sum(random_w)
        se=math.sqrt(1/sum(random_w)); z=effect/se
        p=math.erfc(abs(z)/math.sqrt(2))
        i2=max(0.0,(q-df)/q) if q>0 else 0.0
        concordant=sum(1 for r in rows if (r["effect"]>=0)==(effect>=0))
        pooled.append({"gene":gene,"log2FC":effect,"SE":se,"pvalue":p,"study_n":len(rows),
                       "direction_consistency":concordant,"I2":i2,"tau2":tau2,
                       "study_effects":rows,"method":"DerSimonian-Laird random effects; normal-approximation inference"})
    adjusted=_bh([r["pvalue"] for r in pooled]) if pooled else []
    for row,q in zip(pooled,adjusted):
        row["FDR"]=q
        row["high_priority"]=(q<filters["fdr_threshold"] and abs(row["log2FC"])>=filters["abs_log2fc_min"]
                               and row["direction_consistency"]/row["study_n"]>=filters["meta_direction_concordance_min"]
                               and row["I2"]<=filters["meta_i2_max"])
    pooled.sort(key=lambda r:(r["FDR"],-abs(r["log2FC"]),r["gene"]))
    return {"state":"completed" if pooled else "insufficient_independent_direct_cohorts_or_standard_errors",
            "method":"study-wise direct-cohort effects only; DerSimonian-Laird random effects; BH correction across tested genes; no raw-matrix pooling",
            "filters":filters,
            "cohort_count":len(cohorts),"tested_gene_count":len(pooled),
            "high_priority_count":sum(bool(r["high_priority"]) for r in pooled),"results":pooled,
            "limitations":["Proxy/animal/cell perturbation cohorts are excluded from this synthesis.",
                            "Meta-analysis requires at least two independent direct cohorts with standard errors; otherwise no pooled estimate is produced.",
                            "Normal-approximation inference is used for the random-effects pooled estimate; interpret with small cohort counts cautiously."]}


def prioritize_gene_evidence(report: dict, criteria=None) -> dict:
    """Create evidence tiers; exploratory and proxy signals never enter drug queries."""
    filters=normalize_gene_filter_criteria(criteria or report.get("gene_filter_criteria"))
    direct=defaultdict(list); proxy=defaultdict(list); exploration=defaultdict(list); genetic=defaultdict(list)
    manifest={r["canonical_accession"]:r for r in (report.get("dataset_manifest") or {}).get("records",[])}
    analyses=list(report.get("omics_cohorts") or [])
    if isinstance(report.get("omics"),dict): analyses.append(report["omics"])
    seen=set()
    for analysis in analyses:
        ds=analysis.get("dataset") or {}; acc=_accession(ds); canonical=_canonical_id(acc)
        if canonical in seen: continue
        seen.add(canonical); entry=manifest.get(canonical,{})
        design=analysis.get("design") or {}
        direct_design=(design.get("phenotype_scope")=="direct" or design.get("evidence_scope")=="direct_disease_candidate")
        scope=entry.get("evidence_class") or "unresolved_metadata"
        if scope=="unresolved_metadata" and direct_design: scope="direct_disease_candidate"
        for item in ((analysis.get("analysis") or {}).get("results") or []):
            if not isinstance(item,dict): continue
            gene=str(item.get("gene") or "").strip().upper()
            try: fc=float(item.get("log2FC")); p=float(item.get("pvalue",item.get("p"))); fdr=float(item.get("FDR"))
            except (TypeError,ValueError): continue
            if not gene: continue
            evidence={"gene":gene,"study_id":acc,"log2FC":fc,"pvalue":p,"FDR":fdr,
                      "source":"GEO","scope":scope or "unresolved","modality":entry.get("modality")}
            if scope=="direct_disease_candidate":
                if fdr<filters["fdr_threshold"] and abs(fc)>=filters["abs_log2fc_min"]: direct[gene].append(evidence)
                elif p<filters["exploratory_p_threshold"]: exploration[gene].append(evidence)
            else:
                proxy[gene].append(evidence)
    meta=random_effects_meta(report, filters); meta_by={r["gene"]:r for r in meta["results"]}
    candidates={}
    for gene, rows in direct.items():
        meta_row=meta_by.get(gene)
        replicated=bool(meta_row and meta_row["high_priority"])
        tier="high" if replicated else "medium"
        candidates[gene]={"gene":gene,"tier":tier,"reason":"独立直接疾病队列随机效应Meta通过FDR及方向/效应量规则" if replicated else "单个直接疾病队列通过FDR<0.05且|log2FC|≥0.5；等待独立复现",
                         "sources":rows,"meta_analysis":meta_row}
    # A reproducible random-effects signal can pass the predeclared pooled
    # criteria even when no single cohort crosses FDR on its own. Preserve
    # those genes as high priority, using only direct eligible cohorts already
    # admitted by random_effects_meta().
    for gene, meta_row in meta_by.items():
        if not meta_row.get("high_priority"):
            continue
        rows = [{"gene": gene, "study_id": effect.get("study_id"),
                 "log2FC": effect.get("effect"), "SE": effect.get("SE"),
                 "pvalue": effect.get("pvalue"), "FDR": effect.get("FDR"),
                 "source": "GEO", "scope": "direct_disease_candidate",
                 "evidence_type": "independent_direct_cohort_random_effects_meta"}
                for effect in meta_row.get("study_effects", [])]
        if gene not in candidates:
            candidates[gene] = {"gene": gene, "tier": "high",
                                "reason": "独立直接疾病队列随机效应Meta通过FDR、效应量及方向一致性规则",
                                "sources": rows, "meta_analysis": meta_row}
    multisource=report.get("multi_source_target_evidence") or (
        report.get("target_evidence") if isinstance(report.get("target_evidence"), dict) else {}
    )
    for target_set in ((multisource.get("OpenTargets") or {}).get("target_sets") or []):
        scope = target_set.get("scope")
        for item in target_set.get("targets",[]) or []:
            target=item.get("target") or {}; gene=str(target.get("approvedSymbol") or "").upper()
            try: score=float(item.get("score") or 0)
            except (TypeError,ValueError): score=0
            if gene and score>=filters["opentargets_score_min"] and scope == "exact":
                row={"gene":gene,"tier":"medium","reason":f"Open Targets 精确疾病实体关联分数≥{filters['opentargets_score_min']:g}；数据库关联，不等于因果","sources":[{"source":"Open Targets","scope":"exact","score":score,"study_id":target_set.get("disease_id")}],"meta_analysis":None}
                if gene not in candidates: candidates[gene]=row
                else: candidates[gene]["sources"].extend(row["sources"])
            elif gene and score>=filters["opentargets_score_min"]:
                proxy[gene].append({"gene":gene,"source":"Open Targets","scope":"related_proxy",
                                    "evidence_type":"broader disease ontology association",
                                    "score":score,"study_id":target_set.get("disease_id"),
                                    "disease_name":target_set.get("disease_name")})
    for assoc_set in ((multisource.get("GWAS_Catalog") or {}).get("association_sets") or []):
        scope = assoc_set.get("scope")
        for item in assoc_set.get("associations",[]) or []:
            try: p=float(item.get("p_value"))
            except (TypeError,ValueError): continue
            genes=item.get("mapped_genes") or item.get("mapped_gene") or []
            if isinstance(genes,str): genes=[x.strip() for x in genes.split(",") if x.strip()]
            if p>=filters["gwas_p_threshold"]: continue
            for symbol in genes:
                gene=str(symbol).upper()
                if scope != "exact":
                    proxy[gene].append({"gene":gene,"source":"GWAS Catalog","scope":"related_proxy",
                                        "evidence_type":"broader phenotype genome-wide significant locus; mapped gene is not proven causal",
                                        "pvalue":p,"study_id":item.get("accession_id"),
                                        "trait":assoc_set.get("trait",{}).get("efo_trait")})
                    continue
                # Catalog genome-wide significance is locus-level association,
                # not proof that the mapped nearest gene is causal. Keep it
                # eligible as a medium-priority genetic lead pending fine-mapping.
                row={"gene":gene,"tier":"medium","reason":f"GWAS Catalog 精确疾病关联 P<{filters['gwas_p_threshold']:g}；位点映射基因仍需精细定位/共定位确认","sources":[{"source":"GWAS Catalog","scope":"exact","pvalue":p,"study_id":item.get("accession_id"),"causal_gene_validated":False}],"meta_analysis":None}
                if gene not in candidates: candidates[gene]=row
                else: candidates[gene]["sources"].extend(row["sources"])
    for item in report.get("target_evidence", []) or []:
        if not isinstance(item, dict): continue
        source=str(item.get("source_database") or item.get("source") or "").casefold()
        scope=str(item.get("mapping_scope") or item.get("scope") or ("exact" if "open target" in source else "unresolved")).casefold()
        gene=str(item.get("gene") or item.get("symbol") or item.get("target_symbol") or "").strip().upper()
        try: score=float(item.get("association_score", item.get("score", 0)) or 0)
        except (TypeError,ValueError): score=0
        direct_scope = scope in {"exact", "direct", "direct_phenotype", "exact_clinical_concept", "direct_disease_candidate"}
        is_open_targets = "open target" in source
        is_gwas = "gwas" in source
        is_geo = source.strip() == "geo"
        try: pvalue = float(item.get("pvalue", item.get("p_value")))
        except (TypeError, ValueError): pvalue = math.inf
        try: fdr = float(item.get("FDR", item.get("fdr")))
        except (TypeError, ValueError): fdr = math.inf
        try: effect = abs(float(item.get("log2FC", item.get("log2fc"))))
        except (TypeError, ValueError): effect = 0.0
        source_qualifies = ((is_open_targets and score >= filters["opentargets_score_min"])
                            or (is_gwas and pvalue < filters["gwas_p_threshold"])
                            or (is_geo and fdr < filters["fdr_threshold"] and effect >= filters["abs_log2fc_min"]))
        if gene and direct_scope and source_qualifies:
            evidence={**item,"source":str(item.get("source_database") or item.get("source") or "direct disease evidence"),
                      "scope":scope,"score":score if is_open_targets else item.get("score"),
                      "pvalue":pvalue if is_gwas else item.get("pvalue", item.get("p_value")),
                      "study_id":item.get("source_id")}
            if gene not in candidates:
                reason = (f"Open Targets 精确疾病实体关联分数≥{filters['opentargets_score_min']:g}；数据库关联，不等于因果"
                          if is_open_targets else "明确标记的直接疾病证据；仍需核查来源设计与阈值")
                candidates[gene]={"gene":gene,"tier":"medium","reason":reason,"sources":[evidence],"meta_analysis":None}
            else: candidates[gene]["sources"].append(evidence)
    for item in (((report.get("multi_cohort_cad_proxy") or {}).get("meta_analysis") or {}).get("results") or []):
        gene = str(item.get("gene") or "").upper()
        if gene:
            proxy[gene].append({"gene":gene,"study_id":item.get("source_id"),"source":"broad CAD proxy",
                                "scope":"broad_CAD_proxy","evidence_type":"proxy_meta_analysis",**item})
    eligible=sorted(candidates.values(),key=lambda r:(0 if r["tier"]=="high" else 1,r["gene"]))
    return {"state":"ranked" if eligible else "no_qualified_genes",
            "filters":filters,
            "thresholds":{"single_direct_cohort":f"FDR < {filters['fdr_threshold']:g} and |log2FC| >= {filters['abs_log2fc_min']:g}; medium pending replication",
                          "replicated_direct_cohort":f"random-effects meta FDR < {filters['fdr_threshold']:g}, |pooled log2FC| >= {filters['abs_log2fc_min']:g}, direction concordance >= {filters['meta_direction_concordance_min']:g}, I2 <= {filters['meta_i2_max']:g}; high",
                          "GWAS_Catalog":f"exact phenotype scope and P < {filters['gwas_p_threshold']:g}; medium-priority locus evidence pending fine-mapping/colocalization",
                          "OpenTargets":f"exact phenotype scope and score >= {filters['opentargets_score_min']:g}; medium database evidence",
                          "exploratory":f"nominal P < {filters['exploratory_p_threshold']:g} or effect-size-filtered; proxy-only evidence is separate; neither enters drug lookup"},
            "eligible_genes":eligible,"high_priority_genes":[r for r in eligible if r["tier"]=="high"],
            "medium_priority_genes":[r for r in eligible if r["tier"]=="medium"],
            "exploratory_genes":[{"gene":g,"evidence":v} for g,v in sorted(exploration.items())],
            "proxy_genes":[{"gene":g,"evidence":v} for g,v in sorted(proxy.items())],
            "random_effects_meta":meta,
            "drug_lookup_gene_count":len(eligible),
            "limitations":["Proxy, animal, cell/perturbation and nominal-only genes are preserved separately and are not queried as disease drug targets.",
                            "GWAS Catalog mapped genes are not equivalent to fine-mapped causal genes; colocalization and tissue context remain follow-up steps.",
                            "The available Open Targets score is an evidence aggregation score, not a causal probability."]}
