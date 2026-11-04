"""Local development evidence store; production PostgreSQL schema is separate."""
from __future__ import annotations

import hashlib
import json
import sqlite3
import uuid
from datetime import datetime, timezone
from pathlib import Path
from artifact_store import configured_object_store


SCHEMA = """
PRAGMA foreign_keys=ON;
CREATE TABLE IF NOT EXISTS entity(entity_id TEXT PRIMARY KEY, entity_type TEXT NOT NULL, namespace TEXT, canonical_name TEXT NOT NULL, species TEXT, metadata_json TEXT NOT NULL DEFAULT '{}');
CREATE TABLE IF NOT EXISTS entity_alias(alias_id INTEGER PRIMARY KEY, entity_id TEXT NOT NULL REFERENCES entity(entity_id), source TEXT NOT NULL, source_id TEXT, raw_name TEXT NOT NULL, mapping_status TEXT NOT NULL, mapping_algorithm TEXT, mapping_version TEXT);
CREATE TABLE IF NOT EXISTS study(study_id TEXT PRIMARY KEY, source_database TEXT NOT NULL, source_accession TEXT NOT NULL, title TEXT, assay TEXT, modality TEXT, species TEXT, tissue TEXT, cell_type TEXT, design_json TEXT NOT NULL DEFAULT '{}', source_version TEXT, source_url TEXT);
CREATE TABLE IF NOT EXISTS sample(sample_id TEXT PRIMARY KEY, study_id TEXT NOT NULL REFERENCES study(study_id), group_label TEXT, diagnosis TEXT, species TEXT, tissue TEXT, cell_type TEXT, dose REAL, dose_unit TEXT, duration REAL, duration_unit TEXT, metadata_json TEXT NOT NULL DEFAULT '{}');
CREATE TABLE IF NOT EXISTS intervention(intervention_id TEXT PRIMARY KEY, intervention_type TEXT NOT NULL, entity_id TEXT REFERENCES entity(entity_id), raw_name TEXT NOT NULL, dose REAL, dose_unit TEXT, duration REAL, duration_unit TEXT, metadata_json TEXT NOT NULL DEFAULT '{}');
CREATE TABLE IF NOT EXISTS evidence_edge(evidence_id TEXT PRIMARY KEY, subject_id TEXT NOT NULL REFERENCES entity(entity_id), predicate TEXT NOT NULL, object_id TEXT NOT NULL REFERENCES entity(entity_id), source_database TEXT NOT NULL, source_id TEXT NOT NULL, evidence_type TEXT NOT NULL, evidence_level TEXT, species TEXT, tissue TEXT, cell_type TEXT, dose REAL, time_value REAL, assay TEXT, direction TEXT, effect_size REAL, p_value REAL, adjusted_p REAL, confidence REAL, uncertainty_json TEXT NOT NULL DEFAULT '{}', version TEXT, payload_json TEXT NOT NULL DEFAULT '{}');
CREATE INDEX IF NOT EXISTS edge_subject_idx ON evidence_edge(subject_id,predicate);
CREATE INDEX IF NOT EXISTS edge_source_idx ON evidence_edge(source_database,source_id);
CREATE TABLE IF NOT EXISTS model_run(run_id TEXT PRIMARY KEY, task TEXT NOT NULL, module TEXT NOT NULL, model_name TEXT, model_version TEXT, status TEXT NOT NULL, input_json TEXT NOT NULL DEFAULT '{}', parameters_json TEXT NOT NULL DEFAULT '{}', metrics_json TEXT NOT NULL DEFAULT '{}', assumptions_json TEXT NOT NULL DEFAULT '{}', started_at TEXT, finished_at TEXT, code_version TEXT);
CREATE TABLE IF NOT EXISTS artifact(artifact_id TEXT PRIMARY KEY, run_id TEXT REFERENCES model_run(run_id), artifact_type TEXT NOT NULL, uri TEXT NOT NULL, sha256 TEXT, source_license TEXT, metadata_json TEXT NOT NULL DEFAULT '{}');
"""


def _dump(value):
    return json.dumps(value, ensure_ascii=False, default=str, separators=(",", ":"))


def _entity(db, entity_id, entity_type, name, namespace=None, species=None, metadata=None):
    db.execute("INSERT OR IGNORE INTO entity VALUES(?,?,?,?,?,?)",
               (entity_id, entity_type, namespace, name or entity_id, species, _dump(metadata or {})))


def persist_report(report: dict, database_path: str | Path, artifact_dir: str | Path):
    database_path, artifact_dir = Path(database_path), Path(artifact_dir)
    database_path.parent.mkdir(parents=True, exist_ok=True)
    artifact_dir.mkdir(parents=True, exist_ok=True)
    run_id = "run_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
    now = datetime.now(timezone.utc).isoformat()
    with sqlite3.connect(database_path) as db:
        db.executescript(SCHEMA)
        db.execute("INSERT INTO model_run VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)", (
            run_id, str(report.get("disease", "")), "database_first_research_pipeline", "modular_evidence_workflow",
            "mvp-2026.10", "completed_with_separate_evidence_layers", _dump(report.get("normalized", {})),
            _dump({"database_analysis": "cached public GEO series matrix", "HERB": "official read-only HERB 2.0 API"}),
            _dump({"differential_genes": report.get("omics", {}).get("analysis", {}).get("significant_n"),
                   "HERB_associations": sum(len((report.get("external_tcm_database") or {}).get(key, [])) for key in ("herb_candidates", "formula_candidates", "ingredient_candidates"))}),
            _dump({"causal_claims": False, "clinical_recommendation": False}), now, now, "source-recovered+local-adapters"))

        omics = report.get("omics") or {}
        dataset = omics.get("dataset") or {}
        accession = dataset.get("accession")
        if accession:
            study_id = f"GEO:{accession}"
            db.execute("INSERT OR IGNORE INTO study VALUES(?,?,?,?,?,?,?,?,?,?,?,?)", (
                study_id, "GEO", accession, dataset.get("title"), "microarray", "bulk expression", "Homo sapiens",
                "Peripheral plasma stimulated coronary artery endothelial cells", "coronary artery endothelial cells",
                _dump(omics.get("design", {})), dataset.get("platform"), f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}"))
            design = omics.get("design") or {}
            for sid, title in zip(design.get("sample_ids", []), design.get("sample_titles", [])):
                sid = str(sid).strip('"')
                db.execute("INSERT OR IGNORE INTO sample VALUES(?,?,?,?,?,?,?,?,?,?,?,?)", (
                    sid, study_id, "ACS" if str(title).startswith("ACS_") else "stable angina", 
                    "ACS" if str(title).startswith("ACS_") else "stable angina", "Homo sapiens", "coronary artery", 
                    "coronary artery endothelial cells", None, None, None, None, _dump({"sample_title": title, "is_precomputed_log_ratio": True})))
            disease_id = "DISEASE:stable-angina"
            _entity(db, disease_id, "Disease", str(report.get("disease", "stable angina")), "platform-normalized")
            for gene in (omics.get("analysis") or {}).get("results", []):
                symbol = str(gene.get("gene", "")).strip()
                if not symbol:
                    continue
                gene_id = "HGNC_SYMBOL:" + symbol
                _entity(db, gene_id, "Gene", symbol, "HGNC_SYMBOL", "Homo sapiens")
                effect = float(gene.get("log2FC", 0))
                db.execute("INSERT OR IGNORE INTO evidence_edge VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (
                    f"{run_id}:deg:{symbol}", disease_id, "has_group_differential_expression", gene_id, "GEO",
                    accession, "measured_database_expression", "E1_exploratory", "Homo sapiens",
                    "coronary artery", "coronary artery endothelial cells", None, None, "GPL6480 microarray",
                    "higher_in_ACS" if effect > 0 else "higher_in_stable_angina", effect, gene.get("pvalue"), gene.get("FDR"),
                    None, _dump({"sample_n": [8, 4], "comparison_scope": design.get("comparison_scope"), "limitations": omics.get("limits", [])}),
                    "GSE54977/GPL6480", _dump(gene)))
            for term in omics.get("pathways", [])[:300]:
                pathway_id = "GO_TERM:" + str(term.get("name", "")).casefold().replace(" ", "_")[:160]
                _entity(db, pathway_id, "Pathway", term.get("name", "GO term"), "GPL6480_GO_process")
                db.execute("INSERT OR IGNORE INTO evidence_edge VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (
                    f"{run_id}:go:{hashlib.sha256(pathway_id.encode()).hexdigest()[:16]}", disease_id, "associated_with", pathway_id,
                    "GEO", accession, "exploratory_annotation_enrichment", "E1_exploratory", "Homo sapiens",
                    "coronary artery", "coronary artery endothelial cells", None, None, "GPL6480 GO:Process ORA",
                    term.get("direction"), None, term.get("pvalue"), term.get("adjusted_p"), None,
                    _dump({"limitations": omics.get("limits", [])}), "GPL6480", _dump(term)))

        external = report.get("external_tcm_database") or {}
        matched = external.get("matched_disease") or {}
        disease_key = "HERB_DISEASE:" + str(matched.get("disease_id", "unmapped"))
        _entity(db, disease_key, "Disease", matched.get("name", "HERB disease"), "HERB", "Homo sapiens")
        for category, key, name_field in (("Herb", "herb_candidates", "Herb cn name"), ("Formula", "formula_candidates", "Formula cn name"), ("Ingredient", "ingredient_candidates", "Ingredient name")):
            for item in external.get(key, []):
                target_id = "HERB_ENTITY:" + str(item.get("entity_id", ""))
                _entity(db, target_id, category, item.get(name_field) or item.get("entity_id"), "HERB", None,
                        {"pinyin_or_alias": item.get("Herb pinyin name") or item.get("Ingredient alias name")})
                db.execute("INSERT OR IGNORE INTO evidence_edge VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (
                    f"{run_id}:herb:{item.get('entity_id')}", disease_key, "statistically_associated_with", target_id,
                    "HERB 2.0", str(item.get("entity_id")), item.get("evidence_type", "HERB_indirect_statistical_disease_association"),
                    "E0_database_inferred_prior", None, None, None, None, None, "HERB2 shared-target inference",
                    "unknown", None, item.get("p_value"), item.get("fdr_bh"), None,
                    _dump({"not_efficacy": True, "method": item.get("association_method"), "scope": matched}),
                    "HERB 2.0", _dump(item)))

        # Persist measured GEO herbal perturbation evidence separately from
        # HERB's shared-target disease prior. Keep only FDR-significant edges
        # in SQLite; complete tested-gene tables remain in the run artifacts.
        geo_analysis = report.get("geo_herbal_perturbation_analysis") or {}
        for signature in geo_analysis.get("signatures", []):
            context = signature.get("context") or {}
            accession = str(signature.get("source_accession") or "GSE244687;GSE244682")
            study_id = "GEO:GSE244682"
            db.execute("INSERT OR IGNORE INTO study VALUES(?,?,?,?,?,?,?,?,?,?,?,?)", (
                study_id, "GEO", "GSE244682", "KORE-Map 1.0 tonifying-herb A549 RNA-seq",
                "RNA-seq", "bulk expression perturbation", context.get("species"), context.get("tissue"),
                context.get("cell_type"), _dump(context), "GSE244687/GSE244682",
                "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE244682"))
            herb_id = "HERB_ENTITY:" + str(signature.get("candidate_id", ""))
            _entity(db, herb_id, "Herb", signature.get("name") or signature.get("candidate_name"), "HERB", "Homo sapiens",
                    {"geo_perturbation_source": accession})
            treatment = str(context.get("treatment") or "unspecified")
            for gene in signature.get("genes", []):
                if gene.get("fdr") is None or float(gene["fdr"]) >= 0.05:
                    continue
                symbol = str(gene.get("gene") or "").strip()
                if not symbol:
                    continue
                gene_id = "HGNC_SYMBOL:" + symbol
                _entity(db, gene_id, "Gene", symbol, "HGNC_SYMBOL", context.get("species"))
                effect = float(gene.get("log2FC", 0))
                evidence_id = f"{run_id}:geo-herb:{signature.get('candidate_id')}:{hashlib.sha256(treatment.encode()).hexdigest()[:10]}:{symbol}"
                db.execute("INSERT OR IGNORE INTO evidence_edge VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (
                    evidence_id, herb_id, "differentially_expressed_after_herbal_extract", gene_id,
                    "GEO", "GSE244682", "measured_herbal_perturbation_differential_expression",
                    "E1_exploratory", context.get("species"), context.get("tissue"), context.get("cell_type"),
                    None, None, "RNA-seq; PyDESeq2; FDR<0.05", "up" if effect > 0 else "down", effect,
                    gene.get("pvalue"), gene.get("fdr"), None,
                    _dump({"dose_and_extract": treatment, "time": context.get("time"),
                           "control": context.get("control"), "source_samples": context.get("source_experiments"),
                           "processing": context.get("processing"), "interpretation": "A549 context; not disease efficacy evidence"}),
                    "GSE244687/GSE244682", _dump(gene)))

        target_screen = report.get("target_based_herb_screen") or {}
        for candidate in target_screen.get("candidates", []):
            herb_id = "HERB_ENTITY:" + str(candidate.get("candidate_id", ""))
            _entity(db, herb_id, "Herb", candidate.get("name"), "HERB", None)
            for target in candidate.get("shared_targets", []):
                symbol = str(target.get("gene") or "").strip()
                disease_evidence = target.get("disease_evidence") or {}
                if not symbol:
                    continue
                gene_id = "HGNC_SYMBOL:" + symbol
                _entity(db, gene_id, "Gene", symbol, "HGNC_SYMBOL", "Homo sapiens")
                edge_id = f"{run_id}:target-overlap:{candidate.get('candidate_id')}:{symbol}"
                db.execute("INSERT OR IGNORE INTO evidence_edge VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (
                    edge_id, herb_id, "shares_prioritized_disease_target", gene_id,
                    "HERB 2.0 + GEO multi-cohort CAD", "GSE12288;GSE20680;GSE20681",
                    "exploratory_herb_disease_target_overlap", "E0_exploratory_proxy", "Homo sapiens",
                    "peripheral blood", None, None, None, "HERB herb_target/drug_paper_target ∩ CAD meta candidate genes",
                    "unknown", disease_evidence.get("log2FC"), disease_evidence.get("pvalue"), disease_evidence.get("FDR"),
                    None, _dump({"overlap_count": candidate.get("shared_target_count"),
                                "overlap_p_value": candidate.get("overlap_p_value"),
                                "overlap_fdr_bh": candidate.get("overlap_fdr_bh"),
                                "disease_gene_rank": disease_evidence.get("rank"),
                                "herb_target_evidence": target.get("herb_target_evidence"),
                                "measured_perturbation_at_target": target.get("measured_perturbation_at_target"),
                                "limitations": target_screen.get("limitations", [])}),
                    "target-screen-v1", _dump(target)))

        # Save a portable run snapshot and hash it for audit/provenance.
        output_file = artifact_dir / f"{run_id}.json"
        output_file.write_text(json.dumps(report, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
        digest = hashlib.sha256(output_file.read_bytes()).hexdigest()
        object_ref = configured_object_store(artifact_dir / "objects").put_file(output_file, output_file.name)
        db.execute("INSERT INTO artifact VALUES(?,?,?,?,?,?,?)", (
            run_id + ":report", run_id, "run_report", object_ref["uri"], digest,
            "GEO public data; HERB 2.0 data attribution/source license should be checked for deployment", _dump({"run_id": run_id})))
    return {"run_id": run_id, "database_path": str(database_path.resolve()), "report_path": str(output_file.resolve()), "report_sha256": digest}


def persist_target_screen(report: dict, database_path: str | Path, artifact_path: str | Path):
    """Persist a refreshed target-first screen without duplicating a large report snapshot."""
    database_path, artifact_path = Path(database_path), Path(artifact_path)
    if not artifact_path.is_file():
        raise FileNotFoundError(artifact_path)
    run_id = "target_screen_" + datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
    now = datetime.now(timezone.utc).isoformat()
    screen = report.get("target_based_herb_screen") or {}
    disease_by_gene = {str(item.get("gene")): item for item in screen.get("disease_targets", []) if item.get("gene")}
    digest = hashlib.sha256(artifact_path.read_bytes()).hexdigest()
    with sqlite3.connect(database_path) as db:
        db.executescript(SCHEMA)
        existing = db.execute(
            "SELECT run_id FROM artifact WHERE artifact_type=? AND sha256=? ORDER BY rowid DESC LIMIT 1",
            ("target_based_herb_screen", digest),
        ).fetchone()
        if existing:
            return {"run_id": existing[0], "database_path": str(database_path.resolve()),
                    "artifact_path": str(artifact_path.resolve()), "artifact_sha256": digest,
                    "reused_existing": True}
        db.execute("INSERT INTO model_run VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?)", (
            run_id, str(report.get("disease", "")), "target_first_herb_screen", "explicit_target_overlap",
            "target-screen-v1", "completed_with_low_confidence" if "low_confidence" in str(screen.get("state")) else "completed",
            _dump({"disease": report.get("disease"), "disease_candidate_target_count": screen.get("disease_candidate_target_count")}),
            _dump(screen.get("selection_rule", {})),
            _dump({"candidate_count": screen.get("candidate_count"), "disease_fdr_significant_gene_count": screen.get("disease_fdr_significant_gene_count")}),
            _dump({"efficacy_claim": False, "causal_claim": False, "global_expression_reversal_used": False}),
            now, now, "target_based_screen.py"))
        for candidate in screen.get("candidates", []):
            herb_id = "HERB_ENTITY:" + str(candidate.get("candidate_id", ""))
            _entity(db, herb_id, "Herb", candidate.get("name"), "HERB", None)
            for target in candidate.get("shared_targets", []):
                symbol = str(target.get("gene") or "").strip()
                if not symbol:
                    continue
                gene_id = "HGNC_SYMBOL:" + symbol
                _entity(db, gene_id, "Gene", symbol, "HGNC_SYMBOL", "Homo sapiens")
                disease = target.get("disease_evidence") or disease_by_gene.get(symbol, {})
                source_records = disease.get("sources", [])
                disease_source_names = sorted({str(x.get("source_database")) for x in source_records if x.get("source_database")})
                source_ids = sorted({str(x.get("source_id")) for x in source_records if x.get("source_id")})
                relation_records = target.get("herb_target_evidence", [])
                herb_sources = sorted({str(x.get("relation_source") or x.get("source")) for x in relation_records
                                       if x.get("relation_source") or x.get("source")})
                source_ids = sorted(set(source_ids) | {"HERB2:" + x for x in herb_sources})
                evidence_id = f"{run_id}:{candidate.get('candidate_id')}:{symbol}"
                db.execute("INSERT OR IGNORE INTO evidence_edge VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (
                    evidence_id, herb_id, "shares_prioritized_disease_target", gene_id,
                    " + ".join(sorted(set(disease_source_names + herb_sources))), ";".join(source_ids),
                    "exploratory_herb_disease_target_overlap", "E0_exploratory_proxy", "Homo sapiens",
                    "peripheral blood", None, None, None, "HERB herb_target/drug_paper_target ∩ CAD meta candidate genes",
                    "unknown", disease.get("log2FC"), disease.get("pvalue"), disease.get("FDR"), None,
                    _dump({"run_id": run_id, "overlap_count": candidate.get("shared_target_count"),
                           "weighted_target_coverage": candidate.get("weighted_target_coverage"),
                           "evidence_records": relation_records,
                           "herb_target_evidence": target.get("herb_target_evidence"),
                           "measured_perturbation_at_target": target.get("measured_perturbation_at_target"),
                           "limitations": screen.get("limitations", [])}),
                    "target-screen-v1", _dump(target)))
        for candidate in screen.get("ingredient_candidates", []):
            ingredient_id = "HERB_INGREDIENT:" + str(candidate.get("candidate_id", ""))
            _entity(db, ingredient_id, "Ingredient", candidate.get("name"), "HERB", None)
            for target in candidate.get("shared_targets", []):
                symbol = str(target.get("gene") or "").strip()
                if not symbol:
                    continue
                gene_id = "HGNC_SYMBOL:" + symbol
                _entity(db, gene_id, "Gene", symbol, "HGNC_SYMBOL", "Homo sapiens")
                disease = target.get("disease_evidence") or disease_by_gene.get(symbol, {})
                disease_sources = sorted({str(x.get("source_database")) for x in disease.get("sources", []) if x.get("source_database")})
                relation_records = target.get("herb_target_evidence", [])
                relation_sources = sorted({str(x.get("relation_source") or x.get("source")) for x in relation_records
                                           if x.get("relation_source") or x.get("source")})
                source_ids = sorted({str(x.get("source_id")) for x in disease.get("sources", []) if x.get("source_id")})
                source_ids = sorted(set(source_ids) | {"HERB2:" + x for x in relation_sources})
                evidence_id = f"{run_id}:INGREDIENT:{candidate.get('candidate_id')}:{symbol}"
                db.execute("INSERT OR IGNORE INTO evidence_edge VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)", (
                    evidence_id, ingredient_id, "ingredient_has_prioritized_disease_target", gene_id,
                    " + ".join(sorted(set(disease_sources + relation_sources))), ";".join(source_ids),
                    "exploratory_ingredient_disease_target_overlap", "E0_exploratory_proxy", "Homo sapiens",
                    None, None, None, None, "HERB ingredient_target/drug_paper_target ∩ multi-source disease candidate genes",
                    "unknown", disease.get("log2FC"), disease.get("pvalue"), disease.get("FDR"), None,
                    _dump({"run_id": run_id, "overlap_count": candidate.get("shared_target_count"),
                           "weighted_target_coverage": candidate.get("weighted_target_coverage"),
                           "evidence_records": relation_records, "limitations": screen.get("limitations", [])}),
                    "target-screen-v1", _dump(target)))
        db.execute("INSERT INTO artifact VALUES(?,?,?,?,?,?,?)", (
            run_id + ":target_screen", run_id, "target_based_herb_screen", str(artifact_path.resolve()), digest,
            "HERB 2.0 source attribution/licensing should be checked for deployment", _dump({"run_id": run_id})))
    return {"run_id": run_id, "database_path": str(database_path.resolve()),
            "artifact_path": str(artifact_path.resolve()), "artifact_sha256": digest}
