"""Gene-first drug interaction retrieval with explicit evidence provenance.

This module retrieves database-recorded drug-gene interactions. It does not
infer clinical efficacy or convert a database interaction into a prescription.
"""
from __future__ import annotations

import hashlib
import json
import re
import time
from pathlib import Path
from urllib.error import HTTPError, URLError


DGIDB_GRAPHQL = "https://dgidb.org/api/graphql"
MAX_QUERY_GENES = 40
GENE_RE = re.compile(r"^[A-Z0-9][A-Z0-9._-]{0,19}$")


def _scope_rank(scope: str) -> int:
    value = str(scope or "").casefold()
    if value in {"exact", "direct", "exact_clinical_concept", "direct_phenotype", "direct_disease_candidate"}:
        return 0
    if "proxy" in value or "broad" in value or "comparative" in value:
        return 2
    return 1


def collect_gene_evidence(report: dict) -> dict[str, dict]:
    """Return qualified and proxy genes with their tiers kept explicit.

    Downstream drug lookup must filter to high/medium tiers; proxy evidence is
    retained for review and is never silently promoted to a direct disease gene.
    """
    from gene_evidence_pipeline import prioritize_gene_evidence
    # Reports produced by earlier standalone discovery runs used the shorter
    # ``gene_prioritization`` key. Accept both canonical runtime and archived
    # report shapes so already-computed eligible genes are not silently lost.
    prioritization = (report.get("gene_evidence_prioritization")
                      or report.get("gene_prioritization"))
    raw_evidence_available = bool(report.get("multi_source_target_evidence") or report.get("omics_cohorts")
                                  or report.get("omics") or report.get("target_evidence"))
    if (not isinstance(prioritization, dict)
            or not isinstance(prioritization.get("eligible_genes"), list)
            or (raw_evidence_available and not isinstance(prioritization.get("filters"), dict))):
        prioritization = prioritize_gene_evidence(report)
    genes: dict[str, dict] = {}
    tiered_items = [(item, False) for item in prioritization.get("eligible_genes", [])]
    tiered_items.extend((item, True) for item in prioritization.get("proxy_genes", []))
    for item, is_proxy in tiered_items:
        if not isinstance(item, dict): continue
        symbol = str(item.get("gene") or "").strip().upper()
        if not GENE_RE.fullmatch(symbol): continue
        tier = str(item.get("tier") or ("proxy" if is_proxy else "medium"))
        row = genes.setdefault(symbol, {"gene":symbol,"sources":[],"evidence":[],"tier":tier})
        if tier == "high": row["tier"] = "high"
        evidence_rows = item.get("sources") or item.get("evidence") or []
        row["sources"] = sorted(set(row["sources"] + [str(x.get("source") or "evidence") for x in evidence_rows]))
        for evidence in evidence_rows:
            record = {"source":str(evidence.get("source") or "evidence"),
                      "scope":str(evidence.get("scope") or ("related_proxy" if is_proxy else "direct_phenotype")),
                      "evidence_type":str(evidence.get("evidence_type") or item.get("reason") or "prioritized disease evidence"),
                      "priority_tier":tier,"reason":item.get("reason")}
            record.update({k:v for k,v in evidence.items() if v is not None})
            signature=(record["source"],record["scope"],record["evidence_type"],record.get("study_id"))
            if not any((x.get("source"),x.get("scope"),x.get("evidence_type"),x.get("study_id"))==signature for x in row["evidence"]):
                row["evidence"].append(record)
    return genes


def _priority(row: dict) -> tuple:
    evidence = row.get("evidence", [])
    direct = any(_scope_rank(x.get("scope")) == 0 for x in evidence)
    source_count = len({x.get("source") for x in evidence})
    fdrs = [x.get("FDR") for x in evidence if isinstance(x.get("FDR"), (int, float))]
    pvalues = [x.get("pvalue") for x in evidence if isinstance(x.get("pvalue"), (int, float))]
    return (row.get("tier") != "high", not direct, -source_count, min(fdrs, default=1.0), min(pvalues, default=1.0), row["gene"])


def _unique_values(values: list) -> list:
    """Deduplicate scalar and JSON-like evidence values without hashing dicts."""
    result = []
    seen = set()
    for value in values:
        marker = json.dumps(value, sort_keys=True, ensure_ascii=False, default=str) if isinstance(value, (dict, list)) else str(value)
        if marker not in seen:
            seen.add(marker)
            result.append(value)
    return result


def _error_label(exc: Exception) -> str:
    code = getattr(exc, "code", None)
    if code == 429 or "429" in str(exc):
        return "上游数据库限流（HTTP 429）；稍后可重试。"
    if code in {401, 403} or any(f"{x}" in str(exc) for x in ("401", "403")):
        return f"上游数据库拒绝请求（HTTP {code or '401/403'}）；请检查服务状态或访问策略。"
    if code in {500, 502, 503, 504} or any(f"{x}" in str(exc) for x in ("500", "502", "503", "504")):
        return "上游数据库暂时不可用；其他数据源结果仍保留。"
    if isinstance(exc, (TimeoutError,)) or "timed out" in str(exc).lower() or "timeout" in str(exc).lower():
        return "上游数据库请求超时；已保留其他数据源结果。"
    if isinstance(exc, (URLError, OSError)) or any(x in str(exc).lower() for x in ("dns", "name or service", "connection", "network")):
        return "网络连接失败；请检查服务器出站 HTTPS 和 DNS。"
    if isinstance(exc, json.JSONDecodeError):
        return "上游返回内容无法解析为 JSON。"
    return (str(exc).replace("\n", " ").strip()[:240] or "未分类错误")


def _load_cached(cache_file: Path):
    try:
        obj = json.loads(cache_file.read_text(encoding="utf-8"))
        return obj if isinstance(obj, dict) else None
    except (OSError, json.JSONDecodeError):
        return None


def _retrieve_dgidb_batch(genes: list[str], network, cache_dir: str | Path) -> dict:
    requested = list(dict.fromkeys(str(g).strip().upper() for g in genes if GENE_RE.fullmatch(str(g).strip().upper())))
    requested = requested[:MAX_QUERY_GENES]
    if not requested:
        return {"state": "skipped_no_gene_candidates", "queried_genes": [], "interactions": [], "error": None}
    cache_root = Path(cache_dir)
    cache_root.mkdir(parents=True, exist_ok=True)
    cache_key = hashlib.sha256("\n".join(requested).encode("utf-8")).hexdigest()[:24]
    cache_file = cache_root / f"interactions_{cache_key}.json"
    cached = _load_cached(cache_file)
    if cached is not None:
        cached["cache"] = "hit"
        return cached

    gql = """query DrugInteractionsByGene($genes: [String!]!) {
      genes(names: $genes) { nodes { name longName interactions {
        drug { name conceptId }
        interactionScore
        interactionTypes { type directionality }
        interactionAttributes { name value }
        publications { pmid }
        sources { sourceDbName }
      } } }
    }"""
    payload = json.dumps({"query": gql, "variables": {"genes": requested}}).encode("utf-8")
    headers = {"Content-Type": "application/json", "Accept": "application/json",
               "User-Agent": "ZhiHengResearch/1.0 gene-drug-discovery"}
    errors = []
    response_obj = None
    for attempt in range(3):
        try:
            raw, _ = network(DGIDB_GRAPHQL, payload, headers, timeout=35)
            response_obj = json.loads(raw.decode("utf-8"))
            if response_obj.get("errors"):
                raise RuntimeError("DGIdb GraphQL error: " + json.dumps(response_obj["errors"], ensure_ascii=False)[:400])
            break
        except Exception as exc:
            errors.append(_error_label(exc))
            if attempt < 2 and ("限流" in errors[-1] or "暂时不可用" in errors[-1] or "超时" in errors[-1] or "网络连接失败" in errors[-1]):
                time.sleep(0.5 * (2 ** attempt))
                continue
            return {"state": "failed", "queried_genes": requested, "interactions": [],
                    "error": errors[-1], "attempt_errors": errors}

    nodes = (((response_obj or {}).get("data") or {}).get("genes") or {}).get("nodes") or []
    interactions = []
    returned_symbols = set()
    for node in nodes:
        gene = str(node.get("name") or "").strip().upper()
        if gene not in requested:
            continue
        returned_symbols.add(gene)
        for item in node.get("interactions") or []:
            drug = item.get("drug") or {}
            drug_name = str(drug.get("name") or "").strip()
            if not drug_name:
                continue
            interactions.append({
                "gene": gene, "gene_name": node.get("longName"),
                "drug": drug_name, "drug_id": drug.get("conceptId"),
                "interaction_score": item.get("interactionScore"),
                "interaction_types": item.get("interactionTypes") or [],
                "attributes": item.get("interactionAttributes") or [],
                "pmids": sorted({str(p.get("pmid")) for p in item.get("publications") or [] if p.get("pmid")}),
                "sources": sorted({str(s.get("sourceDbName")) for s in item.get("sources") or [] if s.get("sourceDbName")}),
                "evidence_database": "DGIdb", "evidence_scope": "drug-gene interaction; disease efficacy not established",
                "url": f"https://dgidb.org/results/genes?genes={gene}",
            })
    output = {"state": "completed", "queried_genes": requested,
              "mapped_genes": sorted(returned_symbols), "unmapped_genes": sorted(set(requested) - returned_symbols),
              "interactions": interactions, "retrieved_interaction_count": len(interactions),
              "error": None, "cache": "miss"}
    temporary = cache_file.with_suffix(".tmp")
    temporary.write_text(json.dumps(output, ensure_ascii=False), encoding="utf-8")
    temporary.replace(cache_file)
    return output


def retrieve_dgidb_interactions(genes: list[str], network, cache_dir: str | Path) -> dict:
    """Query every eligible gene in bounded, resumable DGIdb batches."""
    requested = list(dict.fromkeys(str(g).strip().upper() for g in genes
                                   if GENE_RE.fullmatch(str(g).strip().upper())))
    if not requested:
        return {"state": "skipped_no_gene_candidates", "queried_genes": [], "interactions": [],
                "retrieved_interaction_count": 0, "error": None, "cache": "not_applicable",
                "batches": []}
    batches = [requested[i:i + MAX_QUERY_GENES] for i in range(0, len(requested), MAX_QUERY_GENES)]
    results = []
    for batch in batches:
        results.append(_retrieve_dgidb_batch(batch, network, cache_dir))
        if len(results) < len(batches):
            time.sleep(0.1)
    successful = [item for item in results if item.get("state") == "completed"]
    interactions = [row for item in successful for row in item.get("interactions", [])]
    queried = [gene for item in results for gene in item.get("queried_genes", [])]
    mapped = sorted({gene for item in successful for gene in item.get("mapped_genes", [])})
    unmapped = sorted({gene for item in successful for gene in item.get("unmapped_genes", [])})
    failed_batches = [{"queried_genes": item.get("queried_genes", []), "error": item.get("error")}
                      for item in results if item.get("state") != "completed"]
    if not failed_batches:
        state = "completed"
    elif successful:
        state = "partial"
    else:
        state = "failed"
    return {
        "state": state, "queried_genes": queried, "mapped_genes": mapped,
        "unmapped_genes": unmapped, "interactions": interactions,
        "retrieved_interaction_count": len(interactions),
        "error": "; ".join(dict.fromkeys(x["error"] for x in failed_batches if x.get("error"))) or None,
        "cache": "hit" if all(item.get("cache") == "hit" for item in results) else "mixed_or_miss",
        "batch_size": MAX_QUERY_GENES, "batch_count": len(batches),
        "completed_batch_count": len(successful), "failed_batches": failed_batches,
    }


def attach_gene_to_herb_candidates(report: dict, cache_dir: str | Path, retrieve=None) -> dict:
    """Reverse-match qualified disease genes to HERB herbs, formulas and ingredients."""
    from target_based_screen import reverse_lookup_genes, screen_from_reverse_target_links

    eligible = set(reverse_lookup_genes(report))
    if not eligible:
        return {"state": "no_qualified_disease_genes", "genes_queried": [],
                "candidate_count": 0, "candidates": [], "formula_candidate_count": 0,
                "formula_candidates": [], "ingredient_candidate_count": 0,
                "ingredient_candidates": [], "limitations": [
                    "当前疾病证据未达到靶点反查纳入条件；不能从未验证或不显著基因生成中药候选。"]}

    evidence = collect_gene_evidence(report)
    ranked = sorted((evidence[gene] for gene in eligible if gene in evidence), key=_priority)
    genes = [item["gene"] for item in ranked]
    if retrieve is None:
        from herb_target_reverse_lookup import retrieve_target_links
        retrieve = retrieve_target_links
    reverse = retrieve(genes, cache_dir, pause_s=0.08, max_workers=4)
    screen = screen_from_reverse_target_links(report, reverse, signatures=(), catalog_totals={})
    screen["genes_omitted_by_query_limit"] = 0
    screen["genes_queried"] = genes
    screen["candidate_count"] = len(screen.get("candidates") or [])
    screen["formula_candidate_count"] = len(screen.get("formula_candidates") or [])
    screen["ingredient_candidate_count"] = len(screen.get("ingredient_candidates") or [])
    screen.setdefault("limitations", []).append(
        "HERB 共享靶点关系仅用于按基因重合排序，不代表疾病疗效、直接结合、剂量或安全性。")
    if reverse.get("errors"):
        screen["limitations"].append(
            f"HERB 反查有 {len(reverse.get('errors') or [])} 个失败记录；已查询基因与失败项分别保留。")
    return screen


def attach_gene_drug_discovery(report: dict, network, cache_dir: str | Path, herb_retrieve=None) -> dict:
    """Attach transparent gene-prioritized database interactions to a generic report."""
    all_gene_evidence = collect_gene_evidence(report)
    gene_map = {gene: row for gene, row in all_gene_evidence.items() if row.get("tier") in {"high", "medium"}}
    ordered = sorted(gene_map.values(), key=_priority)
    query_genes = [row["gene"] for row in ordered]
    dgidb = (retrieve_dgidb_interactions(query_genes, network, cache_dir) if query_genes else
             {"state": "not_run_no_eligible_genes", "queried_genes": [], "interactions": [],
              "retrieved_interaction_count": 0, "error": None, "cache": "not_applicable"})

    chembl_rows = []
    for item in report.get("drug_candidates", []) or []:
        if not isinstance(item, dict) or not item.get("gene") or not item.get("molecule"):
            continue
        # ChEMBL records may be present in the broad report for unrelated
        # targets. Only attach mechanisms for genes admitted by the same
        # evidence gate used for DGIdb, so the downstream drug set cannot
        # silently expand beyond the prioritized disease genes.
        if str(item.get("gene")).strip().upper() not in gene_map:
            continue
        chembl_rows.append({
            "gene": str(item.get("gene")).upper(), "drug": item.get("molecule"),
            "drug_id": item.get("molecule_chembl_id") or "", "target_id": item.get("target") or "",
            "mechanism": item.get("mechanism"), "action_type": item.get("action_type"),
            "direct_interaction": item.get("direct_interaction"), "pmids": [], "sources": ["ChEMBL"],
            "evidence_database": "ChEMBL", "evidence_scope": "reported drug mechanism; not disease-specific efficacy",
            "url": item.get("url") or "https://www.ebi.ac.uk/chembl/",
        })

    merged: dict[tuple[str, str], dict] = {}
    for item in chembl_rows + (dgidb.get("interactions") or []):
        gene = str(item.get("gene") or "").upper()
        drug = str(item.get("drug") or "").strip()
        key = (gene, str(item.get("drug_id") or drug).casefold())
        row = merged.setdefault(key, {**item, "evidence_databases": [], "interaction_types": [], "sources": [], "pmids": []})
        db = item.get("evidence_database")
        if db and db not in row["evidence_databases"]:
            row["evidence_databases"].append(db)
        for field in ("interaction_types", "sources", "pmids"):
            values = item.get(field) or []
            if not isinstance(values, list):
                values = [values]
            row[field] = _unique_values(row.get(field, []) + values)
        if item.get("mechanism") and not row.get("mechanism"):
            row["mechanism"] = item["mechanism"]
        if item.get("interaction_score") is not None:
            row["interaction_score"] = item["interaction_score"]

    for item in merged.values():
        evidence = gene_map.get(item["gene"], {})
        item["gene_sources"] = evidence.get("sources", [])
        item["gene_evidence"] = evidence.get("evidence", [])
        item["gene_priority_tier"] = evidence.get("tier", "unclassified")
        item["gene_priority"] = _priority(evidence) if evidence else (True, 0, 1.0, 1.0, item["gene"])
        item["evidence_databases"].sort()
        item["sources"].sort()
        item["pmids"].sort()
    candidates = sorted(merged.values(), key=lambda x: (x["gene_priority"], -len(x["evidence_databases"]), x["gene"], str(x["drug"]).casefold()))
    for item in candidates:
        item.pop("gene_priority", None)

    status = ("completed" if dgidb.get("state") == "completed" else
              "partial" if dgidb.get("state") == "partial" or chembl_rows else "failed")
    if not gene_map:
        status = "no_candidate_genes"
    chembl_symbols = sorted({str(item.get("gene") or "").upper() for item in chembl_rows if item.get("gene")})
    chembl_errors = [str(gap)[:240] for gap in report.get("gaps", []) or []
                     if "药物靶点接口" in str(gap) or "ChEMBL" in str(gap)]
    result = {
        "state": status,
        "method": "疾病候选基因 → DGIdb 基因-药物互作 + ChEMBL 靶点机制；按疾病基因证据来源分层展示",
        "candidate_gene_count": len(gene_map), "genes_queried": len(query_genes),
        "genes_not_queried_due_to_limit": 0,
        "dgidb_query_batch_size": MAX_QUERY_GENES,
        "candidate_genes": [{"gene": row["gene"], "tier": row.get("tier", "unclassified"), "sources": row["sources"], "evidence": row["evidence"]} for row in ordered],
        "exploratory_proxy_gene_count": sum(row.get("tier") == "proxy" for row in all_gene_evidence.values()),
        "source_status": {
            "DGIdb": {k: dgidb.get(k) for k in ("state", "queried_genes", "mapped_genes", "unmapped_genes", "retrieved_interaction_count", "error", "cache", "batch_size", "batch_count", "completed_batch_count", "failed_batches")},
            "ChEMBL": {"state": "records_returned" if chembl_rows else ("request_error" if chembl_errors else "no_records_returned_or_not_reported"),
                       "genes_with_returned_records": chembl_symbols,
                       "retrieved_interaction_count": len(chembl_rows), "error": "; ".join(chembl_errors) or None},
        },
        "candidates": candidates,
        "limitations": [
            "基因-药物互作是数据库已收录的关联，不表示该药对当前疾病有效。",
            "疾病基因证据范围、组织/物种和统计阈值随基因保留；代理疾病证据不会标成直接证据。",
            "未验证药物暴露、剂量、安全性或临床结局；需核对来源论文与药品说明。",
        ],
    }
    try:
        result["traditional_medicine"] = attach_gene_to_herb_candidates(
            report, Path(cache_dir).parent / "HERB2" / "target_reverse_lookup", retrieve=herb_retrieve
        )
    except Exception as exc:
        result["traditional_medicine"] = {
            "state": "source_unavailable", "genes_queried": [], "candidate_count": 0,
            "candidates": [], "formula_candidate_count": 0, "formula_candidates": [],
            "ingredient_candidate_count": 0, "ingredient_candidates": [],
            "error": _error_label(exc), "limitations": [
                "HERB 反查失败不会清除 DGIdb 或 ChEMBL 的药物—基因互作结果。"],
        }
    return result
