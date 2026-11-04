import gzip
import json
import os
from pathlib import Path

from target_based_screen import reverse_lookup_genes, screen_from_reverse_target_links
from platform_store import persist_target_screen
from multi_source_target_evidence import retrieve_multisource_target_evidence, build_database_coverage
from herb_target_reverse_lookup import retrieve_target_links, get_catalog_totals


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "runs" / "database_first_stable_angina_20261101"
REPORT_GZ = OUT / "report.json.gz"

with gzip.open(REPORT_GZ, "rt", encoding="utf-8") as handle:
    report = json.load(handle)

signatures = (
    (report.get("measured_herb_perturbation_signatures") or {}).get("signatures", [])
    + (report.get("geo_herbal_perturbation_analysis") or {}).get("signatures", [])
)
normalized = report.get("normalized") or {}
report["multi_source_target_evidence"] = retrieve_multisource_target_evidence(
    report.get("disease", "stable angina"), normalized.get("synonyms", []),
    normalized.get("related_terms", []), ROOT / "runs" / "database_data" / "target_sources"
)
herb_reverse_cache = ROOT / "runs" / "database_data" / "HERB2" / "target_reverse_lookup"
target_genes = reverse_lookup_genes(report)
reverse_lookup = retrieve_target_links(
    target_genes, herb_reverse_cache,
    progress=lambda i, n, gene, state: print(f"HERB target reverse lookup {i}/{n}: {gene} ({state})", flush=True) if i == n or i % 50 == 0 else None
)
report["herb_target_reverse_lookup"] = {
    "state": reverse_lookup.get("state"),
    "requested_genes": reverse_lookup.get("requested_genes"),
    "completed_genes": reverse_lookup.get("completed_genes"),
    "errors": reverse_lookup.get("errors", []),
    "per_gene_cache": "runs/database_data/HERB2/target_reverse_lookup/reverse_target_<gene-hash>.json",
    "detail": "Full target edges are retained in per-gene cache files and in matched candidate evidence; duplicate all-gene payload is omitted from the report.",
}
catalog_totals = get_catalog_totals(herb_reverse_cache)
report["database_source_coverage"] = build_database_coverage(report)
report["database_source_coverage"].append({
    "database": "HERB 2.0 target reverse lookup",
    "status": reverse_lookup["state"],
    "records_or_entities": reverse_lookup["completed_genes"],
    "detail": f"queried {reverse_lookup['requested_genes']} prioritized disease genes; full catalog sizes: "
              + ", ".join(f"{key} {value['total']}" for key, value in catalog_totals.items()),
})
screen = screen_from_reverse_target_links(report, reverse_lookup, signatures, catalog_totals)
screen["source_coverage"] = report["database_source_coverage"]
report["target_based_herb_screen"] = screen
report["expression_reversal_diagnostic"] = {
    "state": "not_used_for_primary_selection",
    "method": "show measured perturbation expression only at disease-herb shared targets",
    "selection_rule": "HERB herb/formula/ingredient candidates are discovered by reverse lookup of prioritized disease target genes; HERB disease association rows are annotations only.",
    "note": "No transcriptome-wide reversal score is used as an independent or competing primary ranking.",
}
report.pop("custom_perturbation_prediction", None)
architecture = report.setdefault("architecture_run", {})
architecture["architecture"] = (
    "single target-first herbal-screening pipeline: multi-source disease candidates -> "
    "explicit HERB herb/formula/ingredient target overlap -> target-specific perturbation evidence"
)
architecture["target_screen_refresh"] = (
    "recomputed from saved GEO/CAD/HERB target/signature caches; no new network "
    "discovery or differential-expression rerun"
)

compact = OUT / "artifacts" / "target_based_herb_screen.json"
compact.write_text(
    json.dumps(
        {
            "disease": report.get("disease"),
            "target_based_herb_screen": screen,
            "expression_reversal_diagnostic": report["expression_reversal_diagnostic"],
        },
        ensure_ascii=False,
        indent=2,
        default=str,
    ),
    encoding="utf-8",
)
architecture["target_screen_store"] = persist_target_screen(
    report, ROOT / "runs" / "research_evidence.sqlite3", compact
)

tmp_gz = OUT / "report.json.gz.tmp"
with gzip.open(tmp_gz, "wt", encoding="utf-8") as handle:
    json.dump(report, handle, ensure_ascii=False, default=str)
os.replace(tmp_gz, REPORT_GZ)

# The gzip report is now a recovery copy. Replace the old expanded report to free
# its space before writing the updated expanded report.
report_json = OUT / "report.json"
if report_json.exists():
    report_json.unlink()
report_json.write_text(
    json.dumps(report, ensure_ascii=False, indent=2, default=str), encoding="utf-8"
)

markdown_path = OUT / "report.md"
markdown = markdown_path.read_text(encoding="utf-8")
start_candidates = [
    markdown.find("## 靶点驱动的药材与方剂筛选"),
    markdown.find("## 靶点驱动的药材筛选（唯一主筛选逻辑）"),
    markdown.find("## 靶点驱动的药材筛选（本轮主排序）"),
]
start_candidates = [position for position in start_candidates if position >= 0]
start = min(start_candidates) if start_candidates else -1
end = markdown.find("## 可复现产物", start)
if start < 0 or end < 0:
    raise RuntimeError("Could not locate the target-screen section in report.md")

lines = [
    "## 靶点驱动的药材与方剂筛选",
    "",
    "**疾病命名核查：**稳定型心绞痛是明确的临床疾病概念（MeSH D060050、SNOMED CT 233819005、UMLS C0340288）。Open Targets 使用的 EFO_0003913 是较宽泛的 angina pectoris，不能区分稳定型与不稳定型心绞痛。因此，精确稳定型心绞痛靶点为 0 表示本轮查询的数据库没有返回该亚型专属靶点集，不代表疾病不存在；宽泛心绞痛/CAD 结果只能作为代理证据单独标记。",
    "",
    "本轮将 Open Targets、GWAS Catalog、GEO 单队列和多队列 Meta 的合格基因合并后反查 HERB 2.0 全库靶点—药材/方剂关系。纳入规则按来源分别执行并写入结构化审计：Open Targets score≥0.1、GWAS P≤5×10⁻⁸、GEO FDR<0.05（无 FDR 时 P<0.05），多队列 CAD 代理 Meta 使用 P<0.1、至少 2 个队列且方向一致数≥2（探索性）。同一基因的所有来源证据保留；精确疾病和代理证据分别标记。",
    "",
    f"本轮合并 {((screen.get('selection_rule') or {}).get('candidate_gene_selection') or {}).get('merged_unique_gene_count', 0)} 个来源基因，其中 {screen.get('genes_qualified_for_reverse_lookup', 0)} 个通过来源特异规则，HERB 成功完成 {screen.get('genes_successfully_reverse_queried', 0)} 个查询、{screen.get('genes_not_successfully_reverse_queried', 0)} 个未完成；精确稳定型心绞痛靶点 {screen.get('exact_disease_target_count', 0)} 个，已检索到的代理证据靶点 {screen.get('proxy_evidence_candidate_target_count', 0)} 个。在成功取得的 HERB 关系中得到 {screen.get('candidate_count', 0)} 味药材、{screen.get('formula_candidate_count', 0)} 个方剂、{screen.get('ingredient_candidate_count', 0)} 个成分交集候选。HERB 目录总量见候选库覆盖表。",
    "",
    "### 数据库检索与实际使用范围",
    "",
    "| 数据库 | 状态 | 记录/实体数 | 本轮用途与范围 |",
    "|---|---|---:|---|",
]
for source in screen.get("source_coverage", []):
    lines.append(
        f"| {source.get('database')} | {source.get('status')} | "
        f"{source.get('records_or_entities', '—')} | {source.get('detail', '')} |"
    )
lines.extend([
    "",
    "| 排名 | 证据范围 | 药材 | 共同靶点 | 交集靶点数 | 多来源疾病支持靶点 | 加权靶点覆盖率 | HERB疾病先验FDR |",
    "|---:|---|---|---|---:|---:|---:|---:|",
])
for index, candidate in enumerate(screen.get("candidates", []), 1):
    genes = ", ".join(target["gene"] for target in candidate.get("shared_targets", []))
    prior_fdr = (candidate.get("disease_association_prior") or {}).get("fdr_bh")
    lines.append(
        f"| {index} | {'直接证据' if candidate.get('evidence_label') == 'direct_disease_target_overlap' else '代理证据'} | {candidate.get('name')} | {genes} | "
        f"{candidate.get('shared_target_count')} | {candidate.get('multi_database_shared_target_count')} | "
        f"{candidate.get('weighted_target_coverage', 0):.4f} | {prior_fdr if prior_fdr is not None else '—'} |"
    )
lines.extend(["", "### 方剂候选", "", "| 排名 | 方剂 | 共同靶点 | 交集靶点数 | 多来源疾病支持靶点 | 加权靶点覆盖率 |", "|---:|---|---|---:|---:|---:|"])
for index, candidate in enumerate(screen.get("formula_candidates", []), 1):
    genes = ", ".join(target["gene"] for target in candidate.get("shared_targets", []))
    lines.append(f"| {index} | {candidate.get('name')} | {genes} | {candidate.get('shared_target_count')} | {candidate.get('multi_database_shared_target_count')} | {candidate.get('weighted_target_coverage', 0):.4f} |")
lines.extend(["", "### 成分靶点候选", "", "成分结果单独展示。成分—靶点重合不能直接替代药材/方剂候选，也不能视为药效证据。", "", "| 排名 | 成分 | HERB ID | 共同靶点 | 交集靶点数 |", "|---:|---|---|---|---:|"])
for index, candidate in enumerate(screen.get("ingredient_candidates", []), 1):
    genes = ", ".join(target["gene"] for target in candidate.get("shared_targets", []))
    lines.append(f"| {index} | {candidate.get('name')} | {candidate.get('candidate_id')} | {genes} | {candidate.get('shared_target_count')} |")
lines.extend(
    [
        "",
        "### 共同靶点证据和表达线索",
        "",
        "| 候选 | 靶点 | 疾病来源 | 疾病 P/FDR | HERB靶点关系来源 | 文献PMID/分级 | 该靶点扰动表达 |",
        "|---|---|---:|---|---|---|",
    ]
)
disease_targets_by_gene = {str(item.get("gene")): item for item in screen.get("disease_targets", []) if item.get("gene")}
for candidate_kind, candidate_rows in (("药材", screen.get("candidates", [])), ("方剂", screen.get("formula_candidates", [])), ("成分", screen.get("ingredient_candidates", []))):
    for candidate in candidate_rows:
      for target in candidate.get("shared_targets", []):
        disease = target.get("disease_evidence") or disease_targets_by_gene.get(str(target.get("gene", "")), {})
        sources = sorted({evidence.get("relation_source", "HERB target relation") for evidence in target.get("herb_target_evidence", [])})
        papers = [f"{ref.get('pmid') or 'NA'} {ref.get('grade') or ''} {ref.get('relationship') or ''}" for evidence in target.get("herb_target_evidence", []) for ref in evidence.get("references", [])]
        expression = "; ".join(
            f"{item.get('accession')} {item.get('condition')} "
            f"log2FC={item.get('log2FC')} FDR={item.get('fdr')}"
            for item in target.get("measured_perturbation_at_target", [])
        ) or "暂无覆盖"
        lines.append(
            f"| {candidate_kind} {candidate.get('name')} | {target.get('gene')} | "
            f"{', '.join(target.get('disease_source_databases', []))} | {disease.get('pvalue')} / {disease.get('FDR')} | "
            f"{'; '.join(sources)} | {'; '.join(papers) or '无直接文献边'} | {expression} |"
        )
lines.extend(
    [
        "",
        "### 结果边界",
        "",
        "HERB疾病关联不再充当候选入口。HERB herb_target 与 paper_herb_target 关系仍需按物种、组织、模型、剂量和方向复核；加权覆盖分数是排序指标，不是概率或疗效估计。当前疾病侧没有 FDR 显著基因，且CAD代理不等同稳定型心绞痛。方剂共享靶点仅按HERB方剂—靶点关系计算，尚未验证成分剂量、方剂制备和暴露。",
        "",
    ]
)
markdown_path.write_text(
    markdown[:start] + "\n".join(lines) + "\n" + markdown[end:], encoding="utf-8"
)

print(
    json.dumps(
        {
            "state": screen.get("state"),
            "candidate_count": screen.get("candidate_count"),
            "formula_candidate_count": screen.get("formula_candidate_count"),
            "ingredient_candidate_count": screen.get("ingredient_candidate_count"),
            "disease_candidate_targets": screen.get("disease_candidate_target_count"),
            "disease_fdr_significant": screen.get("disease_fdr_significant_gene_count"),
            "top": [
                (item.get("name"), item.get("shared_target_count"), [t.get("gene") for t in item.get("shared_targets", [])])
                for item in screen.get("candidates", [])[:10]
            ],
            "report_json_bytes": report_json.stat().st_size,
            "compact_artifact": str(compact.resolve()),
        },
        ensure_ascii=False,
        indent=2,
    )
)
