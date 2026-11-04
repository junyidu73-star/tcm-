"""Run the database-first MVP on a saved disease-discovery report."""
from __future__ import annotations

import gzip
import json
import urllib.request
from datetime import datetime, timezone
from pathlib import Path

from database_matrix_analysis import analyze_gse54977
from cad_multicohort_analysis import analyze_cad_multicohort
from herb2_connector import Herb2Client, match_herb_disease
from herbal_signature_pipeline import analyze_gse56885, retrieve_candidate_signatures
from platform_store import persist_report
from perturbation_data_discovery import discover_candidate_perturbation_studies
from geo_herbal_perturbation import analyze_koremap_herbal_signatures
from scientific_constraints import build_training_manifest
from target_based_screen import reverse_lookup_genes, screen_from_reverse_target_links
from multi_source_target_evidence import retrieve_multisource_target_evidence, build_database_coverage
from herb_target_reverse_lookup import retrieve_target_links, get_catalog_totals


ROOT = Path(__file__).resolve().parent


def latest_cached_stable_angina_report():
    files = sorted((ROOT / "runs").glob("stable_angina_rerun_*.json.gz"), key=lambda p: p.stat().st_mtime, reverse=True)
    if not files:
        raise FileNotFoundError("No saved stable-angina discovery report is available")
    with gzip.open(files[0], "rt", encoding="utf-8") as handle:
        return json.load(handle), files[0]


def run():
    report, base_report = latest_cached_stable_angina_report()
    omics = analyze_gse54977(ROOT / "runs" / "database_data" / "GSE54977",
                             ROOT / "runs" / "analyses" / "GSE54977")
    disease_proxy = analyze_gse56885(ROOT / "runs" / "database_data" / "GSE56885",
                                     ROOT / "runs" / "analyses" / "GSE56885")
    cad_proxy = analyze_cad_multicohort(ROOT / "runs" / "database_data",
                                        ROOT / "runs" / "database_data" / "annotations",
                                        ROOT / "runs" / "analyses" / "CAD_multi_cohort")
    report["omics"] = omics
    report["database_omics_analyses"] = [omics, disease_proxy]
    report["multi_cohort_cad_proxy"] = cad_proxy
    report["pathways"] = [{"accession": "GSE54977", "direction": "ACS_vs_stable_angina", "terms": omics["pathways"]}]
    report.setdefault("attempts", []).append({
        "accession": "GSE54977", "status": "analyzed_from_public_database_matrix",
        "reason": "GEO Series Matrix and GPL6480 annotation reanalyzed locally; scope and sample design are attached to omics.",
    })

    normalized = report.get("normalized") or {}
    herb = match_herb_disease(
        Herb2Client(ROOT / "runs" / "database_data" / "HERB2" / "cache"),
        report.get("disease", "stable angina"),
        normalized.get("disease_en", "Stable Angina"),
        normalized.get("synonyms", []),
    )
    report["external_tcm_database"] = herb
    herb_client = Herb2Client(ROOT / "runs" / "database_data" / "HERB2" / "cache")
    perturbation_inventory = retrieve_candidate_signatures(
        herb_client, herb.get("herb_candidates", []), ROOT / "runs" / "database_data" / "HERB2")
    geo_herbal = analyze_koremap_herbal_signatures(
        ROOT / "runs" / "database_data" / "GEO_herbal_perturbation",
        ROOT / "runs" / "analyses" / "GEO_herbal_perturbation",
        herb.get("herb_candidates", []))
    all_herb_signatures = perturbation_inventory["signatures"] + geo_herbal["signatures"]
    report["measured_herb_perturbation_signatures"] = perturbation_inventory
    report["geo_herbal_perturbation_analysis"] = geo_herbal
    report["multi_source_target_evidence"] = retrieve_multisource_target_evidence(
        report.get("disease", "stable angina"), normalized.get("synonyms", []),
        normalized.get("related_terms", []), ROOT / "runs" / "database_data" / "target_sources")
    reverse_cache = ROOT / "runs" / "database_data" / "HERB2" / "target_reverse_lookup"
    target_genes = reverse_lookup_genes(report)
    herb_reverse = retrieve_target_links(
        target_genes, reverse_cache,
        progress=lambda i, n, gene, state: print(f"HERB target reverse lookup {i}/{n}: {gene} ({state})", flush=True))
    report["herb_target_reverse_lookup"] = herb_reverse
    catalog_totals = get_catalog_totals(reverse_cache)
    report["database_source_coverage"] = build_database_coverage(report)
    report["database_source_coverage"].append({
        "database": "HERB 2.0 target reverse lookup",
        "status": herb_reverse["state"],
        "records_or_entities": herb_reverse["completed_genes"],
        "detail": f"queried {herb_reverse['requested_genes']} prioritized disease genes; full catalog sizes: "
                  + ", ".join(f"{key} {value['total']}" for key, value in catalog_totals.items()),
    })
    target_screen = screen_from_reverse_target_links(report, herb_reverse, all_herb_signatures, catalog_totals)
    target_screen["source_coverage"] = report["database_source_coverage"]
    report["target_based_herb_screen"] = target_screen
    report["expression_reversal_diagnostic"] = {
        "state": "not_used_for_primary_selection",
        "method": "check measured treatment expression only at the target genes shared by disease and herb",
        "selection_rule": "A herb is shortlisted only when its HERB target genes overlap the disease candidate target set.",
        "note": "Global transcriptome-wide reversal ranking is not used as the primary medicine-selection method.",
    }
    def public_get(url):
        request = urllib.request.Request(url, headers={"User-Agent": "ZhiHeng-Research/1.0"})
        with urllib.request.urlopen(request, timeout=45) as response:
            return response.read()

    perturbation_studies = discover_candidate_perturbation_studies(
        public_get, lambda message: print(str(message).encode("ascii", "backslashreplace").decode("ascii"), flush=True),
        herb.get("herb_candidates", []), herb.get("formula_candidates", []))
    report["herbal_perturbation_study_discovery"] = perturbation_studies
    report.pop("custom_perturbation_prediction", None)
    report["literature_constraint_training"] = {
        "state": "training_corpus_stage_only",
        "search_corpus": {"retrieved_abstracts": len(report.get("papers", [])), "screening_audit": report.get("screening_audit", {})},
        "model_fit": False,
        "reason": "The saved corpus contains machine-screened abstracts, not a human-adjudicated dataset of directional mechanism claims. It can seed a curation queue, but cannot be treated as gold labels for a scientific-constraint model.",
        "separate_from_measured_omics_and_herb_associations": True,
        "required_training_fields": ["entity IDs", "claim relation and direction", "species/tissue/cell", "dose/time", "assay", "source sentence and PMID", "support/refute/uncertain label", "independent curator adjudication"],
    }
    report["architecture_run"] = {
        "run_at_utc": datetime.now(timezone.utc).isoformat(),
        "input_report": str(base_report.resolve()),
        "modules": ["cached literature and multi-database discovery", "GEO matrix differential analysis", "literature-guided multi-cohort CAD proxy analysis", "HERB 2.0 herb/formula/ingredient evidence", "Open Targets complete target pagination", "GWAS Catalog trait and association pagination", "HERB herb-target and literature-target retrieval", "GEO KORE-Map measured herbal perturbation differential analysis", "GEO and BioStudies/ArrayExpress candidate perturbation study discovery", "unified multi-source disease-target overlap screening", "shared-target-specific perturbation evidence", "database-by-database retrieval audit", "literature constraint training readiness", "relational evidence persistence"],
        "architecture": "modular evidence edges; SQLite local development store; PostgreSQL production DDL; filesystem artifact store",
    }
    report.setdefault("gaps", []).extend([
        "本次组学结果来自GSE54977的ACS与稳定型心绞痛血浆诱导内皮反应比较；不作为稳定型心绞痛对健康人的药物逆转签名。",
        "药材候选按CAD候选靶点与HERB药材/文献靶点显式交集筛选；表达扰动结果只用于检查这些共享靶点，不再以全转录组逆转分数另行排序。",
        "新增GSE12288、GSE20680和GSE20681三个人血管造影CAD队列并按研究分别分析、随机效应合并；这三项是广义CAD代理，不能替代稳定型心绞痛特异队列。",
        f"已针对HERB药材和方剂候选另行搜索GEO与BioStudies/ArrayExpress；本轮发现 {perturbation_studies['study_record_count']} 条候选研究元数据。元数据命中尚未自动视为扰动签名，需审核样本分组和实际表达矩阵后才能打分。",
        "本地科学约束训练集尚无人工审核的方向性机制标签，因此本次没有训练约束模型。",
    ])

    output_dir = ROOT / "runs" / "database_first_stable_angina_20261101"
    output_dir.mkdir(parents=True, exist_ok=True)
    constraint_manifest = build_training_manifest([], output_dir / "scientific_constraints_claims.json")
    report["literature_constraint_training"]["training_manifest"] = constraint_manifest
    report["literature_constraint_training"]["claim_file"] = str((output_dir / "scientific_constraints_claims.json").resolve())
    store_manifest = persist_report(report, ROOT / "runs" / "research_evidence.sqlite3", output_dir / "artifacts")
    report["architecture_run"]["store"] = store_manifest
    json_path = output_dir / "report.json"
    json_path.write_text(json.dumps(report, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    with gzip.open(output_dir / "report.json.gz", "wt", encoding="utf-8") as handle:
        json.dump(report, handle, ensure_ascii=False, default=str)

    db = report.get("external_tcm_database") or {}
    discovered_records = perturbation_studies.get("records", [])
    unique_perturbation_studies = len({str(x.get("accession")) for x in discovered_records if x.get("accession")})
    perturbation_triage = {}
    for record in discovered_records:
        label = record.get("sample_design_triage") or "sample_metadata_unavailable"
        perturbation_triage[label] = perturbation_triage.get(label, 0) + 1
    lines = [
        "# 稳定型心绞痛：数据库优先平台运行报告",
        "",
        f"运行时间（UTC）：{report['architecture_run']['run_at_utc']}",
        f"基础检索报告：`{base_report.name}`（复用已完成检索结果；本次重新执行矩阵统计和 HERB 公开数据库查询）",
        "",
        "**疾病命名与证据范围：**稳定型心绞痛是明确的临床疾病概念（MeSH D060050、SNOMED CT 233819005、UMLS C0340288）。Open Targets 的 EFO_0003913 是较宽泛的 angina pectoris，不能区分稳定型与不稳定型。精确靶点为 0 表示本次数据库没有返回该亚型专属靶点集；CAD/宽泛心绞痛结果单独作为代理证据，不能当作直接证据。",
        "",
        "## 模块结果",
        "",
        f"- GEO GSE54977：检验 {omics['analysis']['tested_genes']:,} 个基因，FDR<0.05 且 |log2FC|≥0.5 共 {omics['analysis']['significant_n']} 个。方向为 ACS 减稳定型心绞痛。",
        f"- GEO GSE56885：稳定冠心病 PBMC（n={disease_proxy['design']['case_n']}）对健康对照（n={disease_proxy['design']['control_n']}），检验 {disease_proxy['analysis']['tested_genes']:,} 个基因；当前 Welch/BH 流程得到 FDR<0.05 基因 {disease_proxy['analysis']['significant_n']} 个。该队列仅作稳定 CAD 代理，样本量小。",
        f"- GEO 多队列 CAD 代理：{len(cad_proxy['cohorts'])} 个文献常用血液队列，病例 {sum(x['design']['case_n'] for x in cad_proxy['cohorts'])}、对照 {sum(x['design']['control_n'] for x in cad_proxy['cohorts'])}，另有 GSE20680 中间狭窄样本 {cad_proxy['cohorts'][1]['design']['intermediate_n_excluded']} 例未纳入主比较；研究内差异分析后做随机效应 Meta，共同检验 {cad_proxy['meta_analysis']['tested_genes']:,} 个基因，FDR<0.05 {cad_proxy['meta_analysis']['significant_n']} 个。此分支是广义 CAD 代理，不混入稳定型心绞痛药材排序。",
        f"- HERB 2.0：疾病实体映射 {((db.get('matched_disease') or {}).get('name') or '未匹配')}；药材 {len(db.get('herb_candidates', []))}、方剂 {len(db.get('formula_candidates', []))}、成分 {len(db.get('ingredient_candidates', []))} 条共享靶点统计关联。",
        f"- HERB 2.0 人源药材扰动签名：检查 {perturbation_inventory['candidate_entities_checked']} 个 HERB 药材候选，取得可用合并人源签名 {perturbation_inventory['signatures_found']} 个。",
        f"- GEO KORE-Map：解析 SuperSeries {geo_herbal['super_series']} 的 {geo_herbal['sample_records_in_super_series']} 个样本记录和 {geo_herbal['matrix_columns']} 列表达矩阵；按 A549 Batch A、同批 DMSO 对照、溶剂和剂量分层，完成 {geo_herbal['signature_count']} 个实测药材差异签名，涉及 {geo_herbal['candidate_herbs_with_signatures']} 味药材。RSEM expected counts 四舍五入供 PyDESeq2 探索性分析。",
        f"- 靶点驱动筛选：按至少两个队列、方向一致和 Meta 名义 P<0.1 的规则，反向查询 HERB 全库靶点关联，得到 {target_screen.get('candidate_count', 0)} 味药材和 {target_screen.get('formula_candidate_count', 0)} 个方剂候选；原先 24 味 HERB 疾病关联药材只作为先验注释。CAD Meta 无 FDR 显著基因，结果仍为代理性探索候选。",
        f"- 文献科学约束训练：已有 {len(report.get('papers', [])):,} 篇摘要可作为证据检索语料；人工裁定机制标签 {constraint_manifest['curated_claims']} 条，本轮未训练模型。",
        "",
        "## HERB 2.0 共享靶点关联（数据库先验）",
        "",
        "下列排序采用 HERB 2.0 对目标疾病的共享靶点 Fisher 检验及 BH 校正结果。它们不是本平台新模型分数，也不是药效或临床推荐。",
        "",
        "| 类别 | 候选 | HERB ID | P 值 | FDR-BH | 人源实测扰动签名 |",
        "|---|---|---|---:|---:|---|",
    ]
    measured_ids = {str(x.get("candidate_id")) for x in all_herb_signatures}
    signature_failures = {str(x.get("candidate_id")): x.get("reason", "无") for x in perturbation_inventory["failures"]}
    for category, field, name_key in (("药材", "herb_candidates", "Herb cn name"), ("方剂", "formula_candidates", "Formula cn name"), ("成分", "ingredient_candidates", "Ingredient name")):
        ordered = sorted(db.get(field, []), key=lambda row: row.get("fdr_bh") if row.get("fdr_bh") is not None else 1)
        # List every herb association; only formulas and ingredients stay capped.
        for item in (ordered if category == "药材" else ordered[:12]):
            name = item.get(name_key) or item.get("entity_id")
            alias = item.get("Herb pinyin name") or item.get("Formula pinyin name")
            if alias and alias != "NA":
                name = f"{name}（{alias}）"
            entity_id = str(item.get("entity_id"))
            signature_status = ("已取得；只检查交集靶点变化" if entity_id in measured_ids else
                                "无可用签名：" + signature_failures.get(entity_id, "未取得")) if category == "药材" else "未批量检索"
            pvalue = item.get("p_value")
            fdr = item.get("fdr_bh")
            lines.append(f"| {category} | {name} | {entity_id} | {pvalue:.3g} | {fdr:.3g} | {signature_status} |" if pvalue is not None and fdr is not None else f"| {category} | {name} | {entity_id} | — | — | {signature_status} |")
    lines.extend(["", "## GEO 与 BioStudies/ArrayExpress 中药扰动研究元数据检索", "",
                  f"检索 {perturbation_studies['candidate_count']} 个药材/方剂实体，共执行 {perturbation_studies['query_count']} 个 GEO 与 BioStudies/ArrayExpress 名称查询，发现 {perturbation_studies['study_record_count']} 条记录（去重后 {unique_perturbation_studies} 个研究编号）。样本标签初筛：处理/对照线索 {perturbation_triage.get('candidate_and_control_labels_found_review_required', 0)} 条、设计未确认 {perturbation_triage.get('metadata_hit_design_not_confirmed', 0)} 条、样本设计字段缺失 {perturbation_triage.get('sample_metadata_unavailable', 0)} 条。元数据命中仍需核实物种、组织、对照和矩阵，不能直接当作签名或药效。", "",
                  "| 候选实体 | 类型 | 研究编号 | 来源 | 数据类型/物种 | 候选样本/对照线索 | 样本标题线索 | 研究标题 | 检索词 | 状态 |", "|---|---|---|---|---|---|---|---|---|---|"])
    for item in perturbation_studies["records"]:
        lines.append(f"| {item.get('candidate_name')} | {item.get('candidate_type')} | [{item.get('accession')}]({item.get('url')}) | {', '.join(item.get('sources') or [item.get('source','')])} | {item.get('data_type','')} / {item.get('organism','')} | {item.get('candidate_sample_count','—')} / {item.get('control_like_sample_count','—')} | {'；'.join(item.get('candidate_sample_titles',[])).replace('|',' ')} | {str(item.get('title','')).replace('|',' ')} | {', '.join(item.get('matched_terms',[]))} | {item.get('sample_design_triage','元数据待审核')} |")
    lines.extend(["", "## GEO 中药实测表达签名分析", "",
                  "以下是实际下载表达矩阵后，在同一批次内按药材、溶剂、剂量与 DMSO 对照分别运行 PyDESeq2 的结果；与 HERB 共享靶点关联及元数据检索分开。A549 为肺腺癌细胞，不能当作心血管组织证据。", "",
                  "| 药材 | 提取溶剂 | 剂量 | 处理样本数 | 对照样本数 | 检验基因数 | FDR<0.05 | 分析文件 |", "|---|---|---|---:|---:|---:|---:|---|"])
    for group in geo_herbal["group_analyses"]:
        lines.append(f"| {group['candidate_name']} | {group['solvent']} | {group['dose']} | {group['case_n']} | {group['control_n']} | {group['tested_genes']} | {group['significant_n']} | `{Path(group['output_file']).name}` |")
    lines.extend([
        "",
        "## 文献指导的多队列 CAD 代理分析",
        "",
        "参考已发表的 CAD 外周血转录组研究，纳入 GSE12288、GSE20680、GSE20681。每个队列先按自己的分组和设计独立估计效应，再以 DerSimonian–Laird 随机效应模型合并；没有把不同研究原始矩阵直接拼接。GSE20681 按公开的匹配对做配对检验；GSE20680 主分析只比较严重狭窄组与低/无狭窄组，中间组排除，并额外输出两种中间组归类的敏感性结果。",
        "",
        "| 队列 | 分组 | 病例 | 对照 | 设计 | 映射基因 | FDR<0.05 |",
        "|---|---|---:|---:|---|---:|---:|",
    ])
    for item in cad_proxy["cohorts"]:
        d, a = item["design"], item["analysis"]
        lines.append(f"| {item['dataset']['accession']} | {d['case_group']} vs 对照 | {d['case_n']} | {d['control_n']} | {'配对' if d['paired'] else '非配对'} | {a['tested_genes']:,} | {a['significant_n']} |")
    meta = cad_proxy["meta_analysis"]
    lines.extend([
        f"| 随机效应 Meta | 三队列共同基因 | — | — | 逐队列效应合并 | {meta['tested_genes']:,} | {meta['significant_n']} |",
        "",
        "这些队列扩大了血液 CAD 机制探索的样本量，但临床定义是血管造影 CAD，不是稳定型心绞痛患者对健康人的精确对照。多队列结果可支持广义 CAD 的候选机制与复现分析，不能改名为稳定心绞痛直接证据，也不直接生成稳定心绞痛中药推荐。GSE20680 中间狭窄的两种重分类敏感性结果见结构化报告和分析文件。",
    ])
    lines.extend([
        "",
        "## 靶点驱动的药材与方剂筛选（本轮主排序）",
        "",
        f"疾病侧使用 {target_screen.get('disease_candidate_target_count', 0)} 个跨队列方向一致的探索性 CAD 代理基因（Meta 名义 P<0.1，至少两个队列且方向一致）；FDR 显著基因数为 {target_screen.get('disease_fdr_significant_gene_count', 0)}。药材/方剂从 HERB 全库 target-to-entity 关系反向检索，HERB 疾病关联页不再限制候选。",
        "GEO/HERB 实测扰动签名只在已重合的靶点上展示处理后表达变化，不再用全转录组负相关另排一张药材榜。当前药材候选来自广义 CAD 代理，所有结果仅用于提出机制假设。",
        "",
        "| 排名 | 药材 | 共同靶点 | 交集数 | 多来源支持靶点 | 加权覆盖分数 | HERB疾病先验FDR |",
        "|---:|---|---|---:|---:|---:|---:|",
    ])
    for index, candidate in enumerate(target_screen.get("candidates", []), 1):
        genes = ", ".join(x["gene"] for x in candidate["shared_targets"])
        disease_prior_fdr = (candidate.get("disease_association_prior") or {}).get("fdr_bh")
        disease_prior_text = f"{disease_prior_fdr:.3g}" if disease_prior_fdr is not None else "—"
        lines.append(f"| {index} | {candidate.get('name')} | {genes} | {candidate.get('shared_target_count')} | {candidate.get('multi_database_shared_target_count')} | {candidate.get('weighted_target_coverage', 0):.4f} | {disease_prior_text} |")
    lines.extend(["", "### 方剂候选", "", "| 排名 | 方剂 | 共同靶点 | 交集数 | 多来源支持靶点 | 加权覆盖分数 |", "|---:|---|---|---:|---:|---:|"])
    for index, candidate in enumerate(target_screen.get("formula_candidates", []), 1):
        genes = ", ".join(x["gene"] for x in candidate["shared_targets"])
        lines.append(f"| {index} | {candidate.get('name')} | {genes} | {candidate.get('shared_target_count')} | {candidate.get('multi_database_shared_target_count')} | {candidate.get('weighted_target_coverage', 0):.4f} |")
    lines.extend(["", "### 共同靶点证据明细", "",
                  "| 候选 | 共同靶点 | 疾病来源 | log2FC | 疾病 P / FDR | HERB关系来源 | PMID/分级/方向 | 扰动签名线索 |", "|---|---|---|---:|---|---|---|---|"])
    for candidate_kind, candidate_rows in (("药材", target_screen.get("candidates", [])), ("方剂", target_screen.get("formula_candidates", []))):
      for candidate in candidate_rows:
        for target in candidate.get("shared_targets", []):
            disease = target.get("disease_evidence", {})
            sources = sorted({evidence.get("relation_source", "HERB target relation") for evidence in target.get("herb_target_evidence", [])})
            papers = [f"{ref.get('pmid') or 'NA'} {ref.get('grade') or ''} {ref.get('relationship') or ''}" for evidence in target.get("herb_target_evidence", []) for ref in evidence.get("references", [])]
            perturbation_evidence = "; ".join(f"{x.get('accession')} {x.get('condition')} log2FC={x.get('log2FC')} FDR={x.get('fdr')}" for x in target.get("measured_perturbation_at_target", [])) or "无已分析签名覆盖该靶点"
            effect = disease.get("log2FC")
            effect_text = f"{effect:.4g}" if isinstance(effect, (int, float)) else "—"
            pvalue, fdr = disease.get("pvalue"), disease.get("FDR")
            stats_text = f"{pvalue:.3g} / {fdr:.3g}" if isinstance(pvalue, (int, float)) and isinstance(fdr, (int, float)) else f"{pvalue if pvalue is not None else '—'} / {fdr if fdr is not None else '—'}"
            lines.append(f"| {candidate_kind} {candidate.get('name')} | {target.get('gene')} | {', '.join(target.get('disease_source_databases', []))} | {effect_text} | {stats_text} | {'；'.join(sources)} | {'；'.join(papers) or '无文献关系边'} | {perturbation_evidence} |")
    lines.extend(["", "### 解释边界", "",
        "本轮没有 FDR 显著的疾病 Meta 基因，代理疾病也不能等同稳定型心绞痛。加权覆盖是候选排序指标，不是概率或疗效结论；需逐条核查靶点的来源、物种、组织、疾病模型和调控方向。方剂名单基于 HERB 方剂—靶点关联，未整合各方剂的成分剂量与暴露数据。"])
    lines.extend([
        "",
        "## 可复现产物",
        "",
        f"- 完整报告：`{json_path.name}` / `report.json.gz`",
        f"- 关系证据库：`{Path(store_manifest['database_path']).name}`",
        f"- 证据审计快照：`{Path(store_manifest['report_path']).name}`，SHA-256 `{store_manifest['report_sha256']}`",
        "- PostgreSQL 部署表结构：`postgres_schema.sql`",
        "- GSE54977 与 GSE56885 差异表达结果分别见 `runs/analyses/GSE54977/` 和 `runs/analyses/GSE56885/`。",
        "- GSE56885 代理队列的全部基因排序不代表统计学显著发现；HERB 共享靶点先验与实测扰动签名匹配是两个独立证据支路。",
        "",
        "报告中的 HERB 条目需按公开数据库来源、实体和关联方法解释；不构成治疗建议。",
    ])
    md_path = output_dir / "report.md"
    md_path.write_text("\n".join(lines), encoding="utf-8")
    return {"report_json": str(json_path.resolve()), "report_gzip": str((output_dir / "report.json.gz").resolve()),
            "report_markdown": str(md_path.resolve()), **store_manifest,
            "herb_counts": {key: len(db.get(key, [])) for key in ("herb_candidates", "formula_candidates", "ingredient_candidates")},
            "omics": {"tested_genes": omics["analysis"]["tested_genes"], "significant_genes": omics["analysis"]["significant_n"]}}


if __name__ == "__main__":
    print(json.dumps(run(), ensure_ascii=False, indent=2))
