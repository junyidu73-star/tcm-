"""Reliability patches for the recoverable compiled discovery engine.

This module keeps task-level fixes in readable source while the historical
research_engine implementation is being restored from its Python cache.
"""
from __future__ import annotations

import re
import time
import os
import urllib.parse
import urllib.error
from collections import Counter
from pathlib import Path as _Path


_INVALID_SCREENING_REASONS = {
    "模型筛选缺失或无效",
    "筛选引用片段无法在原始摘要核实",
}


class ScreeningValidationError(ValueError):
    pass


def _clean_error(error, config):
    message = str(error)[:500]
    key = str((config or {}).get("api_key") or "")
    if key:
        message = message.replace(key, "[REDACTED]")
    message = re.sub(r"(?i)(api[_-]?key|authorization|bearer)\s*[:=]\s*[^\s,;]+", r"\1=[REDACTED]", message)
    return message


def _transient_error(error):
    code = getattr(error, "code", None)
    if code in {408, 425, 429, 500, 502, 503, 504}:
        return True
    message = str(error).casefold()
    return isinstance(error, (TimeoutError, ConnectionError, OSError, urllib.error.URLError)) or any(
        token in message for token in ("http 408", "http 425", "http 429", "http 500", "http 502", "http 503", "http 504",
                                       "timed out", "timeout", "temporarily unavailable", "connection reset", "name or service not known")
    )


def _install_database_retry_patch(engine_module):
    engine_class = engine_module.ResearchEngine
    if getattr(engine_class.get, "_database_retry_patch", False):
        return
    original_get = engine_class.get

    def get(self, url, payload=None):
        for attempt in range(3):
            try:
                return original_get(self, url, payload)
            except Exception as exc:
                if not _transient_error(exc) or attempt == 2:
                    retries = getattr(self, "job", {}).setdefault("database_retries", [])
                    if attempt:
                        safe_url = urllib.parse.urlsplit(str(url))
                        retries.append({"host": safe_url.hostname or "unknown", "path": safe_url.path[:120],
                                        "attempts": attempt + 1, "error": _clean_error(exc, getattr(self, "config", {}))[:240]})
                    raise
                time.sleep(0.5 * (2 ** attempt))

    get._database_retry_patch = True
    engine_class.get = get
    if not getattr(engine_class.stage, "_activity_heartbeat_patch", False):
        original_stage = engine_class.stage
        def stage(self, text):
            result = original_stage(self, text)
            self.job["last_activity_at"] = time.time()
            self.job["last_activity_kind"] = "stage"
            return result
        stage._activity_heartbeat_patch = True
        engine_class.stage = stage


def _is_screening_request(instruction, data):
    return isinstance(data, dict) and isinstance(data.get("papers"), list) and "筛选" in str(instruction)


def _normalise_quote(text):
    return " ".join(str(text or "").split()).casefold()


def _validate_screening_result(result, papers):
    if not isinstance(result, dict) or not isinstance(result.get("decisions"), list):
        raise ScreeningValidationError("模型筛选结果缺少 decisions 列表")
    expected = {str(row.get("pmid", "")): row for row in papers if row.get("pmid")}
    decisions = {}
    allowed = {"include", "exclude", "uncertain"}
    for row in result["decisions"]:
        if not isinstance(row, dict):
            raise ScreeningValidationError("模型筛选结果包含非对象记录")
        pmid = str(row.get("pmid", ""))
        if not pmid or pmid not in expected or pmid in decisions:
            raise ScreeningValidationError("模型筛选 PMID 缺失、重复或不属于当前批次")
        if row.get("decision") not in allowed:
            raise ScreeningValidationError("模型筛选 decision 值无效")
        abstract = _normalise_quote(expected[pmid].get("abstract", ""))
        quote = _normalise_quote(row.get("quote", ""))
        if abstract and (not quote or quote not in abstract):
            raise ScreeningValidationError("模型筛选引文无法在对应摘要中核实")
        if not isinstance(row.get("reason"), str) or not row["reason"].strip():
            raise ScreeningValidationError("模型筛选缺少理由")
        decisions[pmid] = row
    missing = set(expected) - set(decisions)
    if missing:
        raise ScreeningValidationError("模型筛选缺少 PMID：" + ",".join(sorted(missing)[:12]))
    return result


def _note_failure(engine, papers, error, phase):
    job = getattr(engine, "job", None)
    if not isinstance(job, dict):
        return
    failures = job.setdefault("screening_failures", [])
    failures.append({
        "phase": phase,
        "pmids": [str(row.get("pmid", "")) for row in papers if row.get("pmid")],
        "count": len(papers),
        "error": _clean_error(error, getattr(engine, "config", {})),
    })


def _install_ask_patch(engine_module):
    engine_class = engine_module.ResearchEngine
    if getattr(engine_class.ask, "_reliable_screening_patch", False):
        return
    original_ask = engine_class.ask

    def ask(self, instruction, data, tokens=5000):
        job = getattr(self, "job", {})
        fatal = job.get("model_fatal_error") if isinstance(job, dict) else None
        if fatal:
            raise RuntimeError(fatal)
        if not _is_screening_request(instruction, data):
            for attempt in range(3):
                try:
                    return original_ask(self, instruction, data, tokens)
                except Exception as exc:
                    malformed_json = isinstance(exc, (TypeError, ValueError)) and "json" in str(exc).casefold()
                    if attempt == 2 or (not _transient_error(exc) and not malformed_json):
                        raise RuntimeError(_clean_error(exc, getattr(self, "config", {}))) from None
                    time.sleep(0.5 * (2 ** attempt))

        papers = data["papers"]
        if not papers:
            return {"decisions": [], "herbal_interventions": []}

        def request_batch(rows, depth=0):
            batch_data = {**data, "papers": rows}
            from scientist_workflow import LITERATURE_PROTOCOL
            extra = (
                "\n本批必须为每个输入 PMID 恰好返回一条 decision；不得漏项或重复。"
                "decision 只能是 include、exclude、uncertain；quote 必须是该 PMID 摘要中的原文连续片段。"
                "只输出符合结构的 JSON。"
            ) + "\n研究问题和筛选协议：" + LITERATURE_PROTOCOL
            last_error = None
            max_attempts = 3
            for attempt in range(max_attempts):
                try:
                    result = original_ask(self, instruction + extra, batch_data, tokens)
                    return _validate_screening_result(result, rows)
                except Exception as exc:
                    last_error = exc
                    message = str(exc).lower()
                    fatal_auth = any(token in message for token in ("http 401", "http 402", "http 403", "unauthorized", "insufficient_quota", "余额不足"))
                    if fatal_auth:
                        safe = _clean_error(exc, getattr(self, "config", {}))
                        self.job["model_fatal_error"] = safe
                        _note_failure(self, rows, safe, "provider_auth_or_balance")
                        raise RuntimeError(safe) from None
                    transient = isinstance(exc, (TimeoutError, ConnectionError, OSError, urllib.error.URLError)) or any(token in message for token in ("http 429", "http 500", "http 502", "http 503", "timeout", "timed out", "urlopen", "connection", "temporarily unavailable"))
                    if transient and attempt + 1 < max_attempts:
                        time.sleep(min(2 ** attempt, 4))
                        continue
                    # A truncated or malformed large response often succeeds
                    # when the batch is divided; a transport error does not.
                    if transient:
                        _note_failure(self, rows, last_error, "transient_provider_error")
                        raise RuntimeError(_clean_error(last_error, getattr(self, "config", {}))) from None
                    if isinstance(exc, ScreeningValidationError):
                        break
                    break

            if len(rows) > 1 and depth < 5:
                middle = len(rows) // 2
                left = request_batch(rows[:middle], depth + 1)
                right = request_batch(rows[middle:], depth + 1)
                combined = {"decisions": left["decisions"] + right["decisions"], "herbal_interventions": []}
                for partial in (left, right):
                    combined["herbal_interventions"].extend(partial.get("herbal_interventions", []))
                return combined

            _note_failure(self, rows, last_error or "未知筛选错误", "invalid_model_output")
            raise RuntimeError(_clean_error(last_error or "未知筛选错误", getattr(self, "config", {}))) from None

        return request_batch(papers)

    ask._reliable_screening_patch = True
    engine_class.ask = ask


def _screening_audit(report, failures):
    papers = report.get("papers", []) if isinstance(report, dict) else []
    counts = Counter()
    failed_pmids = set()
    for failure in failures:
        failed_pmids.update(str(pmid) for pmid in failure.get("pmids", []) if pmid)
    record_failures = []
    for paper in papers:
        screening = paper.get("screening") or {}
        reason = str(screening.get("reason", ""))
        decision = screening.get("decision", "uncertain")
        if reason in _INVALID_SCREENING_REASONS or str(paper.get("pmid", "")) in failed_pmids:
            counts["technical_failure"] += 1
            record_failures.append(str(paper.get("pmid", "")))
        else:
            counts[decision if decision in {"include", "exclude", "uncertain"} else "technical_failure"] += 1
    return {
        "total_records": len(papers),
        "valid_decisions": sum(counts[k] for k in ("include", "exclude", "uncertain")),
        "included": counts["include"],
        "excluded": counts["exclude"],
        "uncertain": counts["uncertain"],
        "technical_failures": counts["technical_failure"],
        "failed_pmids": sorted(set(record_failures))[:500],
        "complete": bool(papers) and counts["technical_failure"] == 0 and len(papers) == sum(counts[k] for k in ("include", "exclude", "uncertain")),
    }


def _add_scope_summary(report):
    sources = {"GEO": report.get("datasets", {})}
    sources.update(report.get("databases", {}) or {})
    summary = {}
    for name in ("GEO", "ArrayExpress", "PRIDE"):
        source = sources.get(name, {}) or {}
        tiers = source.get("relevance_summary") or {}
        if not tiers:
            counts = Counter(row.get("relevance", {}).get("tier", "unclassified") for row in source.get("records", []))
            tiers = dict(counts)
        source_ids = {
            str(value)
            for row in source.get("records", [])
            for value in (row.get("accession"), row.get("canonical_study"))
            if value
        }
        summary[name] = {
            "direct_phenotype": int(tiers.get("direct_phenotype", 0)),
            "related_proxy": int(tiers.get("related_proxy", 0)),
            "unrelated_or_unverified": int(tiers.get("unrelated_or_unverified", 0)),
            "analysis_attempts": sum(1 for row in report.get("attempts", []) if str(row.get("accession", "")) in source_ids),
        }
    report["cohort_scope_summary"] = summary
    proxy_total = sum(item["related_proxy"] for item in summary.values())
    if proxy_total:
        report.setdefault("gaps", []).append(
            f"另有 {proxy_total} 条冠心病/相关代理表型研究保留为单独证据；本轮未与目标疾病合并计算，需按样本诊断、组织和对照审核后开展分支分析。"
        )


def _install_dataset_first_analysis_patch(engine_module):
    """Analyze direct GEO cohorts one by one and preserve a status row for every dataset."""
    engine_class = engine_module.ResearchEngine
    if getattr(engine_class.run, "_dataset_first_analysis_patch", False):
        return
    original_run = engine_class.run

    def run(self, disease):
        report = original_run(self, disease)
        if not isinstance(report, dict):
            return report
        from gene_evidence_pipeline import build_dataset_manifest

        analyses = list(report.get("omics_cohorts") or [])
        if isinstance(report.get("omics"), dict):
            primary = report["omics"]
            primary_id = str((primary.get("dataset") or {}).get("accession") or "")
            if not any(str((item.get("dataset") or {}).get("accession") or "") == primary_id for item in analyses):
                analyses.append(primary)
        attempted = {str((item.get("dataset") or {}).get("accession") or "").upper() for item in analyses}
        attempted.update(str(item.get("accession") or "").upper() for item in report.get("attempts", []) or [])
        records = list((report.get("datasets") or {}).get("records", []) or [])
        direct = [r for r in records if isinstance(r, dict) and (r.get("relevance") or {}).get("tier") == "direct_phenotype"]
        # Default to every eligible direct GEO cohort; operators can set an
        # explicit batch ceiling for constrained machines.
        default_limit = max(1, len(direct) + len(analyses))
        try:
            total_limit = max(1, int(os.environ.get("ZHIHENG_MAX_DIRECT_GEO_ANALYSES", str(default_limit))))
        except ValueError:
            total_limit = default_limit
        remaining = max(0, total_limit - len(analyses))
        pending = []
        english = str((report.get("normalized") or {}).get("disease_en") or disease)
        for record in direct:
            accession = str(record.get("accession") or "").upper()
            if not accession or accession in attempted:
                continue
            unsupported = self.unsupported_modality(record)
            if unsupported:
                report.setdefault("attempts", []).append({"accession": accession, "status": "unsupported_modality", "reason": unsupported})
                attempted.add(accession)
                continue
            data_type = str(record.get("data_type") or "").casefold()
            if "expression profiling" not in data_type or not any(x in data_type for x in ("array", "high throughput sequencing", "rna-seq")):
                report.setdefault("attempts", []).append({"accession": accession, "status": "unsupported_data_type", "reason": "数据集元数据未证实是当前 bulk 表达矩阵分析器支持的类型；记录保留待人工核对。"})
                attempted.add(accession)
                continue
            pending.append(record)
        if remaining:
            for record in pending[:remaining]:
                accession = str(record.get("accession") or "")
                self.stage(f"逐队列独立分析：GEO {accession}")
                self.job.setdefault("events", []).append({"time": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "message": f"进入直接疾病候选队列 {accession}；单独审核样本分组、数据尺度和病例/对照设计。"})
                try:
                    cohort = self.analyze_geo(record, english)
                    if not isinstance(cohort, dict):
                        raise ValueError("分析器未返回结构化队列结果")
                    cohort.setdefault("dataset", {}).setdefault("accession", accession)
                    cohort["dataset"].setdefault("evidence_scope", "direct_disease_candidate")
                    cohort.setdefault("design", {})["evidence_scope"] = "direct_disease_candidate"
                    analyses.append(cohort)
                    for pathway in self.enrich(cohort.get("analysis") or []):
                        report.setdefault("pathways", []).append({**pathway, "accession": accession, "evidence_scope": "direct_disease_candidate"})
                    report.setdefault("attempts", []).append({"accession": accession, "status": "independently_analyzed", "evidence_scope": "direct_disease_candidate"})
                except Exception as exc:
                    report.setdefault("attempts", []).append({"accession": accession, "status": "analysis_failed", "evidence_scope": "direct_disease_candidate", "reason": _clean_error(exc, getattr(self, "config", {}))})
                attempted.add(accession.upper())
        for record in pending[remaining:]:
            report.setdefault("attempts", []).append({"accession": record.get("accession"), "status": "deferred_by_run_budget", "reason": f"已保留元数据；本轮直接 GEO 分析上限为 {total_limit} 个，未进入表达矩阵分析。可提高 ZHIHENG_MAX_DIRECT_GEO_ANALYSES 或后续续跑。", "evidence_scope": "direct_disease_candidate"})
        if analyses:
            report["omics_cohorts"] = analyses
            report["omics"] = analyses[0]
        report["dataset_manifest"] = build_dataset_manifest(report)
        report["dataset_analysis_policy"] = {
            "direct_geo_analysis_limit": total_limit,
            "direct_geo_analyzed_this_run": len(analyses),
            "direct_geo_pending": sum(1 for x in report.get("attempts", []) if x.get("status") == "deferred_by_run_budget"),
            "analysis_unit": "one dataset per analysis; no cross-study raw-matrix pooling",
            "non_direct_evidence": "proxy, animal, cell/perturbation and unresolved records remain separately labeled and never enter direct disease meta-analysis",
        }
        return report

    run._dataset_first_analysis_patch = True
    engine_class.run = run


def _install_proxy_analysis_patch(engine_module):
    """Analyze a small, conservative set of related GEO cohorts separately.

    Proxy evidence is deliberately kept out of the direct-disease pathway,
    target, and drug fields. The recovered engine's sample-level design audit
    remains the gate for every matrix.
    """
    engine_class = engine_module.ResearchEngine
    if getattr(engine_class.run, "_proxy_analysis_patch", False):
        return
    original_run = engine_class.run

    def contains_any(text, terms):
        folded = str(text or "").casefold()
        return any(str(term).casefold() in folded for term in terms if str(term).strip())

    def run(self, disease):
        report = original_run(self, disease)
        if not isinstance(report, dict):
            return report
        normalized = report.get("normalized") or {}
        related_terms = [str(term).strip() for term in normalized.get("related_terms", []) if str(term).strip()]
        direct_terms = [str(term).strip() for term in normalized.get("synonyms", []) if str(term).strip()]
        records = (report.get("datasets") or {}).get("records", [])
        proxy_rows = []
        attempted = set()
        excluded_acute = ("acute coronary syndrome", "acute myocardial infarction", "unstable angina", "stemi", "nstemi")
        unsupported = ("mirna", "micro-rna", "microrna", "circrna", "lncrna", "long non-coding", "non-coding rna", "single cell", "single-cell")
        # Bound extra remote matrix/model work per disease task; the report
        # lists every eligible but deferred record so the limit is explicit.
        eligible = []
        for record in records:
            if not isinstance(record, dict) or record.get("relevance", {}).get("tier") != "related_proxy":
                continue
            accession = str(record.get("accession", ""))
            if not accession or accession in attempted:
                continue
            text = " ".join((str(record.get("title", "")), str(record.get("summary", "")), str(record.get("data_type", ""))))
            if contains_any(text, excluded_acute):
                proxy_rows.append({"accession": accession, "status": "excluded_comparison_cohort", "reason": "研究涉及急性冠脉综合征、心梗或不稳定心绞痛，不能作为慢性相关疾病代理队列。", "evidence_scope": "related_proxy"})
                attempted.add(accession)
                continue
            if contains_any(text, unsupported):
                proxy_rows.append({"accession": accession, "status": "unsupported_modality", "reason": "非编码RNA或单细胞数据不适用当前批量mRNA表达分析器。", "evidence_scope": "related_proxy"})
                attempted.add(accession)
                continue
            data_type = str(record.get("data_type", "")).casefold()
            if "expression profiling" not in data_type or not any(x in data_type for x in ("array", "high throughput sequencing", "rna-seq")):
                proxy_rows.append({"accession": accession, "status": "unsupported_data_type", "reason": "元数据未确认是可分析的表达谱数据。", "evidence_scope": "related_proxy"})
                attempted.add(accession)
                continue
            if contains_any(text, direct_terms):
                # A record marked proxy but explicitly describing the target
                # phenotype is not relabelled as a proxy analysis.
                proxy_rows.append({"accession": accession, "status": "scope_conflict", "reason": "队列同时命中目标疾病同义词；需要人工确认分组后再归类。", "evidence_scope": "related_proxy"})
                attempted.add(accession)
                continue
            matched = next((term for term in related_terms if term.casefold() in text.casefold()), None)
            if not matched:
                proxy_rows.append({"accession": accession, "status": "not_specific_enough", "reason": "标题或摘要未明确命中规范化相关疾病词，暂不自动作为代理队列。", "evidence_scope": "related_proxy"})
                attempted.add(accession)
                continue
            eligible.append((record, matched))

        for record, proxy_disease in eligible[:5]:
            accession = str(record.get("accession", ""))
            attempted.add(accession)
            try:
                result = self.analyze_geo(record, proxy_disease)
                if not isinstance(result, dict) or not isinstance(result.get("analysis"), dict):
                    raise ValueError("表达分析未返回可核查的分析对象")
                analysis = result["analysis"]
                pathways = []
                if analysis.get("significant_n", 0) > 0:
                    try:
                        enrichment = self.enrich(analysis)
                        if enrichment:
                            pathways = enrichment if isinstance(enrichment, list) else [enrichment]
                    except Exception as exc:
                        result.setdefault("limits", []).append("代理队列通路富集失败：" + _clean_error(exc, getattr(self, "config", {})))
                result["evidence_scope"] = "related_proxy"
                result["proxy_disease"] = proxy_disease
                result["interpretation"] = "相关疾病代理证据；不能视为目标疾病直接证据，也未用于目标疾病靶点或药物排序。"
                proxy_rows.append({"accession": accession, "status": "analyzed", "proxy_disease": proxy_disease, "result": result, "pathways": pathways, "evidence_scope": "related_proxy"})
            except Exception as exc:
                proxy_rows.append({"accession": accession, "status": "analysis_failed", "reason": _clean_error(exc, getattr(self, "config", {})), "evidence_scope": "related_proxy"})

        deferred = [record for record, _ in eligible[5:]]
        for record in deferred:
            accession = str(record.get("accession", ""))
            attempted.add(accession)
            proxy_rows.append({"accession": accession, "status": "deferred_limit", "reason": "本轮代理队列自动分析上限为5；元数据已保留，可在后续批次分析。", "evidence_scope": "related_proxy"})
        report["proxy_omics_cohorts"] = proxy_rows
        report["proxy_analysis_policy"] = {
            "max_auto_analyses": 5,
            "separate_from_direct_evidence": True,
            "excludes_acute_comparators": True,
            "sample_design_gate": "recovered_engine_analyze_geo",
            "eligible_found": len(eligible),
        }
        report.setdefault("cohort_scope_summary", {})
        summary = report["cohort_scope_summary"].setdefault("GEO", {})
        summary["proxy_analysis_attempts"] = sum(row.get("status") in {"analyzed", "analysis_failed"} for row in proxy_rows)
        summary["proxy_analyzed"] = sum(row.get("status") == "analyzed" for row in proxy_rows)
        return report

    run._proxy_analysis_patch = True
    engine_class.run = run


def _install_herb2_evidence_patch(engine_module):
    """Attach HERB's exact-disease database links as a distinct evidence branch."""
    engine_class = engine_module.ResearchEngine
    if getattr(engine_class.run, "_herb2_evidence_patch", False):
        return
    original_run = engine_class.run

    def run(self, disease):
        report = original_run(self, disease)
        if not isinstance(report, dict):
            return report
        normalized = report.get("normalized") or {}
        try:
            from herb2_connector import Herb2Client, match_herb_disease
            root = _Path(__file__).resolve().parent
            external = match_herb_disease(
                Herb2Client(root / "runs" / "database_data" / "HERB2" / "cache"),
                disease,
                str(normalized.get("disease_en", "")),
                normalized.get("synonyms", []),
            )
            report["external_tcm_database"] = external
            from perturbation_prediction import empty_prediction_state
            omics = report.get("omics") or {}
            reason = "尚未接入符合物种、组织、剂量和时间条件的实测中药扰动签名。"
            report["custom_perturbation_prediction"] = empty_prediction_state(
                reason, (omics.get("dataset") or {}).get("accession"))
            report["literature_constraint_training"] = {
                "state": "awaiting_human_adjudicated_claims",
                "model_fit": False,
                "abstract_count": len(report.get("papers", [])),
                "reason": "摘要检索和模型筛选不等于方向性机制金标准；需要带物种、组织、剂量、时间和实验方法并由研究者裁定的支持/反证记录。",
            }
        except Exception as exc:
            report["external_tcm_database"] = {
                "state": "source_unavailable",
                "error": _clean_error(exc, getattr(self, "config", {})),
                "herb_candidates": [], "formula_candidates": [], "ingredient_candidates": [],
            }
            report.setdefault("gaps", []).append("HERB 2.0疾病关联支路获取失败；未用文献提及名称代替该数据库结果。")
        return report

    run._herb2_evidence_patch = True
    engine_class.run = run


def apply_reliability_patches(engine_module):
    """Patch the recovered module once after loading its compiled code object."""
    _install_ask_patch(engine_module)
    _install_database_retry_patch(engine_module)
    _install_dataset_first_analysis_patch(engine_module)
    _install_proxy_analysis_patch(engine_module)
    _install_herb2_evidence_patch(engine_module)
    if getattr(engine_module.run_research, "_reliability_patch", False):
        return
    original_run_research = engine_module.run_research

    def run_research(job, config, model, network, pubmed):
        from scientist_workflow import LITERATURE_PROTOCOL, review, targeted_literature
        original_model = model
        original_pubmed = pubmed
        job['scientist_model_calls'] = 0
        def model(*args, **kwargs):
            if job.get('cancel_requested') or job['scientist_model_calls'] >= 60:
                raise RuntimeError('科学审查已取消或达到本轮60次模型调用预算。')
            job['scientist_model_calls'] += 1
            return original_model(*args, **kwargs)
        def pubmed(*args, **kwargs):
            if job.get('cancel_requested'):
                raise RuntimeError('任务已取消。')
            return original_pubmed(*args, **kwargs)
        job["literature_protocol"] = LITERATURE_PROTOCOL
        job.setdefault("events", []).append({"message": "已载入疾病机制、靶点因果、药代安全、中药临床等文献筛选协议。"})
        job["state"] = "running"
        job.setdefault("started_at", time.time())
        job.setdefault("last_activity_at", job["started_at"])
        job["stage"] = "科学推理：研究问题与证据缺口"
        initial_review = review("研究规划：疾病名称和适用范围", job.get("disease", ""),
                                [{"disease": job.get("disease", ""), "literature_topics": LITERATURE_PROTOCOL}], config, model)
        try:
            original_run_research(job, config, original_model, network, pubmed)
        except Exception as exc:
            # An uncaught exception in the recovered engine used to terminate
            # the worker while leaving the API job permanently marked running.
            message = _clean_error(exc, config)
            job["state"] = "failed"
            job["stage"] = "研究流程异常退出"
            job["error"] = message
            job["completed_at"] = time.time()
            job.setdefault("events", []).append({
                "time": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
                "message": "研究流程发生未处理异常，任务已标记失败并释放工作线程；错误摘要已脱敏。",
            })
            return
        report = job.get("result")
        if not isinstance(report, dict):
            job["state"] = "failed"
            job["stage"] = "未生成结构化研究报告"
            job["error"] = "研究引擎已结束，但未返回结构化报告。"
            job["completed_at"] = time.time()
            return
        report["gene_filter_criteria"] = dict(job.get("gene_filter_criteria") or {})
        report["scientific_context"] = job.get("scientific_context") or {"disease": job.get("disease", "")}
        report["rule_snapshot"] = job.get("rule_snapshot") or {"id": None, "rules": []}
        report["literature_protocol"] = LITERATURE_PROTOCOL
        report["scientist_reviews"] = {"planning": initial_review}
        job["state"] = "running"
        job["stage"] = "按八类研究问题补充定向文献检索"
        normalized = report.get("normalized") or {}
        report["targeted_literature"] = targeted_literature(normalized.get("disease_en") or job.get("disease", ""), pubmed)
        if any(x["state"] == "failed" for x in report["targeted_literature"]):
            report.setdefault("gaps", []).append("定向文献检索有失败方向，需重试；结果尚未完成纳排筛选。")
        job["stage"] = "整理数据集清单、独立队列结果和基因证据等级"
        job.setdefault("events", []).append({"time": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "message": "数据集检索和单队列分析结束；正在分开汇总直接疾病、代理/实验及遗传证据。"})
        try:
            from gene_evidence_pipeline import build_dataset_manifest, prioritize_gene_evidence
            report["dataset_manifest"] = build_dataset_manifest(report)
            normalized = report.get("normalized") or {}
            from multi_source_target_evidence import retrieve_multisource_target_evidence
            target_cache = _Path(__file__).resolve().parent / "runs" / "database_data" / "target_sources"
            target_terms = list(dict.fromkeys([str(report.get("disease") or ""), str(normalized.get("disease_en") or "")] + list(normalized.get("synonyms") or [])))
            target_terms = [x.strip() for x in target_terms if x.strip()][:5]
            job["stage"] = "查询 Open Targets 与 GWAS Catalog 疾病遗传证据"
            inflight = job.setdefault("in_flight", {})
            inflight["target_databases"] = inflight.get("target_databases", 0) + 1
            try:
                report["multi_source_target_evidence"] = retrieve_multisource_target_evidence(
                    str(report.get("disease") or ""), target_terms[1:], list(normalized.get("related_terms") or [])[:3], target_cache
                )
            finally:
                inflight["target_databases"] = max(0, inflight.get("target_databases", 1) - 1)
                job["last_activity_at"] = time.time()
            report["gene_evidence_prioritization"] = prioritize_gene_evidence(
                report, report["gene_filter_criteria"]
            )
            from decision_rules import apply_gene_rules
            report["baseline_gene_prioritization"] = report["gene_evidence_prioritization"]
            report["gene_evidence_prioritization"] = apply_gene_rules(
                report["baseline_gene_prioritization"], report["rule_snapshot"]["rules"], report["scientific_context"])
            job["stage"] = "科学推理：审查数据集和基因的竞争机制"
            manifest = report.get("dataset_manifest") or {}
            report["scientist_reviews"]["dataset"] = review("数据集适用性与混杂", job.get("disease", ""), manifest.get("records") or [], config, model)
            priority = report["gene_evidence_prioritization"]
            rows = (priority.get("high_priority_genes") or []) + (priority.get("medium_priority_genes") or [])
            report["scientist_reviews"]["gene"] = review("基因筛选：驱动、代偿、损伤及干预方向", job.get("disease", ""), rows, config, model)
        except Exception as exc:
            report["gene_evidence_prioritization"] = {"state": "partial", "eligible_genes": [], "high_priority_genes": [], "medium_priority_genes": [], "exploratory_genes": [], "proxy_genes": [], "error": _clean_error(exc, config)}
            report.setdefault("gaps", []).append("数据集清单/跨队列基因分级或 Open Targets/GWAS Catalog 证据整理失败：" + _clean_error(exc, config))
        job["stage"] = "仅按高/中优先级直接疾病或遗传基因查询药物记录"
        job.setdefault("events", []).append({"time": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "message": f"基因分级完成：{len((report.get('gene_evidence_prioritization') or {}).get('high_priority_genes', []))} 个高优先级、{len((report.get('gene_evidence_prioritization') or {}).get('medium_priority_genes', []))} 个中优先级；探索性和代理基因不进入药物查询。"})
        try:
            from gene_drug_discovery import attach_gene_drug_discovery
            root = _Path(__file__).resolve().parent
            report["gene_drug_discovery"] = attach_gene_drug_discovery(
                report, network, root / "runs" / "database_data" / "DGIdb"
            )
            from decision_rules import apply_drug_rules
            report["gene_drug_discovery"] = apply_drug_rules(report["gene_drug_discovery"],
                report["gene_evidence_prioritization"], report["rule_snapshot"]["rules"], report["scientific_context"])
            job["stage"] = "科学推理：审查药物作用、暴露和反证"
            report["scientist_reviews"]["drug"] = review("药物筛选：作用方向、暴露、安全及反证", job.get("disease", ""), report["gene_drug_discovery"].get("candidates") or [], config, model)
            state = report["gene_drug_discovery"].get("state")
            if state in {"partial", "failed"}:
                report.setdefault("gaps", []).append(
                    "基因-药物查询未能从全部来源完成；请检查 gene_drug_discovery.source_status 中的逐源错误，已有其他数据库结果仍保留。"
                )
            elif state == "no_candidate_genes":
                report.setdefault("gaps", []).append(
                    "本次上游分析未产生可查询的疾病候选基因；未生成药物关联。请先处理疾病映射、队列设计或差异分析缺口。"
                )
        except Exception as exc:
            report["gene_drug_discovery"] = {
                "state": "failed", "method": "按疾病候选基因查询已登记药物-基因关系",
                "candidate_gene_count": 0, "genes_queried": 0, "candidates": [],
                "source_status": {"gene_drug_discovery": {"state": "failed", "error": _clean_error(exc, config)}},
                "limitations": ["该数据库分支失败不会清除已获得的组学、靶点或 ChEMBL 证据。"],
            }
            report.setdefault("gaps", []).append("基因-药物数据库分支失败：" + _clean_error(exc, config))
        failures = job.get("screening_failures", [])
        for name, assessment in report["scientist_reviews"].items():
            if assessment.get("state") in {"failed", "partial"}:
                report.setdefault("gaps", []).append("科学推理环节失败：" + name + "；候选未经过该环节审查。")
        audit = _screening_audit(report, failures)
        report["screening_audit"] = audit
        report["screening_failures"] = failures[:500]
        _add_scope_summary(report)
        if not audit["complete"]:
            job["state"] = "completed_with_gaps"
            report.setdefault("gaps", []).append(
                f"摘要筛选未完整：{audit['valid_decisions']}/{audit['total_records']} 条有有效筛选决定，{audit['technical_failures']} 条发生模型输出或引用校验失败；失败记录已列入 screening_audit.failed_pmids，须重试后才能称筛选完成。"
            )
            if audit["technical_failures"]:
                report["scope"] = str(report.get("scope", "")) + "；摘要筛选存在技术失败，结果为部分完成"
        try:
            from platform_store import persist_report
            root = _Path(__file__).resolve().parent
            store = persist_report(report, root / "runs" / "research_evidence.sqlite3", root / "runs" / "artifacts")
            report["architecture_run"] = {"run_id": store["run_id"], "store": store,
                                           "evidence_storage": "local SQLite development adapter; PostgreSQL schema in postgres_schema.sql",
                                           "artifact_storage": "local content-addressed store unless S3 environment is configured"}
        except Exception as exc:
            report.setdefault("gaps", []).append("证据关系入库失败：" + _clean_error(exc, config))
        job["state"] = "completed_with_gaps" if report.get("gaps") or not audit["complete"] else "completed"
        job["stage"] = "研究报告与基因-药物证据汇总完成"
        job["completed_at"] = time.time()
        job.setdefault("events", []).append({"message": "基因到药物关系已按证据来源汇总；查看逐源状态、候选基因和证据限制。"})
        return

    run_research._reliability_patch = True
    engine_module.run_research = run_research

