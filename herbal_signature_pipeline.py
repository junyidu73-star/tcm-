"""Disease-to-herbal perturbation signature matching with explicit study context."""
from __future__ import annotations

import csv
import gzip
import json
import re
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, ttest_ind

from database_matrix_analysis import _bh, _read_series


def _geo_fields(line: str) -> list[str]:
    return next(csv.reader([line.rstrip("\n")], delimiter="\t", quotechar='"'))


def _symbol_map_from_platform(platform_file: Path, cache_file: Path) -> dict[str, str]:
    if cache_file.is_file():
        frame = pd.read_csv(cache_file, sep="\t", dtype=str)
        return dict(zip(frame["probe"], frame["gene"]))
    symbols = {}
    if platform_file.is_file():
        handle = platform_file.open("rt", encoding="utf-8", errors="replace", newline="")
    else:
        url = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL15207&targ=self&form=text&view=data"
        handle = urllib.request.urlopen(urllib.request.Request(url, headers={"User-Agent": "ZhiHeng-Research/1.0"}), timeout=90)
    with handle:
        for line in handle:
            if isinstance(line, bytes):
                line = line.decode("utf-8", "replace")
            if line.startswith("ID\t"):
                header = _geo_fields(line)
                probe_i = header.index("ID")
                gene_i = header.index("Gene Symbol")
                break
        else:
            raise ValueError("GPL15207 annotation table is missing its ID/Gene Symbol header")
        reader = csv.reader((line.decode("utf-8", "replace") if isinstance(line, bytes) else line for line in handle), delimiter="\t", quotechar='"')
        for row in reader:
            if len(row) <= max(probe_i, gene_i):
                continue
            probe, gene = row[probe_i].strip(), row[gene_i].strip()
            if probe and gene and gene not in {"---", "NA", "null"} and not re.search(r"///|;|,", gene):
                symbols.setdefault(probe, gene)
    pd.DataFrame(symbols.items(), columns=["probe", "gene"]).to_csv(cache_file, sep="\t", index=False)
    return symbols


def analyze_gse56885(data_dir: str | Path, output_dir: str | Path):
    """Analyze stable-CAD versus healthy PBMCs; small, proxy evidence for stable angina."""
    data_dir, output_dir = Path(data_dir), Path(output_dir)
    matrix_path = data_dir / "GSE56885_series_matrix.txt.gz"
    platform_path = data_dir / "GPL15207_data.soft"
    if not matrix_path.is_file():
        url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE56nnn/GSE56885/matrix/GSE56885_series_matrix.txt.gz"
        request = urllib.request.Request(url, headers={"User-Agent": "ZhiHeng-Research/1.0"})
        with urllib.request.urlopen(request, timeout=90) as response, matrix_path.open("wb") as output:
            output.write(response.read())

    metadata, probes = _read_series(matrix_path)
    sample_ids = metadata.get("!Sample_geo_accession", [])
    titles = metadata.get("!Sample_title", [])
    characteristics = []
    with gzip.open(matrix_path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_characteristics_ch1"):
                characteristics.append(_geo_fields(line)[1:])
            elif line.startswith("!series_matrix_table_begin"):
                break
    diagnosis = next((row for row in characteristics if row and all("diagnosis:" in v.casefold() for v in row)), None)
    cell_type = next((row for row in characteristics if row and all("cell type:" in v.casefold() for v in row)), None)
    if not diagnosis or len(diagnosis) != probes.shape[1] or not cell_type or len(cell_type) != probes.shape[1]:
        raise ValueError("GSE56885 diagnosis/cell-type sample metadata could not be aligned")
    case_idx = [i for i, value in enumerate(diagnosis) if "coronary artery disease" in value.casefold()]
    control_idx = [i for i, value in enumerate(diagnosis) if "healthy" in value.casefold()]
    if len(case_idx) != 4 or len(control_idx) != 2 or any("pbmc" not in cell_type[i].casefold() for i in case_idx + control_idx):
        raise ValueError(f"Unexpected GSE56885 design: CAD={len(case_idx)}, healthy={len(control_idx)}")

    symbol_map = _symbol_map_from_platform(platform_path, data_dir / "GPL15207_probe_gene.tsv")
    mapped = probes.loc[probes.index.intersection(symbol_map)].copy()
    mapped.index = [symbol_map[p] for p in mapped.index]
    genes = mapped.groupby(level=0).median()
    case_values, control_values = genes.iloc[:, case_idx].to_numpy(), genes.iloc[:, control_idx].to_numpy()
    tests = ttest_ind(case_values, control_values, axis=1, equal_var=False, nan_policy="omit")
    result = pd.DataFrame({
        "gene": genes.index.astype(str),
        "log2FC": np.nanmean(case_values, axis=1) - np.nanmean(control_values, axis=1),
        "pvalue": tests.pvalue,
        "n_case": np.sum(np.isfinite(case_values), axis=1),
        "n_control": np.sum(np.isfinite(control_values), axis=1),
    }).dropna(subset=["log2FC", "pvalue"])
    result["FDR"] = _bh(result["pvalue"].to_numpy())
    result["significant"] = result["FDR"].lt(0.05)
    result = result.sort_values(["FDR", "pvalue", "gene"]).reset_index(drop=True)
    output_dir.mkdir(parents=True, exist_ok=True)
    out = output_dir / "GSE56885_differential_expression.tsv"
    result.to_csv(out, sep="\t", index=False)
    return {
        "dataset": {"accession": "GSE56885", "platform": "GPL15207", "tissue": "PBMC", "species": "Homo sapiens"},
        "design": {"case_group": "stable coronary artery disease", "case_n": len(case_idx),
                   "control_group": "healthy", "control_n": len(control_idx),
                   "sample_ids": sample_ids, "sample_titles": titles,
                   "comparison_scope": "stable_CAD_vs_healthy_PBMC_proxy_for_stable_angina",
                   "paired": False},
        "analysis": {"method": "Welch t-test on GEO processed array values; BH correction; probe symbols collapsed by median",
                     "tested_genes": int(len(result)), "significant_n": int(result["significant"].sum()),
                     "results": result.to_dict(orient="records")},
        "limits": ["Stable CAD is a phenotype proxy, not a stable-angina-specific cohort.",
                   "Only 4 cases and 2 healthy controls; all findings are exploratory and require independent validation."],
        "output_files": [str(out.resolve())],
    }


def _extract_herb_exp_signature(path: Path, candidate_name: str | None = None):
    data = json.loads(path.read_text(encoding="utf-8"))
    summaries = data.get("Experiment_summary", [])
    if len(summaries) < 2:
        return None
    headers, row = summaries[0], summaries[1]
    meta = dict(zip(headers, row))
    detail = (data.get("Experiment_detail") or {}).get(f"{meta.get('EXP id')} Data Detail", {}).get("data", {})
    table = detail.get("diff_top300_csv")
    if not isinstance(table, list) or len(table) < 2:
        return None
    genes = []
    for gene, logfc, pvalue, fdr in table[1:]:
        try:
            value = float(logfc)
        except (TypeError, ValueError):
            continue
        if gene and str(gene).upper() not in {"NA", "NONE"}:
            genes.append({"gene": str(gene).strip(), "log2FC": value,
                          "pvalue": _float_or_none(pvalue), "fdr": _float_or_none(fdr)})
    extract = lambda key: _text_cell(meta.get(key))
    organism = extract("Organism")
    return {
        "candidate_id": _text_cell(meta.get("Herb id")) or _text_cell(meta.get("Compound id")) or _text_cell(meta.get("Formula id")),
        "name": candidate_name or extract("Herb name") or extract("Ingredient name") or extract("Formula name"),
        "candidate_type": extract("Subject type"),
        "source_database": "HERB 2.0 measured GEO perturbation experiment",
        "source_id": extract("EXP id"),
        "source_accession": extract("GSE id"),
        "genes": genes,
        "context": {"species": organism, "tissue": extract("Tissue"), "cell_type": extract("Cell type"),
                    "experiment_detail": extract("Experiment detail"), "control_condition": extract("Control condition"),
                    "control_samples": extract("Control samples"), "treatment_condition": extract("Experiment subject detail"),
                    "treatment_samples": extract("Treatment samples"), "platform": extract("Platform"),
                    "dose_time": extract("Experiment detail")},
    }


def _float_or_none(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _text_cell(value):
    if isinstance(value, dict):
        return str(value.get("title") or value.get("text") or "")
    return str(value or "")


def load_herb_experiment_signature(experiment_file: str | Path, candidate_name: str | None = None):
    """Read a HERB experiment-detail JSON artifact into a measured DE signature."""
    return _extract_herb_exp_signature(Path(experiment_file), candidate_name)


def _merged_candidate_signature(candidate, detail):
    merged = detail.get("Merge_diff_gene") or {}
    table = merged.get("diff_sig_gene#human")
    if not isinstance(table, list) or len(table) < 2:
        return None
    headers = [str(x.get("title") or "") if isinstance(x, dict) else str(x) for x in table[0]]
    rows = []
    for values in table[1:]:
        if not isinstance(values, list) or len(values) < len(headers):
            continue
        item = dict(zip(headers, values))
        gene = str(item.get("Gene name") or "").strip()
        logfc = _float_or_none(item.get("Log2 (FC_avg)"))
        if gene and logfc is not None and gene.upper() not in {"NA", "NONE"}:
            rows.append({"gene": gene, "log2FC": logfc, "fdr": _float_or_none(item.get("FDR BH"))})
    if len(rows) < 25:
        return None

    experiment_table = merged.get("drug_experiment") or []
    experiment_rows = []
    if experiment_table:
        exp_headers = [str(x) for x in experiment_table[0]]
        for values in experiment_table[1:]:
            if isinstance(values, list) and len(values) >= len(exp_headers):
                exp = dict(zip(exp_headers, values))
                if _text_cell(exp.get("Organism")).casefold() == "homo sapiens":
                    experiment_rows.append({
                        "experiment_id": _text_cell(exp.get("EXP id")),
                        "gse": _text_cell(exp.get("GSE id")),
                        "detail": _text_cell(exp.get("Experiment detail")),
                        "control": _text_cell(exp.get("Control condition")),
                        "treatment": _text_cell(exp.get("Experiment subject detail")),
                        "tissue": _text_cell(exp.get("Experiment detail")),
                        "control_samples": _text_cell(exp.get("Control samples")),
                        "treatment_samples": _text_cell(exp.get("Treatment samples")),
                    })
    tissues = sorted({x["tissue"] for x in experiment_rows if x["tissue"]})
    return {
        "candidate_id": candidate.get("entity_id"),
        "candidate_name": candidate.get("Herb cn name") or candidate.get("Herb en name"),
        "candidate_fdr_bh": candidate.get("fdr_bh"),
        "name": candidate.get("Herb cn name") or candidate.get("Herb en name"),
        "candidate_type": "Herb",
        "source_database": "HERB 2.0 merged measured human treatment-control experiments",
        "source_id": ";".join(x["experiment_id"] for x in experiment_rows if x["experiment_id"]),
        "source_accession": ";".join(sorted({x["gse"] for x in experiment_rows if x["gse"]})),
        "genes": rows,
        "context": {"species": "Homo sapiens", "tissue": "; ".join(tissues) or "unspecified",
                    "cell_type": "pooled/unspecified across HERB experiments",
                    "source_experiments": experiment_rows,
                    "dose_time": "See source_experiments; HERB signature is merged across experiments."},
    }


def retrieve_candidate_signatures(client, candidates, data_dir: str | Path, max_workers=4):
    """Retrieve HERB's measured human treatment signatures for candidate herbs."""
    data_dir = Path(data_dir)
    data_dir.mkdir(parents=True, exist_ok=True)
    signatures, failures = [], []

    def fetch(candidate):
        herb_id = str(candidate.get("entity_id") or "")
        if not herb_id:
            return candidate, None, "missing HERB entity id"
        local = data_dir / f"{herb_id}.json"
        try:
            if local.is_file():
                detail = json.loads(local.read_text(encoding="utf-8"))
            else:
                detail = client.herb_detail(herb_id)
                local.write_text(json.dumps(detail, ensure_ascii=False), encoding="utf-8")
            return candidate, detail, None
        except Exception as exc:
            return candidate, None, str(exc)[:180]

    with ThreadPoolExecutor(max_workers=max_workers) as pool:
        futures = [pool.submit(fetch, candidate) for candidate in candidates]
        for future in as_completed(futures):
            candidate, detail, error = future.result()
            herb_id = candidate.get("entity_id")
            if error:
                failures.append({"candidate_id": herb_id, "reason": f"HERB herb detail unavailable: {error}"})
                continue
            signature = _merged_candidate_signature(candidate, detail)
            if signature:
                signatures.append(signature)
            else:
                failures.append({"candidate_id": herb_id, "reason": "HERB 2.0 has no eligible human differential treatment signature for this herb."})
    signatures.sort(key=lambda x: (x.get("candidate_fdr_bh") if x.get("candidate_fdr_bh") is not None else 1,
                                  x.get("candidate_id") or ""))
    return {"signatures": signatures, "failures": failures,
            "candidate_entities_checked": len(candidates),
            "signatures_found": len(signatures)}


def score_reversal(disease_results, treatment_signatures, min_shared_genes=100, disease_context=None,
                   disease_significant_n=None, min_disease_significant_genes=10,
                   allow_low_confidence=False):
    disease = {str(row["gene"]).upper(): float(row["log2FC"])
               for row in disease_results if row.get("gene") and row.get("log2FC") is not None}
    disease_tissue = str((disease_context or {}).get("tissue") or "").casefold()
    scores, rejected, diagnostics = [], [], []
    for signature in treatment_signatures:
        context = signature.get("context") or {}
        genes = {str(row["gene"]).upper(): float(row["log2FC"])
                 for row in signature.get("genes", []) if row.get("gene") and row.get("log2FC") is not None}
        common = sorted(set(disease) & set(genes))
        if context.get("species") != "Homo sapiens":
            rejected.append({"source_id": signature.get("source_id"), "reason": "非人源扰动签名；保留为跨物种证据，不进入主分数。"})
            continue
        if len(common) < min_shared_genes:
            rejected.append({"source_id": signature.get("source_id"), "reason": "共同可测基因数不足。", "shared_genes": len(common)})
            continue
        rho_value = spearmanr([disease[g] for g in common], [genes[g] for g in common]).statistic
        if not np.isfinite(rho_value):
            rejected.append({"source_id": signature.get("source_id"), "reason": "共同基因效应无变异或相关系数不可计算。", "shared_genes": len(common)})
            continue
        rho = float(rho_value)
        treatment_significant = {
            str(row["gene"]).upper(): float(row["log2FC"])
            for row in signature.get("genes", [])
            if row.get("gene") and row.get("log2FC") is not None
            and row.get("fdr") is not None and float(row["fdr"]) < 0.05
        }
        significant_common = sorted(set(disease) & set(treatment_significant))
        significant_rho = None
        if len(significant_common) >= 3:
            value = spearmanr([disease[g] for g in significant_common],
                              [treatment_significant[g] for g in significant_common]).statistic
            if np.isfinite(value):
                significant_rho = float(value)
        treatment_tissue = str(context.get("tissue") or "").casefold()
        if disease_tissue and disease_tissue in treatment_tissue:
            context_match = "same_context"
        elif any(word in disease_tissue for word in ("blood", "pbmc")) and any(word in treatment_tissue for word in ("blood", "pbmc")):
            context_match = "related_blood_context"
        else:
            context_match = "context_not_matched"
        underpowered = (disease_significant_n is not None
                        and disease_significant_n < min_disease_significant_genes)
        if underpowered and not allow_low_confidence:
            diagnostics.append({
                "candidate_id": signature.get("candidate_id"), "name": signature.get("name"),
                "source_id": signature.get("source_id"), "source_accession": signature.get("source_accession"),
                "raw_spearman_rho": rho, "shared_genes": len(common), "context_match": context_match,
                "treatment_fdr05_shared_genes": len(significant_common),
                "treatment_fdr05_raw_spearman_rho": significant_rho,
                "label": "diagnostic_only_not_ranked",
            })
            continue
        scores.append({
            "candidate_id": signature.get("candidate_id"), "name": signature.get("name"),
            "candidate_name": signature.get("candidate_name"),
            "candidate_type": signature.get("candidate_type"), "source_database": signature.get("source_database"),
            "source_id": signature.get("source_id"), "source_accession": signature.get("source_accession"),
            "connectivity_score": -rho, "raw_spearman_rho": rho, "shared_genes": len(common),
            "context": context,
            "context_match": context_match,
            "evidence_label": ("exploratory_proxy_signature_match_low_confidence" if underpowered
                               else "exploratory_observed_signature_match"),
            "interpretation": ("使用未达到预设显著基因门槛的疾病效应排序计算，仅用于提出低置信度验证假设；"
                               "不是疾病特异结论、药效证明或临床推荐。" if underpowered else
                               "负相关提示表达逆转样模式；这是探索性假设，不是疗效证明。"),
        })
    scores.sort(key=lambda row: row["connectivity_score"], reverse=True)
    underpowered = (disease_significant_n is not None
                    and disease_significant_n < min_disease_significant_genes)
    state = (("exploratory_low_confidence" if underpowered else "scored_exploratory") if scores
             else "insufficient_disease_signature" if underpowered
             else "no_eligible_measured_signatures")
    return {"state": state,
            "method": "negative Spearman correlation over shared measured genes",
            "min_shared_genes": min_shared_genes, "min_disease_significant_genes": min_disease_significant_genes,
            "disease_significant_genes": disease_significant_n, "scores": scores,
            "diagnostic_correlations": diagnostics, "rejected": rejected,
            "p_value": None, "fdr": None,
            "limitations": ["When a disease signature misses the significance gate, any emitted score is a low-confidence hypothesis based on ranked effects, not a validated disease signature.",
                            "HERB experiment differential genes are precomputed by HERB 2.0 and each treatment arm has small sample size.",
                            "This score prioritizes validation experiments; it is not efficacy or clinical recommendation evidence."]}
