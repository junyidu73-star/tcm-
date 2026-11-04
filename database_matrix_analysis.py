"""Reproducible exploratory analysis for a cached GEO series matrix."""
from __future__ import annotations

import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import hypergeom, ttest_ind


def _read_series(path: Path):
    metadata = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("!Sample_title") or line.startswith("!Sample_geo_accession"):
                key, *values = line.rstrip("\n").split("\t")
                metadata[key] = [value.strip().strip('"') for value in values]
            if line.startswith("!series_matrix_table_begin"):
                columns = [value.strip().strip('"') for value in handle.readline().rstrip("\n").split("\t")]
                frame = pd.read_csv(handle, sep="\t", header=None, names=columns, quotechar='"', na_values=["NA", "null"])
                break
    frame = frame.set_index("ID_REF")
    frame.columns = [value.strip('"') for value in frame.columns]
    return metadata, frame.apply(pd.to_numeric, errors="coerce")


def _read_annotation(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ID\t"):
                columns = line.rstrip("\n").split("\t")
                return pd.read_csv(handle, sep="\t", names=columns, dtype=str, na_values=["NA"], quoting=csv.QUOTE_NONE)
    raise ValueError("GPL6480 annotation table header is missing")


def _bh(values):
    values = np.asarray(values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = ranked * len(values) / np.arange(1, len(values) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    output = np.empty_like(adjusted)
    output[order] = np.clip(adjusted, 0, 1)
    return output


def analyze_gse54977(data_dir: str | Path, output_dir: str | Path):
    """ACS-vs-stable-angina comparison of patient-specific endothelial response ratios."""
    data_dir, output_dir = Path(data_dir), Path(output_dir)
    matrix_file = data_dir / "GSE54977_series_matrix.txt.gz"
    annotation_file = data_dir / "GPL6480.annot.gz"
    if not matrix_file.is_file() or not annotation_file.is_file():
        raise FileNotFoundError("GSE54977 processed matrix or GPL6480 annotation is missing")
    metadata, probes = _read_series(matrix_file)
    titles = metadata.get("!Sample_title", [])
    if len(titles) != probes.shape[1] or len(titles) != 12:
        raise ValueError(f"Expected 12 sample-level log-ratio columns; found {probes.shape[1]}")
    sa_idx = [i for i, title in enumerate(titles) if title.casefold().startswith("sa_")]
    acs_idx = [i for i, title in enumerate(titles) if title.casefold().startswith("acs_")]
    if len(sa_idx) != 4 or len(acs_idx) != 8:
        raise ValueError(f"Unexpected groups parsed from GEO titles: SA={len(sa_idx)}, ACS={len(acs_idx)}")

    annotation = _read_annotation(annotation_file).set_index("ID")
    annotation = annotation[~annotation.index.duplicated(keep="first")]
    common = probes.index.intersection(annotation.index)
    joined = annotation.reindex(common).copy()
    joined["_symbol"] = joined.get("Gene symbol", pd.Series(index=joined.index, dtype=object)).fillna("").astype(str).str.strip()
    joined = joined[joined["_symbol"].ne("") & ~joined["_symbol"].str.contains(r"///|;|,", regex=True)]
    joined["_value"] = probes.loc[joined.index].median(axis=1, skipna=True)
    # Collapse duplicate probes using the median, preserving one patient-level value per gene.
    values = probes.loc[joined.index]
    values.index = joined["_symbol"].to_numpy()
    values = values.groupby(level=0).median()

    sa_values, acs_values = values.iloc[:, sa_idx].to_numpy(), values.iloc[:, acs_idx].to_numpy()
    tests = ttest_ind(acs_values, sa_values, axis=1, equal_var=False, nan_policy="omit")
    result = pd.DataFrame({"gene": values.index.astype(str),
                           "log2FC": np.nanmean(acs_values, axis=1) - np.nanmean(sa_values, axis=1),
                           "pvalue": tests.pvalue,
                           "mean_stable_angina": np.nanmean(sa_values, axis=1),
                           "mean_acs": np.nanmean(acs_values, axis=1),
                           "n_stable_angina": np.sum(np.isfinite(sa_values), axis=1),
                           "n_acs": np.sum(np.isfinite(acs_values), axis=1)})
    result = result.dropna(subset=["pvalue", "log2FC"])
    result["FDR"] = _bh(result["pvalue"].to_numpy())
    result["significant"] = (result["FDR"] < 0.05) & (result["log2FC"].abs() >= 0.5)
    result = result.sort_values(["FDR", "pvalue", "gene"]).reset_index(drop=True)

    # Annotation-based exploratory ORA, using all mapped genes tested as the universe.
    gene_to_go = defaultdict(set)
    go_names = {}
    go_columns = joined[["Gene symbol", "GO:Process ID", "GO:Process"]].fillna("")
    for symbol_value, id_value, name_value in go_columns.itertuples(index=False, name=None):
        symbol = str(symbol_value).strip()
        ids = str(id_value).split("///")
        names = str(name_value).split("///")
        for go_id, go_name in zip(ids, names):
            go_id, go_name = go_id.strip(), go_name.strip()
            if go_id.startswith("GO:") and symbol:
                gene_to_go[go_id].add(symbol)
                if go_name:
                    go_names[go_id] = go_name
    tested = set(result["gene"])
    pathways = []
    for direction, condition in (("ACS_higher", result["significant"] & (result["log2FC"] > 0)),
                                 ("stable_angina_higher", result["significant"] & (result["log2FC"] < 0))):
        selected = set(result.loc[condition, "gene"])
        if not selected:
            continue
        raw = []
        for go_id, members in gene_to_go.items():
            members = members & tested
            overlap = selected & members
            if len(members) < 5 or len(overlap) < 2:
                continue
            p = hypergeom.sf(len(overlap) - 1, len(tested), len(members), len(selected))
            raw.append((go_id, go_names.get(go_id, go_id), len(overlap), len(members), p, sorted(overlap)))
        if raw:
            q = _bh([item[4] for item in raw])
            for item, adjusted in zip(raw, q):
                pathways.append({"name": item[1], "source": "GPL6480 GO:Process annotation; exploratory ORA",
                                 "direction": direction, "adjusted_p": float(adjusted), "pvalue": item[4],
                                 "overlap_n": item[2], "term_gene_n": item[3], "genes": item[5]})
    pathways.sort(key=lambda item: (item["adjusted_p"], item["name"]))
    output_dir.mkdir(parents=True, exist_ok=True)
    full_path = output_dir / "GSE54977_differential_expression.tsv"
    result.to_csv(full_path, sep="\t", index=False)
    pathway_path = output_dir / "GSE54977_GO_process_ORA.tsv"
    pd.DataFrame(pathways).to_csv(pathway_path, sep="\t", index=False)
    selected = result[result["significant"]]
    summary_path = output_dir / "GSE54977_analysis_summary.md"
    output_paths = [str(full_path.resolve()), str(pathway_path.resolve()), str(summary_path.resolve())]
    report = {
        "dataset": {"accession": "GSE54977", "platform": "GPL6480", "title": "Different gene expression profile in acute coronary syndrome patients compared to stable angina patients"},
        "design": {"description": "患者血浆刺激的人冠状动脉内皮细胞；每个样本列为冠状窦相对主动脉血浆刺激效应的预计算log比值。比较 ACS 与稳定型心绞痛患者。",
                   "case_group": "急性冠脉综合征（ACS）", "case_n": 8, "control_group": "稳定型心绞痛（SA）", "control_n": 4,
                   "sample_ids": metadata.get("!Sample_geo_accession", []), "sample_titles": titles,
                   "comparison_scope": "ACS_vs_stable_angina_plasma_induced_endothelial_response", "paired": False},
        "annotation": {"platform": "GPL6480 Agilent-014850 Whole Human Genome Microarray", "mapping": "GPL6480 Gene symbol; duplicate probes collapsed by gene median; ambiguous/missing symbols excluded", "tested_genes": len(result)},
        "analysis": {"case_n": 8, "control_n": 4, "paired": False,
                     "method": "Welch two-sample t-test on processed patient-level log ratios; Benjamini-Hochberg FDR",
                     "test": "ACS minus stable angina", "tested_genes": len(result), "significant_n": int(selected.shape[0]),
                     "results": selected.to_dict(orient="records")},
        "pathways": pathways,
        "limits": ["这是 ACS 与稳定型心绞痛患者血浆诱导内皮反应的差异，不是稳定型心绞痛患者组织对健康人的差异。",
                   "样本量小（SA n=4，ACS n=8）；使用预处理后的log比值与探索性Welch检验，结果需独立队列复核。",
                   "GO富集仅依赖GPL6480平台GO注释，属于探索性注释富集，不是疾病因果通路证据。",
                   "本队列没有中药扰动签名；不能据此单独生成中药药效排序或治疗推荐。"],
        "output_files": output_paths,
    }
    (output_dir / "GSE54977_analysis_report.json").write_text(
        __import__("json").dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    top_genes = selected.head(20)
    gene_table = "\n".join(f"| {r.gene} | {r.log2FC:.3f} | {r.FDR:.3g} |" for r in top_genes.itertuples())
    term_table = "\n".join(f"| {p['direction']} | {p['name']} | {p['adjusted_p']:.3g} | {', '.join(p['genes'][:8])} |" for p in pathways[:20])
    summary_path.write_text(
        "# GSE54977 数据库表达矩阵探索分析\n\n"
        "## 队列与比较范围\n\n"
        "GEO GSE54977 / GPL6480。矩阵为患者冠状窦相对主动脉血浆刺激冠状动脉内皮细胞的预计算 log 比值；比较 ACS（n=8）与稳定型心绞痛（n=4）。这不是稳定型心绞痛患者组织对健康人的差异。\n\n"
        "## 方法和结果\n\n"
        f"映射并检验 {len(result):,} 个基因；Welch 双样本 t 检验，BH 校正；FDR<0.05 且 |log2FC|≥0.5 得到 {len(selected)} 个基因。log2FC 为 ACS 减稳定型心绞痛。\n\n"
        "| Gene | log2FC (ACS-SA) | FDR |\n|---|---:|---:|\n" + gene_table + "\n\n"
        "## 平台注释 GO 过程富集（探索性）\n\n"
        "| 方向 | GO 过程 | 校正 P 值 | 重叠基因示例 |\n|---|---|---:|---|\n" + term_table + "\n\n"
        "## 限制与中药预测状态\n\n"
        "样本量小，且为体外血浆刺激的代理反应；当前只有疾病队列表达签名，没有接入中药/成分/方剂的实测扰动签名和经校准的连接分数。因此这次不能报告中药推荐或药物排名。HERB 2.0 等数据库的签名连接分析需真实上传并取得分数后才能作为候选线索。通路富集是平台注释关联，不代表因果机制。\n\n"
        "## 文件\n\n"
        "- 完整差异表达结果：`GSE54977_differential_expression.tsv`\n"
        "- GO 过程富集：`GSE54977_GO_process_ORA.tsv`\n"
        "- 结构化结果：`GSE54977_analysis_report.json`\n", encoding="utf-8")
    return report
