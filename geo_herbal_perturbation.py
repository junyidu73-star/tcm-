"""Analyze a public human herbal-extract RNA-seq perturbation resource in GEO.

The first adapter targets GEO SuperSeries GSE244687 (KORE-Map tonifying-herb
transcriptomes). It analyzes only explicitly matched treatment/vehicle samples
within A549 Batch A, keeping each herb, solvent and dose as a separate contrast.
"""
from __future__ import annotations

import gzip
import json
import re
import shlex
import time
import urllib.request
from pathlib import Path

import numpy as np
import pandas as pd


GEO_SAMPLE_SOFT = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244687/soft/GSE244687_family.soft.gz"
GEO_RCM = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE244nnn/GSE244682/suppl/GSE244682_RCM_tonifying_A549.txt.gz"
NCBI_GENE_INFO = "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"


def _download(url: str, path: Path, retries=4):
    if path.is_file() and path.stat().st_size:
        return path
    path.parent.mkdir(parents=True, exist_ok=True)
    last = None
    for attempt in range(retries):
        try:
            request = urllib.request.Request(url, headers={"User-Agent": "ZhiHeng-Research/1.0", "Connection": "close"})
            with urllib.request.urlopen(request, timeout=90) as response:
                content = response.read()
            if len(content) < 100:
                raise ValueError("Public GEO response was unexpectedly short")
            temporary = path.with_suffix(path.suffix + ".part")
            temporary.write_bytes(content)
            temporary.replace(path)
            return path
        except Exception as exc:
            last = exc
            time.sleep(0.5 * (attempt + 1))
    raise RuntimeError(f"Could not download public resource {url}: {str(last)[:200]}")


def _parse_soft_samples(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        content = handle.read()
    samples = []
    for block in re.split(r"(?m)(?=^\^SAMPLE = )", content):
        if not block.startswith("^SAMPLE = "):
            continue
        item = {"accession": block.splitlines()[0].split("=", 1)[1].strip(), "characteristics": []}
        for line in block.splitlines()[1:]:
            if not line.startswith("!Sample_") or " = " not in line:
                continue
            key, value = line[len("!Sample_"):].split(" = ", 1)
            if key == "characteristics_ch1":
                item["characteristics"].append(value.strip())
            elif key in {"title", "description", "source_name_ch1", "organism_ch1"}:
                item[key] = value.strip()
        if item.get("title") and item.get("description"):
            samples.append(item)
    return samples


def _load_rcm(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        header = shlex.split(handle.readline())
        genes, rows = [], []
        expected = len(header) + 1
        for line in handle:
            fields = shlex.split(line)
            if len(fields) != expected:
                continue
            genes.append(fields[0])
            rows.append(fields[1:])
    values = np.asarray(rows, dtype=np.float64)
    return header, genes, values


def _load_ensembl_to_symbol(path: Path):
    if not path.is_file() or not path.stat().st_size:
        _download(NCBI_GENE_INFO, path)
    result = {}
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 6 or fields[0] != "9606":
                continue
            symbol = fields[2]
            if not symbol or symbol in {"-", "LOC"}:
                continue
            for ensembl_id in re.findall(r"Ensembl:(ENSG\d+)(?:\.\d+)?", fields[5]):
                result.setdefault(ensembl_id, symbol)
    return result


def _condition_from_title(title: str):
    match = re.match(r"^A549_BatchA_(water|ethanol) extract of (.+?) \((High|Medium|Low) dose\)_rep\d+$", title, re.I)
    if match:
        return {"kind": "treatment", "solvent": match.group(1).casefold(),
                "herb": match.group(2).strip(), "dose": match.group(3).casefold()}
    if re.match(r"^A549_BatchA_Vehicle\(DMSO\)_rep\d+$", title, re.I):
        return {"kind": "control", "solvent": "vehicle"}
    return None


def _candidate_for_extract(extract_name: str, herbs: list[dict]):
    normalized = re.sub(r"[^a-z0-9]+", "", extract_name.casefold())
    # Specific, audited name resolution for the two exact plant identities in
    # this public study. Red ginseng remains distinct from unprocessed ginseng.
    explicit = {"panaxginseng": "HERB004609", "glycyrrhizaeradix": "HERB001779"}
    selected_id = explicit.get(normalized)
    if not selected_id:
        return None
    return next((x for x in herbs if str(x.get("entity_id")) == selected_id), None)


def analyze_koremap_herbal_signatures(data_dir: str | Path, output_dir: str | Path,
                                      herbs: list[dict]):
    """Analyze direct measured treatment-control signatures for mapped herbs.

    RSEM expected counts are rounded to integers for PyDESeq2. This is an
    explicit exploratory approximation; the original FASTQ can be reprocessed
    for a confirmatory analysis.
    """
    data_dir, output_dir = Path(data_dir), Path(output_dir)
    cache_path = output_dir / "GEO_herbal_signature_analysis.json"
    candidate_ids = sorted(str(x.get("entity_id")) for x in herbs if x.get("entity_id"))
    if cache_path.is_file():
        try:
            cached = json.loads(cache_path.read_text(encoding="utf-8"))
            if cached.get("candidate_ids") == candidate_ids and all(
                Path(x.get("analysis_file", "")).is_file() for x in cached.get("result", {}).get("signatures", [])
            ):
                return cached["result"]
        except (OSError, ValueError, KeyError):
            pass
    soft_path = _download(GEO_SAMPLE_SOFT, data_dir / "GSE244687_family.soft.gz")
    matrix_path = _download(GEO_RCM, data_dir / "GSE244682_RCM_tonifying_A549.txt.gz")
    gene_info_path = _download(NCBI_GENE_INFO, data_dir / "Homo_sapiens.gene_info.gz")
    samples = _parse_soft_samples(soft_path)
    header, ensembl_ids, matrix = _load_rcm(matrix_path)
    col_index = {name: i for i, name in enumerate(header)}
    observations = []
    for sample in samples:
        code = re.sub(r"^sample\s+", "", sample["description"], flags=re.I).strip()
        if code not in col_index:
            continue
        condition = _condition_from_title(sample["title"])
        if not condition:
            continue
        condition.update({"sample_id": sample["accession"], "title": sample["title"], "code": code,
                          "column": col_index[code]})
        if condition["kind"] == "treatment":
            candidate = _candidate_for_extract(condition["herb"], herbs)
            if not candidate:
                continue
            condition["candidate_id"] = candidate["entity_id"]
            condition["candidate_name"] = candidate.get("Herb cn name") or condition["herb"]
        observations.append(condition)

    groups = {}
    controls = [x for x in observations if x["kind"] == "control"]
    for row in observations:
        if row["kind"] != "treatment":
            continue
        key = (row["candidate_id"], row["herb"], row["solvent"], row["dose"])
        groups.setdefault(key, []).append(row)

    gene_map = _load_ensembl_to_symbol(gene_info_path)
    signatures, group_summaries = [], []
    output_dir.mkdir(parents=True, exist_ok=True)
    try:
        from pydeseq2.dds import DeseqDataSet
        from pydeseq2.default_inference import DefaultInference
        from pydeseq2.ds import DeseqStats
    except ImportError as exc:
        raise RuntimeError("PyDESeq2 is required to analyze the public GEO count matrix") from exc

    # Fit once across all eligible groups, then estimate each pre-specified
    # treatment-vs-vehicle contrast. This avoids refitting the same count matrix
    # for every dose while keeping all treatment arms separate.
    eligible_groups = {key: rows for key, rows in groups.items()
                       if len(rows) >= 3 and len(controls) >= 3}
    group_labels = {key: "treatment_" + str(index) for index, key in enumerate(sorted(eligible_groups))}
    selected = controls + [row for key in sorted(eligible_groups) for row in eligible_groups[key]]
    condition_by_sample = {row["sample_id"]: "control" for row in controls}
    for key, rows in eligible_groups.items():
        for row in rows:
            condition_by_sample[row["sample_id"]] = group_labels[key]
    if selected:
        col_idx = [x["column"] for x in selected]
        raw = np.rint(matrix[:, col_idx]).clip(min=0).astype(np.int64)
        # Standard low-count prefilter: expression >=10 in at least 3 samples.
        keep = (raw >= 10).sum(axis=1) >= 3
        raw = raw[keep]
        gene_ids = np.asarray(ensembl_ids, dtype=object)[keep]
        counts = pd.DataFrame(raw.T, index=[x["sample_id"] for x in selected], columns=gene_ids)
        metadata = pd.DataFrame({"condition": [condition_by_sample[x["sample_id"]] for x in selected]},
                                index=counts.index)
        # Windows sandbox may deny joblib child-process creation; use the
        # deterministic serial backend so the analysis remains executable.
        inference = DefaultInference(n_cpus=1)
        dds = DeseqDataSet(counts=counts, metadata=metadata, design="~condition", refit_cooks=False,
                           inference=inference, quiet=True)
        dds.deseq2()
    for key, treated in sorted(eligible_groups.items()):
        candidate_id, extract_name, solvent, dose = key
        stats = DeseqStats(dds, contrast=["condition", group_labels[key], "control"], cooks_filter=False,
                           independent_filter=True, inference=inference, quiet=True)
        stats.summary()
        results = stats.results_df.copy()
        rows = []
        for ensembl_id, row in results.iterrows():
            symbol = gene_map.get(str(ensembl_id).split(".")[0])
            lfc, pvalue, fdr = row.get("log2FoldChange"), row.get("pvalue"), row.get("padj")
            if not symbol or pd.isna(lfc) or pd.isna(pvalue):
                continue
            rows.append({"gene": symbol, "log2FC": float(lfc),
                         "pvalue": float(pvalue), "fdr": None if pd.isna(fdr) else float(fdr)})
        if not rows:
            continue
        table = pd.DataFrame(rows).groupby("gene", as_index=False).agg(
            log2FC=("log2FC", "median"), pvalue=("pvalue", "min"), fdr=("fdr", "min"))
        result_path = output_dir / f"{candidate_id}_{solvent}_{dose}_GSE244687_DE.tsv"
        table.sort_values(["fdr", "pvalue"], na_position="last").to_csv(result_path, sep="\t", index=False)
        contrast_samples = controls + treated
        sample_titles = [x["title"] for x in contrast_samples]
        signature = {
            "candidate_id": candidate_id, "candidate_name": f"{treated[0]['candidate_name']}（{solvent}提取，{dose}剂量）",
            "name": treated[0]["candidate_name"], "candidate_type": "Herb",
            "source_database": "GEO", "source_id": "GSE244687/GSE244682",
            "source_accession": "GSE244687;GSE244682",
            "genes": table.to_dict(orient="records"),
            "context": {"species": "Homo sapiens", "tissue": "A549 lung adenocarcinoma cell line",
                        "cell_type": "A549", "treatment": f"{extract_name} {solvent} extract, {dose} dose",
                        "time": "24 h", "control": "vehicle DMSO, within A549 Batch A",
                        "source_experiments": sample_titles,
                        "processing": "RSEM expected counts rounded to integers; one PyDESeq2 model across treatment groups, separate treatment-vs-vehicle contrasts"},
            "analysis_file": str(result_path.resolve()),
        }
        signatures.append(signature)
        group_summaries.append({"candidate_id": candidate_id, "candidate_name": treated[0]["candidate_name"],
                                "extract": extract_name, "solvent": solvent, "dose": dose,
                                "case_n": len(treated), "control_n": len(controls),
                                "tested_genes": len(table),
                                "significant_n": int((table["fdr"] < 0.05).sum()),
                                "source_accession": "GSE244687/GSE244682",
                                "output_file": str(result_path.resolve())})
    result = {"state": "analyzed" if signatures else "no_eligible_mapped_signatures",
            "source": "GEO KORE-Map 1.0 tonifying-herb RNA-seq",
            "super_series": "GSE244687", "series_used": "GSE244682",
            "sample_records_in_super_series": len(samples),
            "matrix_columns": len(header), "mapped_treatment_control_samples": len(observations),
            "candidate_herbs_with_signatures": len({x["candidate_id"] for x in signatures}),
            "signature_count": len(signatures), "signatures": signatures,
            "group_analyses": group_summaries,
            "limitations": ["Only A549 Batch A water/ethanol extract versus within-batch DMSO vehicle comparisons are analyzed.",
                            "A549 cancer-cell perturbations are not cardiovascular tissue responses; context is not matched to blood CAD signatures.",
                            "RSEM expected counts were rounded for PyDESeq2; reprocessing raw reads is required for confirmatory use.",
                            "This herbal resource covers only herbs actually present in the GEO sample titles; no signature is imputed for absent candidates."]}
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    cache_path.write_text(json.dumps({"candidate_ids": candidate_ids, "result": result},
                                     ensure_ascii=False, separators=(",", ":")), encoding="utf-8")
    return result
