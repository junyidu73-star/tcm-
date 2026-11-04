"""Study-wise analysis and random-effects synthesis of public blood CAD cohorts.

These cohorts represent angiographic coronary artery disease, not a stable-angina
case-control definition. They are retained as a separate proxy evidence branch.
"""
from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t, ttest_1samp, ttest_ind

from database_matrix_analysis import _bh, _read_series


COHORTS = {
    "GSE12288": {"platform": "GPL96", "grouping": "CADi case/control", "contrast": "angiographic CAD vs angiography-negative control"},
    "GSE20680": {"platform": "GPL4133", "grouping": "stenosis grade", "contrast": "severe CAD (grade 2) vs control (grade 0); grade 1 excluded"},
    "GSE20681": {"platform": "GPL4133", "grouping": "matched case-control pairs", "contrast": "angiographic CAD vs matched control"},
}


def _sample_characteristics(path: Path, n_samples: int):
    fields = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            if line.startswith("!Sample_characteristics_ch1"):
                row = next(csv.reader([line.rstrip("\n")], delimiter="\t", quotechar='"'))
                fields.append(row[1:])
            elif line.startswith("!series_matrix_table_begin"):
                break
    if any(len(row) != n_samples for row in fields):
        raise ValueError("sample characteristics do not align with expression matrix")
    return fields


def _annotation_map(path: Path):
    with gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="") as handle:
        for line in handle:
            if line.startswith("ID\t"):
                header = next(csv.reader([line.rstrip("\n")], delimiter="\t", quotechar='"'))
                break
        else:
            raise ValueError(f"GEO platform annotation has no ID header: {path}")
        folded = {name.strip().casefold(): i for i, name in enumerate(header)}
        probe_i = folded.get("id")
        gene_i = next((folded[x] for x in ("gene symbol", "gene symbols", "symbol") if x in folded), None)
        if probe_i is None or gene_i is None:
            raise ValueError(f"GEO platform annotation has no probe/gene-symbol columns: {path}")
        output = {}
        for row in csv.reader(handle, delimiter="\t", quotechar='"'):
            if len(row) <= max(probe_i, gene_i):
                continue
            probe, gene = row[probe_i].strip(), row[gene_i].strip()
            if probe and gene and gene.casefold() not in {"na", "null", "---", "n/a"} and not re.search(r"///|;|,", gene):
                output.setdefault(probe, gene)
    return output


def _get_characteristic(fields, n, key):
    prefix = key.casefold() + ":"
    for row in fields:
        if row and all(str(value).strip().casefold().startswith(prefix) for value in row):
            return [str(value).split(":", 1)[1].strip() for value in row]
    return None


def _collapse_genes(probes: pd.DataFrame, annotation: dict[str, str]):
    probes = probes.copy()
    probes.index = probes.index.map(str)
    common = probes.index.intersection(annotation)
    if len(common) < 1000:
        raise ValueError(f"only {len(common)} matrix probes mapped to platform symbols")
    values = probes.loc[common].copy()
    values.index = [annotation[x] for x in common]
    return values.groupby(level=0).median()


def _study_effect(values, case_idx, control_idx, paired=False, pair_ids=None):
    case = values.iloc[:, case_idx].to_numpy(dtype=float)
    control = values.iloc[:, control_idx].to_numpy(dtype=float)
    if paired:
        pairs = {}
        for idx in case_idx:
            pairs.setdefault(pair_ids[idx], {})["case"] = idx
        for idx in control_idx:
            pairs.setdefault(pair_ids[idx], {})["control"] = idx
        valid = [pair for pair in pairs.values() if "case" in pair and "control" in pair]
        if len(valid) < 10:
            raise ValueError(f"only {len(valid)} valid matched pairs")
        diffs = np.column_stack([values.iloc[:, pair["case"]].to_numpy(dtype=float) -
                                 values.iloc[:, pair["control"]].to_numpy(dtype=float) for pair in valid])
        n = np.sum(np.isfinite(diffs), axis=1)
        effect = np.nanmean(diffs, axis=1)
        se = np.nanstd(diffs, axis=1, ddof=1) / np.sqrt(n)
        test = ttest_1samp(diffs, 0, axis=1, nan_policy="omit")
        return effect, se, test.pvalue, len(valid), len(valid), "paired t-test on within-pair expression differences"
    tests = ttest_ind(case, control, axis=1, equal_var=False, nan_policy="omit")
    effect = np.nanmean(case, axis=1) - np.nanmean(control, axis=1)
    se = np.sqrt(np.nanvar(case, axis=1, ddof=1) / case.shape[1] +
                 np.nanvar(control, axis=1, ddof=1) / control.shape[1])
    return effect, se, tests.pvalue, case.shape[1], control.shape[1], "Welch t-test on GEO processed values"


def _cohort_analysis(accession, data_dir: Path, annotation_file: Path, contrast_mode=None):
    matrix_path = data_dir / accession / f"{accession}_series_matrix.txt.gz"
    if not matrix_path.is_file():
        raise FileNotFoundError(f"missing public GEO matrix: {matrix_path}")
    metadata, probes = _read_series(matrix_path)
    probes.index = probes.index.map(str)
    sample_ids = metadata.get("!Sample_geo_accession", [])
    if len(sample_ids) != probes.shape[1]:
        raise ValueError("sample accession count differs from matrix columns")
    titles = metadata.get("!Sample_title", [])
    fields = _sample_characteristics(matrix_path, probes.shape[1])
    annotation = _annotation_map(annotation_file)
    values = probes.copy()
    if accession == "GSE12288":
        # The GEO matrix stores non-log fluorescence intensities for this array.
        if float(np.nanquantile(values.to_numpy(dtype=float), 0.9)) > 50:
            values = np.log2(values.clip(lower=0) + 1)
        groups = _get_characteristic(fields, len(sample_ids), "patient group")
        if not groups:
            groups = ["case" if "patient group: case" in title.casefold() else "control" for title in titles]
        case_idx = [i for i, group in enumerate(groups) if group.casefold() == "case"]
        control_idx = [i for i, group in enumerate(groups) if group.casefold() == "control"]
        intermediate_n, pair_ids, paired = 0, None, False
        scope = "broad angiographic CAD proxy; not stable-angina-specific"
    elif accession == "GSE20680":
        states = _get_characteristic(fields, len(sample_ids), "disease state")
        if not states:
            raise ValueError("GSE20680 disease-state metadata is missing")
        if contrast_mode == "include_intermediate_as_cases":
            case_idx = [i for i, state in enumerate(states) if state.casefold().startswith(("case (2)", "intermediate (1)"))]
        elif contrast_mode == "include_intermediate_as_controls":
            case_idx = [i for i, state in enumerate(states) if state.casefold().startswith("case (2)")]
        else:
            case_idx = [i for i, state in enumerate(states) if state.casefold().startswith("case (2)")]
        control_idx = [i for i, state in enumerate(states) if state.casefold().startswith("control (0)")]
        if contrast_mode == "include_intermediate_as_controls":
            control_idx += [i for i, state in enumerate(states) if state.casefold().startswith("intermediate (1)")]
        intermediate_n = sum(state.casefold().startswith("intermediate (1)") for state in states)
        pair_ids, paired = None, False
        scope = "strict angiographic CAD proxy; intermediate stenosis excluded from primary contrast"
    elif accession == "GSE20681":
        states = _get_characteristic(fields, len(sample_ids), "disease state")
        pair_text = _get_characteristic(fields, len(sample_ids), "pair")
        if not states or not pair_text:
            raise ValueError("GSE20681 case/control or pair metadata is missing")
        case_idx = [i for i, state in enumerate(states) if state.casefold().startswith("case (1)")]
        control_idx = [i for i, state in enumerate(states) if state.casefold().startswith("control (0)")]
        pair_ids = [re.sub(r"^pair\s*", "", x.strip(), flags=re.I).casefold() for x in pair_text]
        paired, intermediate_n = True, 0
        scope = "matched angiographic CAD proxy; not stable-angina-specific"
    else:
        raise ValueError(f"unsupported cohort: {accession}")
    if len(case_idx) < 10 or len(control_idx) < 10:
        raise ValueError(f"unexpected small or misparsed groups: case={len(case_idx)}, control={len(control_idx)}")

    genes = _collapse_genes(values, annotation)
    effect, se, pvalue, case_n, control_n, method = _study_effect(genes, case_idx, control_idx, paired, pair_ids)
    result = pd.DataFrame({"gene": genes.index.astype(str), "log2FC": effect, "SE": se, "pvalue": pvalue})
    result = result.replace([np.inf, -np.inf], np.nan).dropna(subset=["log2FC", "SE", "pvalue"])
    result["FDR"] = _bh(result["pvalue"].to_numpy())
    result["significant"] = result["FDR"].lt(0.05)
    result = result.sort_values(["FDR", "pvalue", "gene"]).reset_index(drop=True)
    return {
        "dataset": {"accession": accession, "platform": COHORTS[accession]["platform"], "species": "Homo sapiens", "tissue": "peripheral blood"},
        "design": {"case_group": "angiographic coronary artery disease", "control_group": "angiography-negative or low-stenosis controls",
                   "case_n": case_n, "control_n": control_n, "intermediate_n_excluded": intermediate_n,
                   "paired": paired, "comparison_scope": scope, "case_sample_ids": [sample_ids[i] for i in case_idx],
                   "control_sample_ids": [sample_ids[i] for i in control_idx]},
        "analysis": {"method": method, "tested_genes": int(len(result)), "significant_n": int(result["significant"].sum()),
                     "mapped_probes": int(len(probes.index.intersection(set(annotation)))),
                     "results": result.to_dict(orient="records")},
        "_table": result,
    }


def _random_effects_meta(study_tables: list[pd.DataFrame]):
    merged = None
    for i, table in enumerate(study_tables):
        part = table[["gene", "log2FC", "SE"]].rename(columns={"log2FC": f"b{i}", "SE": f"se{i}"})
        merged = part if merged is None else merged.merge(part, on="gene", how="inner")
    rows = []
    for row in merged.itertuples(index=False, name=None):
        gene, *values = row
        effects = np.asarray(values[0::2], dtype=float)
        ses = np.asarray(values[1::2], dtype=float)
        good = np.isfinite(effects) & np.isfinite(ses) & (ses > 0)
        effects, ses = effects[good], ses[good]
        if len(effects) < 2:
            continue
        fixed_w = 1 / ses**2
        fixed = float(np.sum(fixed_w * effects) / np.sum(fixed_w))
        q = float(np.sum(fixed_w * (effects - fixed)**2))
        df = len(effects) - 1
        c = float(np.sum(fixed_w) - np.sum(fixed_w**2) / np.sum(fixed_w))
        tau2 = max(0.0, (q - df) / c) if c > 0 else 0.0
        random_w = 1 / (ses**2 + tau2)
        pooled = float(np.sum(random_w * effects) / np.sum(random_w))
        hk_q = float(np.sum(random_w * (effects - pooled)**2) / df)
        pooled_se = float(np.sqrt(max(1.0, hk_q) / np.sum(random_w)))
        statistic = pooled / pooled_se
        rows.append({"gene": gene, "log2FC": pooled, "SE": pooled_se,
                     "pvalue": float(2 * t.sf(abs(statistic), df)), "study_n": len(effects),
                     "direction_consistency": int(np.sum(np.sign(effects) == np.sign(pooled))),
                     "tau2": tau2, "I2": max(0.0, (q - df) / q) if q > 0 else 0.0})
    result = pd.DataFrame(rows)
    result["FDR"] = _bh(result["pvalue"].to_numpy())
    result["significant"] = result["FDR"].lt(0.05)
    return result.sort_values(["FDR", "pvalue", "gene"]).reset_index(drop=True)


def analyze_cad_multicohort(data_dir: str | Path, annotation_dir: str | Path, output_dir: str | Path):
    """Analyze three larger blood CAD cohorts separately and synthesize effects.

    GSE20680/81 and GSE12288 are broader angiographic CAD proxies. They must not
    be mislabeled or directly substituted for stable-angina case-control data.
    """
    data_dir, annotation_dir, output_dir = Path(data_dir), Path(annotation_dir), Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    studies, tables = [], []
    for accession in COHORTS:
        result = _cohort_analysis(accession, data_dir, annotation_dir / f"{COHORTS[accession]['platform']}.annot.gz")
        table = result.pop("_table")
        path = output_dir / f"{accession}_differential_expression.tsv"
        table.to_csv(path, sep="\t", index=False)
        result["output_file"] = str(path.resolve())
        result["design"].pop("case_sample_ids", None)
        result["design"].pop("control_sample_ids", None)
        studies.append(result)
        tables.append(table)
    sensitivity = []
    for mode, label in (("include_intermediate_as_cases", "intermediate stenosis grouped with CAD"),
                        ("include_intermediate_as_controls", "intermediate stenosis grouped with controls")):
        result = _cohort_analysis("GSE20680", data_dir, annotation_dir / "GPL4133.annot.gz", contrast_mode=mode)
        table = result.pop("_table")
        path = output_dir / f"GSE20680_sensitivity_{mode}.tsv"
        table.to_csv(path, sep="\t", index=False)
        sensitivity.append({"definition": label, "case_n": result["design"]["case_n"],
                            "control_n": result["design"]["control_n"],
                            "tested_genes": result["analysis"]["tested_genes"],
                            "significant_n": result["analysis"]["significant_n"],
                            "output_file": str(path.resolve())})
    shared = set.intersection(*(set(table["gene"]) for table in tables))
    meta = _random_effects_meta([table[table["gene"].isin(shared)] for table in tables])
    meta_path = output_dir / "CAD_proxy_random_effects_meta.tsv"
    meta.to_csv(meta_path, sep="\t", index=False)
    return {
        "evidence_scope": "related_proxy_not_stable_angina_specific",
        "selection_basis": "three human peripheral-blood angiographic CAD cohorts repeatedly used in CAD transcriptome studies; all study matrices analyzed separately",
        "cohorts": studies,
        "case_definition_sensitivity": sensitivity,
        "meta_analysis": {"method": "DerSimonian-Laird random-effects meta-analysis with modified Hartung-Knapp t inference; Benjamini-Hochberg correction",
                          "cohort_n": len(studies), "shared_genes": len(shared), "tested_genes": len(meta),
                          "significant_n": int(meta["significant"].sum()),
                          "direction_concordant_n": int(((meta["direction_consistency"] >= 2) & meta["significant"]).sum()),
                          "results": meta.to_dict(orient="records"), "output_file": str(meta_path.resolve())},
        "limits": ["All cohorts represent broad angiographic CAD, not a stable-angina-specific patient-versus-healthy comparison.",
                   "GSE20680 primary contrast uses strict grade-2 cases versus grade-0 controls; grade-1 intermediate stenosis is excluded.",
                   "Study-wise estimates are meta-analyzed; raw matrices are not pooled across different cohorts and platforms. Public accession IDs were checked for duplicate samples, but participant-level overlap cannot be fully ruled out where GEO does not expose patient identifiers.",
                   "The meta-analysis is a discovery proxy and must stay separate from direct stable-angina evidence and clinical recommendations."],
    }
