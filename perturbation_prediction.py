"""A separate connectivity model for disease versus measured treatment signatures."""
from __future__ import annotations

from collections import defaultdict

import numpy as np
from scipy.stats import spearmanr


def score_treatment_signatures(disease_results, treatment_signatures, min_shared_genes=25):
    """Rank observed perturbations by reversal of a disease expression signature.

    Every treatment record must include gene, log2FC and experiment context.
    This function does not infer missing treatment profiles or impute unknown
    herbs. It returns no score when context or overlap is insufficient.
    """
    disease = {str(row.get("gene", "")).upper(): float(row["log2FC"])
               for row in disease_results if row.get("gene") and row.get("log2FC") is not None}
    scored, rejected = [], []
    for signature in treatment_signatures:
        context = signature.get("context") or {}
        profile = signature.get("genes") or []
        if context.get("species") != "Homo sapiens":
            rejected.append({"candidate_id": signature.get("candidate_id"), "reason": "Species must be harmonized before scoring."})
            continue
        treatment = {str(row.get("gene", "")).upper(): float(row["log2FC"])
                     for row in profile if row.get("gene") and row.get("log2FC") is not None}
        common = sorted(set(disease) & set(treatment))
        if len(common) < min_shared_genes:
            rejected.append({"candidate_id": signature.get("candidate_id"), "reason": "Insufficient shared measured genes.", "shared_genes": len(common)})
            continue
        rho = float(spearmanr([disease[g] for g in common], [treatment[g] for g in common]).statistic)
        if not np.isfinite(rho):
            rejected.append({"candidate_id": signature.get("candidate_id"), "reason": "Rank correlation was undefined.", "shared_genes": len(common)})
            continue
        scored.append({
            "candidate_id": signature.get("candidate_id"),
            "name": signature.get("name"),
            "candidate_type": signature.get("candidate_type"),
            "source_database": signature.get("source_database"),
            "source_id": signature.get("source_id"),
            "connectivity_score": -rho,
            "raw_spearman_rho": rho,
            "shared_genes": len(common),
            "context": context,
            "interpretation": "Negative correlation suggests a reversal-like expression pattern; this is a hypothesis, not efficacy evidence.",
        })
    scored.sort(key=lambda row: row["connectivity_score"], reverse=True)
    return {"state": "scored" if scored else "no_eligible_measured_signatures", "method": "negative Spearman correlation over shared measured genes",
            "min_shared_genes": min_shared_genes, "scores": scored, "rejected": rejected,
            "limitations": ["A reversal score is context-specific and does not prove clinical benefit.",
                            "The input treatment profiles must be actual measured perturbation signatures, not literature names or predicted targets."]}


def empty_prediction_state(reason, disease_accession=None):
    return {"state": "not_scored", "module": "independent transcriptomic connectivity model",
            "input_dataset": disease_accession, "scores": [], "reason": reason,
            "requires_observed_perturbation_signatures": True}
