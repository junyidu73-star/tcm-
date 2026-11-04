"""Evidence-to-constraint intake and guarded training for literature claims."""
from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
import json


@dataclass
class MechanismClaim:
    claim_id: str
    subject_id: str
    predicate: str
    object_id: str
    direction: str
    species: str
    tissue: str
    cell_type: str
    dose: str
    duration: str
    assay: str
    source_id: str
    source_sentence: str
    adjudication: str = "unreviewed"
    reviewer: str = ""
    contradiction_group: str = ""


ALLOWED_ADJUDICATIONS = {"support", "refute", "uncertain", "unreviewed"}


def validate_claim(claim: dict):
    required = ("claim_id", "subject_id", "predicate", "object_id", "source_id", "source_sentence", "species", "assay", "adjudication")
    missing = [key for key in required if not str(claim.get(key, "")).strip()]
    if missing:
        return False, "missing required provenance/context fields: " + ", ".join(missing)
    if claim.get("adjudication") not in ALLOWED_ADJUDICATIONS:
        return False, "invalid adjudication label"
    if claim.get("adjudication") != "unreviewed" and not str(claim.get("reviewer", "")).strip():
        return False, "reviewed labels require a named curator"
    if claim.get("direction") not in {"increase", "decrease", "activate", "inhibit", "bind", "no_effect", "unknown"}:
        return False, "invalid/unknown direction must be recorded explicitly"
    return True, "eligible for curation; not automatically eligible for model training"


def build_training_manifest(claims, output_path: str | Path):
    """Export claims and report whether a source-grouped model fit is defensible."""
    rows = [asdict(claim) if isinstance(claim, MechanismClaim) else dict(claim) for claim in claims]
    valid = [row for row in rows if validate_claim(row)[0] and row.get("adjudication") in {"support", "refute", "uncertain"}]
    labels = {label: sum(row["adjudication"] == label for row in valid) for label in ("support", "refute", "uncertain")}
    source_groups = {row.get("source_id") for row in valid}
    fit_ready = len(valid) >= 100 and len(source_groups) >= 10 and all(labels[label] >= 20 for label in ("support", "refute"))
    manifest = {
        "status": "training_ready" if fit_ready else "awaiting_human_adjudicated_claims",
        "claims_total": len(rows), "curated_claims": len(valid), "label_counts": labels,
        "independent_source_groups": len(source_groups), "fit_performed": False,
        "readiness_rule": ">=100 adjudicated claims from >=10 independent sources; >=20 support and >=20 refute labels; split by source to prevent citation leakage",
        "features_must_include_context": ["species", "tissue", "cell_type", "dose", "duration", "assay"],
        "output": "scientific_constraints_claims.json",
    }
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps({"manifest": manifest, "claims": rows}, ensure_ascii=False, indent=2), encoding="utf-8")
    return manifest
