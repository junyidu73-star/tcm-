"""Read-only client for the public HERB 2.0 database; no third-party code."""
from __future__ import annotations

import gzip
import hashlib
import json
import time
import urllib.error
import urllib.request
from pathlib import Path


HERB_APIS = ("http://herb.ac.cn/chedi/api/", "http://47.92.70.12/chedi/api/")


class Herb2Error(RuntimeError):
    pass


class Herb2Client:
    def __init__(self, cache_dir: str | Path, timeout: int = 45):
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self.timeout = timeout

    def _post(self, payload: dict, cache_name: str):
        path = self.cache_dir / cache_name
        if path.is_file():
            try:
                return json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                path.unlink(missing_ok=True)
        body = json.dumps(payload, ensure_ascii=False).encode("utf-8")
        error = None
        for endpoint in HERB_APIS:
            request = urllib.request.Request(
                endpoint, data=body,
                headers={"Content-Type": "application/json", "Accept-Encoding": "gzip", "User-Agent": "ZhiHeng-Research/1.0"},
            )
            for attempt in range(2):
                try:
                    with urllib.request.urlopen(request, timeout=self.timeout) as response:
                        raw = response.read()
                        if response.headers.get("Content-Encoding", "").lower() == "gzip" or raw[:2] == b"\x1f\x8b":
                            raw = gzip.decompress(raw)
                    result = json.loads(raw.decode("utf-8", "replace"))
                    if not isinstance(result, dict) or "message" in result:
                        raise Herb2Error("HERB 2.0 returned an error response")
                    path.write_text(json.dumps(result, ensure_ascii=False), encoding="utf-8")
                    return result
                except (urllib.error.URLError, TimeoutError, OSError, json.JSONDecodeError, Herb2Error) as exc:
                    error = exc
                    if attempt < 1:
                        time.sleep(0.75)
        raise Herb2Error(f"HERB 2.0 request failed: {str(error)[:240]}") from None

    def search_disease(self, keyword: str):
        key = "disease_search_" + hashlib.sha256(keyword.casefold().encode("utf-8")).hexdigest()[:16] + ".json"
        return self._post({"func_name": "search_api", "keyword": keyword, "label": "Disease"}, key)

    def disease_detail(self, key_id: str):
        safe_id = "".join(ch for ch in key_id if ch.isalnum() or ch in "_-.")
        if not safe_id:
            raise ValueError("Invalid HERB disease ID")
        return self._post({"func_name": "detail_api", "key_id": safe_id, "label": "Disease"}, f"disease_{safe_id}.json")

    def herb_detail(self, key_id: str):
        safe_id = "".join(ch for ch in key_id if ch.isalnum() or ch in "_-." )
        if not safe_id:
            raise ValueError("Invalid HERB herb ID")
        return self._post({"func_name": "detail_api", "key_id": safe_id, "label": "Herb"}, f"herb_{safe_id}.json")

def _cell_text(value):
    if isinstance(value, dict):
        return str(value.get("title") or value.get("text") or "")
    return str(value or "")


def _rows_as_dicts(table, id_header, name_fields):
    if not isinstance(table, list) or not table:
        return []
    headers = [str(value) for value in table[0]]
    result = []
    for row in table[1:]:
        if not isinstance(row, list) or len(row) != len(headers):
            continue
        item = dict(zip(headers, row))
        key = _cell_text(item.get(id_header))
        if not key:
            continue
        names = {field: _cell_text(item.get(field)) for field in name_fields}
        result.append({
            "entity_id": key,
            "entity_type": id_header.split()[0].lower(),
            **names,
            "p_value": _as_float(item.get("P value")),
            "fdr_bh": _as_float(item.get("FDR BH")),
            "source_database": "HERB 2.0",
            "source_record_id": key,
            "evidence_type": "HERB_indirect_statistical_disease_association",
            "association_method": "shared-target Fisher exact test with Benjamini-Hochberg correction (HERB 2.0) — database prior, not treatment efficacy",
        })
    return result


def _as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def match_herb_disease(client: Herb2Client, disease: str, disease_en: str = "", synonyms=()):
    """Resolve a disease term and return separate HERB evidence layers."""
    terms = []
    for term in (disease_en, disease, *(synonyms or ())):
        term = str(term or "").strip()
        if term and term.casefold() not in {value.casefold() for value in terms}:
            terms.append(term)
    attempted = []
    matches = []
    for term in terms[:5]:
        response = client.search_disease(term)
        table = response.get("res_data", [])
        attempted.append({"keyword": term, "records": max(0, len(table) - 1) if isinstance(table, list) else 0})
        if isinstance(table, list) and table:
            headers = [str(value) for value in table[0]]
            for row in table[1:]:
                if not isinstance(row, list) or len(row) < len(headers):
                    continue
                record = dict(zip(headers, row))
                disease_id = _cell_text(record.get("Disease ID"))
                name = _cell_text(record.get("Disease Name"))
                alias = _cell_text(record.get("Disease Alias Name"))
                exact = any(term.casefold() in {name.casefold(), alias.casefold()} for term in terms)
                if disease_id and exact:
                    matches.append({"disease_id": disease_id, "name": name, "alias": alias, "matched_keyword": term})
            if matches:
                break
    # Deduplicate and keep only exact normalized-name or alias matches.
    unique = {}
    for item in matches:
        unique[item["disease_id"]] = item
    if not unique:
        return {"state": "no_exact_herb_disease_mapping", "queries": attempted, "matched_disease": None,
                "herb_candidates": [], "formula_candidates": [], "ingredient_candidates": [], "limits": ["HERB disease synonym lookup returned no exact disease entity; no broad disease transfer performed."]}
    chosen = next(iter(unique.values()))
    detail = client.disease_detail(chosen["disease_id"])
    return {
        "state": "matched",
        "queries": attempted,
        "matched_disease": chosen,
        "herb_candidates": _rows_as_dicts(detail.get("herb_disease", []), "Herb id", ["Herb pinyin name", "Herb cn name", "Herb en name", "Herb latin name"]),
        "formula_candidates": _rows_as_dicts(detail.get("formula_disease", []), "Formula id", ["Formula pinyin name", "Formula cn name", "Formula en name"]),
        "ingredient_candidates": _rows_as_dicts(detail.get("ingredient_disease", []), "Ingredient id", ["Ingredient name", "Ingredient alias name"]),
        "transcriptomic_signature_status": "available" if isinstance(detail.get("Merge_diff_gene", {}).get("diff_sig_gene#human"), list) else "no_significant_disease_signature_in_HERB",
        "disease_experiments": detail.get("Merge_diff_gene", {}).get("disease_experiment", []),
        "limits": [
            "HERB 2.0 association rows are statistical links inferred through shared targets; they are database priors, not evidence of clinical efficacy.",
            "HERB 2.0's own disease transcriptomic signature is absent when diff_sig_gene#human is not a result table.",
            "These associations are kept separate from this platform's custom transcriptomic perturbation score.",
        ],
    }
