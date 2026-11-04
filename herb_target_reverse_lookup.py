"""Retrieve HERB 2.0 herb/formula target links from prioritized disease genes.

This uses HERB's public JSON API directly. Responses are cached per target so a
partial run can resume without re-downloading completed targets.
"""
from __future__ import annotations

import hashlib
import json
import re
import time
import urllib.error
import urllib.request
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

APIS = (
    "http://herb.ac.cn/chedi/api/",
    "http://47.92.70.12/chedi/api/",
)
UA = "ZhiHeng-Research/1.0 (HERB 2.0 target reverse lookup)"


def _post(payload: dict, timeout: int = 40):
    errors = []
    for endpoint in APIS:
        request = urllib.request.Request(
            endpoint,
            data=json.dumps(payload, ensure_ascii=False).encode("utf-8"),
            headers={"Accept": "application/json, text/plain, */*",
                     "Content-Type": "application/json", "User-Agent": UA},
        )
        try:
            with urllib.request.urlopen(request, timeout=timeout) as response:
                raw = response.read()
            result = json.loads(raw.decode("utf-8", "replace"))
            if result is None:
                raise ValueError(f"HERB returned an empty response from {endpoint}")
            return result
        except Exception as exc:
            errors.append(f"{endpoint}: {str(exc)[:180]}")
    raise RuntimeError("; ".join(errors))


def _decode_text(value):
    if isinstance(value, dict):
        value = value.get("title") or value.get("text") or ""
    if not isinstance(value, str):
        return str(value or "")
    decoded = re.sub(r"\\u([0-9a-fA-F]{4})", lambda m: chr(int(m.group(1), 16)), value)
    try:
        return decoded.encode("utf-16", "surrogatepass").decode("utf-16")
    except UnicodeError:
        return decoded


def _table_rows(table):
    if not isinstance(table, list) or not table:
        return []
    headers = [_decode_text(value) for value in table[0]]
    rows = []
    for row in table[1:]:
        if isinstance(row, list) and len(row) == len(headers):
            rows.append({key: value for key, value in zip(headers, row)})
    return rows


def _id(value):
    return _decode_text(value)


def _parse_target_relation(row, relation, gene):
    entity_id = _id(row.get("Herb id") or row.get("Formula id") or row.get("Ingredient id"))
    if not entity_id:
        return None
    if entity_id.startswith("HBFO") or "Formula id" in row:
        entity_type, prefix = "Formula", "Formula"
    elif entity_id.startswith("HBIN") or "Ingredient id" in row:
        entity_type, prefix = "Ingredient", "Ingredient"
    else:
        entity_type, prefix = "Herb", "Herb"
    names = [_decode_text(row.get(f"{prefix} {field}")) for field in ("name", "cn name", "Chinese name", "Chinese Name", "pinyin name", "PinYin Name", "en name", "English name", "English Name", "latin name", "Latin Name")]
    names = list(dict.fromkeys(name for name in names if name and name != "NA"))
    refs = []
    reference = row.get("Reference")
    if isinstance(reference, list):
        refs = _table_rows(reference)
    if not refs:
        pmid = _id(row.get("PubMed id") or row.get("PubMed ID") or row.get("Pubmed id"))
        paper_id = _id(row.get("Paper id") or row.get("Reference id"))
        if pmid or paper_id:
            refs = [{"pmid": pmid, "reference_id": paper_id}]
    normalized_fields = {re.sub(r"[^a-z0-9]", "", str(key).casefold()): value for key, value in row.items()}
    return {
        "entity_id": entity_id, "entity_type": entity_type,
        "name": names[0] if names else entity_id, "aliases": names[1:],
        "gene": gene, "relation_source": relation,
        "p_value": normalized_fields.get("pvalue"), "fdr_bh": normalized_fields.get("fdrbh"),
        "references": [{
            "reference_id": _id(ref.get("Reference ID") or ref.get("reference_id")),
            "pmid": _id(ref.get("PubMed ID") or ref.get("pmid")), "title": _decode_text(ref.get("Reference title") or ref.get("title")),
            "relationship": _decode_text(ref.get("Relationship") or ref.get("relationship")),
            "grade": _decode_text(ref.get("Grade") or ref.get("grade")),
            "supporting_sentence": _decode_text(ref.get("Supporting sentences") or ref.get("supporting_sentence")),
        } for ref in refs],
    }


def _upgrade_cached_target_edges(gene, links, cache_file, cache_dir):
    """Reparse cached official HERB detail tables with the current field mapping."""
    if links.get("relation_parser_version", 0) >= 3:
        return links
    parsed_by_kind = {"herb": [], "formula": [], "ingredient": []}
    for target_id in links.get("target_ids", []):
        detail_path = cache_dir / f"target_detail_{target_id}.json"
        if not detail_path.is_file():
            continue
        try:
            detail = json.loads(detail_path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            continue
        for key, source, destination in (
            ("herb_target", "HERB herb_target", "herb"),
            ("paper_herb_target", "HERB paper_herb_target", "herb"),
            ("formula_target", "HERB formula_target", "formula"),
            ("paper_formula_target", "HERB paper_formula_target", "formula"),
            ("ingredient_target", "HERB ingredient_target", "ingredient"),
            ("paper_ingredient_target", "HERB paper_ingredient_target", "ingredient"),
        ):
            for row in _table_rows(detail.get(key)):
                parsed = _parse_target_relation(row, source, gene)
                if parsed:
                    parsed_by_kind[destination].append(parsed)
    for kind, parsed in parsed_by_kind.items():
        if parsed:
            links[kind] = parsed
        else:
            links.setdefault(kind, [])
    links["relation_parser_version"] = 3
    tmp_cache = cache_file.with_suffix(".json.tmp")
    tmp_cache.write_text(json.dumps(links, ensure_ascii=False), encoding="utf-8")
    tmp_cache.replace(cache_file)
    return links


def retrieve_target_links(gene_symbols, cache_dir, progress=None, pause_s=0.08, max_workers=8):
    """Reverse-map gene symbols to all HERB herb/formula target associations.

    Searches every supplied symbol and retrieves all target detail rows. Search
    results are exact-symbol matched; ambiguity is retained in the audit data.
    """
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    by_gene, errors = {}, []
    symbols = list(dict.fromkeys(str(g).strip().upper() for g in gene_symbols if str(g).strip()))

    def gene_cache_path(gene):
        digest = hashlib.sha256(gene.encode()).hexdigest()[:16]
        return cache_dir / f"reverse_target_{digest}.json"

    pending = []
    for gene in symbols:
        cache_file = gene_cache_path(gene)
        if cache_file.is_file():
            try:
                cached = json.loads(cache_file.read_text(encoding="utf-8"))
                by_gene[gene] = _upgrade_cached_target_edges(gene, cached, cache_file, cache_dir)
                continue
            except (OSError, json.JSONDecodeError):
                cache_file.unlink(missing_ok=True)
        pending.append(gene)

    prefetched_search = None
    if pending:
        try:
            prefetched_search = _post({"func_name": "search_api", "keyword": pending[0], "label": "Target"})
        except Exception as exc:
            message = str(exc)[:300]
            for index, gene in enumerate(pending, 1):
                by_gene[gene] = {"gene": gene, "target_ids": [], "herb": [], "formula": [],
                                 "complete": False, "error": message}
                if progress:
                    progress(len(by_gene), len(symbols), gene, "source_unavailable")
            errors.append({"source": "HERB 2.0 target API", "error": message,
                           "affected_genes": len(pending), "state": "source_unavailable"})
            return {"state": "source_unavailable", "requested_genes": len(symbols),
                    "completed_genes": sum(bool(x.get("complete")) for x in by_gene.values()),
                    "errors": errors, "by_gene": by_gene}

    def retrieve_one(gene, initial_search=None):
        cache_file = gene_cache_path(gene)
        try:
            search = initial_search if initial_search is not None else _post({"func_name": "search_api", "keyword": gene, "label": "Target"})
            table = search.get("res_data", [])
            if not isinstance(table, list) or not table:
                # HERB legitimately returns an empty table for symbols absent
                # from its target catalog. Cache this as a successful no-match.
                links = {"gene": gene, "target_ids": [], "herb": [], "formula": [], "complete": True,
                         "ingredient": [], "relation_parser_version": 3, "match_state": "no_target_record"}
                tmp_cache = cache_file.with_suffix(".json.tmp")
                tmp_cache.write_text(json.dumps(links, ensure_ascii=False), encoding="utf-8")
                tmp_cache.replace(cache_file)
                return gene, links, None, "no_match"
            headers = [_decode_text(x) for x in table[0]]
            symbol_index = next((headers.index(key) for key in ("Gene Symbol", "Gene name", "Gene Name") if key in headers), 1)
            id_index = next((headers.index(key) for key in ("Target ID", "Target id", "Target Id") if key in headers), 0)
            matches = []
            for row in table[1:]:
                if not isinstance(row, list) or len(row) <= max(symbol_index, id_index):
                    continue
                if _decode_text(row[symbol_index]).upper() == gene:
                    match_id = _id(row[id_index])
                    if match_id:
                        matches.append(match_id)
            links = {"gene": gene, "target_ids": matches, "herb": [], "formula": [], "ingredient": [],
                     "relation_parser_version": 3, "complete": True}
            for target_id in matches:
                detail_cache = cache_dir / f"target_detail_{target_id}.json"
                detail = None
                if detail_cache.is_file():
                    try:
                        detail = json.loads(detail_cache.read_text(encoding="utf-8"))
                    except (OSError, json.JSONDecodeError):
                        detail_cache.unlink(missing_ok=True)
                if detail is None:
                    last_error = None
                    for attempt in range(3):
                        try:
                            detail = _post({"func_name": "detail_api", "key_id": target_id, "label": "Target", "v": target_id})
                            if not isinstance(detail, dict):
                                raise ValueError("invalid HERB target detail")
                            detail_cache.write_text(json.dumps(detail, ensure_ascii=False), encoding="utf-8")
                            break
                        except (urllib.error.URLError, TimeoutError, OSError, ValueError, json.JSONDecodeError) as exc:
                            last_error = str(exc)
                            if attempt < 2:
                                time.sleep(0.5 * (attempt + 1))
                    if detail is None:
                        raise RuntimeError(f"HERB detail failed for {target_id}: {last_error}")
                    if pause_s:
                        time.sleep(pause_s)
                for source_key, source_name, destination in (
                    ("herb_target", "HERB herb_target", "herb"),
                    ("paper_herb_target", "HERB paper_herb_target", "herb"),
                    ("formula_target", "HERB formula_target", "formula"),
                    ("paper_formula_target", "HERB paper_formula_target", "formula"),
                    ("ingredient_target", "HERB ingredient_target", "ingredient"),
                    ("paper_ingredient_target", "HERB paper_ingredient_target", "ingredient"),
                ):
                    for row in _table_rows(detail.get(source_key)):
                        parsed = _parse_target_relation(row, source_name, gene)
                        if parsed:
                            links[destination].append(parsed)
            tmp_cache = cache_file.with_suffix(".json.tmp")
            tmp_cache.write_text(json.dumps(links, ensure_ascii=False), encoding="utf-8")
            tmp_cache.replace(cache_file)
            return gene, links, None, "retrieved"
        except Exception as exc:
            error = {"gene": gene, "error": str(exc)[:300]}
            failed = {"gene": gene, "target_ids": [], "herb": [], "formula": [], "complete": False, "error": error["error"]}
            return gene, failed, error, "error"

    completed = 0
    worker_count = max(1, min(int(max_workers), 12))
    with ThreadPoolExecutor(max_workers=worker_count) as pool:
        futures = [pool.submit(retrieve_one, gene, prefetched_search if gene == pending[0] else None) for gene in pending]
        for future in as_completed(futures):
            gene, result, error, state = future.result()
            by_gene[gene] = result
            if error:
                errors.append(error)
            completed += 1
            if progress:
                progress(len(by_gene), len(symbols), gene, state)
    return {"state": "complete" if not errors else "partial", "requested_genes": len(symbols),
            "completed_genes": sum(bool(x.get("complete")) for x in by_gene.values()),
            "errors": errors, "by_gene": by_gene}


def get_catalog_totals(cache_dir):
    """Read catalog sizes from HERB's browse endpoint without truncating rows."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    path = cache_dir / "catalog_totals.json"
    if path.is_file():
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            path.unlink(missing_ok=True)
    totals = {}
    for label in ("Herb", "Formula", "Ingredient", "Target"):
        data = _post({"func_name": "browse_api", "label": label, "page": 1, "page_size": 1})
        totals[label] = {"total": int(data.get("total_num") or 0), "complete": "total_num" in data}
    path.write_text(json.dumps(totals, ensure_ascii=False), encoding="utf-8")
    return totals

