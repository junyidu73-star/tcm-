"""Discover public perturbation-expression studies for HERB herb/formula candidates.

This module performs metadata discovery only. A hit is not treated as a usable
signature until its samples, intervention/control, organism, tissue and matrix
have been reviewed and analyzed.
"""
from __future__ import annotations

import json
import re
import time
import urllib.parse


def _json(get, url):
    return json.loads(get(url))


def _terms_for(entity, kind):
    keys = (("Herb en name", "Herb latin name", "Herb pinyin name", "Herb cn name")
            if kind == "herb" else
            ("Formula en name", "Formula pinyin name", "Formula cn name"))
    terms = []
    for key in keys:
        value = str(entity.get(key) or "").strip()
        for alias in re.split(r"[;/|]+", value):
            alias = alias.strip()
            if alias and alias.casefold() not in {"na", "n/a", "none", "null"} and alias.casefold() not in {x.casefold() for x in terms}:
                terms.append(alias)
    # Avoid broad Chinese two-character queries where no scientific alias exists.
    return terms


def _geo_search(get, term, pause):
    params = {"db": "gds", "term": f'"{term}"[All Fields] AND gse[ETYP]',
              "retmode": "json", "retmax": "10000", "retstart": "0"}
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?" + urllib.parse.urlencode(params)
    time.sleep(pause)
    result = _json(get, url).get("esearchresult", {})
    ids = list(result.get("idlist", []))
    total = int(result.get("count", len(ids)))
    complete = total <= len(ids)
    rows = []
    for start in range(0, len(ids), 500):
        batch = ids[start:start + 500]
        summary_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?" + urllib.parse.urlencode(
            {"db": "gds", "id": ",".join(batch), "retmode": "json"})
        time.sleep(pause)
        summaries = _json(get, summary_url).get("result", {})
        for uid in summaries.get("uids", []):
            item = summaries.get(uid, {})
            accession = item.get("accession") or item.get("accessionversion")
            if accession and str(accession).upper().startswith("GSE"):
                sample_rows = item.get("samples") or []
                tokens = [x for x in re.findall(r"[a-z0-9]+", term.casefold()) if len(x) > 2]
                matched_samples = [sample for sample in sample_rows
                                   if tokens and all(token in str(sample.get("title", "")).casefold() for token in tokens)]
                controls = [sample for sample in sample_rows if re.search(
                    r"vehicle|control|untreated|DMSO|PBS", str(sample.get("title", "")), re.I)]
                rows.append({"accession": str(accession).split(".")[0].upper(),
                             "title": item.get("title", ""), "summary": item.get("summary", ""),
                             "platform": item.get("gpl", []),
                             "organism": item.get("taxon", ""), "data_type": item.get("gdstype", ""),
                             "sample_count": len(sample_rows),
                             "candidate_sample_count": len(matched_samples),
                             "control_like_sample_count": len(controls),
                             "candidate_sample_titles": [sample.get("title", "") for sample in matched_samples],
                             "candidate_sample_accessions": [sample.get("accession", "") for sample in matched_samples],
                             "sample_design_triage": ("candidate_and_control_labels_found_review_required"
                                                       if len(matched_samples) >= 2 and len(controls) >= 2
                                                       else "metadata_hit_design_not_confirmed"),
                             "url": f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
                             "source": "GEO/E-utilities", "metadata_only": True})
    return {"total": total, "retrieved": len(rows), "complete": complete, "records": rows}


def _biostudies_search(get, term, pause):
    cursor, seen, rows, pages = None, set(), [], 0
    while True:
        # Exact phrase search avoids tens of thousands of noisy hits for common
        # transliterations while still searching the full study record.
        phrase = '"' + term.replace('"', " ").strip() + '"'
        params = {"query": phrase, "collection": "arrayexpress",
                  "pageSize": "100", "pagination": "cursor"}
        if cursor:
            params["cursor"] = cursor
        url = "https://www.ebi.ac.uk/biostudies/api/v1/search?" + urllib.parse.urlencode(params)
        time.sleep(pause)
        obj = _json(get, url)
        hits = obj.get("hits")
        if not isinstance(hits, list):
            raise ValueError("BioStudies search response has no hits array")
        pages += 1
        for hit in hits:
            accession = str(hit.get("accession") or "")
            if accession:
                rows.append({"accession": accession,
                             "canonical_study": re.sub(r"^E-GEOD-", "GSE", accession, flags=re.I),
                             "title": hit.get("title", ""), "summary": hit.get("content", ""),
                             "url": "https://www.ebi.ac.uk/biostudies/arrayexpress/studies/" + urllib.parse.quote(accession),
                             "source": "BioStudies/ArrayExpress", "metadata_only": True})
        next_cursor = obj.get("nextCursor")
        if not next_cursor:
            return {"pages": pages, "complete": True, "records": rows}
        if next_cursor in seen:
            raise ValueError("BioStudies pagination cursor repeated")
        seen.add(next_cursor)
        cursor = next_cursor


def discover_candidate_perturbation_studies(get, stage, herbs, formulas=(), pause=0.36):
    """Search GEO and ArrayExpress/BioStudies for every named herb/formula.

    Results remain an evidence-discovery queue; this function deliberately does
    not auto-score metadata hits or assume that a title match is an experiment.
    """
    entities = [("herb", item) for item in herbs] + [("formula", item) for item in formulas]
    found, query_audit, errors = {}, [], []
    geo_cache, biostudies_cache = {}, {}
    for kind, entity in entities:
        entity_id = str(entity.get("entity_id") or "")
        label = entity.get("Herb cn name") or entity.get("Formula cn name") or entity_id
        for term in _terms_for(entity, kind):
            query_audit.append({"entity_id": entity_id, "entity_type": kind, "label": label,
                                "term": term, "source": "GEO/E-utilities"})
            try:
                stage(f"扰动表达数据检索：GEO · {label} · {term}")
                if term not in geo_cache:
                    geo_cache[term] = _geo_search(get, term, pause)
                result = geo_cache[term]
                query_audit[-1].update({k: result[k] for k in ("total", "retrieved", "complete")})
                for record in result["records"]:
                    key = (entity_id, record["accession"])
                    row = found.setdefault(key, {**record, "candidate_id": entity_id,
                                                 "candidate_name": label, "candidate_type": kind,
                                                 "matched_terms": [], "sources": []})
                    row["sample_count"] = max(row.get("sample_count", 0), record.get("sample_count", 0))
                    row["candidate_sample_count"] = max(row.get("candidate_sample_count", 0), record.get("candidate_sample_count", 0))
                    row["control_like_sample_count"] = max(row.get("control_like_sample_count", 0), record.get("control_like_sample_count", 0))
                    row["candidate_sample_titles"] = list(dict.fromkeys(
                        row.get("candidate_sample_titles", []) + record.get("candidate_sample_titles", [])))
                    row["candidate_sample_accessions"] = list(dict.fromkeys(
                        row.get("candidate_sample_accessions", []) + record.get("candidate_sample_accessions", [])))
                    if record.get("sample_design_triage") == "candidate_and_control_labels_found_review_required":
                        row["sample_design_triage"] = record["sample_design_triage"]
                    if term not in row["matched_terms"]:
                        row["matched_terms"].append(term)
                    if record["source"] not in row["sources"]:
                        row["sources"].append(record["source"])
            except Exception as exc:
                errors.append({"candidate_id": entity_id, "source": "GEO/E-utilities",
                               "term": term, "reason": str(exc)[:240]})
            query_audit.append({"entity_id": entity_id, "entity_type": kind, "label": label,
                                "term": term, "source": "BioStudies/ArrayExpress"})
            try:
                stage(f"扰动表达数据检索：BioStudies/ArrayExpress · {label} · {term}")
                if term not in biostudies_cache:
                    biostudies_cache[term] = _biostudies_search(get, term, pause)
                result = biostudies_cache[term]
                query_audit[-1].update({"pages": result["pages"], "complete": result["complete"]})
                for record in result["records"]:
                    key = (entity_id, record.get("canonical_study") or record["accession"])
                    row = found.setdefault(key, {**record, "candidate_id": entity_id,
                                                 "candidate_name": label, "candidate_type": kind,
                                                 "matched_terms": [], "sources": []})
                    if term not in row["matched_terms"]:
                        row["matched_terms"].append(term)
                    if record["source"] not in row["sources"]:
                        row["sources"].append(record["source"])
            except Exception as exc:
                errors.append({"candidate_id": entity_id, "source": "BioStudies/ArrayExpress",
                               "term": term, "reason": str(exc)[:240]})
    records = sorted(found.values(), key=lambda x: (x["candidate_type"], x["candidate_name"], x["accession"]))
    source_query_counts = {}
    for query in query_audit:
        source_query_counts[query.get("source", "unknown")] = source_query_counts.get(query.get("source", "unknown"), 0) + 1
    triage_counts = {}
    for record in records:
        triage = record.get("sample_design_triage") or "sample_metadata_unavailable"
        triage_counts[triage] = triage_counts.get(triage, 0) + 1
    return {"state": "metadata_discovered" if records else "no_metadata_hits",
            "candidate_count": len(entities), "query_count": len(query_audit),
            "study_record_count": len(records),
            "unique_study_count": len({str(x.get("accession")) for x in records if x.get("accession")}),
            "source_query_counts": source_query_counts, "sample_design_triage_counts": triage_counts,
            "records": records,
            "queries": query_audit, "errors": errors,
            "interpretation": "命中是候选研究元数据，不代表已确认药材干预、对照设计或表达签名；需继续检查样本、物种、组织、剂量、时间和矩阵。",
            "sources": ["GEO/E-utilities", "BioStudies/ArrayExpress"]}
