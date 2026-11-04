"""Paginated, provenance-preserving disease targets from Open Targets and GWAS Catalog."""
from __future__ import annotations

import hashlib
import json
import re
import time
import urllib.parse
import urllib.request
from pathlib import Path

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
GWAS_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2"
UA = "ZhiHengResearch/1.0 (public target evidence aggregation)"


def disease_ontology_crosswalk(disease):
    """Avoid disease-specific hardcoding; mappings come from live source queries."""
    return {
        "state": "dynamic_mapping_required",
        "requested_disease": str(disease or "").strip(),
        "entities": [],
        "interpretation": "No static disease crosswalk is applied. Exact and broader mappings are returned by queried databases and retain their scope labels.",
    }


def _get_json(url, network=None):
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": UA})
    if network is not None:
        raw, _ = network(url, headers=dict(req.header_items()), timeout=45)
        return json.loads(raw.decode("utf-8"))
    with urllib.request.urlopen(req, timeout=45) as response:
        return json.loads(response.read().decode("utf-8"))


def _post_ot(query, network=None):
    body = json.dumps({"query": query}).encode("utf-8")
    req = urllib.request.Request(OT_URL, data=body, headers={
        "Accept": "application/json", "Content-Type": "application/json", "User-Agent": UA,
    })
    if network is not None:
        raw, _ = network(OT_URL, body, dict(req.header_items()), timeout=45)
        obj = json.loads(raw.decode("utf-8"))
    else:
        with urllib.request.urlopen(req, timeout=45) as response:
            obj = json.loads(response.read().decode("utf-8"))
    if obj.get("errors"):
        raise RuntimeError(str(obj["errors"])[:400])
    return obj.get("data") or {}


def _cache(cache_dir, prefix, key):
    digest = hashlib.sha256(key.casefold().encode()).hexdigest()[:20]
    return cache_dir / f"{prefix}_{digest}.json"


def _ontology_synonyms(curie, cache_dir, network=None):
    """Return authoritative synonyms for an Open Targets disease concept."""
    curie = str(curie or "").strip()
    match = re.fullmatch(r"([A-Za-z][A-Za-z0-9]*)[_:]([A-Za-z0-9.]+)", curie)
    if not match:
        return {"label": None, "synonyms": [], "source": "OLS4", "complete": False}
    prefix, local_id = match.groups()
    ontology = prefix.casefold()
    if prefix.upper() == "EFO":
        iri = f"http://www.ebi.ac.uk/efo/{prefix.upper()}_{local_id}"
    else:
        iri = f"http://purl.obolibrary.org/obo/{prefix.upper()}_{local_id}"
    path = _cache(Path(cache_dir), "ols_term", curie)
    if path.exists():
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError):
            pass
    url = "https://www.ebi.ac.uk/ols4/api/ontologies/" + urllib.parse.quote(ontology, safe="") + "/terms?" + urllib.parse.urlencode({"iri": iri})
    headers = {"Accept": "application/json", "User-Agent": UA}
    if network is not None:
        raw, _ = network(url, headers=headers, timeout=20)
    else:
        req = urllib.request.Request(url, headers=headers)
        with urllib.request.urlopen(req, timeout=20) as response:
            raw = response.read()
    obj = json.loads(raw.decode("utf-8"))
    rows = (obj.get("_embedded") or {}).get("terms") or []
    term = next((row for row in rows if row.get("short_form", "").casefold() == curie.replace(":", "_").casefold()), rows[0] if rows else {})
    output = {"label": term.get("label"), "synonyms": term.get("synonyms") or [], "source": "OLS4",
              "ontology": prefix.upper(), "curie": curie, "complete": bool(term)}
    path.write_text(json.dumps(output, ensure_ascii=False), encoding="utf-8")
    return output


def _ot_search(term, cache_dir, network=None):
    path = _cache(cache_dir, "ot_search", term)
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    query = f'''query {{ search(queryString:{json.dumps(term)}, entityNames:["disease"], page:{{index:0,size:100}}) {{ hits {{ id name entity }} }} }}'''
    hits = _post_ot(query, network=network).get("search", {}).get("hits", [])
    output = {"query": term, "hits": hits, "retrieved": len(hits), "complete": True}
    path.write_text(json.dumps(output, ensure_ascii=False), encoding="utf-8")
    return output


def _ot_targets(efo_id, cache_dir, network=None):
    path = _cache(cache_dir, "ot_targets", efo_id)
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    page, size, rows, total = 0, 500, [], None
    while total is None or len(rows) < total:
        query = f'''query {{ disease(efoId:"{efo_id}") {{ id name associatedTargets(page:{{index:{page},size:{size}}}) {{ count rows {{ score target {{ id approvedSymbol approvedName }} datatypeScores {{ id score }} }} }} }} }}'''
        disease = _post_ot(query, network=network).get("disease")
        if not disease:
            break
        assoc = disease.get("associatedTargets") or {}
        total = int(assoc.get("count") or 0)
        batch = assoc.get("rows") or []
        if not batch:
            break
        rows.extend(batch)
        page += 1
        if page > 1000:
            raise RuntimeError("Open Targets pagination safety limit exceeded")
    output = {"disease_id": efo_id, "disease_name": disease.get("name") if disease else None,
              "count": total or 0, "retrieved": len(rows), "complete": len(rows) >= (total or 0), "targets": rows}
    path.write_text(json.dumps(output, ensure_ascii=False), encoding="utf-8")
    return output


def _gwas_page(resource, params, network=None):
    return _get_json(f"{GWAS_URL}/{resource}?{urllib.parse.urlencode(params)}", network=network)


def _gwas_traits(term, cache_dir, network=None):
    path = _cache(cache_dir, "gwas_traits", term)
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    page, size, rows, total = 0, 100, [], None
    while total is None or len(rows) < total:
        obj = _gwas_page("efo-traits", {"efo_trait": term, "page": page, "size": size}, network=network)
        total = int((obj.get("page") or {}).get("totalElements", 0))
        batch = (obj.get("_embedded") or {}).get("efo_traits") or []
        rows.extend(batch)
        if not batch:
            break
        page += 1
    output = {"query": term, "total": total or 0, "retrieved": len(rows), "complete": len(rows) >= (total or 0), "traits": rows}
    path.write_text(json.dumps(output, ensure_ascii=False), encoding="utf-8")
    return output


def _gwas_associations(efo_id, cache_dir, network=None):
    path = _cache(cache_dir, "gwas_assoc_efo_id_v2", efo_id)
    if path.exists():
        return json.loads(path.read_text(encoding="utf-8"))
    page, size, rows, total = 0, 100, [], None
    while total is None or len(rows) < total:
        obj = _gwas_page("associations", {"efo_id": efo_id, "show_child_traits": "false", "page": page, "size": size}, network=network)
        total = int((obj.get("page") or {}).get("totalElements", 0))
        batch = (obj.get("_embedded") or {}).get("associations") or []
        rows.extend(batch)
        if not batch:
            break
        page += 1
        time.sleep(0.05)
    output = {"efo_id": efo_id, "total": total or 0, "retrieved": len(rows), "complete": len(rows) >= (total or 0), "associations": rows}
    path.write_text(json.dumps(output, ensure_ascii=False), encoding="utf-8")
    return output


def retrieve_multisource_target_evidence(disease, synonyms=(), related_terms=(), cache_dir="runs/database_data/target_sources", network=None):
    """Query exact synonyms and explicit proxies. Proxy evidence stays separately labeled."""
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    exact = list(dict.fromkeys(str(x).strip() for x in (disease, *(synonyms or ())) if str(x).strip()))
    proxies = list(dict.fromkeys(str(x).strip() for x in (related_terms or ()) if str(x).strip()))
    exact_names = {re.sub(r"\s+", " ", x).strip().casefold() for x in exact}
    proxy_names = {re.sub(r"\s+", " ", x).strip().casefold() for x in proxies}
    ot_searches, concepts = [], {}
    for term in dict.fromkeys(exact + proxies):
        try:
            result = _ot_search(term, cache_dir, network=network)
            ot_searches.append(result)
            is_exact_query = re.sub(r"\s+", " ", term).strip().casefold() in exact_names
            for rank, hit in enumerate(result.get("hits", [])[:100]):
                label = re.sub(r"\s+", " ", str(hit.get("name") or "")).strip().casefold()
                scope = "exact" if label in exact_names else ("proxy" if label in proxy_names else None)
                scope_basis = "label_exact_match" if label in exact_names else ("label_proxy_match" if label in proxy_names else None)
                if scope != "exact" and is_exact_query and rank < 10:
                    try:
                        ontology_term = _ontology_synonyms(str(hit.get("id") or ""), cache_dir, network=network)
                        synonym_names = {re.sub(r"\s+", " ", str(value)).strip().casefold() for value in ontology_term.get("synonyms", [])}
                        if exact_names.intersection(synonym_names):
                            scope, scope_basis = "exact", "OLS4_exact_synonym"
                    except Exception:
                        pass
                if scope is None:
                    continue
                concept = concepts.setdefault(str(hit["id"]), {**hit, "scope": scope, "scope_basis": scope_basis, "matched_queries": [], "matched_labels": []})
                concept["matched_queries"].append(term)
                concept["matched_labels"].append({"query": term, "label": hit.get("name"), "scope": scope, "basis": scope_basis})
                if scope == "exact":
                    concept["scope"] = scope
                    concept["scope_basis"] = scope_basis
        except Exception as exc:
            ot_searches.append({"query": term, "hits": [], "complete": False, "error": str(exc)[:240]})
    ot_sets = []
    for concept in concepts.values():
        try:
            ot_sets.append({**_ot_targets(concept["id"], cache_dir, network=network), "scope": concept["scope"], "matched_queries": concept["matched_queries"]})
        except Exception as exc:
            ot_sets.append({"disease_id": concept["id"], "disease_name": concept.get("name"), "scope": concept["scope"], "complete": False, "error": str(exc)[:240], "targets": []})

    trait_searches, gwas_sets = [], []
    for term in dict.fromkeys(exact + proxies):
        try:
            result = _gwas_traits(term, cache_dir, network=network)
            trait_searches.append(result)
            for trait in result.get("traits", []):
                trait_name = str(trait.get("efo_trait") or "").casefold()
                scope = "exact" if trait_name in exact_names else "proxy"
                scope_basis = "label_exact_match" if scope == "exact" else "related_phenotype"
                if scope != "exact":
                    try:
                        ontology_term = _ontology_synonyms(str(trait.get("efo_id") or ""), cache_dir, network=network)
                        synonym_names = {re.sub(r"\s+", " ", str(value)).strip().casefold() for value in ontology_term.get("synonyms", [])}
                        if exact_names.intersection(synonym_names):
                            scope, scope_basis = "exact", "OLS4_exact_synonym"
                    except Exception:
                        pass
                try:
                    data = _gwas_associations(str(trait.get("efo_id") or ""), cache_dir, network=network)
                    gwas_sets.append({**data, "trait": trait, "scope": scope, "scope_basis": scope_basis})
                except Exception as exc:
                    gwas_sets.append({"trait": trait, "scope": scope, "complete": False, "error": str(exc)[:240], "associations": []})
        except Exception as exc:
            trait_searches.append({"query": term, "traits": [], "complete": False, "error": str(exc)[:240]})
    return {
        "state": "retrieved" if any(x.get("complete") for x in ot_searches + trait_searches) else "source_errors",
        "requested_disease": disease, "exact_terms": exact, "proxy_terms": proxies,
        "ontology_crosswalk": disease_ontology_crosswalk(disease),
        "OpenTargets": {"searches": ot_searches, "disease_concepts": list(concepts.values()), "target_sets": ot_sets,
                        "retrieved_target_count": sum(x.get("retrieved", 0) for x in ot_sets),
                        "complete": all(x.get("complete", False) for x in ot_sets) if ot_sets else bool(ot_searches)},
        "GWAS_Catalog": {"trait_searches": trait_searches, "association_sets": gwas_sets,
                         "retrieved_association_count": sum(x.get("retrieved", 0) for x in gwas_sets),
                         "complete": all(x.get("complete", False) for x in trait_searches + gwas_sets) if trait_searches else False},
        "interpretation": "Exact phenotype records and broader proxy concepts stay separate. GWAS Catalog contains curated top associations, not all genome-wide summary statistics.",
    }


def build_database_coverage(report):
    """Audit what was actually searched, analyzed, and used in a run."""
    retrieval = report.get("retrieval", [])
    rows = []
    def add(name, status, records=None, detail=""):
        rows.append({"database": name, "status": status, "records_or_entities": records, "detail": detail})
    for db_name in ("PubMed", "Europe PMC"):
        hits = [x for x in retrieval if x.get("database") == db_name]
        add(db_name, "queried" if hits else "not_queried", sum(x.get("retrieved", 0) for x in hits),
            f"{sum(x.get('total', 0) for x in hits)} total hits across {len(hits)} search branches; completeness is branch-specific")
    datasets = report.get("datasets") or {}
    geo = datasets.get("GEO", {}) if isinstance(datasets, dict) else {}
    geo_records = geo.get("records", []) if isinstance(geo, dict) else []
    cad_proxy = report.get("multi_cohort_cad_proxy") or {}
    cad_cohorts = cad_proxy.get("cohorts", [])
    geo_count = max(len(geo_records), len(cad_cohorts), len(report.get("database_omics_analyses", [])))
    add("GEO", "analyzed_and_searched" if report.get("database_omics_analyses") else "metadata_search_only",
        geo_count, f"{len(cad_cohorts)} disease proxy cohorts entered multi-cohort synthesis; {len(report.get('database_omics_analyses', []))} additional saved omics analyses; not all GEO records were matrix-analyzed")
    databases = report.get("databases") or {}
    for name in ("ArrayExpress", "PRIDE"):
        item = databases.get(name, {})
        add("BioStudies/ArrayExpress" if name == "ArrayExpress" else "PRIDE/ProteomeXchange",
            item.get("state", "not_queried"), item.get("retrieved", len(item.get("records", []))),
            "metadata/file inventory only; no blanket proteomics or expression-matrix analysis" if name == "PRIDE" else "metadata and file inventory; duplicate GEO accessions are not independent cohorts")
    herb = report.get("external_tcm_database") or {}
    add("HERB 2.0", herb.get("state", "not_queried"), sum(len(herb.get(k, [])) for k in ("herb_candidates", "formula_candidates", "ingredient_candidates")),
        f"{len(herb.get('herb_candidates', []))} herbs, {len(herb.get('formula_candidates', []))} formulas, {len(herb.get('ingredient_candidates', []))} ingredients")
    multi = report.get("multi_source_target_evidence") or {}
    for name in ("OpenTargets", "GWAS_Catalog"):
        source = multi.get(name) or {}
        exact = sum(x.get("scope") == "exact" for x in (source.get("target_sets", []) if name == "OpenTargets" else source.get("association_sets", [])))
        proxy = sum(x.get("scope") == "proxy" for x in (source.get("target_sets", []) if name == "OpenTargets" else source.get("association_sets", [])))
        count = source.get("retrieved_target_count", source.get("retrieved_association_count", 0))
        add("Open Targets" if name == "OpenTargets" else "GWAS Catalog", "queried" if source else "not_queried", count,
            f"exact result sets={exact}; broader proxy result sets={proxy}; complete={source.get('complete')}")
    add("ChEMBL", "not_queried", 0, "connector is not yet part of this run's herb-ingredient target ranking")
    add("LINCS/CMap", "not_loaded", 0, "full perturbation signature set is not integrated in this run")
    add("PRIDE quantitative proteomics", "metadata_only", databases.get("PRIDE", {}).get("retrieved", 0), "protein quantification tables were not downloaded and analyzed")
    return rows
