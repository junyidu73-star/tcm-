"""Public study discovery and file inventories; no downloaded file is claimed analyzed."""
import json
import re
import time
import urllib.parse


def _json(get, url):
    return json.loads(get(url))


def search_biostudies(get, stage, terms):
    records, queries = {}, []
    for term in terms:
        cursor = None
        seen = set()
        item = {"query": term, "search_expression": f'title:"{term.replace(chr(34), " ")}"', "pages": 0, "complete": False, "error": None}
        queries.append(item)
        while True:
            stage(f"BioStudies 数据检索：{term}，第 {item['pages']+1} 页")
            # Searching entire record contents expands common disease words
            # (for example "stable" cell lines) and returns large noisy sets.
            # Restrict the text search to the study title for reproducible,
            # clinically relevant dataset discovery.
            params = {"query": item["search_expression"], "collection": "arrayexpress", "pageSize": 100, "pagination": "cursor"}
            if cursor: params["cursor"] = cursor
            try:
                obj = _json(get, "https://www.ebi.ac.uk/biostudies/api/v1/search?" + urllib.parse.urlencode(params))
                hits = obj.get("hits")
                if not isinstance(hits, list): raise ValueError("缺少研究记录列表")
                item["pages"] += 1
                item["total"] = obj.get("totalHits")
                item["total_exact"] = obj.get("isTotalHitsExact")
                for r in hits:
                    acc = r.get("accession", "")
                    if not acc: continue
                    entry = records.setdefault(acc, {"accession":acc,"title":r.get("title", ""),"summary":r.get("content", ""),
                        "url":"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/"+urllib.parse.quote(acc),
                        "canonical_study":re.sub(r"^E-GEOD-", "GSE", acc),"matched_queries":[],
                        "status":"研究元数据；样本和文件需进一步审核"})
                    entry["matched_queries"].append(term)
                if "nextCursor" not in obj:
                    raise ValueError("接口未返回 nextCursor，无法确认分页完成")
                next_cursor = obj.get("nextCursor")
                if not next_cursor:
                    item["complete"] = True
                    break
                if next_cursor in seen: raise ValueError("分页游标重复")
                seen.add(next_cursor)
                cursor = next_cursor
                time.sleep(.36)
            except Exception as exc:
                item["error"] = str(exc)[:300]
                break
    return {"state":"retrieved_metadata" if all(q["complete"] for q in queries) else "partial_metadata",
            "complete":all(q["complete"] for q in queries),"queries":queries,"retrieved":len(records),"records":list(records.values())}


def search_pride(get, stage, terms):
    records, queries = {}, []
    for term in terms:
        page = 0
        fingerprints = set()
        item = {"query":term,"pages":0,"complete":False,"error":None}
        queries.append(item)
        while True:
            stage(f"PRIDE 数据检索：{term}，第 {page+1} 页")
            url = "https://www.ebi.ac.uk/pride/ws/archive/v2/search/projects?" + urllib.parse.urlencode({"keyword":term,"pageSize":100,"page":page})
            try:
                obj = _json(get, url)
                rows = obj if isinstance(obj, list) else obj.get("_embedded", {}).get("projects") if isinstance(obj, dict) else None
                if not isinstance(rows, list): raise ValueError("项目列表格式不受支持")
                item["pages"] += 1
                fingerprint = tuple(str(r.get("accession", "")) for r in rows)
                if fingerprint and fingerprint in fingerprints: raise ValueError("接口重复返回相同页面")
                fingerprints.add(fingerprint)
                for r in rows:
                    acc = r.get("accession", "")
                    if not acc: continue
                    entry = records.setdefault(acc, {"accession":acc,"title":r.get("title", ""),"summary":r.get("projectDescription", ""),
                        "url":"https://www.ebi.ac.uk/pride/archive/projects/"+urllib.parse.quote(acc),"matched_queries":[],
                        "status":"研究元数据；蛋白定量表尚未分析"})
                    entry["matched_queries"].append(term)
                if len(rows) < 100:
                    item["complete"] = True
                    break
                page += 1
                time.sleep(.36)
            except Exception as exc:
                item["error"] = str(exc)[:300]
                break
    return {"state":"retrieved_metadata" if all(q["complete"] for q in queries) else "partial_metadata",
            "complete":all(q["complete"] for q in queries),"queries":queries,"retrieved":len(records),"records":list(records.values())}


def geo_inventory(get, accession):
    if not re.fullmatch(r"GSE\d+", accession): raise ValueError("无效的 GEO 编号")
    digits = accession[3:]
    bucket = "GSE" + (digits[:-3] if len(digits) > 3 else "") + "nnn"
    root = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{bucket}/{accession}/"
    inventory = {"files":[],"errors":[],"status":"inventory_only","sample_metadata_url":f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={accession}",
                 "raw_sequence_search_url":"https://www.ncbi.nlm.nih.gov/sra/?term="+urllib.parse.quote(accession)}
    for directory in ("matrix/", "suppl/"):
        try:
            html = get(root + directory).decode("utf-8", "replace")
            for href in re.findall(r'href=["\']([^"\']+)["\']', html, re.I):
                name = urllib.parse.unquote(href)
                if name in ("../", "./") or any(c in name for c in ("/", "\\", "?", "#")): continue
                inventory["files"].append({"name":name,"category":"processed_matrix" if directory=="matrix/" else "supplement",
                                           "url":root+directory+urllib.parse.quote(name),"access":"listed_not_downloaded"})
        except Exception as exc:
            inventory["errors"].append({"directory":directory,"reason":str(exc)[:200]})
    return inventory


def pride_inventory(get, accession):
    if not re.fullmatch(r"PXD\d+", accession): raise ValueError("无效的 PRIDE 编号")
    files, errors, page, seen = [], [], 0, set()
    while True:
        url = f"https://www.ebi.ac.uk/pride/ws/archive/v2/projects/{accession}/files?" + urllib.parse.urlencode({"page":page,"pageSize":100})
        try:
            obj = _json(get, url)
            rows = obj if isinstance(obj, list) else obj.get("_embedded", {}).get("files") if isinstance(obj, dict) else None
            if not isinstance(rows, list): raise ValueError("文件列表格式不受支持")
            fingerprint = tuple(str(r.get("accession") or r.get("fileName") or "") for r in rows)
            if fingerprint and fingerprint in seen: raise ValueError("文件分页重复")
            seen.add(fingerprint)
            for r in rows:
                location = r.get("publicFileLocations") or []
                # PRIDE still returns FTP as the public location. Its archive
                # mirror also serves HTTPS, which browsers can open safely.
                download_url = next((
                    str(loc.get("value", "")).replace("ftp://ftp.pride.ebi.ac.uk/", "https://ftp.pride.ebi.ac.uk/")
                    for loc in location if isinstance(loc, dict)
                    and str(loc.get("value", "")).startswith("ftp://ftp.pride.ebi.ac.uk/")
                ), None)
                files.append({"name":r.get("fileName"),"accession":r.get("accession"),"size":r.get("fileSizeBytes",r.get("fileSize")),"url":download_url,
                              "checksum":r.get("checksum"),"category":r.get("fileCategory"),"locations":location,"access":"listed_not_downloaded"})
            if len(rows)<100: break
            page += 1
            time.sleep(.36)
        except Exception as exc:
            errors.append({"page":page,"reason":str(exc)[:200]})
            break
    return {"status":"inventory_only","files":files,"errors":errors,
            "sample_and_files_url":f"https://www.ebi.ac.uk/pride/archive/projects/{accession}/files"}


def biostudies_inventory(get, accession):
    if not re.fullmatch(r"[A-Za-z0-9-]+", accession): raise ValueError("无效的 BioStudies 编号")
    inventory = {"status":"inventory_only","files":[],"errors":[],"sample_and_files_url":"https://www.ebi.ac.uk/biostudies/arrayexpress/studies/"+accession,
                 "raw_sequence_search_url":"https://www.ebi.ac.uk/ena/browser/text-search?query="+urllib.parse.quote(accession)}
    try:
        obj = _json(get, "https://www.ebi.ac.uk/biostudies/api/v1/studies/"+urllib.parse.quote(accession))
        seen_files=set()
        def walk(node, in_files=False):
            if isinstance(node, list):
                for child in node: walk(child, in_files)
            elif isinstance(node, dict):
                if in_files:
                    path = node.get("path") or node.get("name")
                    if isinstance(path, str) and path and not path.startswith("/") and ".." not in path.split("/"):
                        if path not in seen_files:
                            seen_files.add(path)
                            inventory["files"].append({"name":path,"url":"https://www.ebi.ac.uk/biostudies/files/"+accession+"/"+urllib.parse.quote(path, safe="/"),"access":"listed_not_downloaded"})
                for key, value in node.items():
                    if isinstance(value, (dict, list)): walk(value, key == "files")
        walk(obj)
    except Exception as exc:
        inventory["errors"].append({"reason":str(exc)[:200]})
    return inventory
