"""Transparent research utilities used by the integrated workspace pages.

These routines operate only on user supplied, already extracted study data.
They never infer that a search hit is eligible or that an association is causal.
"""
from __future__ import annotations

import csv
import io
import math
import statistics
from collections import defaultdict


def _number(row, key, line):
    try:
        value = float(row[key])
    except (KeyError, TypeError, ValueError):
        raise ValueError(f"第 {line} 行缺少有效数值列 {key}。") from None
    if not math.isfinite(value):
        raise ValueError(f"第 {line} 行 {key} 不是有限数值。")
    return value


def _read_csv(text):
    if not isinstance(text, str) or not text.strip():
        raise ValueError("请上传含表头的 CSV 数据。")
    if len(text.encode("utf-8")) > 900_000:
        raise ValueError("文件超过 900 KB，请先精简或拆分。")
    reader = csv.DictReader(io.StringIO(text.lstrip("\ufeff")))
    rows = [{str(k).lstrip("\ufeff"): v for k, v in row.items() if k is not None} for row in reader]
    if not rows or len(rows) > 5000:
        raise ValueError("CSV 需要 1–5000 行数据。")
    return rows


def _study_effect(row, line, kind):
    if "yi" in row and ("sei" in row or "vi" in row):
        yi = _number(row, "yi", line)
        vi = _number(row, "vi", line) if row.get("vi", "").strip() else _number(row, "sei", line) ** 2
        if vi <= 0:
            raise ValueError(f"第 {line} 行方差必须大于 0。")
        return yi, vi, "user_effect"
    if kind == "continuous":
        nt, mt, st = (_number(row, k, line) for k in ("n_t", "mean_t", "sd_t"))
        nc, mc, sc = (_number(row, k, line) for k in ("n_c", "mean_c", "sd_c"))
        if nt < 2 or nc < 2 or st <= 0 or sc <= 0:
            raise ValueError(f"第 {line} 行连续结局样本量需≥2且标准差需>0。")
        df = nt + nc - 2
        pooled = math.sqrt(((nt - 1) * st * st + (nc - 1) * sc * sc) / df)
        d = (mt - mc) / pooled
        correction = 1 - 3 / (4 * df - 1)
        yi = correction * d
        vi = (nt + nc) / (nt * nc) + yi * yi / (2 * df)
        return yi, vi, "hedges_g"
    if kind == "binary":
        et, nt, ec, nc = (_number(row, k, line) for k in ("event_t", "total_t", "event_c", "total_c"))
        if min(et, ec) < 0 or min(nt, nc) <= 0 or et > nt or ec > nc:
            raise ValueError(f"第 {line} 行二分类事件数/总数无效。")
        cells = [et, nt - et, ec, nc - ec]
        if any(v == 0 for v in cells):
            et, nt, ec, nc = et + .5, nt + 1, ec + .5, nc + 1
        yi = math.log((et / nt) / (ec / nc))
        vi = 1 / et - 1 / nt + 1 / ec - 1 / nc
        return yi, vi, "log_risk_ratio"
    raise ValueError("CSV 需含 yi 与 sei/vi，或选择连续/二分类结局并提供对应列。")


def _normal_p(z):
    return math.erfc(abs(z) / math.sqrt(2))


def _chi_square_sf(value, df):
    if df < 1 or value <= 0:
        return 1.0
    shape, x = df / 2, value / 2
    # Regularized upper incomplete gamma for the integer/half-integer shapes
    # used by subgroup Q tests.
    if df % 2 == 0:
        term = total = 1.0
        for j in range(1, df // 2):
            term *= x / j
            total += term
        return min(1.0, math.exp(-x) * total)
    q = math.erfc(math.sqrt(x))
    current = 0.5
    while current < shape:
        q = q + (x ** current) * math.exp(-x) / math.gamma(current + 1)
        current += 1.0
    return min(1.0, max(0.0, q))


def _pooled(items):
    k = len(items)
    if k == 0:
        return None
    fixed_w = [1 / x["vi"] for x in items]
    sw = sum(fixed_w)
    fixed = sum(w * x["yi"] for w, x in zip(fixed_w, items)) / sw
    q = sum(w * (x["yi"] - fixed) ** 2 for w, x in zip(fixed_w, items))
    c = sw - sum(w * w for w in fixed_w) / sw if sw else 0
    tau2 = max(0.0, (q - (k - 1)) / c) if k > 1 and c > 0 else 0.0
    rw = [1 / (x["vi"] + tau2) for x in items]
    rsw = sum(rw)
    mu = sum(w * x["yi"] for w, x in zip(rw, items)) / rsw
    se = math.sqrt(1 / rsw)
    i2 = max(0.0, (q - (k - 1)) / q) * 100 if q > 0 and k > 1 else 0.0
    # A t critical value would be preferable for small k; report normal limits
    # explicitly and keep the small-sample caution visible in the UI.
    ci = [mu - 1.96 * se, mu + 1.96 * se]
    pi_se = math.sqrt(se * se + tau2)
    return {"k": k, "estimate": mu, "se": se, "ci95": ci, "tau2_dl": tau2,
            "q": q, "q_df": max(0, k - 1), "i2_percent": i2,
            "prediction_interval_normal_approx": [mu - 1.96 * pi_se, mu + 1.96 * pi_se],
            "fixed_effect_estimate": fixed, "fixed_effect_se": math.sqrt(1 / sw)}


def _egger(items):
    # Egger regression: standardized effect on precision; weighted by inverse
    # variance. This is a diagnostic for small-study effects, not proof of bias.
    points = [(x["yi"] / math.sqrt(x["vi"]), 1 / math.sqrt(x["vi"]), 1 / x["vi"]) for x in items]
    sw = sum(p[2] for p in points)
    sx = sum(p[2] * p[1] for p in points) / sw
    sy = sum(p[2] * p[0] for p in points) / sw
    sxx = sum(p[2] * (p[1] - sx) ** 2 for p in points)
    if sxx <= 0:
        return {"available": False, "reason": "精度无变异，无法拟合回归。"}
    slope = sum(p[2] * (p[1] - sx) * (p[0] - sy) for p in points) / sxx
    intercept = sy - slope * sx
    residual = sum(p[2] * (p[0] - intercept - slope * p[1]) ** 2 for p in points)
    se = math.sqrt((residual / max(1, len(points) - 2)) * (1 / sw + sx * sx / sxx))
    z = intercept / se if se else float("inf")
    return {"available": True, "intercept": intercept, "se": se, "z": z,
            "p_normal_approx": _normal_p(z),
            "caution": "小研究效应/漏斗图不对称并非发表偏倚的特异证据。"}


def run_meta_analysis(payload):
    rows = _read_csv(payload.get("csv", ""))
    kind = str(payload.get("outcome_type", "effect"))
    group_col = str(payload.get("group_column", "")).strip()
    risk_col = str(payload.get("risk_column", "")).strip()
    items = []
    effect_types = set()
    for idx, row in enumerate(rows, 2):
        yi, vi, effect_type = _study_effect(row, idx, kind)
        study = (row.get("study") or row.get("study_id") or f"研究 {idx - 1}").strip()
        items.append({"study": study[:200], "yi": yi, "vi": vi,
                      "subgroup": row.get(group_col, "未分组").strip() if group_col else "",
                      "risk": row.get(risk_col, "未记录").strip() if risk_col else "未记录"})
        effect_types.add(effect_type)
    if len(items) < 2:
        raise ValueError("至少需要 2 项独立研究才能计算合并效应。")
    if len(effect_types) > 1:
        raise ValueError("同一分析中不能混用不同效应量尺度。")
    pooled = _pooled(items)
    subgroups = []
    if group_col:
        groups = defaultdict(list)
        for item in items:
            groups[item["subgroup"] or "未分组"].append(item)
        if len(groups) < 2:
            raise ValueError("所选分组列只有一个取值，不能进行亚组分析。")
        for name, group_items in sorted(groups.items()):
            subgroups.append({"name": name, "result": _pooled(group_items)})
        within_q = sum(x["result"]["q"] for x in subgroups if x["result"])
        q_between = max(0.0, pooled["q"] - within_q)
        subgroup_test = {"q_between": q_between, "df": len(groups) - 1,
                         "p_chi_square_approx": _chi_square_sf(q_between, len(groups) - 1)}
    else:
        subgroup_test = None
    loo = []
    for index, item in enumerate(items):
        remaining = items[:index] + items[index + 1:]
        calc = _pooled(remaining)
        loo.append({"omitted_study": item["study"], "estimate": calc["estimate"], "i2_percent": calc["i2_percent"]})
    high_risk = [x for x in items if x["risk"].casefold() in {"high", "高", "高风险", "serious"}]
    low_risk_sensitivity = _pooled([x for x in items if x not in high_risk]) if high_risk and len(items) - len(high_risk) >= 2 else None
    k = len(items)
    return {"state": "completed_with_caveats", "method": "DerSimonian–Laird 随机效应；正态近似区间；未校正多重亚组探索",
            "effect_type": next(iter(effect_types)), "overall": pooled, "studies": items,
            "subgroups": subgroups, "subgroup_test": subgroup_test, "leave_one_out": loo,
            "high_risk_exclusion_sensitivity": low_risk_sensitivity,
            "egger": _egger(items) if k >= 10 else {"available": False, "reason": f"仅 {k} 项研究；常规漏斗图不对称检验通常至少需要约10项，当前不计算。"},
            "funnel": [{"study": x["study"], "effect": x["yi"], "se": math.sqrt(x["vi"])} for x in items],
            "heterogeneity_guidance": ("异质性较高；先检查研究设计、组织/模型、剂量、时间和结局定义等预先指定的调节因素。平台不自动遍历分组以追求较低 I²。" if pooled["i2_percent"] >= 50 else "当前 I² 未达到常用的中高异质性提示界值；仍应结合 τ²、研究差异和临床/生物学背景判断。"),
            "warnings": ["异质性不能单独决定合理分组；亚组应由临床/生物学依据预先定义。",
                         "发表偏倚检验效能有限；漏斗图不对称也可能来自真实异质性或方法差异。",
                         "当前效应区间采用正态近似；研究数少时应使用合适的小样本方法并由统计人员复核。",
                         "该分析只使用上传表格中的汇总数据，不自动读取论文全文或验证数据提取正确性。"]}


def run_result_diagnostic(payload):
    rows = _read_csv(payload.get("csv", ""))
    domain = str(payload.get("domain", "")).strip()[:80]
    prediction = str(payload.get("prediction", "")).strip()[:2000]
    columns = list(rows[0].keys())
    numeric = {}
    for column in columns:
        values = []
        for row in rows:
            try:
                value = float(row.get(column, ""))
                if math.isfinite(value):
                    values.append(value)
            except (TypeError, ValueError):
                pass
        if values:
            numeric[column] = {"n": len(values), "missing": len(rows) - len(values),
                               "mean": statistics.mean(values),
                               "sd": statistics.stdev(values) if len(values) > 1 else None,
                               "min": min(values), "max": max(values)}
    if not numeric:
        raise ValueError("没有可分析的数值列。请上传包含测量值的 CSV。")
    warnings = []
    for column, stats in numeric.items():
        if stats["missing"]:
            warnings.append(f"{column} 有 {stats['missing']} 个缺失值；先查明缺失机制，勿默认填零。")
        if stats["n"] < 3:
            warnings.append(f"{column} 有效观测少于3个，不能据此做稳定的模型修正。")
    if not any(name in columns for name in ("group", "condition", "treatment", "组别", "处理组")):
        warnings.append("未识别到组别列；目前无法重算处理效应或病例-对照差异。")
    if not any(name in columns for name in ("sample_id", "sample", "样本", "subject_id", "patient_id")):
        warnings.append("未识别到样本/受试者 ID；无法检查重复测量、配对或伪重复。")
    if not prediction:
        warnings.append("尚未提供原预测、预期方向或模型版本，不能量化预测偏差。")
    categories = ["实验设计/对照与重复", "数据质量、批次和缺失", "测量尺度与预处理", "样本/疾病/组织范围映射", "模型适用域与不确定性", "文献证据和科学约束"]
    return {"state": "exploratory_diagnostic", "domain": domain, "row_count": len(rows),
            "columns": columns, "numeric_summary": numeric, "review_categories": categories,
            "warnings": warnings or ["未发现基础格式问题；这不表示原假设成立或数据已通过统计审查。"],
            "next_steps": ["先核对原始文件、单位、实验批次、随机化/盲法、对照和生物学重复。",
                           "按研究设计重算主要终点并报告效应量和区间，不只比较显著性。",
                           "检索并人工核实同物种、组织、干预、剂量、时间和结局的独立证据。",
                           "将新规律写成带来源、适用域和反例的约束；在留出队列/实验上比较旧模型与新版本。",
                           "只有外部验证改善且无关键亚组退化时，才生成候选算法版本；保留旧版和回滚路径。"],
            "model_update": "未自动训练或覆盖任何模型。当前结果是数据质量与分析流程检查，不是模型诊断结论。"}
