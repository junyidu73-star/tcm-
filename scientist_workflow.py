"""Evidence-bounded hypothesis review for the disease-to-drug workflow."""
import json

LITERATURE_PROTOCOL = (
    "按研究用途检索并筛选：疾病起始、持续驱动和复发；保护性代偿、病程与不可逆损伤；"
    "组织/细胞特异性和混杂；遗传、扰动、细胞内靶点结合和救援；"
    "靶点作用方向与分子结合、功能及构效；游离浓度、组织暴露、吸收代谢清除；"
    "酶/转运体相互作用、脱靶、长期安全及年龄肝肾功能；制剂质量和生产；"
    "中药临床试验、证候、剂量疗程及系统评价。每篇按研究对象、疾病、组织、"
    "干预、对照、时间、结局和设计判定用途；无足够信息标 uncertain。"
    "原始实验、综述、数据库关联和计算预测分开；不要把相关性当因果，"
    "不要臆造引用、效果、剂量或缺失条件；保留反证、阴性结果和排除理由。"
)

TOPICS = {
    "起始与持续驱动": '(pathogenesis OR "disease progression" OR recurrence)',
    "代偿损伤与病程": '(compensatory OR protective OR irreversible OR longitudinal)',
    "靶点因果与救援": '("gene perturbation" OR knockout OR "rescue experiment" OR "target engagement")',
    "作用方向与构效": '(pharmacology OR inhibition OR activation OR "structure activity relationship")',
    "暴露和清除": '(pharmacokinetics OR "unbound concentration" OR distribution OR clearance)',
    "安全和相互作用": '(toxicity OR "drug interactions" OR transporter OR "off target")',
    "中药临床和证候": '("traditional Chinese medicine" OR herbal OR syndrome) AND (clinical OR trial OR meta-analysis)',
    "提取制剂和质量": '(herbal OR "natural product") AND (extraction OR formulation OR quality)',
}

def targeted_literature(disease, pubmed):
    """Independent bounded discovery queries, retaining errors and actual coverage."""
    term = str(disease).replace('"', ' ').strip()[:150]
    results = []
    for topic, clause in TOPICS.items():
        query = f'"{term}" AND {clause}'
        try:
            result = pubmed(query, limit=None, with_meta=True)
            results.append({"topic": topic, "query": query, "state": "retrieved", "result": result,
                            "limit": None, "scope": "按接口分页检索；完整性以返回元数据为准，尚未完成纳排筛选"})
        except Exception:
            results.append({"topic": topic, "query": query, "state": "failed", "result": None,
                            "error": "检索失败，请重试。"})
    return results


def review(stage, disease, records, config, model):
    """Request bounded, auditable hypotheses; unvalidated output cannot change tiers."""
    entries = []
    for index, item in enumerate(records[:20]):
        entries.append({"id": f"R{index + 1}", "record": item})
    if not entries:
        return {"state": "skipped_no_evidence", "stage": stage, "reviews": []}
    instruction = (
        "模拟科学研究中的生成假设、寻找反证、比较解释和修订步骤。"
        "仅基于输入记录，不能编造来源或外部检索结果。每条判断引用输入id，"
        "缺失证据写未知；驱动、代偿、损伤和混杂要作为竞争假设。"
        "药物须区分结合、功能、方向、游离暴露和临床疗效。"
        "返回JSON对象 {\"reviews\":[{\"ids\":[\"R1\"],\"hypothesis\":\"...\","
        "\"counterevidence_or_unknown\":\"...\",\"next_test\":\"...\"}]}。"
        "只返回JSON，不得要求改变统计阈值或将关联升级为因果。"
    )
    try:
        result = model(config, [{"role": "user", "content": json.dumps(
            {"disease": disease, "stage": stage, "records": entries}, ensure_ascii=False, default=str)}],
            system=instruction, max_tokens=1800)
        value = json.loads(result)
        reviews = value.get("reviews")
        if not isinstance(reviews, list) or not reviews:
            raise ValueError("missing reviews")
        valid = {item["id"] for item in entries}
        for row in reviews:
            if (not isinstance(row, dict) or not isinstance(row.get("ids"), list)
                    or not row["ids"] or not set(row["ids"]).issubset(valid)
                    or any(not isinstance(row.get(k), str) or not row[k].strip()
                           for k in ("hypothesis", "counterevidence_or_unknown", "next_test"))):
                raise ValueError("invalid or ungrounded review")
        if {ident for row in reviews for ident in row['ids']} != valid:
            raise ValueError('review does not cover all supplied records')
        return {"state": "completed", "stage": stage, "input_records": entries,
                "reviews": reviews, "scope": "研究假设，未独立核验；不改变证据等级"}
    except Exception:
        return {"state": "failed", "stage": stage, "reviews": [],
                "error": "模型调用或结构与引用校验失败，请检查配置后重试。"}

_one_pass = review

def review(stage, disease, records, config, model):
    """Batched proposer -> critic -> revision; all decisions remain hypotheses."""
    if not records:
        return {'state':'skipped_no_evidence','stage':stage,'reviews':[],'pending_count':0}
    combined=[]; inputs=[]; traces=[]; processed=0
    for start in range(0,len(records),10):
        batch=records[start:start+10]; passes=[]
        for role in ('提出竞争机制与可检验假设','寻找反证、条件错配和替代解释','依据证据修订并安排区分实验'):
            def role_model(cfg, messages, system, max_tokens):
                context=json.dumps(passes,ensure_ascii=False,default=str)
                return model(cfg,messages,system=system+'\n本轮职责：'+role+'。前轮输出仅为待检验假设：'+context,max_tokens=max_tokens)
            result=_one_pass(stage,disease,batch,config,role_model)
            passes.append({'role':role,'result':result})
            if result['state']!='completed':
                return {'state':'partial' if processed else 'failed','stage':stage,'reviews':combined,
                        'input_records':inputs,'passes':traces+passes,'processed_count':processed,
                        'pending_count':len(records)-processed,'error':'预算、取消、调用或证据校验中断；未审查部分未作判断。'}
        prefix=f'B{start//10+1}-'
        inputs.extend({**entry,'id':prefix+entry['id']} for entry in result['input_records'])
        combined.extend({**row,'ids':[prefix+i for i in row['ids']]} for row in result['reviews'])
        traces.append({'batch_start':start,'passes':passes}); processed+=len(batch)
    return {'state':'completed','stage':stage,'reviews':combined,'input_records':inputs,
            'passes':traces,'processed_count':processed,'pending_count':0,
            'scope':'三角色科学假设迭代；规则生效另需来源及适用条件核验'}
