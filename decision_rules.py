"""Versioned, project-owned scientific rules and deterministic decision gates."""
import hashlib
import json
import math
import sqlite3
import time
from contextlib import contextmanager
from pathlib import Path

DB = Path(__file__).resolve().parent / 'runs' / 'decision_rules.sqlite3'
KINDS = {'gene_role': {'driver', 'protective', 'injury_marker'},
         'target_direction': {'activate', 'inhibit'},
         'drug_direction': {'activate', 'inhibit'},
         'safety': {'contraindicated', 'risk'}, 'exposure': {'insufficient', 'supported'}}

@contextmanager
def connect():
    DB.parent.mkdir(parents=True, exist_ok=True)
    db = sqlite3.connect(DB, timeout=15)
    db.row_factory = sqlite3.Row
    try:
        with db:
            db.execute('CREATE TABLE IF NOT EXISTS revisions (id TEXT PRIMARY KEY, owner TEXT NOT NULL, created REAL, payload TEXT NOT NULL)')
            yield db
    finally:
        db.close()

def validate(rule):
    if not isinstance(rule, dict):
        raise ValueError('规则必须是对象。')
    kind = rule.get('kind')
    if kind not in KINDS or rule.get('value') not in KINDS[kind]:
        raise ValueError('规则类别或结论无效。')
    if kind == 'drug_direction' and not str(rule.get('target_gene') or '').strip():
        raise ValueError('药物作用方向必须明确对应target_gene。')
    for name in ('subject', 'disease', 'species', 'tissue', 'stage', 'source_id', 'quote', 'source_text'):
        if not isinstance(rule.get(name), str) or not rule[name].strip() or len(rule[name]) > 30000:
            raise ValueError('规则缺少条件或来源：' + name)
    if len(rule['quote']) < 12 or rule['quote'] not in rule['source_text']:
        raise ValueError('引文须为提供的原文中至少12个字符的连续片段。')
    if kind == 'exposure':
        if rule.get('comparison_basis') != 'unbound_site_vs_unbound_effective':
            raise ValueError('暴露比较须明确为作用部位游离浓度与游离有效浓度；不能用Kd或总血药浓度替代。')
        for key in ('dose', 'time', 'assay'):
            if not rule.get(key):
                raise ValueError('暴露规则缺少实验条件：' + key)
        for key in ('site_unbound_nM', 'effective_unbound_nM'):
            value = rule.get(key)
            if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
                raise ValueError('浓度必须为同一nM尺度的正有限数。')
        value = 'insufficient' if rule['site_unbound_nM'] < rule['effective_unbound_nM'] else 'supported'
        if rule['value'] != value:
            raise ValueError('暴露结论与浓度比较不符。')
    return dict(rule)

def save_revision(owner, rules, parent=None, reason=''):
    if not isinstance(rules, list) or len(rules) > 200:
        raise ValueError('每个版本最多200条经核验规则。')
    rows = [validate(r) for r in rules]
    if parent:
        get_revision(owner, parent)
    payload = {'rules': rows, 'parent': parent, 'reason': str(reason)[:1000],
               'reviewer_user_id': str(owner), 'created': time.time(), 'status': 'human_approved'}
    encoded = json.dumps(payload, ensure_ascii=False, sort_keys=True)
    ident = hashlib.sha256(encoded.encode()).hexdigest()[:24]
    with connect() as db:
        db.execute('INSERT INTO revisions VALUES(?,?,?,?)', (ident, str(owner), payload['created'], encoded))
    return {'id': ident, **payload}

def get_revision(owner, ident):
    with connect() as db:
        row = db.execute('SELECT payload FROM revisions WHERE id=? AND owner=?', (ident, str(owner))).fetchone()
    if not row:
        raise ValueError('规则版本不存在或不属于当前账号。')
    return {'id': ident, **json.loads(row['payload'])}

def list_revisions(owner):
    with connect() as db:
        rows = db.execute('SELECT id,created,payload FROM revisions WHERE owner=? ORDER BY created DESC', (str(owner),)).fetchall()
    return [{'id': r['id'], 'created': r['created'], 'reason': json.loads(r['payload'])['reason'],
             'rule_count': len(json.loads(r['payload'])['rules'])} for r in rows]

def matches(rule, subject, context):
    if rule['subject'].casefold() != str(subject).casefold():
        return False
    # Unknown conditions never act as wildcards; deliberately broad rules need an explicit "any".
    for field in ('disease', 'species', 'tissue', 'stage'):
        if rule[field] == 'any':
            continue
        value = str(context.get(field) or '').strip()
        if not value or value.casefold() != rule[field].casefold():
            return False
    return True

def applicable(rules, subject, context):
    return [r for r in rules if matches(r, subject, context)]

def decide_gene(row, rules, context):
    out = dict(row)
    evidence = applicable(rules, row.get('gene'), context)
    roles = {r['value'] for r in evidence if r['kind'] == 'gene_role'}
    directions = {r['value'] for r in evidence if r['kind'] == 'target_direction'}
    reasons = []
    hold = False
    if len(roles) > 1 or len(directions) > 1:
        hold = True; reasons.append('适用证据冲突，暂缓进入干预查询，须区分竞争机制。')
    if roles == {'injury_marker'}:
        hold = True; reasons.append('当前适用证据仅支持损伤标志，尚不支持作为干预靶点。')
    out['scientific_decision'] = {'state': 'held' if hold else 'eligible',
        'role': next(iter(roles)) if len(roles) == 1 else 'unknown',
        'direction': next(iter(directions)) if len(directions) == 1 else 'unknown',
        'reasons': reasons, 'source_ids': sorted({r['source_id'] for r in evidence}),
        'unknowns': [x for x, values in [('mechanism_role', roles), ('intervention_direction', directions)] if not values]}
    return out

def apply_gene_rules(priority, rules, context):
    rows = [decide_gene(row, rules, context) for row in priority.get('eligible_genes', [])]
    eligible = [r for r in rows if r['scientific_decision']['state'] == 'eligible']
    return {**priority, 'eligible_genes': eligible,
        'high_priority_genes': [r for r in eligible if r.get('tier') == 'high'],
        'medium_priority_genes': [r for r in eligible if r.get('tier') == 'medium'],
        'scientifically_held_genes': [r for r in rows if r['scientific_decision']['state'] == 'held'],
        'drug_lookup_gene_count': len(eligible), 'scientific_context': context}

def apply_drug_rules(result, priority, rules, context):
    directions = {r['gene']: r.get('scientific_decision', {}).get('direction', 'unknown') for r in priority.get('eligible_genes', [])}
    rows = []
    for item in result.get('candidates', []):
        row = dict(item)
        evidence = applicable(rules, row.get('drug'), context)
        actions = {r['value'] for r in evidence if r['kind'] == 'drug_direction' and str(r.get('target_gene', '')).casefold() == str(row.get('gene', '')).casefold()}
        required = directions.get(row.get('gene'), 'unknown')
        reasons = []
        if row.get('gene') not in directions:
            reasons.append('对应基因未通过当前科学约束，暂缓该基因药物候选。')
        if len(actions) > 1:
            reasons.append('适用药物作用方向证据冲突。')
        elif required != 'unknown' and actions and required not in actions:
            reasons.append('药物已核验作用方向与疾病靶点干预方向相反。')
        if any(r['kind'] == 'safety' and r['value'] == 'contraindicated' for r in evidence):
            reasons.append('适用条件下存在经核验禁忌证据。')
        exposure = {r['value'] for r in evidence if r['kind'] == 'exposure'}
        if 'insufficient' in exposure:
            reasons.append('指定剂量时间下作用部位游离暴露低于匹配的有效浓度；暂停转化优先级。')
        unknowns = [name for name, present in [('direction', bool(actions) and required != 'unknown'),
                    ('exposure', bool(exposure)), ('safety', any(r['kind'] == 'safety' for r in evidence))] if not present]
        state = 'held' if reasons else ('requires_evidence' if unknowns else 'priority_for_validation')
        row['scientific_decision'] = {'state': state, 'reasons': reasons, 'unknowns': unknowns,
            'required_direction': required, 'source_ids': sorted({r['source_id'] for r in evidence})}
        rows.append(row)
    rows.sort(key=lambda r: {'priority_for_validation': 0, 'requires_evidence': 1, 'held': 2}[r['scientific_decision']['state']])
    output = {**result, 'candidates': rows, 'scientific_screening': {'method': '条件化证据规则；原始记录全部保留',
        'counts': {state: sum(r['scientific_decision']['state'] == state for r in rows)
                   for state in ('priority_for_validation', 'requires_evidence', 'held')}}}
    if isinstance(result.get('traditional_medicine'), dict):
        medicine = dict(result['traditional_medicine'])
        for field in ('candidates', 'formula_candidates', 'ingredient_candidates'):
            updated=[]
            for item in medicine.get(field, []):
                evidence=applicable(rules, item.get('name') or item.get('candidate_id'), context)
                reasons=[r['kind']+':'+r['value'] for r in evidence if r['kind'] in {'safety','exposure'} and r['value'] in {'contraindicated','insufficient'}]
                updated.append({**item, 'scientific_decision': {'state':'held' if reasons else 'requires_evidence',
                    'reasons':reasons,'unknowns':['组合效应','方向、制剂及目标组织暴露需独立核验'],
                    'source_ids':sorted({r['source_id'] for r in evidence})}})
            medicine[field]=updated
        output['traditional_medicine']=medicine
    return output
