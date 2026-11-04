"""Small neural ranking prototype with explicit monotonic constraints and group holdout.

Training labels must be experimentally adjudicated; imported rows are never invented.
The model is an offline research artifact until independent validation and deployment.
"""
import hashlib
import json
import math
import threading
from pathlib import Path
import numpy as np

FEATURES = ['direct_replication', 'genetic_functional_support', 'perturbation_rescue',
            'direction_match', 'exposure_match', 'safety_risk']
TRAIN_LOCK = threading.Lock()

def train(payload):
    rows = payload.get('rows')
    if payload.get('task') not in {'gene', 'drug'} or not isinstance(rows, list) or not 40 <= len(rows) <= 5000:
        raise ValueError('请指定gene或drug任务，并提供40—5000条核验标签及独立研究组。')
    groups = []
    X, y = [], []
    for r in rows:
        if not isinstance(r, dict) or type(r.get('label')) is not int or r['label'] not in (0, 1):
            raise ValueError('label必须为人工核验的0或1。')
        if not r.get('source_group') or not r.get('source_id') or r.get('adjudicated') is not True:
            raise ValueError('标签必须包含source_group、source_id及adjudicated=true。')
        values = r.get('features') or {}
        if not isinstance(values, dict):
            raise ValueError('features必须是特征对象。')
        if any(k not in FEATURES for k in values):
            raise ValueError('未知模型特征。')
        vector, missing = [], []
        for key in FEATURES:
            value = values.get(key)
            if value is not None and (isinstance(value, bool) or not isinstance(value, (float, int)) or not math.isfinite(value) or not 0 <= value <= 1):
                raise ValueError('特征须为0—1有限数；未知用null。')
            missing.append(float(value is None))
            vector.append(0.5 if value is None else (1-value if key == 'safety_risk' else value))
        X.append(vector + missing); y.append(r['label']); groups.append(str(r['source_group']))
    unique = sorted(set(groups), key=lambda x: hashlib.sha256(x.encode()).hexdigest())
    if len(unique) < 8:
        raise ValueError('至少需要8个独立研究组，不能随机拆分同一来源的重复记录。')
    # A source id must never span study groups; otherwise the split can leak citations.
    source_groups = {}
    for r in rows:
        source_groups.setdefault(str(r['source_id']), set()).add(str(r['source_group']))
    if any(len(x) > 1 for x in source_groups.values()):
        raise ValueError('同一source_id被分配到不同研究组，存在数据泄漏。')
    held = set(unique[:max(2, len(unique)//4)])
    mask = np.array([g not in held for g in groups]); X=np.array(X); y=np.array(y, dtype=float)
    if min(mask.sum(), (~mask).sum()) < 10 or len(set(y[mask])) < 2 or len(set(y[~mask])) < 2:
        raise ValueError('固定研究组留出后，训练和验证集均须至少10条且包含正负标签。')
    def sigmoid(z): return 1/(1+np.exp(-np.clip(z, -30, 30)))
    rng = np.random.default_rng(42)
    w1 = np.abs(rng.normal(0,.08,(12,8))); b1=np.zeros(8); w2=np.abs(rng.normal(0,.08,8)); b2=0.
    tx, ty = X[mask], y[mask]
    for _ in range(400):
        h = np.tanh(tx@w1+b1); p=sigmoid(h@w2+b2); d=(p-ty)/len(ty)
        dh = d[:,None]*w2*(1-h*h)
        g1=tx.T@dh+.01*w1; g2=h.T@d+.01*w2
        w1-=.15*g1; b1-=.15*dh.sum(0); w2-=.15*g2; b2-=.15*d.sum()
        # Monotonic rules are architectural constraints, not model-generated labels.
        w1[:6]=np.maximum(0,w1[:6]); w2=np.maximum(0,w2)
    def metrics(p, truth):
        p=np.clip(p,1e-7,1-1e-7)
        return {'brier': float(np.mean((p-truth)**2)),
                'log_loss': float(-np.mean(truth*np.log(p)+(1-truth)*np.log(1-p))),
                'accuracy_at_0_5': float(np.mean((p>=.5)==truth))}
    # Logistic regression is fitted using the identical training-only split.
    bw=np.zeros(12); bb=0.
    for _ in range(400):
        d=(sigmoid(tx@bw+bb)-ty)/len(ty)
        bw-=.15*(tx.T@d+.01*bw); bb-=.15*d.sum()
    probability=sigmoid(np.tanh(X[~mask]@w1+b1)@w2+b2)
    neural=metrics(probability,y[~mask]); baseline=metrics(sigmoid(X[~mask]@bw+bb),y[~mask])
    return {'task': payload['task'], 'features': FEATURES, 'seed':42,
        'model_type':'monotonic_tanh_neural_network', 'fit_performed':True,
        'status':'offline_candidate_requires_independent_validation',
        'train_count':int(mask.sum()), 'validation_count':int((~mask).sum()),
        'held_out_groups':sorted(held), 'neural_validation':neural, 'logistic_baseline_validation':baseline,
        'beats_baseline_brier':neural['brier'] < baseline['brier'],
        'assumptions':['证据/方向/暴露支持增大不得降低预测值，安全风险增大不得提高预测值；此为待验证建模假设。',
                       '未知特征填0.5并增加缺失指示，不作为阴性。',
                       '数据集、标签及规则来源需人工审查；本验证集不能循环用于选择规则。'],
        'weights':{'w1':w1.tolist(),'b1':b1.tolist(),'w2':w2.tolist(),'b2':float(b2)},
        'data_hash':hashlib.sha256(json.dumps(rows,sort_keys=True,ensure_ascii=False).encode()).hexdigest()}

def save_artifact(owner, model):
    ident=hashlib.sha256(json.dumps(model,sort_keys=True).encode()).hexdigest()[:24]
    owner_key=hashlib.sha256(str(owner).encode()).hexdigest()[:24]
    folder=Path(__file__).resolve().parent/'runs'/'scientific_models'/owner_key
    folder.mkdir(parents=True,exist_ok=True)
    path=folder/(ident+'.json')
    path.write_text(json.dumps(model,ensure_ascii=False,indent=2),encoding='utf-8')
    return {'id':ident,'model_card':{k:v for k,v in model.items() if k!='weights'}}
