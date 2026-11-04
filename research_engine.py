"""Temporary recovery loader for the compiled research engine backup."""
from pathlib import Path as _Path
import importlib.util as _importlib_util
import marshal as _marshal

_backup = _Path(__file__).with_name("research_engine.recovered.pyc")
_payload = _backup.read_bytes()
if _payload[:4] != _importlib_util.MAGIC_NUMBER:
    raise RuntimeError("恢复的研究引擎缓存与当前Python版本不匹配。")
_compiled_engine = _marshal.loads(_payload[16:])
exec(_compiled_engine, globals())

from research_reliability import apply_reliability_patches as _apply_reliability_patches
_apply_reliability_patches(__import__(__name__))

# PyDESeq2 0.5.2 fails with pandas 3.x; the compatible release is 0.5.4+.
_compiled_count_analysis = differential_expression_counts
def differential_expression_counts(*args, **kwargs):
    from importlib.metadata import version as _package_version
    import pandas as _pd
    import numpy as _np
    _matrix = args[0] if args else kwargs.get("matrix")
    _values = _np.asarray(_matrix)
    if _values.ndim != 2 or _np.any(~_np.isfinite(_values)) or _np.any(_values < 0) or _np.any(_np.abs(_values-_np.rint(_values)) > 1e-8):
        return _compiled_count_analysis(*args, **kwargs)
    _version = tuple(int(part) for part in _package_version("pydeseq2").split(".")[:3])
    _pandas_major = int(_pd.__version__.split(".", 1)[0])
    if _pandas_major >= 3 and _version < (0, 5, 4):
        raise ImportError("当前PyDESeq2与pandas版本不兼容：pandas 3需PyDESeq2 0.5.4或更新版本，请重新安装requirements.txt。")
    return _compiled_count_analysis(*args, **kwargs)

# Add per-gene standard errors so independently analyzed cohorts can be
# synthesized without pooling their raw matrices.
_compiled_differential_expression = differential_expression
def differential_expression(matrix, genes, case, control, pairs=None):
    import numpy as _np
    result = _compiled_differential_expression(matrix, genes, case, control, pairs=pairs)
    values = _np.asarray(matrix, dtype=float)
    if values.ndim != 2 or not isinstance(result, dict):
        return result
    case_idx = _np.asarray(case, dtype=int)
    control_idx = _np.asarray(control, dtype=int)
    if pairs:
        pair_values = _np.asarray(pairs, dtype=int).reshape(len(pairs), 2)
        differences = values[:, pair_values[:, 0]] - values[:, pair_values[:, 1]]
        se_values = _np.std(differences, axis=1, ddof=1) / _np.sqrt(differences.shape[1]) if differences.shape[1] > 1 else _np.full(values.shape[0], _np.nan)
    else:
        case_values, control_values = values[:, case_idx], values[:, control_idx]
        se_values = _np.sqrt(_np.var(case_values, axis=1, ddof=1) / case_values.shape[1] +
                             _np.var(control_values, axis=1, ddof=1) / control_values.shape[1])
    positions = {str(gene).upper(): i for i, gene in enumerate(genes)}
    enriched = []
    for row in result.get("results", []):
        item = dict(row)
        pos = positions.get(str(item.get("gene", "")).upper())
        if pos is not None and _np.isfinite(se_values[pos]):
            item["SE"] = float(se_values[pos])
        if "p" in item and "pvalue" not in item:
            item["pvalue"] = item["p"]
        enriched.append(item)
    return {**result, "results": enriched}
