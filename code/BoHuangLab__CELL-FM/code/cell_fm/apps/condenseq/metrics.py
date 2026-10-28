"""Curve metrics for the condensate titration analysis.

Port of cell_diff/tasks/analysis/utils.py from CELL-Diff2-Dev, restricted to the
three quantities this demo reports:

    AUC   area under the moving-averaged condensate-probability curve
    NR    the same curve with everything past its maximum held flat ("no reentrant")
    AAC   NR - AUC, i.e. the area *above* the curve that reentrant dissolution opens up

All three integrate against log10(intensity) by default and are normalised by the
width of the integration domain, so they live in [0, 1] and are comparable across
sequences.
"""

import numpy as np
import pandas as pd

# NumPy 2 renamed trapz to trapezoid; the old spelling still resolves there but raises a
# DeprecationWarning, and NumPy 1 has only trapz. Pick whichever exists, once, at import.
# The alternative was a notebook aliasing np.trapz onto the numpy module before importing
# this one, which works right up until something else imports numpy first.
_trapezoid = getattr(np, "trapezoid", None) or np.trapz


def moving_average(y, window: int) -> np.ndarray:
    """Centred rolling mean, shrinking the window at the edges rather than padding."""
    y = np.asarray(y, dtype=float)
    window = int(max(1, min(window, len(y))))
    return (
        pd.Series(y)
        .rolling(window=window, center=True, min_periods=1)
        .mean()
        .to_numpy()
    )


def integral(xs: np.ndarray, ys: np.ndarray) -> float:
    xs = np.asarray(xs, dtype=float)
    ys = np.asarray(ys, dtype=float)

    if xs.ndim != 1 or ys.ndim != 1:
        raise ValueError("xs and ys must be 1D arrays")
    if xs.shape[0] != ys.shape[0]:
        raise ValueError("xs and ys must have the same length")
    if xs.shape[0] < 2:
        return 0.0

    return float(_trapezoid(ys, xs))


def _prepare(df: pd.DataFrame, window_size: int, use_log_scale: bool):
    """Return (x, curve) sorted by x, with x log10-scaled when requested.

    Sorting happens once, before the moving average is taken, so the smoothing
    always runs along increasing intensity.
    """
    intensity_levels = df["protein_intensity_level"].to_numpy(dtype=float)
    predicted_classes = df["predicted_class"].to_numpy(dtype=float)

    order = np.argsort(intensity_levels)
    intensity_levels = intensity_levels[order]
    predicted_classes = predicted_classes[order]

    curve = moving_average(predicted_classes, window_size)

    if use_log_scale:
        # clip at 1.0 so that intensities below 1 do not send log10 negative/-inf
        x = np.log10(np.clip(intensity_levels, a_min=1.0, a_max=None))
    else:
        x = intensity_levels

    return x, curve


def _normalise(area: float, x: np.ndarray, normalize: bool) -> float:
    if not normalize:
        return area
    total_area = (x[-1] - x[0]) * 1.0
    return area / total_area


def compute_area_under_curve(
    df: pd.DataFrame, window_size: int, normalize: bool = True, use_log_scale: bool = True
) -> float:
    """AUC: how much of the intensity range this sequence spends in the condensed state."""
    x, curve = _prepare(df, window_size, use_log_scale)
    return _normalise(integral(x, curve), x, normalize)


def no_reentrant_curve(curve: np.ndarray) -> np.ndarray:
    """The curve it would have been without reentrant dissolution.

    Everything past the maximum is held at the maximum, so the difference against
    the real curve is exactly the dissolution that happens at high concentration.
    """
    max_index = int(np.argmax(curve))
    out = curve.copy()
    out[max_index:] = float(np.max(curve))
    return out


def compute_area_if_no_reentrant(
    df: pd.DataFrame, window_size: int, normalize: bool = True, use_log_scale: bool = True
) -> float:
    """Area under the curve after flattening everything past its maximum."""
    x, curve = _prepare(df, window_size, use_log_scale)
    return _normalise(integral(x, no_reentrant_curve(curve)), x, normalize)


def compute_area_if_no_reentrant_minus_auc(
    df: pd.DataFrame, window_size: int, normalize: bool = True, use_log_scale: bool = True
) -> float:
    """AAC, the area *above* the curve: no-reentrant area minus the real area.

    Zero when the curve never comes back down; larger the more the sequence
    dissolves again at high concentration.
    """
    x, curve = _prepare(df, window_size, use_log_scale)
    return _normalise(integral(x, no_reentrant_curve(curve) - curve), x, normalize)


# the paper's name for compute_area_if_no_reentrant_minus_auc
compute_area_above_curve = compute_area_if_no_reentrant_minus_auc


def compute_concentration(
    df: pd.DataFrame, window_size: int, threshold_prob: float = 0.8
) -> float:
    """c_sat: the first intensity at which the smoothed curve reaches threshold_prob.

    inf when the sequence never condenses over the scanned range.
    """
    intensity_levels = df["protein_intensity_level"].to_numpy(dtype=float)
    predicted_classes = df["predicted_class"].to_numpy(dtype=float)

    order = np.argsort(intensity_levels)
    intensity_levels = intensity_levels[order]
    curve = moving_average(predicted_classes[order], window_size)

    hits = np.flatnonzero(curve >= threshold_prob)
    return float(intensity_levels[hits[0]]) if hits.size else float("inf")
