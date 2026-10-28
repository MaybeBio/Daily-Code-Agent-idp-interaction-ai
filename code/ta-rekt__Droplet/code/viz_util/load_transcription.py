"""
load_transcription.py

Load transcription time-series data from multiple independent runs,
potentially spread across different date folders, into k×m arrays.

Array shape: k (parameter values) × m (independent runs).
Each entry [k, m] is a numpy array of varying length (one run).

Usage (from notebook, run from Droplet/ directory):
    import importlib, viz.load_transcription as lt
    importlib.reload(lt)

    # Define groups: each entry is a list of (data_dir, run_indices) pairs
    # sharing the same parameter value.
    groups = lt.make_groups([
        # param_value -> list of (data_dir, run_indices)
        (0.003, [("data/20260407/proximal_promoter_transc_d_50_c_0.003", range(1, 17)),
                 ("data/20260413/proximal_promoter_transc_d_50_c_0.003", range(1, 17))]),
        (0.004, [("data/20260407/proximal_promoter_transc_d_50_c_0.004", range(1, 17)),
                 ("data/20260413/proximal_promoter_transc_d_50_c_0.004", range(1, 17))]),
        # ... add more parameter values as needed
    ])

    on_times, off_times, bound_times, unbound_times, pattern_distances, \
        transcription_activity, binding_activity, local_concentration, run_ids = lt.load(groups, n1=10, n=5000)

    # run_ids[k][m] == (data_dir, alpha_index) — same ordering as every observable.
    # Boolean index array: pattern_distances < 5, same shape as on_times
    close_mask = lt.make_mask(pattern_distances, threshold=5.0)
"""

import json
import numpy as np
from collections import OrderedDict


def make_groups(param_dir_pairs):
    """
    Build a list of (param_value, [(data_dir, run_indices), ...]) tuples.

    param_dir_pairs: iterable of
        (param_value, [(data_dir, run_indices), ...])
    """
    return list(param_dir_pairs)


def load(groups, n1=10, n=5000, keys=None):
    """
    Load data from all groups.

    groups: list of (param_value, [(data_dir, run_indices), ...])
    n1: number of on/off event entries to skip at the beginning of each run
    n: number of pattern_distance / activity entries to skip at the beginning

    Returns one array per observable. Each array has shape (k, m) where:
        k = number of parameter values
        m = total number of independent runs for that parameter value
    Each entry arr[k][m] is a 1-D numpy array.

    Run order is preserved across all observables: arr_a[k][m] and
    arr_b[k][m] always come from the same file and alpha_index.

    run_ids[k][m] is a (data_dir, alpha_index) tuple identifying each run.
    """
    if keys is None:
        keys = [    
            "on_times",
            "off_times",
            "bound_times",
            "unbound_times",
            "pattern_distances",
            "transcription_activity",
            "binding_activity",
            "local_concentration",
        ]
    # skips[key] = how many leading elements to drop
    skips = {k: n1 for k in keys[:4]}
    skips.update({k: n for k in keys[4:]})

    result = {k: [] for k in keys}
    run_ids = []   # run_ids[k] = list of (data_dir, alpha_index) in order

    for _param_val, dir_run_pairs in groups:
        # Collect all (data_dir, alpha_index) pairs for this parameter value
        # in a fixed, explicit order so every observable sees the same ordering.
        ordered_runs = [
            (data_dir, alpha_index)
            for data_dir, run_indices in dir_run_pairs
            for alpha_index in run_indices
        ]

        param_runs = {k: [] for k in keys}
        for data_dir, alpha_index in ordered_runs:
            path = "{}/output_cpu{}.json".format(data_dir, alpha_index)
            with open(path, "r") as f:
                res_ = json.load(f)
            res = res_[str(alpha_index)][str(1)]
            for k in keys:
                skip = skips[k]
                raw = res.get(k, [])
                param_runs[k].append(np.array(raw[skip:]))

        for k in keys:
            result[k].append(param_runs[k])
        run_ids.append(ordered_runs)

    return tuple(result[k] for k in keys) + (run_ids,)


def load_concat(groups, n1=10, n=5000):
    """
    Like load(), but traces sharing the same run index across dirs are
    concatenated along the time axis rather than treated as separate runs.

    groups: same format as load() — list of
        (param_value, [(data_dir, run_indices), ...])

    Runs are matched by index: all (data_dir, i) pairs with the same i are
    concatenated in the order they appear in dir_run_pairs.
    Each arr[k][m] is a 1-D numpy array of combined length.
    """
    keys = [
        "on_times", "off_times", "bound_times", "unbound_times",
        "pattern_distances", "transcription_activity", "binding_activity",
        "local_concentration",
    ]
    skips = {k: n1 for k in ["on_times", "off_times", "bound_times", "unbound_times"]}
    skips.update({k: n for k in ["pattern_distances", "transcription_activity",
                                  "binding_activity", "local_concentration"]})

    result = {k: [] for k in keys}
    run_ids = []

    for _param_val, dir_run_pairs in groups:
        # Build an ordered dict: run_index -> [(data_dir, run_index), ...]
        
        index_dirs = OrderedDict()
        for data_dir, run_indices in dir_run_pairs:
            for alpha_index in run_indices:
                index_dirs.setdefault(alpha_index, []).append(data_dir)

        param_runs = {k: [] for k in keys}
        for alpha_index, dirs in index_dirs.items():
            segments = {k: [] for k in keys}
            for data_dir in dirs:
                path = "{}/output_cpu{}.json".format(data_dir, alpha_index)
                with open(path, "r") as f:
                    res = json.load(f)[str(alpha_index)]["1"]
                for k in keys:
                    segments[k].append(np.array(res.get(k, [])[skips[k]:]))
            for k in keys:
                param_runs[k].append(np.concatenate(segments[k]))
        for k in keys:
            result[k].append(param_runs[k])
        run_ids.append([(dirs, i) for i, dirs in index_dirs.items()])

    return tuple(result[k] for k in keys) + (run_ids,)


def make_mask(pattern_distances, threshold=5.0):
    """
    Build a boolean mask array with the same shape as pattern_distances.

    Each entry mask[i][j] is a boolean numpy array: True where
    pattern_distances[i][j] < threshold.
    """
    mask = []
    for row in pattern_distances:
        mask.append([arr < threshold for arr in row])
    return mask
