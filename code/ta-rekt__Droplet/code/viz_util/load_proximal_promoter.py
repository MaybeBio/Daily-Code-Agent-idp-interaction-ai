"""
load_proximal_promoter.py

Configuration + loader for the proximal_promoter_transc_d_50 experiment.
Merges runs from 20260407 and 20260413 (16 runs each) at each concentration.

Paste this into a notebook cell and run it. Produces:
    on_times, off_times, bound_times, unbound_times,
    pattern_distances, transcription_activity, binding_activity,
    local_concentration   -- each a list of length k (parameters),
                             where each entry is a list of m numpy arrays.
    close_mask            -- boolean version of pattern_distances < 5,
                             same shape.
    params                -- list of parameter values (concentrations).
"""

import importlib
import viz.load_transcription as lt
importlib.reload(lt)

# Concentration values shared across both date folders
_concs = [
    0.003, 0.004, 0.005, 0.006, 0.007, 0.008,
    0.009, 0.010, 0.011, 0.012, 0.013, 0.014,
    0.015, 0.016, 0.017, 0.018,
]

_dates = ["data/20260407", "data/20260413"]
_exp   = "proximal_promoter_transc_d_50_c_{:.3f}"
_runs  = range(1, 17)   # output_cpu1.json ... output_cpu16.json

groups = lt.make_groups([
    (c, [("{}/{}".format(d, _exp.format(c)), _runs) for d in _dates])
    for c in _concs
])

(on_times, off_times, bound_times, unbound_times,
 pattern_distances, transcription_activity, binding_activity,
 local_concentration, run_ids) = lt.load(groups, n1=10, n=5000)

close_mask = lt.make_mask(pattern_distances, threshold=5.0)

params = _concs

print("Loaded {} parameter values, {} runs each.".format(
    len(params), len(on_times[0])))
