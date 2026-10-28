"""The compute plan (design §8), from the pilot timings. Fits steps per second against chain length (power
law) on the pilot proteins' measured speed, projects every protein's single-thread wall time, and applies
the subset rule: proteins projected over 12 hours run 505 frames; if the projected wall time on 6 workers
still exceeds 24 hours, proteins above the longest length that fits are left out. Writes
results/simulations/plan.json, which 05_run.py and 07_analyse.py follow."""

import json
import sys

import pandas as pd
from common import OUT, RUNS
from simlib import compute_plan, fit_speed, makespan

PILOT = sys.argv[1:] or ["SASDJH6", "SASDJT7", "SASDLK9"]
WORKERS, BUDGET_H = 6, 24.0
inp = pd.read_csv(OUT / "inputs.csv")
prog = {
    lab: json.loads((RUNS / lab / "progress.json").read_text())
    for lab in inp["label"]
    if (RUNS / lab / "progress.json").exists()
}
pil = inp[inp["label"].isin(PILOT)].assign(sps=lambda d: [prog[x]["steps_per_second"] for x in d["label"]])
a, b = fit_speed(pil["length"], pil["sps"])
inp["steps_done"] = [prog[x]["steps"] if x in prog else 0 for x in inp["label"]]
plan = compute_plan(inp, a, b, WORKERS, BUDGET_H)
full = (inp["n_frames"] * inp["n_save"]) / (a * inp["length"] ** -b) / 3600
out = {
    "pilot": [
        {
            "label": r.label,
            "length": int(r.length),
            "frames": prog[r.label]["frames"],
            "steps": prog[r.label]["steps"],
            "wall_hours": round(prog[r.label]["wall_seconds"] / 3600, 3),
            "steps_per_second": r.sps,
            "projected_hours_full_run": round(r.n_frames * r.n_save / r.sps / 3600, 2),
        }
        for r in pil.itertuples()
    ],
    "speed_fit": {"a": a, "b": b, "form": "steps/s = a * N^-b"},
    "workers": WORKERS,
    "budget_hours": BUDGET_H,
    "projected_makespan_all_full_hours": round(makespan(full.tolist(), WORKERS), 1),
    "projected_core_hours_all_full": round(float(full.sum()), 1),
    "projected_makespan_plan_hours": round(
        makespan(plan.loc[plan["included"], "projected_hours"].tolist(), WORKERS), 1
    ),
    "n_included": int(plan["included"].sum()),
    "n_reduced_frames": int((plan["included"] & (plan["frames_planned"] < plan["n_frames"])).sum()),
    "excluded": plan.loc[~plan["included"], "label"].tolist(),
    "max_length_included": int(plan.loc[plan["included"], "length"].max()),
    "proteins": plan[
        [
            "label",
            "length",
            "included",
            "frames_planned",
            "discard_planned",
            "projected_hours_full",
            "projected_hours",
        ]
    ]
    .round(3)
    .to_dict(orient="records"),
}
(OUT / "plan.json").write_text(json.dumps(out, indent=1) + "\n")
print(json.dumps({k: v for k, v in out.items() if k != "proteins"}, indent=1))
