"""./noodle analyze <run name>  --  pool the results of every finished task into tables.

Every job measures each of its tasks as soon as it finishes (observables.py; every column
is defined in OBSERVABLES.md): the averages go into the task's row of its job file. This
command gathers them, so it takes seconds and can be rerun at any time while a run is going:

    analysis/replicates.csv  one row per analysed task (replicate): every observable, with an
                             error bar from block averaging of its own frames
    analysis/report.json     how many tasks are analysed, the range of each observable per
                             physics set, and the rows worth a look (drift, too few samples)
    manifest.csv             its results columns: replicates, then every observable averaged
                             over the replicates (the columns written at launch stay as they are)

Error bars (<name>_err) in the manifest: the standard error of the mean over replicates
(their standard deviation / sqrt(replicates)) when there are 2 or more; with a single
replicate, its own block-average error.

The diagnostics, in the manifest: tau_rg_tu and tau_ree_tu are averaged over replicates,
independent_samples is their total (the independent Rg values behind the average, per
chain) and drift counts the replicates that were still drifting.

A completed task its job could not analyse (analysis_error in the job file) is analysed
here instead, from its trajectory, into analysis/reanalysed.csv.
"""

import json
import math
import statistics
from datetime import datetime

import yaml

from runs import (DIAGNOSTICS, OBSERVABLES, read_csv, read_tasks, refresh_registry, run_from_command_line,
                  run_lengths, task_states, trajectory_path, write_csv)

FEW_SAMPLES = 20  # fewer independent samples than this (all replicates together): too short to trust
SHOWN = 5         # rows named per warning on screen (report.json names them all)


# ── Gathering the per-task results ────────────────────────────────────────────
def reanalyse(run_dir, tasks, physics, to_do):
    """Analyse here the completed tasks in `to_do` that their job could not analyse, and
    add them to analysis/reanalysed.csv. Returns {task id: error}
    of the ones that failed again."""
    print(f"  {len(to_do)} finished task(s) were not analysed by their job: analysing them here. "
          "Loading HOOMD first can take a minute.", flush=True)
    from observables import analyse_task  # numpy, scipy and HOOMD: only when needed
    from simulate import system_geometry

    averages_path = run_dir / "analysis" / "reanalysed.csv"
    averages = read_csv(averages_path) if averages_path.exists() else []
    errors = {}
    for task_id in to_do:
        task = tasks[task_id]
        p, sequence = physics[task["physics"]], task["sequence"]
        trajectory = trajectory_path(run_dir, task)
        if not trajectory.exists():
            errors[task_id] = "no trajectory here (moved to storage?): copy it back into trajectories/"
            continue
        print(f"  analysing task {task_id} ({task['sequence_name']}, {task['physics']}) ...", flush=True)
        try:
            frame_time = run_lengths(p, len(sequence))[2] * p["dt"]
            result = analyse_task(trajectory, sequence, p, system_geometry(sequence, p), frame_time)
        except Exception as error:
            errors[task_id] = f"{type(error).__name__}: {error}"
            continue
        averages.append(dict(task_id=task_id, analysed_at=datetime.now().isoformat(timespec="seconds"), **result))
    if averages:
        write_csv(averages, averages_path)
    return errors


def gather(run_dir, tasks, states, physics):
    """(results, errors): every analysed task's averages {task id: {column: value}}, and
    {task id: error} of completed tasks that could not be analysed."""
    completed = [i for i, state in states.items() if state["status"] == "completed"]
    results = {i: states[i] for i in completed if states[i].get("analysed_at")}

    reanalysed_path = run_dir / "analysis" / "reanalysed.csv"
    done_here = {row["task_id"] for row in read_csv(reanalysed_path)} if reanalysed_path.exists() else set()
    to_do = [i for i in completed if i not in results and i not in done_here]
    errors = reanalyse(run_dir, tasks, physics, to_do) if to_do else {}
    if reanalysed_path.exists():
        results.update({row["task_id"]: row for row in read_csv(reanalysed_path) if row["task_id"] not in results})

    return results, errors


# ── Pooling the replicates ────────────────────────────────────────────────────
def replicate_rows(tasks, results):
    """analysis/replicates.csv: one row per analysed task, in task order."""
    rows = []
    for task_id in sorted(results):
        task, result = tasks[task_id], results[task_id]
        run_time = result.get("run_time_s")
        rows.append(dict(task_id=task_id, name=task["sequence_name"], physics=task["physics"],
                         replicate=task["replicate"], seed=task["seed"], N=task["length"], chains=result.get("chains"),
                         frames=result.get("frames"), run_time_h=round(run_time / 3600, 3) if run_time else None,
                         **{c: result.get(c) for name in OBSERVABLES for c in (name, f"{name}_err")},
                         **{c: result.get(c) for c in DIAGNOSTICS}))
    return rows


def pooled(replicates):
    """{column: value} of one protein and physics set from the rows of its replicates."""
    row = {"replicates": len(replicates)}
    for name in OBSERVABLES:
        values = [r[name] for r in replicates if r[name] is not None]
        if len(values) == 0:
            row[name] = row[f"{name}_err"] = None
        elif len(values) == 1:
            row[name], row[f"{name}_err"] = values[0], next(r[f"{name}_err"] for r in replicates if r[name] is not None)
        else:
            row[name], row[f"{name}_err"] = statistics.fmean(values), statistics.stdev(values) / math.sqrt(len(values))
    def values(column):
        return [r[column] for r in replicates if r[column] is not None and not math.isnan(r[column])]
    row["tau_rg_tu"] = statistics.fmean(values("tau_rg_tu")) if values("tau_rg_tu") else None
    row["tau_ree_tu"] = statistics.fmean(values("tau_ree_tu")) if values("tau_ree_tu") else None
    row["independent_samples"] = sum(values("independent_samples")) if values("independent_samples") else None
    row["drift"] = int(sum(values("drift")))
    return row


def add_reduced_density(manifest):
    """reduced_density of each brush row: grafting density x pi <Rg^2>, with the <Rg^2> of
    the same protein as a single chain at the same temperature and salt in this run (a bulk
    set if there is one, else a single grafted chain). Around 1 the chains start to touch:
    below it they are "mushrooms", well above it a brush. Empty without such a single chain."""
    single = {}
    for geometry in ("grafted", "bulk"):  # bulk last, so it wins
        for row in manifest:
            if row.get("chains") == 1 and row["geometry"] == geometry and row.get("rg_sq_nm2") is not None:
                single[row["name"], row["temperature_K"], row["ionic_concentration_M"]] = row["rg_sq_nm2"]
    for row in manifest:
        key = (row["name"], row["temperature_K"], row["ionic_concentration_M"])
        brush = (row.get("chains") or 1) > 1
        row["reduced_density"] = row["grafting_density"] * math.pi * single[key] if brush and key in single else None


def update_manifest(run_dir, replicates):
    """Rewrite manifest.csv: the columns written at launch, then the pooled results."""
    path = run_dir / "manifest.csv"
    manifest = read_csv(path)
    launch_columns = list(manifest[0])
    if "replicates" in launch_columns:  # results of an earlier ./noodle analyze: replaced
        launch_columns = launch_columns[:launch_columns.index("replicates")]
    groups = {}
    for r in replicates:
        groups.setdefault((r["name"], r["physics"]), []).append(r)
    rows = []
    for old in manifest:
        row = {c: old[c] for c in launch_columns}
        group = groups.get((row["name"], row["physics"]))
        row.update(pooled(group) if group else {"replicates": 0})
        rows.append(row)
    add_reduced_density(rows)
    columns = (launch_columns + ["replicates"]
               + [c for name in OBSERVABLES for c in (name, f"{name}_err")]
               + ["reduced_density"] + DIAGNOSTICS)
    write_csv(rows, path, columns)
    return rows


# ── The report ────────────────────────────────────────────────────────────────
def clean_number(value):
    """A number for JSON: 4 significant digits, NaN as null."""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return None
    return float(f"{value:.4g}")


def build_report(run_dir, tasks, states, results, errors, manifest):
    statuses = [state["status"] for state in states.values()]
    sets = {}
    for row in manifest:
        sets.setdefault(row["physics"], []).append(row)
    spread = {}
    for name, rows in sets.items():
        done = [r for r in rows if r["replicates"]]
        observables = {}
        for column in OBSERVABLES + ["reduced_density"]:
            values = [r[column] for r in done if r.get(column) is not None and not math.isnan(r[column])]
            if values:
                observables[column] = dict(min=clean_number(min(values)), median=clean_number(statistics.median(values)),
                                           max=clean_number(max(values)))
        spread[name] = dict(rows=len(rows), analysed_rows=len(done), observables=observables)
    label = lambda r: f"{r['name']} ({r['physics']})"
    return dict(
        run=run_dir.name, written_at=datetime.now().isoformat(timespec="seconds"),
        tasks=dict(total=len(tasks), completed=statuses.count("completed"), analysed=len(results),
                   failed=statuses.count("failed"), could_not_analyse=len(errors)),
        physics_sets=spread,
        look_at=dict(
            drift=[label(r) for r in manifest if r.get("drift")],
            few_independent_samples=[label(r) for r in manifest if r["replicates"]
                                     and (r.get("independent_samples") or 0) < FEW_SAMPLES],
            could_not_analyse=[dict(task_id=i, name=tasks[i]["sequence_name"], physics=tasks[i]["physics"], error=e)
                               for i, e in sorted(errors.items())]))


def print_report(report):
    t = report["tasks"]
    print(f"{t['analysed']} of {t['total']} tasks analysed ({t['completed']} completed"
          + (f", {t['failed']} failed" if t["failed"] else "") + ")")
    headline = {"rg_nm": "Rg", "com_height_nm": "centre-of-mass height", "brush_height_nm": "brush height"}
    for name, s in report["physics_sets"].items():
        parts = [f"{label} {o['min']:g} to {o['max']:g} nm (median {o['median']:g})"
                 for column, label in headline.items() if (o := s["observables"].get(column))]
        print(f"  {name:20s} {s['analysed_rows']} of {s['rows']} proteins   " + " · ".join(parts))
    look = report["look_at"]
    for key, text in (("drift", "were still drifting (column drift): equilibrate longer"),
                      ("few_independent_samples", f"have fewer than {FEW_SAMPLES} independent samples "
                                                  "(column independent_samples): run longer, or more replicates")):
        if look[key]:
            print(f"  [!] {len(look[key])} protein x physics row(s) {text}, e.g. {', '.join(look[key][:SHOWN])}")
    for e in look["could_not_analyse"][:SHOWN]:
        print(f"  [!] task {e['task_id']} ({e['name']}, {e['physics']}) could not be analysed: {e['error']}")


# ── The command ───────────────────────────────────────────────────────────────
def analyze(run_dir):
    """Pool every analysed task into analysis/ and manifest.csv. Returns the ids of the
    tasks whose results are saved (./noodle clean may delete their trajectories)."""
    tasks = {task["task_id"]: task for task in read_tasks(run_dir)}
    print(f"Gathering the results of {run_dir.name} ({len(tasks)} tasks) ...", flush=True)
    physics = yaml.safe_load((run_dir / "physics.yaml").read_text())
    states = task_states(run_dir)
    results, errors = gather(run_dir, tasks, states, physics)

    replicates = replicate_rows(tasks, results)
    write_csv(replicates, run_dir / "analysis" / "replicates.csv")
    manifest = update_manifest(run_dir, replicates)
    report = build_report(run_dir, tasks, states, results, errors, manifest)
    (run_dir / "analysis" / "report.json").write_text(json.dumps(report, indent=2) + "\n")

    refresh_registry(run_dir, states)
    print_report(report)
    print(f"  -> manifest.csv, analysis/replicates.csv, analysis/report.json")
    return set(results)


def main():
    analyze(run_from_command_line("noodle analyze <run>"))


if __name__ == "__main__":
    main()
