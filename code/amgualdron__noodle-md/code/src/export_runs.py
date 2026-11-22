"""./noodle export <run> [<run> ...]  --  the results of one or more runs, together, in exports/<name>/.

It first brings each run's analysis up to date (as ./noodle analyze), then writes:
    manifest.csv     every run's manifest.csv, one after another
    replicates.csv   every run's analysis/replicates.csv
    runs.csv         the runs' rows of runs/runs.csv (numbers, dates, counts, core-hours)
    physics/<run>.yaml   the physics sets each run used, fully written out
The tables start with a column `run`, so rows from different runs stay apart. Copy the
folder to your own computer (scp, rsync, Globus) and work from there; the trajectories
stay in the runs (or on the NAS).

<name> is the run's name for a single run. For several runs it asks, suggesting what
their names have in common. It warns about runs that are not fully analysed, and about a
physics set name that means different settings in different runs.
"""

import contextlib
import io
import json
import os
import re
import shutil
import sys

import yaml

from analyze_run import analyze
from runs import REGISTRY_COLUMNS, ROOT, find_run, list_runs, read_csv, registry, registry_row, run_here, write_csv

EXPORTS = ROOT / "exports"


def ask(question, default):
    try:
        answer = input(f"{question} [{default}] ").strip()
    except EOFError:
        answer = ""
    return answer or default


def runs_from_command_line():
    names = sys.argv[1:]
    if not names:
        here = run_here()
        if here:
            return [here]
        print("Usage: noodle export <run> [<run> ...]   (names or numbers; nothing inside a run's folder)\n")
        list_runs()
        sys.exit(1)
    run_dirs = []
    for name in names:
        run_dir = find_run(name)
        if run_dir is None:
            print(f"There is no run called '{name}'.\n")
            list_runs()
            sys.exit(1)
        if not (run_dir / "tasks.csv").exists():
            sys.exit(f"'{run_dir.name}' has not been launched yet: there is nothing to export.")
        if run_dir not in run_dirs:
            run_dirs.append(run_dir)
    return run_dirs


def export_name(run_dirs):
    """The export folder's name: the run's own, or (several runs) what you choose."""
    if len(run_dirs) == 1:
        return run_dirs[0].name
    common = os.path.commonprefix([run_dir.name for run_dir in run_dirs]).rstrip("_-.")
    name = ask(f"Name of the export folder (exports/<name>/) for these {len(run_dirs)} runs?",
               common or "_".join(run_dir.name for run_dir in run_dirs[:3]))
    if not re.fullmatch(r"[A-Za-z0-9._-]+", name):
        sys.exit("A name can only use letters, digits, '.', '_' and '-'.")
    return name


def physics_conflicts(run_dirs):
    """Physics set names that mean different settings in different runs."""
    seen, conflicts = {}, []
    for run_dir in run_dirs:
        for name, settings in yaml.safe_load((run_dir / "physics.yaml").read_text()).items():
            if name in seen and seen[name][1] != settings and name not in conflicts:
                conflicts.append(name)
                print(f"  [!] physics set '{name}' differs between {seen[name][0]} and {run_dir.name}: "
                      f"compare physics/{seen[name][0]}.yaml and physics/{run_dir.name}.yaml "
                      "(the condition columns of manifest.csv tell the rows apart)")
            seen.setdefault(name, (run_dir.name, settings))
    return conflicts


def main():
    run_dirs = runs_from_command_line()
    name = export_name(run_dirs)
    folder = EXPORTS / name
    if folder.exists():
        if ask(f"exports/{name}/ already exists. Replace it?", "y/N").lower() not in ("y", "yes"):
            print("Nothing was exported.")
            return

    tables = {"manifest.csv": [], "replicates.csv": []}
    for run_dir in run_dirs:
        with contextlib.redirect_stdout(io.StringIO()):  # its own report would be long; one line instead
            analyze(run_dir)
        counts = json.loads((run_dir / "analysis" / "report.json").read_text())["tasks"]
        line = f"  {run_dir.name}: {counts['analysed']} of {counts['total']} tasks analysed"
        if counts["analysed"] < counts["total"]:
            line += f"  [!] not finished ({counts['completed']} completed, {counts['failed']} failed)"
        print(line)
        for table, path in (("manifest.csv", run_dir / "manifest.csv"),
                            ("replicates.csv", run_dir / "analysis" / "replicates.csv")):
            if path.exists():
                tables[table] += [dict(run=run_dir.name, **row) for row in read_csv(path)]
    physics_conflicts(run_dirs)

    tmp = EXPORTS / f".{name}.tmp"  # written in full first, so a cancelled export leaves nothing half-done
    shutil.rmtree(tmp, ignore_errors=True)
    for table, rows in tables.items():
        write_csv(rows, tmp / table, ["run"])
    rows = registry()
    write_csv([registry_row(rows, run_dir.name) for run_dir in run_dirs], tmp / "runs.csv", REGISTRY_COLUMNS)
    (tmp / "physics").mkdir()
    for run_dir in run_dirs:
        shutil.copy(run_dir / "physics.yaml", tmp / "physics" / f"{run_dir.name}.yaml")
    shutil.rmtree(folder, ignore_errors=True)
    tmp.rename(folder)

    size = sum(p.stat().st_size for p in folder.rglob("*") if p.is_file())
    print(f"\nExported {len(run_dirs)} run(s) to exports/{name}/ ({size / 1e6:.1f} MB): manifest.csv "
          f"({len(tables['manifest.csv'])} rows), replicates.csv ({len(tables['replicates.csv'])}), "
          "runs.csv, physics/")
    print(f"To your own computer, e.g.:  rsync -av <you>@<cluster>:{folder}/ ./{name}/")


if __name__ == "__main__":
    main()
