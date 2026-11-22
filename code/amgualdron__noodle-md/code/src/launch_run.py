"""./noodle launch <run name>  --  start (or continue) every unfinished task of a run.

The first time, it reads run.yaml + your CSV, checks the settings, shows the plan (how
many tasks, how many steps, roughly how many core-hours, how much disk) and asks before
doing anything. If you say yes, it freezes the run:
    tasks.csv      one row per task: protein x physics set x replicate, with its seed
    physics.yaml   the physics sets this run uses, fully written out
    manifest.csv   one row per protein and physics set: conditions and sequence descriptors
and starts the tasks. Later launches show what is done, what is still in the queue and
what will be (re)started, then start only the tasks that are neither completed nor in
the queue. Tasks that were stopped resume from their last checkpoint.

Tasks are grouped into jobs: a job runs on one core and works through its tasks one
after another. A run is split into at most max_jobs jobs (run.yaml) of about equal work,
grouped by chain length: the longest proteins get a job each, the shortest are pooled
many to a job, and a run with fewer tasks than max_jobs gets one job per task. The
grouping depends only on the tasks (their chain length and physics set), never on speed
estimates. The jobs with the longest proteins are numbered (and so started) first. On a cluster (when sbatch exists) the jobs are
submitted as a SLURM job array, one array element per job; on your own computer they
run one after another, right here.
"""

import csv
import difflib
import math
import random
import re
import shutil
import signal
import subprocess
import sys
from datetime import datetime
from statistics import median

import yaml

from runs import (RESIDUE_STEPS_PER_SECOND, ROOT, SEQUENCES, about, side_by_side, brush_size, checkpoint_path, compress,
                  conditions, file_size, job_numbers, read_submissions, read_tasks,
                  refresh_registry, run_from_command_line, run_lengths, slurm_queue, task_states,
                  update_registry, write_csv, write_job)

# The 20 standard amino acids: the same letters as the residue table in forcefield.py
# (repeated here so this script does not have to load HOOMD).
AMINO_ACIDS = set("ARNDCQEGHILKMFPSTWYV")
STOPPED_EXIT_CODE = 3  # simulate.py's exit code when a job stopped because it was asked to
MAX_JOBS = 1000        # max_jobs when run.yaml doesn't set it (one SLURM job array on most clusters)

RUN_SETTINGS = ["description", "sequences_csv", "name_column", "sequence_column", "rows", "max_jobs",
                "max_residues", "physics", "replicates", "time_limit_days", "max_running", "extra_sbatch"]
OPTIONAL = ["description", "max_residues", "max_jobs"]  # older run.yaml files may lack the last two
# The sets in physics.yaml that list every setting a bulk / grafted simulation needs
REFERENCE_SETS = {"bulk": "baseline", "grafted": "grafted"}
CHOICES = {"geometry": ["bulk", "grafted"], "graft_terminus": ["N", "C"],
           "wall_potential": ["lj93", "lj126"], "wall_attraction_scale": ["hydropathy", "uniform"]}

BYTES_PER_RESIDUE = 12      # trajectory size per residue per frame: one position (3 x 4 bytes)
BYTES_PER_FRAME = 100       # plus the step number and the file's index, per frame


# ── Checking the settings ─────────────────────────────────────────────────────
def did_you_mean(word, options):
    close = difflib.get_close_matches(word, options, n=1)
    return f" (did you mean '{close[0]}'?)" if close else ""


def stop_if(problems):
    if problems:
        print("Nothing was launched. Please fix:")
        for problem in problems:
            print(f"  - {problem}")
        sys.exit(1)


def check_run_settings(settings):
    """Clear mistakes in run.yaml: unknown or missing settings, impossible values."""
    problems = [f"run.yaml: unknown setting '{key}'{did_you_mean(key, RUN_SETTINGS)}"
                for key in settings if key not in RUN_SETTINGS]
    problems += [f"run.yaml: the setting '{key}' is missing"
                 for key in RUN_SETTINGS if key not in settings and key not in OPTIONAL]
    stop_if(problems)
    max_residues = settings.get("max_residues") or 0  # missing or empty: no limit
    if not isinstance(max_residues, int) or max_residues < 0:
        problems.append("run.yaml: max_residues must be a whole number (0 = no limit)")
    if not isinstance(settings["replicates"], int) or settings["replicates"] < 1:
        problems.append("run.yaml: replicates must be a whole number, 1 or more")
    if not isinstance(settings["time_limit_days"], (int, float)) or settings["time_limit_days"] <= 0:
        problems.append("run.yaml: time_limit_days must be a number of days, e.g. 3 or 0.5")
    max_jobs = settings.get("max_jobs", MAX_JOBS)
    if not isinstance(max_jobs, int) or isinstance(max_jobs, bool) or max_jobs < 0:
        problems.append("run.yaml: max_jobs must be a whole number, e.g. 1000 (0 = one job per task)")
    if not isinstance(settings["extra_sbatch"], list):
        problems.append('run.yaml: extra_sbatch must be a list, e.g. [] or ["--partition=normal"]')
    stop_if(problems)


def check_physics_set(name, physics, all_sets):
    """Clear mistakes in a physics set: misspelled or missing settings, impossible choices."""
    geometry = physics.get("geometry")
    if geometry not in CHOICES["geometry"]:
        return [f"physics set '{name}': geometry must be bulk or grafted, not '{geometry}'"]
    every_setting = resolve_physics(all_sets, REFERENCE_SETS["grafted"])
    needed = resolve_physics(all_sets, REFERENCE_SETS[geometry])
    problems = [f"physics set '{name}': unknown setting '{key}'{did_you_mean(key, every_setting)}"
                for key in physics if key not in every_setting]
    problems += [f"physics set '{name}': the setting '{key}' is missing"
                 for key in needed if key not in physics]
    for key, options in CHOICES.items():
        if key in physics and physics[key] not in options:
            problems.append(f"physics set '{name}': {key} must be {' or '.join(options)}, not '{physics[key]}'")
    if geometry == "grafted":
        density, chains = physics.get("grafting_density"), physics.get("chains")
        # above ~1 chain/nm^2 the grafting points are closer than two residues are wide
        if not isinstance(density, (int, float)) or not 0 <= density <= 1:
            problems.append(f"physics set '{name}': grafting_density must be a number of chains per nm^2 "
                            f"from 0 (single chain) to 1, not '{density}'")
        if not isinstance(chains, int) or chains < 1 or math.isqrt(chains) ** 2 != chains:
            problems.append(f"physics set '{name}': chains must be a square number (4, 9, 16, 25, ...), "
                            f"not '{chains}'")
    return problems


# ── Building the task list (first launch only) ────────────────────────────────
def select_rows(rows, setting):
    """The rows the `rows` setting picks, in file order (counted from 1):
    all | 7 | 101-200 | a list of these, e.g. [1-50, 75, 200-210]."""
    picked = set()
    for item in setting if isinstance(setting, list) else [setting]:
        text = str(item).strip()
        if text == "all":
            return rows
        match = re.fullmatch(r"(\d+)(?:\s*-\s*(\d+))?", text)
        if not match:
            sys.exit(f"rows: '{item}' is not understood. Use all, a row like 7, a range like 101-200, "
                     "or a list of rows and ranges like [1-50, 75, 200-210].")
        first = int(match.group(1))
        last = int(match.group(2) or first)
        if not 1 <= first <= last <= len(rows):
            sys.exit(f"rows: {text} is outside the file, which has {len(rows)} rows.")
        picked.update(range(first, last + 1))
    return [rows[i - 1] for i in sorted(picked)]


def rows_text(setting):
    """The rows setting as written in run.yaml: 1-50 or [1-50, 75]."""
    return f"[{', '.join(map(str, setting))}]" if isinstance(setting, list) else str(setting)


def sequence_file(run_dir, settings):
    """The CSV named by sequences_csv: in the run folder if it is there, else in sequences/.

    The shared sets in sequences/ are read in place (no copy needed): the sequences a run
    uses are written into its tasks.csv at the first launch, so changing or replacing a
    file in sequences/ later never changes a run that has been launched.
    """
    for folder in (run_dir, SEQUENCES):
        if (folder / settings["sequences_csv"]).is_file():
            return folder / settings["sequences_csv"]
    shared = ", ".join(sorted(p.name for p in SEQUENCES.glob("*.csv")))
    sys.exit(f"Can't find the sequence file '{settings['sequences_csv']}' in runs/{run_dir.name}/ "
             f"or in sequences/ ({shared or 'empty'}). Copy it into runs/{run_dir.name}/ "
             "or fix sequences_csv in run.yaml.")


def read_sequences(run_dir, settings):
    """({name: sequence}, [names skipped], [names too long], {name: row}) for the selected
    rows of the run's CSV. Too long: more residues than max_residues (if that is set and
    not 0). row: the protein's row in the CSV, counted from 1 like the rows setting."""
    csv_path = sequence_file(run_dir, settings)
    with open(csv_path, newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    name_col, seq_col = settings["name_column"], settings["sequence_column"]
    for column in (name_col, seq_col):
        if column not in rows[0]:
            sys.exit(f"{csv_path.name} has no column '{column}'. Its columns are: {list(rows[0])}")

    sequences, skipped, too_long, row_numbers = {}, [], [], {}
    max_residues = settings.get("max_residues") or 0
    for number, row in enumerate(rows, 1):
        row["_row"] = number
    for row in select_rows(rows, settings["rows"]):
        # names end up in file names, so keep them to safe characters
        name = re.sub(r"[^A-Za-z0-9._-]+", "_", row[name_col].strip())
        sequence = row[seq_col].strip().upper()
        if name in sequences:
            sys.exit(f"The name '{name}' appears twice in {csv_path.name}; names must be unique.")
        if not sequence or set(sequence) - AMINO_ACIDS:
            skipped.append(name)
            continue
        if max_residues and len(sequence) > max_residues:
            too_long.append(name)
            continue
        sequences[name] = sequence
        row_numbers[name] = row["_row"]
    return sequences, skipped, too_long, row_numbers


def physics_sets(run_dir):
    """Every physics set a run can use: those of physics.yaml, plus any the run folder
    defines for itself in physics_sets.yaml (same format; they may extend the shared sets),
    so sets that only one run needs don't crowd physics.yaml."""
    all_sets = yaml.safe_load((ROOT / "physics.yaml").read_text())
    own = run_dir / "physics_sets.yaml"
    if own.exists():
        local = yaml.safe_load(own.read_text()) or {}
        clash = sorted(set(local) & set(all_sets))
        if clash:
            sys.exit(f"runs/{run_dir.name}/physics_sets.yaml redefines {', '.join(clash)}, already in "
                     "physics.yaml: give its sets other names.")
        all_sets.update(local)
    return all_sets


def resolve_physics(all_sets, name):
    """A physics set with everything it `extends` filled in underneath it."""
    if name not in all_sets:
        sys.exit(f"physics.yaml has no set called '{name}'. Its sets are: {list(all_sets)}")
    settings = dict(all_sets[name])
    parent = settings.pop("extends", None)
    return {**resolve_physics(all_sets, parent), **settings} if parent else settings


def plan_tasks(run_dir, settings):
    """Everything the first launch would create, checked, but not written yet.

    Returns (sequences, skipped, too_long, row_numbers, physics, tasks): the proteins, the
    names of the rows skipped for unknown letters and for being longer than max_residues,
    each protein's row in the CSV, the physics sets {name: settings}, and the task list."""
    sequences, skipped, too_long, row_numbers = read_sequences(run_dir, settings)
    if not sequences:
        sys.exit("No usable sequences were selected.")
    physics_names = settings["physics"] if isinstance(settings["physics"], list) else [settings["physics"]]
    all_sets = physics_sets(run_dir)
    physics = {name: resolve_physics(all_sets, name) for name in physics_names}

    problems = []
    for name, settings_of_set in physics.items():
        problems += check_physics_set(name, settings_of_set, all_sets)
    stop_if(problems)
    for name, settings_of_set in physics.items():
        try:
            run_lengths(settings_of_set, min(map(len, sequences.values())))
        except ValueError as error:
            problems.append(f"physics set '{name}': {error}")
    stop_if(problems)

    tasks = []
    for name, sequence in sequences.items():
        for physics_name in physics_names:
            # One seed per replicate, all different: two replicates with the same seed would
            # be the same simulation. HOOMD seeds are 16-bit numbers.
            seeds = random.sample(range(1, 65536), settings["replicates"])
            for replicate, seed in enumerate(seeds, 1):
                tasks.append(dict(task_id=len(tasks) + 1, sequence_name=name, physics=physics_name,
                                  replicate=replicate, seed=seed, length=len(sequence), sequence=sequence))
    return sequences, skipped, too_long, row_numbers, physics, tasks


def write_tasks(run_dir, settings, sequences, row_numbers, physics, tasks):
    """Freeze the run: write its physics sets, its task list and the start of its manifest."""
    (run_dir / "physics.yaml").write_text(
        "# The physics sets this run uses, copied from physics.yaml when it was launched.\n"
        + yaml.safe_dump(physics, sort_keys=False))
    with open(run_dir / "tasks.csv", "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(tasks[0]))
        writer.writeheader()
        writer.writerows(tasks)
    write_manifest(run_dir, settings, sequences, row_numbers, physics)


def write_manifest(run_dir, settings, sequences, row_numbers, physics):
    """manifest.csv: one row per protein and physics set, with everything known before
    the simulations (name, sequence, where it came from, the physics conditions, the
    sequence descriptors). ./noodle analyze adds the simulation averages to the same rows."""
    from descriptors import descriptors

    print(f"Computing the sequence descriptors of {len(sequences)} proteins for manifest.csv ...", flush=True)
    described = {name: descriptors(sequence) for name, sequence in sequences.items()}
    rows = []
    for name, sequence in sequences.items():
        for physics_name, p in physics.items():
            rows.append(dict(name=name, physics=physics_name, sequence=sequence,
                             source_csv=settings["sequences_csv"], source_row=row_numbers[name],
                             **conditions(p), **described[name]))
    write_csv(rows, run_dir / "manifest.csv")


# ── The plan shown before the first launch ────────────────────────────────────
def describe_physics(p):
    """One line saying what a physics set simulates."""
    if p["geometry"] == "bulk":
        where = "free chain in solution"
    else:
        attraction = ("repulsive" if p["wall_attraction"] == 0 else
                      f"attractive {p['wall_attraction']:g} ({p['wall_attraction_scale']})")
        charge = "neutral" if p["wall_charge_density"] == 0 else f"charge {p['wall_charge_density']:+g} e/nm^2"
        where = f"grafted by the {p['graft_terminus']}-terminus to a {attraction}, {charge} {p['wall_potential']} surface"
        chains, width = brush_size(p)
        if chains > 1:
            where = f"brush: {chains} chains at {p['grafting_density']:g}/nm^2 (box {width:.1f} nm), " + where
    return f"{where}; {p['temperature_K']:g} K, {1000 * p['ionic_concentration_M']:g} mM salt"


def surface_potential(p):
    """e psi / kT at a charged surface: 4 pi l_B sigma / kappa (linearised Debye-Huckel,
    salt water on one side, as in forcefield.charged_surface)."""
    from forcefield import debye_huckel

    bjerrum, kappa = debye_huckel(p["temperature_K"], p["ionic_concentration_M"])
    return 4 * math.pi * bjerrum * p["wall_charge_density"] / kappa


def in_words(number):
    """1_800_000 -> '1.8 million'."""
    for size, word in ((1e9, "billion"), (1e6, "million"), (1e3, "thousand")):
        if number >= size:
            return f"{number / size:.3g} {word}"
    return str(number)


def task_cost(physics, n_residues):
    """Rough run time of one task on one core, in seconds: residues x steps / core speed."""
    chains, _ = brush_size(physics)
    return chains * n_residues * sum(run_lengths(physics, n_residues)[:2]) / RESIDUE_STEPS_PER_SECOND


def trajectory_bytes(physics, n_residues):
    """Size of one task's trajectory: every residue of every chain, in every frame."""
    chains, _ = brush_size(physics)
    return int(physics["frames"]) * (BYTES_PER_RESIDUE * chains * n_residues + BYTES_PER_FRAME)


def show_plan(run_dir, settings, sequences, skipped, too_long, physics, tasks, groups, on_cluster):
    lengths = sorted(map(len, sequences.values()))
    print(f"New run: {run_dir.name}" + (f'   "{settings["description"]}"' if settings.get("description") else ""))
    print()
    sizes = (f"{lengths[0]} residues" if lengths[0] == lengths[-1] else
             f"{lengths[0]} to {lengths[-1]} residues (median {median(lengths):.0f})")
    csv_path = sequence_file(run_dir, settings).relative_to(ROOT)
    print(f"Sequences  {csv_path}, rows {rows_text(settings['rows'])}: {len(sequences)} proteins, {sizes}")
    if skipped:
        print(f"           skipped {len(skipped)} with letters other than the 20 standard amino acids: "
              + ", ".join(skipped[:5]) + (", ..." if len(skipped) > 5 else ""))
    if too_long:
        print(f"           skipped {len(too_long)} longer than max_residues ({settings['max_residues']} residues)")

    print("Physics")
    for name, p in physics.items():
        print(f"  {name:20s} {describe_physics(p)}")
        if p["geometry"] == "grafted" and p["wall_charge_density"] != 0:
            psi = surface_potential(p)
            print(f"  {'':20s} surface potential {psi:+.2f} kT/e"
                  + ("" if abs(psi) <= 1 else ": [!] above ~1 kT/e the linearised (Debye-Huckel) "
                     "surface overestimates its pull; consider a smaller wall_charge_density or more salt"))
        if name.startswith("test_"):
            print(f"  {'':20s} [!] a quick-test set: tiny run length, results are NOT physically meaningful")
    print(f"Tasks      {len(sequences)} proteins x {len(physics)} physics set(s) x "
          f"{settings['replicates']} replicate(s) = {len(tasks)} tasks")

    print("Length     steps per task (equilibration + production), for the shortest and longest protein")
    for name, p in physics.items():
        shortest, longest = (sum(run_lengths(p, n)[:2]) for n in (lengths[0], lengths[-1]))
        steps = (in_words(shortest) if shortest == longest else
                 f"{in_words(shortest)} (N={lengths[0]}) to {in_words(longest)} (N={lengths[-1]})")
        print(f"  {name:20s} {steps}, saving {p['frames']} frames")

    # Only a rough cost: the speed of the cluster's cores is not known here. ./noodle
    # status gives real estimates, from the speeds the run's own tasks measure.
    costs = [task_cost(physics[t["physics"]], t["length"]) for t in tasks]
    print(f"Cost       ~{sum(costs) / 3600:,.0f} core-hours, the longest task {about(max(costs))} "
          f"(rough: {RESIDUE_STEPS_PER_SECOND / 1e6:g} million residue-steps per second per core)")
    limit = settings["time_limit_days"] * 86400
    if on_cluster and max(costs) > limit:
        over = sum(c > limit for c in costs)
        print(f"           [!] {over} task(s) need more than time_limit_days: each stops at the limit and "
              f"continues at the next ./noodle launch ({math.ceil(max(costs) / limit)} launches for the "
              "longest). Consider max_residues, or a longer time limit if your cluster allows it")

    n_bytes = sum(trajectory_bytes(physics[t["physics"]], t["length"]) for t in tasks)
    print(f"Disk       ~{file_size(n_bytes)} of trajectories")

    if on_cluster:
        at_once = f"at most {settings['max_running']} at a time" if settings["max_running"] else "no limit on how many at a time"
        extra = f", extra sbatch options: {' '.join(settings['extra_sbatch'])}" if settings["extra_sbatch"] else ""
        print(f"Cluster    {describe_jobs(groups, tasks, settings)}")
        print(f"           time limit {settings['time_limit_days']:g} day(s) per job, {at_once}{extra}")
        # (a single task longer than the time limit is counted above; here only shared jobs)
        cost = {t["task_id"]: task_cost(physics[t["physics"]], t["length"]) for t in tasks}
        wait = side_by_side([sum(cost[i] for i in g) for g in groups], settings["max_running"])
        print(f"           {about(wait)} until all are done if "
              + (f"{settings['max_running']} jobs run at a time" if settings["max_running"] and
                 settings["max_running"] < len(groups) else "they all run at once")
              + ", plus any wait in the SLURM queue (rough, as the cost)")
        over = [g for g in groups if len(g) > 1 and sum(cost[i] for i in g) > settings["time_limit_days"] * 86400]
        if over:
            print(f"           [!] {len(over)} job(s) of several tasks may need more than the time limit (rough "
                  "estimate). They stop at the limit and continue at the next launch; a higher max_jobs "
                  "splits the work finer")
    else:
        print("Where      no SLURM here: the tasks run in this terminal, one after another "
              "(Ctrl-C stops, ./noodle launch continues)")
    print()
    print("Once launched, the run is frozen: only the cluster settings (max_jobs, time_limit_days, max_running, "
          "extra_sbatch) can still change.")


def confirm(question):
    """Ask a yes/no question; anything but y or yes means no."""
    try:
        return input(f"{question} [y/N] ").strip().lower() in ("y", "yes")
    except EOFError:  # no keyboard (e.g. input from a pipe): treat as no
        print()
        return False


# ── Grouping tasks into jobs ──────────────────────────────────────────────────
def task_work(task, state, physics):
    """What a task still has to compute: residues x steps left. Known exactly from its chain
    length and physics set (and, for a task that ran before, how far it got), so grouping by
    it never depends on how fast a computer is."""
    total = state.get("total_steps") or sum(run_lengths(physics, task["length"])[:2])
    return brush_size(physics)[0] * task["length"] * (total - (state.get("step") or 0))


def group_into_jobs(work, max_jobs):
    """[[task ids of job 1], [task ids of job 2], ...] from {task id: work}: at most max_jobs
    jobs (0 = no limit) of about equal work.

    The tasks are sorted from the most work to the least (for one physics set: from the
    longest chain to the shortest) and cut into consecutive groups. Each job gets about the
    work still left divided by the jobs still available; a task with more than that gets a
    job of its own. So the longest proteins run alone, the shortest are pooled many to a
    job, and every job holds proteins of similar length. With no more tasks than max_jobs,
    every task gets its own job. The same tasks always give the same jobs."""
    order = sorted(work, key=lambda i: (-work[i], i))
    if not max_jobs or len(order) <= max_jobs:
        return [[i] for i in order]
    groups, left, start = [], sum(work.values()), 0
    while start < len(order):
        share = left / (max_jobs - len(groups))
        group, total = [order[start]], work[order[start]]
        # take the next task while that brings the job closer to its share
        while start + len(group) < len(order) and total + work[order[start + len(group)]] / 2 <= share:
            total += work[order[start + len(group)]]
            group.append(order[start + len(group)])
        if len(groups) == max_jobs - 1:  # the last job available takes whatever is left
            group = order[start:]
        groups.append(group)
        left -= sum(work[i] for i in group)
        start += len(group)
    return groups


def describe_jobs(groups, tasks, settings):
    """E.g. '12000 tasks in 1000 SLURM jobs (max_jobs), grouped by chain length: from 1 task
    per job (N = 380-599) to 41 (N = 20-22)'."""
    length = {t["task_id"]: t["length"] for t in tasks}
    n_tasks = sum(map(len, groups))
    text = f"{n_tasks} task(s) in {len(groups)} SLURM job(s)"
    if len(groups) == n_tasks:
        return text + " (one task per job)"

    def lengths(group):
        low, high = min(length[i] for i in group), max(length[i] for i in group)
        return f"N = {low}" if low == high else f"N = {low}-{high}"
    alone = [g for g in groups if len(g) == 1]
    first = (f"1 task per job ({lengths(sum(alone, []))})" if alone
             else f"{len(groups[0])} per job ({lengths(groups[0])})")
    return (text + f" (max_jobs: {settings.get('max_jobs', MAX_JOBS)}), grouped by chain length: from {first} "
            f"to {len(groups[-1])} tasks per job ({lengths(groups[-1])})")


def write_jobs(run_dir, groups, tasks, states):
    """Write a job file for every group of tasks (numbered after the run's existing jobs)
    and return the new job numbers. A task's earlier progress (attempts, run time, steps)
    is carried over into its new row."""
    by_id = {task["task_id"]: task for task in tasks}
    first = max(job_numbers(run_dir), default=0) + 1
    carried = ["chains", "phase", "step", "total_steps", "steps_per_second", "run_time_s", "attempts",
               "host", "started_at", "last_update"]
    for job, task_ids in enumerate(groups, first):
        rows = []
        for task_id in task_ids:
            task, state = by_id[task_id], states.get(task_id, {})
            rows.append(dict(task_id=task_id, sequence_name=task["sequence_name"], physics=task["physics"],
                             replicate=task["replicate"], N=task["length"], status="submitted",
                             **{key: state[key] for key in carried if state.get(key) is not None}))
        write_job(run_dir, job, rows)
    return list(range(first, first + len(groups)))


# ── Deciding what to start (later launches) ───────────────────────────────────
def plan_relaunch(run_dir, on_cluster):
    """Sort the tasks into done / still in the queue / failed / to start, say so, and
    return the ids of the tasks to start (failed ones only if you agree)."""
    tasks = read_tasks(run_dir)
    names = {t["task_id"]: t["sequence_name"] for t in tasks}
    states = task_states(run_dir)
    done = [i for i in names if states.get(i, {}).get("status") == "completed"]
    if len(done) == len(tasks):
        print(f"All {len(tasks)} tasks of '{run_dir.name}' are completed.")
        return []

    unfinished = [i for i in names if i not in done]
    in_queue = []
    if on_cluster:  # in a job that is still waiting or running: left alone
        queue, submitted = slurm_queue(run_dir.name) or {}, read_submissions(run_dir)
        in_queue = [i for i in unfinished if i in states and submitted.get(states[i]["job"]) in queue]
    failed = [i for i in unfinished if i not in in_queue and states.get(i, {}).get("status") == "failed"]
    to_start = [i for i in unfinished if i not in in_queue and i not in failed]
    resuming = [i for i in to_start if checkpoint_path(run_dir, i).exists()]

    parts = [f"{len(done)} done"]
    if in_queue:
        parts.append(f"{len(in_queue)} in the queue (left alone)")
    if failed:
        parts.append(f"{len(failed)} failed")
    parts.append(f"{len(to_start)} to start ({len(resuming)} resume from a checkpoint, "
                 f"{len(to_start) - len(resuming)} from the beginning)")
    print(f"{len(tasks)} tasks: " + " · ".join(parts))

    if failed:
        print(f"\n{len(failed)} task(s) failed with an error last time:")
        for i in failed[:10]:
            print(f"  task {i:>4}  {names[i]}: {states[i]['error']}")
        if len(failed) > 10:
            print(f"  ... and {len(failed) - 10} more (./noodle status {run_dir.name} lists them all)")
        print(f"They will most likely fail again unless the cause has been fixed (full messages in "
              f"runs/{run_dir.name}/jobs/job_<number>.out; the job of each task is in log.csv).")
        if confirm("Start the failed tasks again anyway?"):
            to_start = sorted(to_start + failed)
    return to_start


# ── Starting the jobs ─────────────────────────────────────────────────────────
def max_array_size():
    out = subprocess.run(["scontrol", "show", "config"], capture_output=True, text=True).stdout
    match = re.search(r"MaxArraySize\s*=\s*(\d+)", out)
    return int(match.group(1)) if match else 1001


def submit_to_slurm(run_dir, settings, jobs):
    """Submit the jobs (numbers) as SLURM job arrays and record them in jobs/submitted.csv."""
    name = run_dir.name
    time_limit_minutes = max(1, round(settings["time_limit_days"] * 24 * 60))  # 0 would mean no limit
    slurm_ids = {}  # job number -> SLURM job id

    # Array indices must stay below the cluster's MaxArraySize, so many jobs are sent as
    # several arrays, each with an offset: job number = offset + array index.
    chunk = max_array_size() - 1
    for offset in range(0, max(jobs), chunk):
        numbers = [j for j in jobs if offset < j <= offset + chunk]
        if not numbers:
            continue
        array = compress([j - offset for j in numbers])
        if settings["max_running"]:
            array += f"%{settings['max_running']}"
        # Each job writes its own jobs/job_0007.out (see slurm_job.sh); SLURM's --output only
        # gets SLURM's own messages (e.g. about the time limit), all in jobs/slurm.out.
        command = ["sbatch", "--parsable", f"--job-name={name}", f"--array={array}",
                   f"--time={time_limit_minutes}", f"--chdir={ROOT}", "--open-mode=append",
                   f"--output={run_dir / 'jobs' / 'slurm.out'}", *settings["extra_sbatch"],
                   str(ROOT / "src" / "slurm_job.sh"), f"runs/{name}", str(offset)]
        result = subprocess.run(command, capture_output=True, text=True)
        if result.returncode != 0:
            sys.exit(f"sbatch refused the job:\n{result.stderr}"
                     + (f"\n({len(slurm_ids)} job(s) had already been submitted.)" if slurm_ids else ""))
        array_id = result.stdout.strip().split(";")[0]
        slurm_ids.update({j: f"{array_id}_{j - offset}" for j in numbers})

        # Record which SLURM job runs which job, so ./noodle status can find it in the queue.
        submissions = run_dir / "jobs" / "submitted.csv"
        is_new_file = not submissions.exists()
        with open(submissions, "a", newline="") as f:
            writer = csv.writer(f)
            if is_new_file:
                writer.writerow(["job", "slurm_job", "submitted_at"])
            submitted_at = datetime.now().isoformat(timespec="seconds")
            writer.writerows([j, slurm_ids[j], submitted_at] for j in numbers)
        print(f"Submitted {len(numbers)} job(s) as SLURM job array {array_id}.")

    # Check that the jobs really are in the queue
    queue = slurm_queue(name) or {}
    found = sum(slurm_id in queue for slurm_id in slurm_ids.values())
    if found == len(slurm_ids):
        print(f"Confirmed: {found} of {len(slurm_ids)} job(s) are in the SLURM queue.")
    else:
        print(f"Only {found} of {len(slurm_ids)} job(s) are in the SLURM queue; the others may have "
              f"finished or failed already: see ./noodle status {name}")
    print(f"Check on them any time with:  ./noodle status {name}")


def run_jobs_here(run_dir, jobs):
    print(f"No SLURM here: running {len(jobs)} job(s) one after another. "
          f"Ctrl-C stops; ./noodle launch {run_dir.name} again resumes.\n")
    # Ctrl-C reaches the job too, which then saves a checkpoint and stops cleanly; this
    # program ignores it and simply starts no further jobs.
    signal.signal(signal.SIGINT, signal.SIG_IGN)
    for job in jobs:
        result = subprocess.run([sys.executable, str(ROOT / "src" / "simulate.py"), str(run_dir), str(job)])
        print()
        if result.returncode == STOPPED_EXIT_CODE:
            print(f"Stopped. ./noodle launch {run_dir.name} continues where it left off.")
            return
        if result.returncode != 0:
            print(f"Job {job} ended with an error (above).")
    states = task_states(run_dir)
    ran = [task_id for task_id, state in states.items() if state["job"] in jobs]
    failed = [i for i in ran if states[i]["status"] == "failed"]
    print(f"Finished: {len(ran) - len(failed)} of {len(ran)} task(s) completed"
          + (f", {len(failed)} failed (tasks {compress(sorted(failed))})." if failed else "."))


def work_left(run_dir, task_ids, tasks, physics):
    """{task id: residues x steps it still has to compute}."""
    states = task_states(run_dir) if (run_dir / "jobs").exists() else {}
    by_id = {task["task_id"]: task for task in tasks}
    return {i: task_work(by_id[i], states.get(i, {}), physics[by_id[i]["physics"]]) for i in task_ids}


def main():
    run_dir = run_from_command_line("noodle launch <run>")
    settings = yaml.safe_load((run_dir / "run.yaml").read_text())
    check_run_settings(settings)
    on_cluster = shutil.which("sbatch") is not None

    if not (run_dir / "tasks.csv").exists():
        sequences, skipped, too_long, row_numbers, physics, tasks = plan_tasks(run_dir, settings)
        to_start = [t["task_id"] for t in tasks]
        groups = group_into_jobs(work_left(run_dir, to_start, tasks, physics), settings.get("max_jobs", MAX_JOBS))
        show_plan(run_dir, settings, sequences, skipped, too_long, physics, tasks, groups, on_cluster)
        if not confirm(f"Launch these {len(tasks)} tasks?"):
            print("Nothing was launched.")
            return
        write_tasks(run_dir, settings, sequences, row_numbers, physics, tasks)
        update_registry(run_dir, launched=datetime.now().isoformat(timespec="seconds"))
    else:
        if (run_dir / "run.yaml").stat().st_mtime > (run_dir / "tasks.csv").stat().st_mtime:
            print("Note: run.yaml was edited after this run was launched. Only changes to the cluster settings "
                  "(max_jobs, time_limit_days, max_running, extra_sbatch) take effect; for anything else make "
                  "a new run.\n")
        to_start = plan_relaunch(run_dir, on_cluster)
        if not to_start:
            return
        tasks = read_tasks(run_dir)
        physics = yaml.safe_load((run_dir / "physics.yaml").read_text())
        groups = group_into_jobs(work_left(run_dir, to_start, tasks, physics), settings.get("max_jobs", MAX_JOBS))
        print(f"Starting {describe_jobs(groups, tasks, settings)}.")

    jobs = write_jobs(run_dir, groups, read_tasks(run_dir), task_states(run_dir))
    refresh_registry(run_dir)
    if on_cluster:
        submit_to_slurm(run_dir, settings, jobs)
    else:
        run_jobs_here(run_dir, jobs)


if __name__ == "__main__":
    main()
