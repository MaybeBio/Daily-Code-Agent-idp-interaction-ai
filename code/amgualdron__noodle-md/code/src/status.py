"""./noodle status [run]  --  how a run is going: done, running, waiting, stopped, failed.

With no run (outside a run folder), it lists every run with its number and how far along
it is, from runs/runs.csv, which it brings up to date.

The report is printed, each list showing at most SHOW_ON_SCREEN tasks. Every task, with
its state, progress, job and error, is also written to runs/<run name>/log.csv (open it
in a spreadsheet to sort and filter).

A task's state comes from two places:
    its row in its job file (jobs/job_0007.csv)  how far it got, written by the job itself
    the SLURM queue                              whether its job is waiting or running now
On a computer without SLURM, a task counts as running if its row was updated recently
(a running task updates it at least once an hour).
"""

import os
import sys
from collections import Counter, defaultdict
from datetime import datetime

import yaml

from runs import (RUNS, about, side_by_side, compress, duration, list_runs, measured_speeds, read_submissions, read_tasks,
                  refresh_registry, run_from_command_line, run_here, seconds_left, slurm_queue, task_states,
                  write_csv)

SHOW_ON_SCREEN = 10  # tasks per list on screen (log.csv lists them all)
LIST_DONE_UP_TO = 20  # runs with this many tasks or fewer also list every finished task
STALE_HOURS = 0.5     # a running task that has not reported for this long (3 reports) is stuck or gone

# In the order the report shows them: what needs you first
STATES = ("failed", "stopped", "running", "waiting", "not launched", "done")
COLORS = {"failed": "31", "stopped": "33", "running": "36", "waiting": "2", "not launched": "2",
          "done": "32", "bold": "1", "dim": "2"}  # ANSI codes: red, yellow, cyan, dim, green


def hours_since(timestamp):
    return (datetime.now() - datetime.fromisoformat(timestamp)).total_seconds() / 3600


def task_state(log, slurm_job, queue, job_active=False):
    """The state of one task, with what the report needs to describe it:
        state    one of STATES
        percent  how far it got (None if it never started)
        silent   hours since it last reported, if a running task is overdue (else None)

    log        the task's current row from its job file ({} if it was never in a job)
    slurm_job  the SLURM id of that job (None if it was never submitted)
    queue      {SLURM id: SLURM state} on a cluster; None on a computer without SLURM
    job_active on a computer without SLURM: another task of the same job reported recently
    """
    info = dict(state=None, percent=None, silent=None)
    if log.get("total_steps"):
        info["percent"] = 100 * (log.get("step") or 0) / log["total_steps"]
    status = log.get("status")
    in_queue = queue is not None and slurm_job in queue

    if status == "completed":
        return dict(info, state="done", percent=100)
    if status == "failed":
        return dict(info, state="failed")
    if not log:
        return dict(info, state="not launched")
    if status == "submitted":
        # its job is waiting, or is running the tasks before it; else the job ended first
        return dict(info, state="waiting" if in_queue or job_active else "stopped")
    if status == "stopped":
        return dict(info, state="stopped")  # it stopped itself (time limit, scancel, Ctrl-C)

    # status running: is it really?
    if queue is None:  # own computer: has it reported recently?
        alive = hours_since(log["last_update"]) < STALE_HOURS
    else:              # cluster: is its job still in the queue?
        alive = in_queue
        if alive and hours_since(log["last_update"]) > STALE_HOURS:
            info["silent"] = hours_since(log["last_update"])  # in the queue but not reporting
    return dict(info, state="running" if alive else "stopped")


def add_time_left(rows, physics, speeds):
    """Give every unfinished, not failed task r["left"]: the run time it still needs, in
    seconds (runs.seconds_left). None until some task of the run has measured its speed."""
    for r in rows:
        r["left"] = None
        if r["state"] not in ("done", "failed") and speeds:
            r["left"] = seconds_left(r["task"], r["log"], physics[r["task"]["physics"]], speeds)


def waiting_time(rows, max_running, on_cluster):
    """Roughly how long until the run's waiting and running tasks are all done, once their
    jobs run: the tasks of a job run one after another, and jobs run side by side (at most
    max_running at once; one at a time without SLURM). Time spent waiting in the SLURM
    queue is not included. None if some of these tasks have no estimate."""
    jobs = defaultdict(float)  # job number -> seconds its unfinished tasks still need
    for r in rows:
        if r["state"] in ("running", "waiting"):
            if r["left"] is None:
                return None
            jobs[r["log"]["job"]] += r["left"]
    if not jobs:
        return None
    if not on_cluster:
        return sum(jobs.values())
    return side_by_side([jobs[job] for job in sorted(jobs)], max_running)


# ── The report ────────────────────────────────────────────────────────────────
class Report:
    """Builds the report's lines, in color for the screen (or plain, when printed into a file)."""

    def __init__(self, rows, color, limit):
        self.lines, self.color, self.limit = [], color, limit
        tasks = [r["task"] for r in rows]
        self.id_width = max(len(str(t["task_id"])) for t in tasks)
        self.name_width = min(24, max(len(t["sequence_name"]) for t in tasks))
        self.physics_width = max(len(t["physics"]) for t in tasks)

    def paint(self, text, color):
        return f"\033[{COLORS[color]}m{text}\033[0m" if self.color else text

    def add(self, *lines):
        self.lines += lines

    def task(self, r, details):
        t = r["task"]
        name = t["sequence_name"]
        if len(name) > self.name_width:
            name = name[:self.name_width - 1] + "…"
        self.add(f"  {t['task_id']:>{self.id_width}}  {name:{self.name_width}}  N={t['length']:<4}  "
                 f"{t['physics']:{self.physics_width}}  {details}")

    def tasks(self, rows, details):
        """One line per task, cut off on screen after self.limit lines (unless only a
        couple more would follow)."""
        shown = rows if self.limit is None or len(rows) <= self.limit + 2 else rows[:self.limit]
        for r in shown:
            self.task(r, details(r))
        if len(rows) > len(shown):
            self.add(self.paint(f"  … and {len(rows) - len(shown)} more (all are listed in log.csv)", "dim"))

    def heading(self, state, count, text=""):
        self.add("", self.paint(f"{state.upper()} {count}", state) + (f"  {text}" if text else ""))


def progress_bar(counts, total, width=30):
    """█ done, ▓ running, ░ everything else."""
    done = round(width * counts["done"] / total)
    running = round(width * (counts["done"] + counts["running"]) / total) - done
    return "█" * done + "▓" * running + "░" * (width - done - running)


def build_report(run_dir, settings, rows, speeds, color, limit, on_cluster):
    name = run_dir.name
    by_state = {s: [r for r in rows if r["state"] == s] for s in STATES}
    counts = {s: len(by_state[s]) for s in STATES}
    physics_sets = list(dict.fromkeys(r["task"]["physics"] for r in rows))
    launched = datetime.fromtimestamp((run_dir / "tasks.csv").stat().st_mtime)
    report = Report(rows, color, limit)

    # ── Summary ──
    if counts["done"] == len(rows):  # all done: from launch to the last task finishing
        finished = max((datetime.fromisoformat(r["log"]["completed_at"]) for r in rows
                        if "completed_at" in r["log"]), default=datetime.now())
        elapsed = (report.paint(f"took {duration((finished - launched).total_seconds())}", "bold")
                   + report.paint(f"   launched {launched:%Y-%m-%d %H:%M} · last task done {finished:%Y-%m-%d %H:%M}", "dim"))
    else:
        elapsed = (report.paint(f"{duration((datetime.now() - launched).total_seconds())} since launch", "bold")
                   + report.paint(f"   launched {launched:%Y-%m-%d %H:%M} · checked {datetime.now():%Y-%m-%d %H:%M}", "dim"))
    report.add(report.paint(name, "bold") + f"  {settings.get('description') or '(no description)'}",
               elapsed,
               "",
               f"{progress_bar(counts, len(rows))}  {counts['done']} of {len(rows)} tasks done "
               f"({100 * counts['done'] / len(rows):.0f}%)",
               "  " + " · ".join(report.paint(f"{counts[s]} {s}", s) for s in STATES if counts[s]))

    to_go = [r for r in rows if r["state"] not in ("done", "failed")]
    if to_go:
        if speeds:
            left = [r["left"] for r in to_go]
            wait = waiting_time(rows, settings.get("max_running"), on_cluster)
            active = sum(r["state"] in ("running", "waiting") for r in rows)
            stopped = len(to_go) - active
            computing = (f"{about(sum(left))} of computing in all"
                         + (f" ({stopped} stopped or not launched: noodle launch starts them)" if stopped else "")
                         + f"; from the speeds measured by {len(speeds)} task{'s' * (len(speeds) > 1)} of this run")
            if wait is not None:
                report.add("", f"Time left  {about(wait)} until the {active} running and waiting task"
                               f"{'s are' if active > 1 else ' is'} done"
                               + (", plus any wait in the SLURM queue" if on_cluster else ""),
                           report.paint("           " + computing, "dim"))
            else:
                report.add("", "Time left  " + computing)
        else:
            report.add("", "Time left  not known yet: each task measures its speed a minute after it starts")

    if len(physics_sets) > 1:  # the same counts for each physics set, with done as done/total
        shown = [s for s in STATES if counts[s]]
        table = {}
        for physics in physics_sets:
            in_set = [r for r in rows if r["task"]["physics"] == physics]
            per_state = Counter(r["state"] for r in in_set)
            table[physics] = [f"{per_state[s]}/{len(in_set)}" if s == "done" else str(per_state[s] or "·")
                              for s in shown]
        widths = [max(7, len(s), *(len(cells[i]) for cells in table.values())) for i, s in enumerate(shown)]
        report.add("", "  " + " " * report.physics_width + "".join(f"  {s:>{w}}" for s, w in zip(shown, widths)))
        for physics, cells in table.items():
            report.add(f"  {physics:{report.physics_width}}" + "".join(f"  {c:>{w}}" for c, w in zip(cells, widths)))

    # ── Failed: grouped by error, since one cause often fails many tasks ──
    if counts["failed"]:
        report.heading("failed", counts["failed"], f"fix the cause first · full messages in runs/{name}/jobs/job_<number>.out")
        by_error = defaultdict(list)
        for r in by_state["failed"]:
            by_error[r["log"].get("error") or "(no error recorded)"].append(r["task"]["task_id"])
        groups = sorted(by_error.items(), key=lambda g: -len(g[1]))
        for error, ids in groups if limit is None else groups[:limit]:
            jobs = sorted({r["log"]["job"] for r in by_state["failed"] if r["task"]["task_id"] in ids})
            report.add(f"  {len(ids)} × {error}",
                       report.paint(f"      task{'s' * (len(ids) > 1)} {compress(sorted(ids)).replace(',', ', ')}"
                                    f"  (job{'s' * (len(jobs) > 1)} {compress(jobs).replace(',', ', ')})", "dim"))
        if limit is not None and len(groups) > limit:
            report.add(report.paint(f"  … and {len(groups) - limit} other errors (all are listed in log.csv)", "dim"))

    # ── Stopped ──
    if counts["stopped"]:
        report.heading("stopped", counts["stopped"], f"./noodle launch {name} resumes them from their last checkpoint")
        report.tasks(by_state["stopped"], lambda r: "not started yet" if r["percent"] is None else
                     f"stopped at {r['percent']:3.0f}% after {duration(r['log'].get('run_time_s') or 0)} of run time"
                     + (f" ({r['log']['stop_reason']})" if r["log"].get("stop_reason") else "")
                     + (f", needs {about(r['left'])} more" if r["left"] is not None else ""))

    # ── Running ──
    if counts["running"]:
        report.heading("running", counts["running"])

        def running_details(r):
            if r["percent"] is None:
                return "starting"
            text = (f"{r['log'].get('phase') or 'starting':13s} {r['percent']:3.0f}%  "
                    f"{duration(r['log'].get('run_time_s') or 0)} so far")
            if r["left"] is not None:
                text += f", {about(r['left'])} left"
            if r["silent"]:
                text += report.paint(f"  no report for {duration(r['silent'] * 3600)}: stuck?", "stopped")
            return text
        report.tasks(by_state["running"], running_details)

    # ── Waiting and not launched: just how many ──
    if counts["waiting"]:
        resuming = sum(r["percent"] is not None for r in by_state["waiting"])
        report.heading("waiting", counts["waiting"], "in the SLURM queue, or next in a job that is running"
                       + (f" ({resuming} will resume from a checkpoint)" if resuming else ""))
    if counts["not launched"]:
        report.heading("not launched", counts["not launched"], f"./noodle launch {name} starts them")

    # ── Done ──
    if counts["done"]:
        done = sorted(by_state["done"], key=lambda r: r["log"]["run_time_s"])
        fastest, slowest = done[0], done[-1]
        total = sum(r["log"]["run_time_s"] for r in done)
        report.heading("done", counts["done"],
                       f"run time {duration(fastest['log']['run_time_s'])} (N={fastest['task']['length']}) "
                       f"to {duration(slowest['log']['run_time_s'])} (N={slowest['task']['length']}), "
                       f"{duration(total)} in total")
        if len(rows) <= LIST_DONE_UP_TO:
            report.tasks(done, lambda r: f"{r['log']['frames']} frames in {duration(r['log']['run_time_s'])}")

    # ── What to do next ──
    report.add("")
    if counts["done"] == len(rows):
        report.add(report.paint("All tasks are done.", "done") + f" Next: ./noodle analyze {name}")
    elif counts["stopped"] or counts["not launched"]:
        todo = [f"resumes the {counts['stopped']} stopped" if counts["stopped"] else "",
                f"starts the {counts['not launched']} not launched" if counts["not launched"] else ""]
        report.add(f"Next: ./noodle launch {name} " + " and ".join(t for t in todo if t)
                   + f" task{'s' * (counts['stopped'] + counts['not launched'] > 1)}"
                   + (" (and asks before retrying the failed ones)." if counts["failed"] else "."))
    elif counts["failed"] and not (counts["running"] or counts["waiting"]):
        report.add(f"Next: fix the errors above, then ./noodle launch {name} retries the failed tasks.")
    else:
        report.add(f"Nothing to do: check back later. ./noodle analyze {name} already works on the finished tasks."
                   if counts["done"] else "Nothing to do: check back later.")
    return report.lines


# ── log.csv: every task, as a table ───────────────────────────────────────────
LOG_COLUMNS = ["task_id", "sequence_name", "physics", "replicate", "N", "state", "percent", "phase", "job",
               "slurm_job", "steps_per_second", "run_time_h", "time_left_h", "attempts", "host", "started_at",
               "last_update", "completed_at", "stop_reason", "error"]


def write_log_csv(run_dir, rows, submitted):
    table = []
    for r in rows:
        task, log = r["task"], r["log"]
        table.append(dict(
            task_id=task["task_id"], sequence_name=task["sequence_name"], physics=task["physics"],
            replicate=task["replicate"], N=task["length"], state=r["state"],
            percent=None if r["percent"] is None else round(r["percent"], 1), phase=log.get("phase"),
            job=log.get("job"), slurm_job=submitted.get(log.get("job")), steps_per_second=log.get("steps_per_second"),
            run_time_h=round(log["run_time_s"] / 3600, 3) if log.get("run_time_s") else None,
            time_left_h=round(r["left"] / 3600, 3) if r["left"] is not None else None,
            attempts=log.get("attempts"), host=log.get("host"), started_at=log.get("started_at"),
            last_update=log.get("last_update"), completed_at=log.get("completed_at"),
            stop_reason=log.get("stop_reason"), error=log.get("error")))
    write_csv(table, run_dir / "log.csv", LOG_COLUMNS)


# ── Every run at a glance ─────────────────────────────────────────────────────
def list_all_runs():
    for run_dir in sorted(RUNS.glob("*")) if RUNS.exists() else []:
        if (run_dir / "run.yaml").exists():
            refresh_registry(run_dir)
    list_runs()
    print("\nDetails of one run: noodle status <run name or number>   (all numbers: runs/runs.csv)")


def main():
    if len(sys.argv) == 1 and run_here() is None:
        list_all_runs()
        return
    run_dir = run_from_command_line("noodle status <run>")
    settings = yaml.safe_load((run_dir / "run.yaml").read_text())
    if not (run_dir / "tasks.csv").exists():
        print(f"'{run_dir.name}' has not been launched yet: ./noodle launch {run_dir.name}")
        return

    queue = slurm_queue(run_dir.name)
    submitted = read_submissions(run_dir)
    states = task_states(run_dir)
    tasks = read_tasks(run_dir)
    # without SLURM, a job is going while one of its tasks has reported recently
    active_jobs = {s["job"] for s in states.values()
                   if s["status"] == "running" and hours_since(s["last_update"]) < STALE_HOURS}
    rows = []
    for task in tasks:
        log = states.get(task["task_id"], {})
        rows.append(dict(task=task, log=log, **task_state(log, submitted.get(log.get("job")), queue,
                                                          log.get("job") in active_jobs)))
    speeds = measured_speeds(tasks, states)
    add_time_left(rows, yaml.safe_load((run_dir / "physics.yaml").read_text()), speeds)

    color = sys.stdout.isatty() and "NO_COLOR" not in os.environ
    print("\n".join(build_report(run_dir, settings, rows, speeds, color, SHOW_ON_SCREEN, queue is not None)))
    write_log_csv(run_dir, rows, submitted)
    refresh_registry(run_dir, states)


if __name__ == "__main__":
    main()
