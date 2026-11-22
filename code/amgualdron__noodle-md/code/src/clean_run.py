"""./noodle clean <run name>  --  free the disk space of finished tasks, keeping their results.

Trajectories are big (megabytes per task) and only needed to compute new things, so
once they are copied to long-term storage (the group NAS, by Globus) they can go:
  1. the results are pooled first (as ./noodle analyze). Every completed task was
     analysed by its job (or else is analysed now), and its results stay in jobs/ and
     analysis/, so ./noodle analyze keeps working without the trajectories;
  2. the SHA-256 checksum of each trajectory is added to trajectories.sha256, a record
     of exactly what was deleted that any copy can be checked against;
  3. it shows what it would delete and keep, and deletes only once you type the run's
     name: the trajectories of completed tasks that have been analysed, and any
     checkpoint such a task left behind.

Never deleted: unfinished tasks' trajectories and checkpoints (they are needed to
resume them), and everything small: run.yaml, tasks.csv, physics.yaml, manifest.csv,
jobs/, log.csv, analysis/. It can be run again as more tasks finish.
"""

import hashlib
import sys

from analyze_run import analyze
from runs import file_size, read_tasks, run_from_command_line, task_states, trajectory_path

CHECKSUMS = "trajectories.sha256"  # in the run folder; lists names relative to trajectories/


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as f:
        while chunk := f.read(1 << 20):
            digest.update(chunk)
    return digest.hexdigest()


def add_checksums(run_dir, trajectories):
    """Add the trajectories to trajectories.sha256, in the format `sha256sum -c` reads."""
    path = run_dir / CHECKSUMS
    known = {}
    if path.exists():
        for line in path.read_text().splitlines():
            if line.strip():
                digest, name = line.split(maxsplit=1)
                known[name.lstrip("*")] = digest
    for done, trajectory in enumerate(trajectories, 1):
        known[trajectory.name] = sha256(trajectory)
        if sys.stdout.isatty() and (done % 20 == 0 or done == len(trajectories)):
            print(f"\r  checksums of {done} of {len(trajectories)} trajectories", end="\n" if done == len(trajectories) else "")
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text("".join(f"{digest}  {name}\n" for name, digest in sorted(known.items())))
    tmp.replace(path)


def total_size(paths):
    return file_size(sum(p.stat().st_size for p in paths if p.is_file()))


def main():
    run_dir = run_from_command_line("noodle clean <run>")
    name = run_dir.name
    if not (run_dir / "tasks.csv").exists():
        sys.exit(f"'{name}' has not been launched yet: there is nothing to clean.")

    print("Analysing first (a trajectory is only deleted once its results are saved)")
    analysed = analyze(run_dir)
    print()

    trajectories, leftovers, unfinished, unanalysed = [], [], [], []
    states = task_states(run_dir)
    for task in read_tasks(run_dir):
        task_id = task["task_id"]
        log = states.get(task_id, {})
        stem = f"task_{task_id:04d}"
        own_files = [p for folder in ("trajectories", "checkpoints") for p in (run_dir / folder).glob(f"{stem}[_.]*")]
        if log.get("status") != "completed":
            unfinished += own_files
        elif task_id not in analysed:
            unanalysed += own_files
        else:
            trajectory = trajectory_path(run_dir, task)
            trajectories += [trajectory] if trajectory.exists() else []
            leftovers += [p for p in own_files if p != trajectory]  # checkpoints, half-written .tmp files

    also_kept = []
    if unfinished:
        also_kept.append(f"the files of unfinished tasks ({total_size(unfinished)}): needed to resume them")
    if unanalysed:
        also_kept.append(f"the files of completed tasks that could not be analysed ({total_size(unanalysed)})")
    if not trajectories and not leftovers:
        print(f"Nothing to clean: no completed task of '{name}' still has its trajectory here.")
        for line in also_kept:
            print(f"Kept       {line}")
        return

    if trajectories:
        add_checksums(run_dir, trajectories)
        print(f"  -> {CHECKSUMS}\n")
    freed = total_size(trajectories + leftovers)
    print(f"Delete     the trajectories of {len(trajectories)} completed, analysed task(s)"
          + (f", and {len(leftovers)} leftover file(s) (checkpoints, unfinished writes)" if leftovers else "") + f": {freed}")
    print(f"Keep       analysis/ ({total_size((run_dir / 'analysis').rglob('*'))}), {CHECKSUMS}, jobs/, "
          f"manifest.csv, tasks.csv, physics.yaml, run.yaml")
    for line in also_kept:
        print(f"           {line}")
    print()
    print(f"This can't be undone. First copy runs/{name}/trajectories/ to the NAS (Globus checks each")
    print(f"file as it transfers it). To check a copy yourself, put {CHECKSUMS} in its folder and run")
    print(f"    sha256sum -c {CHECKSUMS}")
    try:
        answer = input(f"Type the run name ({name}) to delete them, anything else to cancel: ")
    except EOFError:
        answer = ""
    if answer.strip() != name:
        print("Cancelled: nothing was deleted.")
        return

    for path in trajectories + leftovers:
        path.unlink(missing_ok=True)
    print(f"Deleted {len(trajectories) + len(leftovers)} file(s), {freed} freed.")
    print(f"./noodle analyze {name} keeps working (the results are in jobs/ and analysis/); to compute "
          f"something new from the trajectories, copy them back into runs/{name}/trajectories/.")


if __name__ == "__main__":
    main()
