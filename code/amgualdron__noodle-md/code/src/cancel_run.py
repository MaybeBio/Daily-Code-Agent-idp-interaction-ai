"""./noodle cancel <run>  --  take a run's jobs out of the SLURM queue.

It counts the run's jobs that are running or waiting, asks, and then cancels them with
scancel. A running task saves a checkpoint before it stops (as at the time limit), so
nothing is lost: ./noodle launch starts the run's unfinished tasks again later, from
their checkpoints.

On your own computer (no SLURM) a run goes in the terminal it was launched from:
Ctrl-C there stops it the same way.
"""

import subprocess
import sys
from collections import Counter

from runs import read_submissions, run_from_command_line, slurm_queue


def main():
    run_dir = run_from_command_line("noodle cancel <run>")
    name = run_dir.name
    queue = slurm_queue(name)
    if queue is None:
        sys.exit("No SLURM here: a run on this computer goes in the terminal it was launched from, "
                 "and Ctrl-C there stops it.")

    # only this run's jobs: the SLURM ids it submitted that are still in the queue
    slurm_ids = [i for i in read_submissions(run_dir).values() if i in queue]
    if not slurm_ids:
        print(f"'{name}' has no jobs in the SLURM queue: nothing to cancel.")
        return
    states = Counter(queue[i] for i in slurm_ids)
    print(f"'{name}' has {len(slurm_ids)} job(s) in the queue: "
          + ", ".join(f"{count} {state.lower()}" for state, count in states.most_common()))
    try:
        answer = input("Cancel them all? [y/N] ")
    except EOFError:
        answer = ""
    if answer.strip().lower() not in ("y", "yes"):
        print("Nothing was cancelled.")
        return

    # Every job array this run submitted belongs to it alone, so cancelling the arrays
    # (1234 rather than 1234_1, 1234_2, ...) cancels exactly its jobs, in one short command.
    arrays = sorted({i.split("_")[0] for i in slurm_ids})
    result = subprocess.run(["scancel", *arrays], capture_output=True, text=True)
    if result.returncode != 0:
        sys.exit(f"scancel refused:\n{result.stderr}")
    print(f"Cancelled. Running tasks save a checkpoint before they stop (a minute or so).\n"
          f"noodle status {name} shows where each one got to; noodle launch {name} continues the run.")


if __name__ == "__main__":
    main()
