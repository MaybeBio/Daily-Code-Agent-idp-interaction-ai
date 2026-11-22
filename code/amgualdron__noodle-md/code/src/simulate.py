"""Run one job of a run: its tasks one after another. Each task builds its chain (or brush),
equilibrates it, then records a trajectory.

    env/bin/python src/simulate.py runs/<run name> <job number>

You normally never call this yourself: ./noodle launch does it (through SLURM, or right
here on your own computer), once per job.

Everything a job needs comes from its run folder:
    jobs/job_0007.csv  which tasks this job runs, in order (written by ./noodle launch)
    tasks.csv          each task's sequence, physics set and random seed
    physics.yaml       the physics settings, frozen when the run was launched
and everything it produces goes back there:
    jobs/job_0007.csv                  each task's state and progress, updated as it goes
                                       (read by ./noodle status)
    trajectories/task_0001_<name>.gsd  a task's production frames
    checkpoints/task_0001.gsd          a task's full state, saved every 10 minutes and when it is
                                       asked to stop; deleted when the task is done

A task that failed is recorded as failed and the job goes on with the next one. A task
that was stopped resumes from its last checkpoint the next time it runs.

HOOMD in brief: a `hoomd.Simulation` holds the `state` (box, particles, bonds) and the
`operations`: an integrator that moves the particles under the forces from
forcefield.py, and writers that save frames to a file. `sim.run(n)` advances n steps.
"""

import os
import signal
import socket
import sys
import time
import traceback
from datetime import datetime
from pathlib import Path

import gsd.fl
import gsd.hoomd
import hoomd
import hoomd.md
import numpy as np
import yaml

from forcefield import (AA_PARAMS, AMINO_ACIDS, residue_pair_forces, thermal_energy,
                        wall_contact_distance, wall_forces)
from observables import analyse_task
from runs import (brush_size, checkpoint_path, read_job, read_tasks, run_lengths, trajectory_path,
                  write_job)

REPORT_EVERY_MINUTES = 10  # how often to print progress, update the job file and save a checkpoint
FIRST_REPORT_SECONDS = 60  # one early report, so you can see the speed right away
CHUNK_SECONDS = 60         # the simulation runs in chunks of about this long; between chunks it
                           # checks the clock and whether it was asked to stop (not a physics setting)
FIRST_CHUNK_STEPS = 1000   # the first chunk, before the speed is known (and the smallest chunk)
STOPPED_EXIT_CODE = 3      # exit code of a job that stopped because it was asked to


# ── Stopping cleanly ──────────────────────────────────────────────────────────
# A job can be asked to stop by a signal:
#   SIGUSR1  SLURM, 10 minutes before the time limit (#SBATCH --signal in slurm_job.sh)
#   SIGTERM  scancel (./noodle cancel)
#   SIGINT   Ctrl-C, when running on your own computer
# HOOMD can't be interrupted in the middle of sim.run(), so the signal only notes the
# request here; after the current chunk (about CHUNK_SECONDS) the running task saves a
# checkpoint and the job stops, losing nothing. ./noodle launch resumes it from there.
STOP_REASONS = {signal.SIGUSR1: "time limit", signal.SIGTERM: "cancelled", signal.SIGINT: "Ctrl-C"}
stop_requests = []  # the reasons received so far


def request_stop(signal_number, frame):
    if stop_requests and signal_number == signal.SIGINT:
        raise KeyboardInterrupt  # a second Ctrl-C: quit right away (the last checkpoint, at most 10 minutes old, stays)
    reason = STOP_REASONS[signal_number]
    stop_requests.append(reason)
    print(f"\n[{now()}] Asked to stop ({reason}): saving a checkpoint after the current chunk"
          + (" (Ctrl-C again to quit right away)" if reason == "Ctrl-C" else ""), flush=True)


# ── Geometry and starting configuration ───────────────────────────────────────
def system_geometry(sequence, physics):
    """The simulation box and, for grafted chains, where the surface, ceiling and anchors are.

    bulk:    one chain in a large periodic cube.
    grafted: the surface sits near the bottom of the box. A purely repulsive ceiling
             sits above it, higher than a fully stretched chain can reach (contour
             length plus wall_headroom), so beads can never wrap through the periodic
             boundary in z. Above the ceiling there is one interaction cutoff of empty
             space, so beads near the ceiling never interact with beads near the surface
             through that boundary either.
             One chain is grafted at the centre of the surface; a brush (grafting_density
             > 0) has its chains grafted at the centres of a square lattice of cells,
             `spacing` apart (runs.brush_size gives their number and the box width).
    Chain k holds particles k*N ... k*N + N-1, in sequence order.
    """
    if physics["geometry"] == "bulk":
        return dict(box=[physics["box_size"]] * 3)
    if physics["geometry"] != "grafted":
        raise ValueError(f"unknown geometry '{physics['geometry']}' (use bulk or grafted)")

    if physics["graft_terminus"] not in ("N", "C"):
        raise ValueError("graft_terminus must be N or C")
    n = len(sequence)
    chains, width = brush_size(physics)
    side = round(chains ** 0.5)
    contour = (n - 1) * physics["harmonic_r0"]
    cutoff = max(physics["lj_r_cut"], physics["yukawa_r_cut"], physics["wall_r_cut"])
    height = contour + physics["wall_headroom"]  # surface to ceiling
    box_z = height + cutoff + 1.0
    z_wall = -box_z / 2 + 0.5
    # the grafted residue of each chain: its first one (N-terminus) or its last one (C-terminus)
    end = 0 if physics["graft_terminus"] == "N" else n - 1
    return dict(box=[width, width, box_z], z_wall=z_wall, z_ceiling=z_wall + height,
                chains=chains, spacing=width / side, anchors=[k * n + end for k in range(chains)])


def grafted_start_positions(sequence, physics, geometry, rng):
    """Random, non-overlapping chains growing upward from their anchors.

    Each anchor sits at its contact distance from the surface, where the wall exerts no
    force on it. Each following bead is one bond length away in a random upward
    direction, retried until it is at least 0.5 nm from every earlier bead of its chain.
    In a brush, each chain grows inside its own column of the lattice, kept 0.25 nm away
    from the column's sides, so beads of different chains also start at least 0.5 nm
    apart: the start needs no extra step to remove overlaps, at any density.
    """
    n, r0 = len(sequence), physics["harmonic_r0"]
    anchor_index = geometry["anchors"][0]  # 0 or n - 1: the same residue in every chain
    anchor_sigma = AA_PARAMS[sequence[anchor_index]]["sigma"]
    z_anchor = geometry["z_wall"] + wall_contact_distance(anchor_sigma, physics)
    side, spacing = round(geometry["chains"] ** 0.5), geometry["spacing"]
    half_width = spacing / 2 - 0.25 if geometry["chains"] > 1 else np.inf  # the column's reach

    chains = []
    for k in range(geometry["chains"]):
        # lattice cell centres, row by row, in a box centred on 0
        centre = np.array([(k % side + 0.5) * spacing, (k // side + 0.5) * spacing]) - spacing * side / 2
        grown = np.zeros((n, 3))
        grown[0] = [*centre, z_anchor]
        for i in range(1, n):
            for _ in range(10_000):
                direction = rng.normal(size=3)
                direction /= np.linalg.norm(direction)
                direction[2] = abs(direction[2])
                candidate = grown[i - 1] + r0 * direction
                if (np.all(np.abs(candidate[:2] - centre) <= half_width)
                        and np.all(np.linalg.norm(grown[:i - 1] - candidate, axis=1) >= 0.5)):
                    break
            else:
                raise RuntimeError("could not grow a non-overlapping starting chain")
            grown[i] = candidate
        # beads were grown from the anchor outward; for C-terminal grafting that is backwards
        chains.append(grown if anchor_index == 0 else grown[::-1])
    return np.concatenate(chains)


def starting_snapshot(sequence, physics, geometry, seed):
    """The starting configuration: one bead per residue, bonded in a chain (or several
    identical chains, for a brush).

    A HOOMD snapshot is a plain description of the system: the box, and per particle
    its type, mass, charge and position, plus the list of bonds.
    """
    n = len(sequence)
    chains = geometry.get("chains", 1)
    if physics["geometry"] == "bulk":
        # started as a straight line along x; equilibration makes it a random coil
        r0 = physics["harmonic_r0"]
        positions = [[i * r0 - n * r0 / 2, 0, 0] for i in range(n)]
    else:
        positions = grafted_start_positions(sequence, physics, geometry, np.random.default_rng(seed))

    snapshot = hoomd.Snapshot()
    snapshot.configuration.box = geometry["box"] + [0, 0, 0]  # Lx, Ly, Lz, no tilt
    snapshot.particles.N = n * chains
    snapshot.particles.types = AMINO_ACIDS
    snapshot.particles.typeid[:] = [AMINO_ACIDS.index(aa) for aa in sequence] * chains
    snapshot.particles.mass[:] = [AA_PARAMS[aa]["mass"] for aa in sequence] * chains
    snapshot.particles.charge[:] = [AA_PARAMS[aa]["charge"] for aa in sequence] * chains
    snapshot.particles.diameter[:] = [AA_PARAMS[aa]["sigma"] for aa in sequence] * chains  # for viewers only
    snapshot.particles.position[:] = positions
    snapshot.bonds.N = (n - 1) * chains
    snapshot.bonds.types = ["AA-bond"]
    snapshot.bonds.typeid[:] = 0
    snapshot.bonds.group[:] = [[k * n + i, k * n + i + 1] for k in range(chains) for i in range(n - 1)]
    return snapshot


# ── Simulation and forces ─────────────────────────────────────────────────────
def create_simulation(seed):
    """An empty HOOMD simulation, on the GPU if this HOOMD was built for one.

    The seed drives the thermostat's random kicks."""
    device = hoomd.device.GPU() if hoomd.version.gpu_enabled else hoomd.device.CPU()
    return hoomd.Simulation(device=device, seed=seed)


def add_forces(sim, physics, geometry):
    """Give a simulation (whose particles already exist) its forces and its integrator."""
    # Neighbour list: a tree is efficient for one chain in a big, mostly empty box; a cell
    # list for a brush, whose box is filled (~1.3x faster there, measured). The buffer
    # makes it rebuild only when some bead has moved more than 0.1 nm.
    if geometry.get("chains", 1) > 1:
        nlist = hoomd.md.nlist.Cell(buffer=0.2, exclusions=["bond"])
    else:
        nlist = hoomd.md.nlist.Tree(buffer=0.2, exclusions=["bond"])
    forces = residue_pair_forces(physics, nlist)
    moving = hoomd.filter.All()
    if physics["geometry"] == "grafted":
        forces += wall_forces(physics, geometry["z_wall"], geometry["z_ceiling"])
        # The anchors are left out of the integrator: they feel forces but never move.
        moving = hoomd.filter.SetDifference(hoomd.filter.All(), hoomd.filter.Tags(geometry["anchors"]))
    # Langevin thermostat: friction plus random kicks that keep the temperature.
    # Friction gamma = 1 for every residue (HOOMD's default), i.e. a velocity
    # relaxation time m/gamma of ~100 time units (~50 ps) for a typical residue.
    langevin = hoomd.md.methods.Langevin(filter=moving, kT=thermal_energy(physics["temperature_K"]))
    sim.operations.integrator = hoomd.md.Integrator(dt=physics["dt"], methods=[langevin], forces=forces)


# ── Checkpoints ───────────────────────────────────────────────────────────────
def save_checkpoint(sim, path):
    """Save the full state (positions and velocities, in double precision).

    Written to a temporary file first, so a task killed while saving still has its
    previous checkpoint."""
    tmp = path.with_suffix(".tmp")
    hoomd.write.GSD.write(state=sim.state, filename=str(tmp), mode="wb", precision="double")
    os.replace(tmp, path)


def drop_frames_after(trajectory, step):
    """When resuming: keep only frames up to the checkpoint; the rest is recomputed.

    The frames are copied chunk by chunk with the low-level gsd.fl module, so the new
    file is in exactly the format HOOMD wrote and HOOMD can keep appending to it.
    """
    tmp = trajectory.with_suffix(".tmp")
    with gsd.fl.open(trajectory, "r") as old, gsd.fl.open(
            tmp, "w", application=old.application, schema=old.schema,
            schema_version=old.schema_version) as new:
        names = old.find_matching_chunk_names("")
        for frame in range(old.nframes):
            if old.read_chunk(frame, "configuration/step")[0] > step:
                break
            for chunk in names:
                if old.chunk_exists(frame, chunk):
                    new.write_chunk(chunk, old.read_chunk(frame, chunk))
            new.end_frame()
    os.replace(tmp, trajectory)


def now():
    return datetime.now().isoformat(timespec="seconds")


# ── One task ──────────────────────────────────────────────────────────────────
def run_task(run_dir, task, physics, row, update):
    """Run one task to the end, or until the job is asked to stop.

    row: the task's row in the job file; update(row, field=value, ...) changes it and
    writes the job file. Returns "completed" or "stopped"; an error is raised."""
    task_id, name, sequence, seed = task["task_id"], task["sequence_name"], task["sequence"], task["seed"]
    trajectory, checkpoint = trajectory_path(run_dir, task), checkpoint_path(run_dir, task_id)
    trajectory.parent.mkdir(exist_ok=True)
    checkpoint.parent.mkdir(exist_ok=True)

    equilibration, production, save_every = run_lengths(physics, len(sequence))
    total = equilibration + production
    geometry = system_geometry(sequence, physics)

    # ── System: resume from the last checkpoint, or start fresh ──
    sim = create_simulation(seed)
    if checkpoint.exists():
        sim.create_state_from_gsd(filename=str(checkpoint))
        if trajectory.exists():
            drop_frames_after(trajectory, sim.timestep)
        print(f"Resuming from checkpoint at step {sim.timestep:,} of {total:,}")
    else:
        sim.create_state_from_snapshot(starting_snapshot(sequence, physics, geometry, seed))
        trajectory.unlink(missing_ok=True)  # leftovers from an attempt that never checkpointed

    add_forces(sim, physics, geometry)

    # Frames every save_every steps, counted from the start of production
    writer = hoomd.write.GSD(filename=str(trajectory), mode="ab",
                             trigger=hoomd.trigger.Periodic(save_every, phase=equilibration % save_every))

    chains = geometry.get("chains", 1)
    phase = "equilibration" if sim.timestep < equilibration else "production"
    update(row, status="running", phase=phase, N=len(sequence), chains=chains, step=sim.timestep,
           total_steps=total, attempts=(row.get("attempts") or 0) + 1, host=socket.gethostname(),
           started_at=row.get("started_at") or now(), last_update=now(), stop_reason=None, error=None)
    brush = (f" x {chains} chains (brush, {physics['grafting_density']:g} chains/nm^2, "
             f"box {geometry['box'][0]:.1f} nm wide)" if chains > 1 else "")
    print(f"Task {task_id}: {name}, {len(sequence)} residues{brush}, physics '{task['physics']}' "
          f"({physics['geometry']}), seed {seed}")
    print(f"  {equilibration:,} equilibration + {production:,} production steps, "
          f"a frame every {save_every:,} steps -> {trajectory.relative_to(run_dir)}", flush=True)

    # ── Run, saving progress every REPORT_EVERY_MINUTES ──
    # time.monotonic() is a stopwatch: unlike the wall clock it never jumps
    previous_run_time = row.get("run_time_s") or 0  # from earlier attempts of this task
    session_start = time.monotonic()
    last_report_time, last_report_step = session_start, sim.timestep
    next_report_time = session_start + FIRST_REPORT_SECONDS

    chunk_steps = FIRST_CHUNK_STEPS
    while sim.timestep < total and not stop_requests:
        if sim.timestep >= equilibration and writer not in sim.operations.writers:
            sim.operations.writers.append(writer)  # production starts: begin saving frames
        end_of_phase = equilibration if sim.timestep < equilibration else total
        steps = min(chunk_steps, end_of_phase - sim.timestep)
        chunk_start = time.monotonic()
        sim.run(steps)
        if steps == chunk_steps:  # a full chunk: size the next one to take about CHUNK_SECONDS
            elapsed = max(time.monotonic() - chunk_start, 1e-3)
            chunk_steps = max(FIRST_CHUNK_STEPS, int(steps / elapsed * CHUNK_SECONDS))

        if time.monotonic() >= next_report_time:
            writer.flush()
            save_checkpoint(sim, checkpoint)
            clock = time.monotonic()
            speed = (sim.timestep - last_report_step) / (clock - last_report_time)
            phase = "equilibration" if sim.timestep < equilibration else "production"
            update(row, step=sim.timestep, phase=phase, steps_per_second=round(speed, 1),
                   run_time_s=round(previous_run_time + clock - session_start), last_update=now())
            print(f"[{now()}] {phase:13s} {100 * sim.timestep / total:5.1f}%  "
                  f"{speed:,.0f} steps/s  ~{(total - sim.timestep) / speed / 3600:.1f} h left", flush=True)
            last_report_time, last_report_step = clock, sim.timestep
            next_report_time = clock + REPORT_EVERY_MINUTES * 60
    writer.flush()
    run_time = round(previous_run_time + time.monotonic() - session_start)

    # ── Asked to stop: save everything, so ./noodle launch resumes from here ──
    if sim.timestep < total:
        save_checkpoint(sim, checkpoint)
        phase = "equilibration" if sim.timestep < equilibration else "production"
        update(row, status="stopped", stop_reason=stop_requests[0], step=sim.timestep, phase=phase,
               run_time_s=run_time, last_update=now())
        print(f"Stopped ({stop_requests[0]}) at {100 * sim.timestep / total:.1f}%, checkpoint saved: "
              "./noodle launch resumes it from here")
        return "stopped"

    # ── Done ──
    with gsd.hoomd.open(trajectory, "r") as f:
        frames = len(f)
    update(row, status="completed", phase="done", step=sim.timestep, frames=frames, run_time_s=run_time,
           completed_at=now(), last_update=now())
    checkpoint.unlink(missing_ok=True)
    print(f"Done: {frames} frames in {run_time / 3600:.2f} h of run time")
    return "completed"


def analyse(run_dir, job, task, physics, row, update):
    """Measure what the finished task's trajectory shows (observables.py) and write the
    averages into its row of the job file."""
    sequence = task["sequence"]
    save_every = run_lengths(physics, len(sequence))[2]
    started = time.monotonic()
    averages = analyse_task(trajectory_path(run_dir, task), sequence, physics,
                            system_geometry(sequence, physics), save_every * physics["dt"])
    update(row, analysed_at=now(), analysis_error=None, **averages)
    print(f"Analysed in {time.monotonic() - started:.0f} s: Rg {averages['rg_nm']:.2f} nm, "
          f"{averages['independent_samples']:.0f} independent samples", flush=True)


# ── The job: its tasks, one after another ─────────────────────────────────────
def run_job(run_dir, job):
    """Run every unfinished task of the job, in the order of its job file. Returns
    "stopped" if the job was asked to stop, else "done"."""
    tasks = {task["task_id"]: task for task in read_tasks(run_dir)}
    physics = yaml.safe_load((run_dir / "physics.yaml").read_text())
    rows = read_job(run_dir, job)

    def update(row, **fields):
        row.update(fields)
        write_job(run_dir, job, rows)

    print(f"Job {job}: {len(rows)} task(s) on {socket.gethostname()}, HOOMD {hoomd.version.version}, "
          f"started {now()}", flush=True)
    for n, row in enumerate(rows, 1):
        if stop_requests:
            break
        if row["status"] == "completed":
            continue
        task = tasks[row["task_id"]]
        print(f"\n── task {task['task_id']} ({n} of {len(rows)} in job {job}) ──", flush=True)
        try:
            if run_task(run_dir, task, physics[task["physics"]], row, update) == "stopped":
                break
        except KeyboardInterrupt:  # Ctrl-C twice: no checkpoint now, the last one (at most 10 minutes old) stays
            update(row, status="stopped", stop_reason="Ctrl-C twice (resumes from the last checkpoint)",
                   last_update=now())
            raise
        except Exception as error:
            # Recorded for ./noodle status and ./noodle launch; the full message is printed
            # (into this job's .out file), and the job goes on with its next task.
            traceback.print_exc()
            update(row, status="failed", error=f"{type(error).__name__}: {error}", last_update=now())
            continue
        try:
            analyse(run_dir, job, task, physics[task["physics"]], row, update)
        except Exception as error:  # the simulation itself is fine: ./noodle analyze tries again
            traceback.print_exc()
            update(row, analysis_error=f"{type(error).__name__}: {error}")
    return "stopped" if stop_requests else "done"


def main():
    for signal_number in STOP_REASONS:
        signal.signal(signal_number, request_stop)
    try:
        result = run_job(Path(sys.argv[1]), int(sys.argv[2]))
    except KeyboardInterrupt:
        result = "stopped"
    if result == "stopped":
        sys.exit(STOPPED_EXIT_CODE)


if __name__ == "__main__":
    main()
