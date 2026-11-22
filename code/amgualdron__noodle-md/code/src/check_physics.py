"""Checks that the physics is what the docs say it is. Run by setup.sh; takes ~1 minute.

    env/bin/python src/check_physics.py

Prints PASS/FAIL per check and exits with an error if anything failed.
"""

import contextlib
import csv
import io
import math
import subprocess
import sys
import tempfile
from pathlib import Path

import gsd.hoomd
import hoomd
import numpy as np
import yaml

from analyze_run import analyze
from forcefield import (AA_PARAMS, AMINO_ACIDS, debye_huckel, residue_pair_forces, thermal_energy,
                        wall_contact_distance, wall_forces)
from launch_run import resolve_physics, write_jobs, write_manifest
from runs import ROOT, brush_size, read_csv, read_tasks, task_states, trajectory_path
from simulate import system_geometry

ALL_SETS = yaml.safe_load((ROOT / "physics.yaml").read_text())
failures = []


def check(description, ok):
    print(f"  {'PASS' if ok else 'FAIL'}  {description}")
    if not ok:
        failures.append(description)


def wall_energy_and_force(physics, residue, height):
    """Energy and z-force on one residue at `height` above the surface (at z = -15)."""
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    snap = hoomd.Snapshot()
    snap.configuration.box = [20, 20, 40, 0, 0, 0]
    snap.particles.N = 1
    snap.particles.types = AMINO_ACIDS
    snap.particles.typeid[0] = AMINO_ACIDS.index(residue)
    snap.particles.position[0] = [0, 0, -15 + height]
    sim.create_state_from_snapshot(snap)
    forces = wall_forces(physics, z_wall=-15, z_ceiling=15)
    sim.operations.integrator = hoomd.md.Integrator(dt=0.001, forces=forces)
    sim.run(0)
    return sum(f.energy for f in forces), sum(f.forces[0][2] for f in forces)


def pair_energy(physics, a, b, r):
    """Energy of residues a and b at distance r, from the forces of forcefield.py."""
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    snap = hoomd.Snapshot()
    snap.configuration.box = [20, 20, 20, 0, 0, 0]
    snap.particles.N = 2
    snap.particles.types = AMINO_ACIDS
    snap.particles.typeid[:] = [AMINO_ACIDS.index(a), AMINO_ACIDS.index(b)]
    snap.particles.charge[:] = [AA_PARAMS[a]["charge"], AA_PARAMS[b]["charge"]]
    snap.particles.position[:] = [[0, 0, 0], [r, 0, 0]]
    snap.bonds.types = ["AA-bond"]
    sim.create_state_from_snapshot(snap)
    forces = residue_pair_forces(physics, hoomd.md.nlist.Cell(buffer=0.2, exclusions=["bond"]))
    sim.operations.integrator = hoomd.md.Integrator(dt=0.001, forces=forces)
    sim.run(0)
    return sum(f.energy for f in forces)


def hps_formula(physics, a, b, r):
    """The HPS pair energy written out (forcefield.residue_pair_forces gives the formulas)."""
    eps = physics["epsilon_lj"]
    sigma = (AA_PARAMS[a]["sigma"] + AA_PARAMS[b]["sigma"]) / 2
    lam = (AA_PARAMS[a]["hydropathy"] + AA_PARAMS[b]["hydropathy"]) / 2
    lj = 4 * eps * ((sigma / r) ** 12 - (sigma / r) ** 6)
    if r < 2 ** (1 / 6) * sigma:
        u = lj + (1 - lam) * eps
    else:
        u = lam * lj if r < physics["lj_r_cut"] else 0.0
    bjerrum, kappa = debye_huckel(physics["temperature_K"], physics["ionic_concentration_M"])
    qq = AA_PARAMS[a]["charge"] * AA_PARAMS[b]["charge"]
    rc = physics["yukawa_r_cut"]
    if r < rc:
        u += qq * bjerrum * thermal_energy(physics["temperature_K"]) * (np.exp(-kappa * r) / r - np.exp(-kappa * rc) / rc)
    return u


def main():
    grafted = resolve_physics(ALL_SETS, "grafted")
    baseline = resolve_physics(ALL_SETS, "baseline")

    print("Residue pairs (HPS: Ashbaugh-Hatch + Debye-Huckel)")
    worst = 0.0
    for a, b in (("K", "E"), ("K", "K"), ("R", "D"), ("H", "E"), ("A", "L"), ("R", "R"), ("K", "A"), ("W", "Y")):
        sigma = (AA_PARAMS[a]["sigma"] + AA_PARAMS[b]["sigma"]) / 2
        for r in (0.9 * sigma, sigma, 1.2 * sigma, 1.5, 2.5, 3.4):
            worst = max(worst, abs(pair_energy(baseline, a, b, r) - hps_formula(baseline, a, b, r)))
    check(f"8 residue pairs at 6 distances match the formulas (largest difference {worst:.1e} kcal/mol)",
          worst < 1e-9)
    check("beyond lj_r_cut, a neutral pair does not interact; a charged pair still does",
          pair_energy(baseline, "K", "A", 2.5) == 0 and pair_energy(baseline, "K", "E", 2.5) < 0)

    print("Brush geometry")
    brush = resolve_physics(ALL_SETS, "brush")
    chains, width = brush_size(brush)
    check(f"brush: {chains} chains in a box {width:.2f} nm wide: {brush['grafting_density']} chains/nm^2",
          chains == brush["chains"] and abs(chains / width**2 - brush["grafting_density"]) < 1e-12)
    densest = dict(brush, grafting_density=1.0)
    chains, width = brush_size(densest)
    check(f"densest allowed brush (1 chain/nm^2) keeps its {brush['chains']} chains: box {width:.2f} nm, at least "
          f"twice the interaction range", chains == brush["chains"] and width >= 2 * (brush["yukawa_r_cut"] + 0.5))
    small_dense = dict(brush, chains=16, grafting_density=0.5)
    chains, width = brush_size(small_dense)
    check(f"a 4 x 4 brush at 0.5 chains/nm^2 grows to {chains} chains, so its box ({width:.2f} nm) is at least "
          f"twice the interaction range", chains > 16 and width >= 2 * (brush["yukawa_r_cut"] + 0.5)
          and math.isqrt(chains) ** 2 == chains)
    check("single grafted chain: 1 chain in a box box_size wide", brush_size(grafted) == (1, grafted["box_size"]))

    print("Electrostatics")
    bjerrum, kappa = debye_huckel(298, 0.150)
    check(f"Bjerrum length at 298 K is {bjerrum:.3f} nm (water: ~0.71)", 0.69 < bjerrum < 0.73)
    check(f"Debye length at 150 mM is {1 / kappa:.3f} nm (expected ~0.79)", 0.76 < 1 / kappa < 0.81)

    print("Neutral, purely repulsive surface (the default)")
    for potential in ("lj93", "lj126"):
        physics = dict(grafted, wall_potential=potential)
        for residue in ("G", "W"):
            contact = wall_contact_distance(AA_PARAMS[residue]["sigma"], physics)
            heights = np.linspace(0.8 * contact, contact, 20)
            energies, forces = zip(*(wall_energy_and_force(physics, residue, h) for h in heights))
            u_far, f_far = wall_energy_and_force(physics, residue, contact * 1.01)
            check(f"{potential} {residue}: zero energy and force beyond contact ({contact:.3f} nm)",
                  u_far == 0 and f_far == 0)
            check(f"{potential} {residue}: energy only rises toward the surface, force pushes away",
                  all(np.diff(energies) <= 0) and all(f >= 0 for f in forces))
            check(f"{potential} {residue}: energy goes smoothly to zero at contact",
                  energies[-1] < 1e-6 and energies[-2] < 1e-3)
        _, f_below = wall_energy_and_force(physics, "G", -0.05)
        check(f"{potential}: a bead that slipped past the surface is pushed back up", f_below > 0)

    print("Attractive surface (wall_attraction > 0)")
    for scale in ("uniform", "hydropathy"):
        physics = dict(grafted, wall_attraction=0.5, wall_attraction_scale=scale)
        lam = 0.5 if scale == "uniform" else 0.5 * AA_PARAMS["F"]["hydropathy"]
        contact = wall_contact_distance(AA_PARAMS["F"]["sigma"], physics)
        u_min, f_min = wall_energy_and_force(physics, "F", contact)
        expected = -physics["wall_epsilon"] * lam
        check(f"{scale}: well depth at contact is {u_min:.4f} kcal/mol (expected {expected:.4f})",
              abs(u_min - expected) < 0.02 * abs(expected))
        check(f"{scale}: no force at the bottom of the well", abs(f_min) < 1e-3)
    physics = dict(grafted, wall_attraction=0.5, wall_attraction_scale="hydropathy")
    contact = wall_contact_distance(AA_PARAMS["R"]["sigma"], physics)
    u_r, _ = wall_energy_and_force(physics, "R", contact + 0.3)
    check("hydropathy: arginine (hydropathy lambda = 0) feels no attraction", u_r == 0)

    print("Charged surface (Debye-Huckel)")
    physics = dict(grafted, wall_charge_density=-0.1)
    bjerrum, kappa = debye_huckel(physics["temperature_K"], physics["ionic_concentration_M"])
    kT = thermal_energy(physics["temperature_K"])
    for h in (1.0, 2.0):
        u, f = wall_energy_and_force(physics, "K", h)
        q_lysine = AA_PARAMS["K"]["charge"]
        analytic = 4 * np.pi * bjerrum * kT * q_lysine * (-0.1) / kappa * np.exp(-kappa * h)
        u_up, _ = wall_energy_and_force(physics, "K", h + 1e-4)
        u_down, _ = wall_energy_and_force(physics, "K", h - 1e-4)
        check(f"lysine at {h} nm: U = {u:.5f}, analytic {analytic:.5f} kcal/mol", abs(u - analytic) < 1e-9)
        check(f"lysine at {h} nm: force = -dU/dh", abs(f + (u_up - u_down) / 2e-4) < 1e-5)
    check("negative surface attracts lysine and repels glutamate",
          wall_energy_and_force(physics, "K", 1.0)[1] < 0 < wall_energy_and_force(physics, "E", 1.0)[1])

    print("Short simulations (bulk; grafted with attraction and charge; grafted_charged; brushes)")
    with tempfile.TemporaryDirectory() as tmp:
        run_dir = Path(tmp)
        sets = {"bulk": resolve_physics(ALL_SETS, "test_bulk"),
                "grafted": dict(resolve_physics(ALL_SETS, "test_grafted"),
                                wall_attraction=0.3, wall_charge_density=-0.05)}
        sets["grafted_C_terminus"] = dict(sets["grafted"], graft_terminus="C")
        sets["grafted_charged"] = resolve_physics(ALL_SETS, "test_grafted_charged")
        sets["brush"] = resolve_physics(ALL_SETS, "test_brush")
        sets["brush_dense_C_terminus"] = dict(sets["brush"], grafting_density=0.4, graft_terminus="C",
                                              wall_attraction=0.3, wall_charge_density=-0.05)
        (run_dir / "physics.yaml").write_text(yaml.safe_dump(sets))
        sequence = "MSEQKKPLGEGFDNRYSTWAVHIC"
        names = list(sets)
        with open(run_dir / "tasks.csv", "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["task_id", "sequence_name", "physics", "replicate", "seed", "length", "sequence"])
            for i, name in enumerate(names, 1):
                writer.writerow([i, "test", name, 1, 42, len(sequence), sequence])
        with contextlib.redirect_stdout(io.StringIO()):
            write_manifest(run_dir, {"sequences_csv": "(check)"}, {"test": sequence}, {"test": 0}, sets)
        # all of them as one job, as ./noodle launch groups short tasks
        tasks = read_tasks(run_dir)
        write_jobs(run_dir, [[task["task_id"] for task in tasks]], tasks, {})
        result = subprocess.run([sys.executable, str(ROOT / "src" / "simulate.py"), str(run_dir), "1"],
                                capture_output=True, text=True)
        check("one job runs every task", result.returncode == 0)
        states = task_states(run_dir)
        for task in tasks:
            name, log = task["physics"], states[task["task_id"]]
            check(f"{name}: task completes", log["status"] == "completed")
            if log["status"] != "completed":
                print(result.stdout[-3000:], result.stderr[-2000:])
                continue
            with gsd.hoomd.open(trajectory_path(run_dir, task)) as traj:
                positions = np.array([frame.particles.position for frame in traj])
            check(f"{name}: {len(positions)} frames, all positions finite", np.isfinite(positions).all())
            if name.startswith(("grafted", "brush")):
                geometry = system_geometry(sequence, sets[name])
                anchors = positions[:, geometry["anchors"], :]
                check(f"{name}: {geometry['chains']} chain(s) of {len(sequence)} residues, anchors never moved",
                      positions.shape[1] == geometry["chains"] * len(sequence) and np.all(anchors == anchors[0]))
                check(f"{name}: every bead stayed between surface and ceiling",
                      positions[:, :, 2].min() > geometry["z_wall"] and positions[:, :, 2].max() < geometry["z_ceiling"])

        print("Analysis of those simulations")
        check("each task was analysed by the job right after it finished",
              all(states[task["task_id"]].get("analysed_at") for task in tasks))
        with contextlib.redirect_stdout(io.StringIO()):
            analyze(run_dir)
        manifest = {row["physics"]: row for row in read_csv(run_dir / "manifest.csv")}
        check("one manifest row per simulation, Rg of every chain between 0.5 and 5 nm",
              set(manifest) == set(names) and all(0.5 < (row["rg_nm"] or 0) < 5 for row in manifest.values()))
        check("heights only for grafted chains, brush columns only for brushes",
              all((row["com_height_nm"] is None) == (row["geometry"] == "bulk")
                  and (row["brush_height_nm"] is None) == (row["chains"] == 1) for row in manifest.values()))
        for name in ("brush", "brush_dense_C_terminus"):
            row = manifest[name]
            check(f"{name}: brush height, concentration and lateral heterogeneity are positive, "
                  "interchain contacts are not negative",
                  all((row[c] or 0) > 0 for c in ("brush_height_nm", "concentration_mg_ml", "lateral_heterogeneity"))
                  and row["interchain_contacts_per_residue"] is not None
                  and row["interchain_contacts_per_residue"] >= 0)
        check("brush: reduced density = grafting density x pi <Rg^2> of the free chain",
              math.isclose(manifest["brush"]["reduced_density"],
                           sets["brush"]["grafting_density"] * math.pi * manifest["bulk"]["rg_sq_nm2"], rel_tol=1e-4))

    print()
    if failures:
        print(f"{len(failures)} check(s) FAILED.")
        sys.exit(1)
    print("All physics checks passed.")


if __name__ == "__main__":
    main()
