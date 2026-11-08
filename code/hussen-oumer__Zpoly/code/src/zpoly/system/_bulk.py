"""Dry bulk system for density and Tg"""

# rq-4a8e3ff0 rq-d8a94b51

from __future__ import annotations

import os
from pathlib import Path

import yaml

from ..forcefield import parametrize
from ..forcefield._dft import dft_charges_from_values
from ..structure import build_polymer, resolve_units
from ..zwitterions import Library, default_library
from ._assemble import (
    AMU_G,
    NAME,
    box_angstrom,
    chain_mass,
    chain_piece,
    counterions_for,
    density_g_cm3,
    ion_piece,
    solvent_parts,
    total_charge,
    write_forcefield,
    write_gro,
    write_system_top,
)
from ._errors import SystemBuildError
from ._pack import pack
from ._protocol import Step, mdp_text, run_script, slurm_script

DEFAULT_TIMES = {"nvt_hot": 1.0, "compress": 2.0, "anneal": 5.0, "cool": 2.0, "deform": 1.0}
PACK_DENSITY = 0.35  # g/cm3, density mode: packed at this density, then GROMACS shrinks the box (deform)


def temperature_ladder(t_high: float, t_low: float, t_step: float) -> list[float]:
    if not (t_high > t_low and t_step > 0):
        raise ValueError("Need t_high > t_low and t_step > 0.")
    temps, t = [], t_high
    while t >= t_low - 1e-9:
        temps.append(round(t, 6))
        t -= t_step
    return temps


def _tname(t: float) -> str:
    return f"cool_{int(t)}K" if float(t).is_integer() else f"cool_{t:g}K"


def build_chains(units, chains, initiator, ends, seed, lib, cap="CH3"):
    """`chains` copies, each with its own seed"""
    structures = [build_polymer(units, initiator=initiator, ends=ends, seed=seed + k, library=lib, cap=cap)
                  for k in range(chains)]
    names0 = [(a.GetProp("res"), a.GetProp("name")) for a in structures[0].mol.GetAtoms()]
    for s in structures[1:]:
        if [(a.GetProp("res"), a.GetProp("name")) for a in s.mol.GetAtoms()] != names0:
            raise SystemBuildError("Internal error: chain copies differ in their atoms.")
    return structures


# rq-d8a94b51
def bulk_system(monomer, n: int | None = None, *, chains: int, box_nm=None, density: float | None = None,
                initiator: str | None = None,
                ends: str | None = None, spacer: int | None = None, alkyl: int | None = None, dft=None,
                dft_parts: str = "whole", charge_scale: float = 1.0, cap: str = "CH3",
                builtin_charges: bool = True, solvent: str = "tip3p",
                ff: str = "gaff2", packer: str = "packmol", t_high: float = 600.0,
                t_low: float = 200.0, t_step: float = 20.0, times_ns: dict | None = None,
                pressure_bar: float = 1.0, compress_bar: float = 1000.0, seed: int = 0,
                solvent_ff: str | os.PathLike | None = None, out: str | os.PathLike = "bulk",
                library: Library | None = None) -> Path:
    """Build a dry bulk system and its Tg protocol."""
    if not isinstance(chains, int) or chains < 1:
        raise ValueError(f"chains must be an integer >= 1, got {chains!r}")
    if (box_nm is None) == (density is None):
        raise ValueError("Give either density (g/cm3) or box_nm")
    if density is not None and (isinstance(density, bool) or not isinstance(density, (int, float)) or density <= 0):
        raise ValueError(f"density must be > 0 g/cm3, got {density!r}")
    box = box_angstrom(box_nm) if box_nm is not None else None
    times = dict(DEFAULT_TIMES)
    unknown = set(times_ns or {}) - set(times)
    if unknown:
        raise ValueError(f"Unknown times_ns key(s) {', '.join(sorted(unknown))}; use {', '.join(times)}")
    times.update(times_ns or {})
    temps = temperature_ladder(t_high, t_low, t_step)
    if solvent == "tip4p2005" and abs(charge_scale - 0.85) > 1e-9:
        raise ValueError("SOLVENT = tip4p2005 uses the Madrid-2019 ions (charges +-0.85 e): "
                         "set CHARGE_SCALE = 0.85 so the system is neutral")
    lib = library if library is not None else default_library()
    units, aliases = resolve_units(monomer, n, spacer, alkyl, lib)

    # rq-758b974f: emc copies one chain
    structures = build_chains(units, 1 if packer == "emc" else chains, initiator, ends, seed, lib, cap)
    structures = structures * chains if packer == "emc" else structures
    first = structures[0]
    if isinstance(dft, dict):
        dft = dft_charges_from_values(dft, first.units, first.initiator, aliases)
    topology = parametrize(first, ff, name=NAME, dft_charges=dft, dft_parts=dft_parts, charge_scale=charge_scale,
                           builtin_charges=builtin_charges)
    if topology.missing_parameters:
        raise SystemBuildError("The topology has missing parameters.")

    ion, n_ions = counterions_for(first, chains)
    parts = solvent_parts(topology, [ion] if ion else [], solvent_ff, solvent=solvent) if ion else None
    mass = chain_mass(topology) * chains + (parts.masses[ion] * n_ions if ion else 0.0)
    pack_box = box
    if density is not None:  # rq-7c268316: cube edge from the density; packed looser, then deformed to it
        edge = (mass * AMU_G / density * 1e24) ** (1 / 3)
        box = (edge, edge, edge)
        loose = (mass * AMU_G / min(density, PACK_DENSITY) * 1e24) ** (1 / 3)
        pack_box = (loose, loose, loose)
    start_density = density_g_cm3(mass, box[0] * box[1] * box[2])
    pieces = ([chain_piece(first, emc=True)] * chains if packer == "emc" else [chain_piece(s) for s in structures]
              ) + ([ion_piece(parts.molecules[ion])] * n_ions if ion else [])
    placed = pack(pieces, pack_box, packer, seed)
    q = total_charge(topology, chains, parts, {ion: n_ions} if ion else {})
    if abs(q) > 1e-3:
        raise SystemBuildError(f"The system's total charge is {q:+.4f}, not 0.")

    folder = Path(out)
    folder.mkdir(parents=True, exist_ok=True)
    write_forcefield(folder, topology, parts)
    title = f"{chains} x {'-'.join(dict.fromkeys(z.name for z in units))} ({len(units)} units) bulk"
    write_system_top(folder, title, chains, [(ion, n_ions)] if ion else [], bool(ion), False, charge_scale)
    extra = [(parts.molecules[ion], xyz) for xyz in placed[chains:]] if ion else []
    natoms = write_gro(folder / "system.gro", title, first, placed[:chains], extra, pack_box)

    steps = [Step("01_em", "em")]
    if pack_box != box:  # heat (velocities), then shrink the cube to the requested density at t_high
        steps += [Step("01_heat", "nvt", t_high, 0.01, generate_velocities=True),
                  Step("01_deform", "nvt", t_high, times["deform"],
                       deform_nm_ps=(box[0] - pack_box[0]) / 10 / (times["deform"] * 1000))]
    steps += [Step("02_nvt_hot", "nvt", t_high, times["nvt_hot"], generate_velocities=True),
             Step("03_npt_compress", "npt", t_high, times["compress"], compress_bar),
             Step("04_npt_anneal", "npt", t_high, times["anneal"], pressure_bar)]
    steps += [Step(_tname(t), "npt", t, times["cool"], pressure_bar) for t in temps]
    for step in steps:
        (folder / f"{step.name}.mdp").write_text(mdp_text(step), encoding="ascii")
    (folder / "run.sh").write_text(run_script(steps), encoding="ascii")
    total_ns = sum(s.time_ns for s in steps)
    (folder / "run_slurm.sh").write_text(slurm_script(folder.name or "bulk", total_ns), encoding="ascii")
    for script in ("run.sh", "run_slurm.sh"):
        (folder / script).chmod(0o755)
    protocol = {
        "kind": "bulk", "title": title, "force_field": topology.force_field, "charge_scale": charge_scale,
        "solvent": solvent,
        "chains": chains,
        "units": [z.name for z in units], "initiator": first.chain.initiator, "counterion": ion,
        "counterions": n_ions, "atoms": natoms, "box_nm": [round(b / 10, 5) for b in box],
        "start_density_g_cm3": round(start_density, 4), "packer": packer, "seed": seed,
        "density_g_cm3": density, "packing_box_nm": [round(b / 10, 5) for b in pack_box],
        "dft_blocks": list(topology.dft_blocks),
        "steps": [{"name": s.name, "kind": s.kind, "temperature": s.temperature, "time_ns": s.time_ns,
                   "pressure_bar": s.pressure_bar} for s in steps],
        "cooling": {"steps": [{"name": _tname(t), "temperature": t, "time_ns": times["cool"]} for t in temps],
                    "average_last_ns": min(1.0, times["cool"])},
        "total_ns": total_ns,
    }
    (folder / "protocol.yaml").write_text(yaml.safe_dump(protocol, sort_keys=False), encoding="utf-8")
    dims = " x ".join(f"{b / 10:g}" for b in box)
    print(f"{folder}: {chains} chains + {n_ions} {ion or 'counterions'}, {natoms} atoms, box {dims} nm "
          f"(starting density {start_density:.3f} g/cm3), "
          f"{len(temps)} cooling steps ({total_ns:g} ns in total). Run {folder / 'run.sh'}, then analyze_tg.")
    return folder
