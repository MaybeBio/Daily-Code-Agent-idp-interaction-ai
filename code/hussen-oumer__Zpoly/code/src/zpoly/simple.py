"""The simple way"""

# rq-b39965b3 rq-735d2b77 rq-9cc74219

from __future__ import annotations

import os
import sys
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from .forcefield import Topology, export_dft_inputs, parametrize, write_gromacs
from .forcefield._dft import dft_charges_from_values
from .structure import PolymerStructure, build_polymer, resolve_units, write_pdb
from .structure._units import variant
from .zwitterions import Library, Zwitterion, default_library


@dataclass(frozen=True, eq=False)
class BuildResult:
    pdb: Path
    top: Path | None
    ff_dir: Path | None
    structure: PolymerStructure
    topology: Topology | None
    net_charge: int
    counterion: str | None
    dft_blocks: tuple[str, ...]


# rq-82b577fa
def dft_inputs(monomer: Zwitterion | str, *, initiator: str | None = None, folder: str | os.PathLike = "dft",
               spacer: int | None = None, alkyl: int | None = None, parts: str = "whole",
               library: Library | None = None) -> list[Path]:
    """Write the DFT input molecules."""
    import yaml

    lib = library if library is not None else default_library()
    unit = variant(monomer, spacer, alkyl, lib)
    manifest = export_dft_inputs(unit, folder, initiator=initiator, library=lib, parts=parts)
    paths = []
    for mol in yaml.safe_load(manifest.read_text(encoding="utf-8"))["molecules"]:
        path = Path(folder) / mol["files"][0]
        paths.append(path)
        print(f"{path}: {len(mol['atoms'])} atoms, charge {mol['net_charge']}, multiplicity {mol['multiplicity']}")
    return paths


# rq-567ac589
def build(monomer, n: int | None = None, *, initiator: str | None = None, ends: str | None = None,
          spacer: int | None = None, alkyl: int | None = None, dft=None, ff: str = "gaff2",
          out: str | os.PathLike = "polymer", seed: int = 0, library: Library | None = None,
          dft_parts: str = "whole", charge_scale: float = 1.0, cap: str = "CH3",
          builtin_charges: bool = True) -> BuildResult:
    """Build a chain; write .pdb (+ .top and .ff/)."""
    lib = library if library is not None else default_library()
    units, aliases = resolve_units(monomer, n, spacer, alkyl, lib)
    structure = build_polymer(units, initiator=initiator, ends=ends, seed=seed, library=lib, cap=cap)
    out = Path(out)
    stem = out.name[:-4] if out.name.endswith(".pdb") else out.name
    folder = out.parent
    chain = structure.chain
    if ff.lower() == "none":
        pdb = write_pdb(structure, folder / f"{stem}.pdb")
        top = ff_dir = topology = None
        dft_blocks: tuple[str, ...] = ()
    else:
        if isinstance(dft, Mapping):
            dft = dft_charges_from_values(dft, structure.units, structure.initiator, aliases)
        topology = parametrize(structure, ff, name=stem, dft_charges=dft, dft_parts=dft_parts,
                               charge_scale=charge_scale, builtin_charges=builtin_charges)
        if topology.dft_missing:
            print(f"note: no DFT charges for {', '.join(topology.dft_missing)}; using {topology.force_field} charges",
                  file=sys.stderr)
        top = write_gromacs(topology, folder, stem)
        pdb, ff_dir, dft_blocks = folder / f"{stem}.pdb", folder / f"{stem}.ff", topology.dft_blocks
    counter = f", counterion {chain.counterion}" if chain.counterion else ""
    files = f"{pdb}" + (f", {top}, {ff_dir}/" if top else "")
    print(f"{len(units)} units, {structure.mol.GetNumAtoms()} atoms, net charge {chain.net_charge:+d}{counter}"
          f" -> {files}")
    return BuildResult(pdb, top, ff_dir, structure, topology, chain.net_charge, chain.counterion, dft_blocks)
