"""One input file of `KEY = value` lines drives every job (rqm/input_file.md)."""

# rq-8ecf5b05 rq-44600338

from __future__ import annotations

import difflib
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import yaml

from .forcefield import (
    OplsaaBackend,
    charge_template,
    export_dft_inputs,
    export_dft_molecule,
    ionic_parameter_template,
    parametrize,
    register_backend,
    write_gromacs,
)
from .structure import (
    alternating,
    block,
    build_polymer,
    homopolymer,
    pattern,
    random_copolymer,
    resolve_units,
    write_pdb,
)
from .system import analyze_slab, analyze_tg, bulk_system, slab_system
from .zwitterions import default_library, load_library

JOBS = ("build", "dft", "bulk", "slab", "analyze", "charges", "params")
POLYMER = ("build", "dft", "bulk", "slab", "charges", "params")
WITH_FF = ("build", "bulk", "slab", "charges", "params")
SYSTEM = ("bulk", "slab")
_UNSAFE = re.compile(r"[^A-Za-z0-9._-]")


class InputError(ValueError):
    def __init__(self, message: str, line: int | None = None, key: str | None = None):
        super().__init__(message)
        self.line = line
        self.key = key


@dataclass(frozen=True)
class Key:
    name: str
    jobs: tuple[str, ...]
    kind: str  # str int float floats list choice
    short: str
    value: str
    default: object = None
    example: str = ""
    choices: tuple[str, ...] = ()
    count: tuple[int, ...] = ()  # allowed numbers of values for "floats"


# rq-f2d8deab
KEYS = {k.name: k for k in (
    Key("JOB", JOBS, "choice", "build | dft | bulk | slab | analyze | ...", " | ".join(JOBS), "build",
        "JOB = slab", JOBS),
    Key("MONOMER", POLYMER, "str", "monomer of a homopolymer (see: list)", "a monomer name", None,
        "MONOMER = DMAPS"),
    Key("SEQUENCE", POLYMER, "list", "explicit chain, monomer by monomer", "monomer names in chain order",
        None, "SEQUENCE = DMAPS DMAPS CBMA DMAPS"),
    Key("COPOLYMER", POLYMER, "list", "alternating / block / random / pattern",
        "alternating A B | block A:10 B:20 | random A:0.5 B:0.5 | pattern AAB A=X B=Y", None,
        "COPOLYMER = random DMAPS:0.5 CBMA:0.5"),
    Key("N", POLYMER, "int", "units per chain", "integer >= 1", None, "N = 12"),
    Key("INITIATOR", POLYMER, "str", "end groups (see: list initiators)", "an initiator name", None,
        "INITIATOR = KPS"),
    Key("ENDS", POLYMER, "choice", "head | tail | both | none", "both | head | tail | none",
        None, "ENDS = head", ("both", "head", "tail", "none")),
    Key("CAP", POLYMER, "choice", "other chain ends: CH3 | H", "CH3 | H", "CH3", "CAP = H",
        ("ch3", "h")),
    Key("BUILTIN_CHARGES", WITH_FF, "choice", "yes | no (DMAPS, A3361, A3367, KPS)",
        "yes | no: the group's RESP charges for DMAPS, A3361, A3367 and KPS", "yes", "BUILTIN_CHARGES = no",
        ("yes", "no")),
    Key("SPACER", POLYMER, "int", "CH2 groups between cation and anion", "integer >= 1", None, "SPACER = 4"),
    Key("ALKYL", POLYMER, "int", "N-alkyl length (1 methyl, 2 ethyl...)", "integer >= 1", None, "ALKYL = 2"),
    Key("CONFORMATION", ("build",), "choice", "starting shape of the chain", "coil | extended", "coil",
        "CONFORMATION = extended", ("coil", "extended")),
    Key("SEED", POLYMER, "int", "random seed", "integer", 0, "SEED = 1"),
    Key("FF", WITH_FF, "str", "oplsaa | gaff2 | none (PDB only)", "gaff2 | oplsaa | none (PDB only)",
        "gaff2", "FF = oplsaa"),
    Key("DFT", ("build", "dft", "bulk", "slab"), "str", "folder of your DFT charges",
        "folder; for JOB = dft the folder the molecules are written to", None, "DFT = dft"),
    Key("DFT_PARTS", ("build", "dft", "bulk", "slab"), "choice", "whole | ions (+ cation, anion, ester)",
        "whole | ions (cation/anion fragments for the ionic atoms)", "whole", "DFT_PARTS = ions",
        ("whole", "ions")),
    Key("CHARGE_SCALE", ("build", "bulk", "slab"), "float", "0.85 = ionic groups and ions x 0.85; 1 = none",
        "number in (0, 1]; scales ionic groups, counterions and salt", 1.0, "CHARGE_SCALE = 0.8"),
    Key("MOLECULE", ("dft",), "str", "zwitterionic liquid for JOB = dft", "a small-molecule name (see: list)",
        None, "MOLECULE = TMAO"),
    Key("CHARGES", ("build",), "str", "your per-atom charge file", "a file from JOB = charges", None,
        "CHARGES = DMAPS_n3_charges.yaml"),
    Key("PARAMS", ("build",), "str", "your bonded parameters (ionic part)", "a file from JOB = params", None,
        "PARAMS = DMAPS_n3_params.yaml"),
    Key("OPLSAA_FF", WITH_FF, "str", "your oplsaa.ff folder (FF = oplsaa)", "folder", None,
        "OPLSAA_FF = /opt/gromacs/share/gromacs/top/oplsaa.ff"),
    Key("OPLS_RULES", WITH_FF, "str", "extra OPLS-AA typing rules", "YAML file", None,
        "OPLS_RULES = my_rules.yaml"),
    Key("CHAINS", SYSTEM, "int", "number of chains (you choose)", "integer >= 1", None, "CHAINS = 50"),
    Key("BOX", SYSTEM, "floats", "box in nm: L or Lx Ly Lz (you choose)",
        "one edge (cube) or Lx Ly Lz, in nm; a slab needs Lx = Ly and Lz >= 3 Lx", None, "BOX = 8 8 24",
        count=(1, 3)),
    Key("SLAB", ("slab",), "float", "starting slab thickness, nm", "number (nm), 0 < SLAB <= Lz", None,
        "SLAB = 6"),
    Key("DENSITY", SYSTEM, "float", "start density, g/cm3 (the box follows)",
        "number (g/cm3): bulk: cube edge (mass/DENSITY)^(1/3); slab: cube Lx = Ly = slab, Lz = ASPECT x Lx",
        None, "DENSITY = 1.0"),
    Key("ASPECT", ("slab",), "float", "Lz / Lx with DENSITY (at least 3)", "number >= 3", 3.0, "ASPECT = 4"),
    Key("SALT", ("slab",), "list", "NaCl 0.15 | KCl 0.1 | none", "NaCl | KCl [molar] | none",
        ("NaCl", 0.15), "SALT = KCl 0.1"),
    Key("TEMPERATURES", ("slab",), "floats", "production temperatures, K", "numbers (K)",
        (280.0, 300.0, 320.0, 340.0, 360.0), "TEMPERATURES = 290 310 330"),
    Key("TEQUIL", ("slab",), "float", "equilibration temperature, K", "number (K)", 300.0, "TEQUIL = 300"),
    Key("TG", ("bulk",), "floats", "cooling: T_high T_low step, K", "three numbers (K)", (600.0, 200.0, 20.0),
        "TG = 600 200 20", count=(3,)),
    Key("TIMES", SYSTEM, "list", "step times, ns (name value ...)",
        "bulk: nvt_hot compress anneal cool; slab: nvt npt_z production", None, "TIMES = cool 2 anneal 5"),
    Key("SOLVENT", SYSTEM, "choice", "tip3p | tip4p2005 (+ Madrid-2019 ions)",
        "tip3p | tip4p2005 (TIP4P/2005 water, Madrid-2019 ions; needs CHARGE_SCALE = 0.85)", "tip3p",
        "SOLVENT = tip4p2005", ("tip3p", "tip4p2005")),
    Key("PACKER", SYSTEM, "choice", "emc | builtin | packmol", "emc (EMC-style Monte Carlo) | builtin | packmol",
        "packmol", "PACKER = emc", ("packmol", "builtin", "emc")),
    Key("SOLVENT_FF", SYSTEM, "str", "GROMACS .ff folder for water/ions", "folder", None,
        "SOLVENT_FF = /opt/gromacs/share/gromacs/top/amber99sb-ildn.ff"),
    Key("OUTPUT", tuple(j for j in JOBS if j != "dft"), "str",
        "output name; for analyze the folder", "name or path (default: from the polymer)", None,
        "OUTPUT = dmaps12"),
    Key("LIBRARY", JOBS, "list", "extra monomer/initiator YAML files", "file names", None,
        "LIBRARY = my_monomers.yaml"),
)}


def suggest(name: str) -> str:
    close = difflib.get_close_matches(name.upper(), KEYS, n=3, cutoff=0.6)
    return f" (did you mean {', '.join(close)}?)" if close else ""


def _convert(key: Key, raw: str, where: str):
    tokens = [t for t in re.split(r"[,\s]+", raw.strip()) if t]
    if not tokens:
        raise InputError(f"{where}: {key.name} has no value; expected {key.value}", key=key.name)
    try:
        if key.kind == "int":
            if len(tokens) != 1:
                raise ValueError
            return int(tokens[0])
        if key.kind == "float":
            if len(tokens) != 1:
                raise ValueError
            return float(tokens[0])
        if key.kind == "floats":
            values = tuple(float(t) for t in tokens)
            if key.count and len(values) not in key.count:
                raise ValueError
            return values
    except ValueError:
        need = {"int": "an integer", "float": "a number"}.get(key.kind, key.value)
        raise InputError(f"{where}: {key.name} needs {need}, got {raw.strip()!r}", key=key.name) from None
    if key.kind == "choice":
        value = raw.strip().lower()
        if value not in key.choices:
            raise InputError(f"{where}: {key.name} must be one of {', '.join(key.choices)}, got {raw.strip()!r}",
                             key=key.name)
        return value
    if key.kind == "list":
        return tuple(tokens)
    return raw.strip()


# rq-b18197c0 rq-cfe286dc
def parse_input(text: str, source: str = "<input>") -> dict:
    """Parse and check an input file"""
    found: dict[str, tuple[int, object]] = {}
    for lineno, line in enumerate(text.splitlines(), 1):
        body = re.split(r"[#!]", line, maxsplit=1)[0].strip()
        if not body:
            continue
        where = f"{source}, line {lineno}"
        name, sep, raw = body.partition("=")
        name = name.strip().upper()
        if not sep or not name:
            raise InputError(f"{where}: expected KEY = value, got {body!r}", line=lineno)
        if name not in KEYS:
            raise InputError(f"{where}: unknown keyword {name}{suggest(name)}", line=lineno, key=name)
        if name in found:
            raise InputError(f"{where}: {name} is already given on line {found[name][0]}", line=lineno, key=name)
        found[name] = (lineno, _convert(KEYS[name], raw, where))

    job = found.get("JOB", (0, "build"))[1]
    for name, (lineno, _) in found.items():
        if job not in KEYS[name].jobs:
            raise InputError(f"{source}, line {lineno}: {name} is for JOB = {' / '.join(KEYS[name].jobs)}, "
                             f"not {job}", line=lineno, key=name)
    opts = {k.name.lower(): k.default for k in KEYS.values() if job in k.jobs}
    opts.update({name.lower(): value for name, (_, value) in found.items()})
    opts["job"] = job
    _check(opts, job, source)
    return opts


def _check(opts: dict, job: str, source: str) -> None:
    scale = opts.get("charge_scale")
    if scale is not None and not 0 < scale <= 1:
        raise InputError(f"{source}: CHARGE_SCALE must be in (0, 1], got {scale:g}", key="CHARGE_SCALE")
    if job == "dft" and opts.get("molecule") is not None:
        given = [k.upper() for k in ("monomer", "sequence", "copolymer", "n", "initiator") if opts.get(k) is not None]
        if given:
            raise InputError(f"{source}: MOLECULE is a single molecule; remove {', '.join(given)}", key="MOLECULE")
    elif job in POLYMER:
        given = [k.upper() for k in ("monomer", "sequence", "copolymer") if opts.get(k) is not None]
        if len(given) != 1:
            raise InputError(f"{source}: give exactly one of MONOMER, SEQUENCE, COPOLYMER"
                             + (f" (found {', '.join(given)})" if given else ""))
        needs_n = opts.get("monomer") is not None or (
            opts.get("copolymer") is not None and opts["copolymer"][0].lower() != "block")
        if needs_n and opts.get("n") is None:
            raise InputError(f"{source}: N (units per chain) is required with {given[0]}", key="N")
        if not needs_n and opts.get("n") is not None:
            raise InputError(f"{source}: N is not used with {given[0]}"
                             + (" block" if given[0] == "COPOLYMER" else "") + "; the chain length follows from it",
                             key="N")
    by_density = opts.get("density") is not None
    if job in SYSTEM and by_density and (opts.get("box") or opts.get("slab")):
        raise InputError(f"{source}: give either DENSITY or BOX (and SLAB), not both", key="DENSITY")
    required = {"bulk": ("chains",) if by_density else ("chains", "box"),
                "slab": ("chains",) if by_density else ("chains", "box", "slab"), "analyze": ("output",)}
    missing = [k.upper() for k in required.get(job, ()) if opts.get(k) is None]
    if missing:
        raise InputError(f"{source}: JOB = {job} needs {', '.join(missing)} (you choose these)")
    ff = (opts.get("ff") or "").lower()
    if ff == "none" and job in ("bulk", "slab", "params"):
        raise InputError(f"{source}: JOB = {job} needs a force field (FF = gaff2 or oplsaa)", key="FF")
    if ff == "none" and job == "build":
        bad = [k.upper() for k in ("dft", "charges", "params") if opts.get(k) is not None]
        if bad:
            raise InputError(f"{source}: {', '.join(bad)} need a force field (FF = gaff2 or oplsaa)", key="FF")
    if (opts.get("oplsaa_ff") or opts.get("opls_rules")) and ff != "oplsaa":
        raise InputError(f"{source}: OPLSAA_FF and OPLS_RULES need FF = oplsaa", key="FF")


def _pairs(items, kind, what: str):
    out = []
    for item in items:
        name, sep, value = item.rpartition(":")
        if not sep or not name:
            raise InputError(f"COPOLYMER: expected MONOMER:{what}, got {item!r}", key="COPOLYMER")
        try:
            out.append((name, kind(value)))
        except ValueError:
            raise InputError(f"COPOLYMER: {value!r} in {item!r} is not a valid {what.lower()}",
                             key="COPOLYMER") from None
    return out


def sequence_and_stem(opts: dict) -> tuple[list[str], str]:
    """The chain as monomer names, and the default output name"""
    n, seed = opts.get("n"), opts.get("seed", 0)
    if opts.get("monomer") is not None:
        seq, stem = homopolymer(opts["monomer"], n), f"{opts['monomer']}_n{n}"
    elif opts.get("sequence") is not None:
        seq = list(opts["sequence"])
        stem = f"seq_n{len(seq)}"
    else:
        form, *args = opts["copolymer"]
        form = form.lower()
        if form == "alternating" and len(args) == 2:
            seq, stem = alternating(args[0], args[1], n), f"{args[0]}-alt-{args[1]}_n{n}"
        elif form == "block" and args:
            blocks = _pairs(args, int, "COUNT")
            seq, stem = block(blocks), "-".join(f"{m}{c}" for m, c in blocks)
        elif form == "random" and args:
            comp = dict(_pairs(args, float, "FRACTION"))
            seq, stem = random_copolymer(comp, n, seed=seed), "-ran-".join(comp) + f"_n{n}_s{seed}"
        elif form == "pattern" and len(args) >= 2:
            mapping = {}
            for item in args[1:]:
                sym, sep, mono = item.partition("=")
                if not sep or len(sym) != 1 or not mono:
                    raise InputError(f"COPOLYMER: expected SYMBOL=MONOMER, got {item!r}", key="COPOLYMER")
                mapping[sym] = mono
            seq, stem = pattern(args[0], mapping, n), f"{args[0]}_n{n}"
        else:
            raise InputError(f"COPOLYMER must be {KEYS['COPOLYMER'].value}, got {' '.join(opts['copolymer'])!r}",
                             key="COPOLYMER")
    stem += "".join(f"_{tag}{opts[k]}" for tag, k in (("sp", "spacer"), ("alk", "alkyl")) if opts.get(k) is not None)
    return seq, stem


def _library(opts):
    return load_library(*opts["library"]) if opts.get("library") else default_library()


def _stem_with_initiator(opts, stem, lib):
    if opts.get("initiator") and opts.get("ends") != "none":
        stem += f"_{lib.initiator(opts['initiator']).name}"
    return _UNSAFE.sub("_", stem)


def _structure(opts, lib):
    seq, stem = sequence_and_stem(opts)
    units, _ = resolve_units(seq, None, opts.get("spacer"), opts.get("alkyl"), lib)
    structure = build_polymer(units, initiator=opts.get("initiator"), ends=opts.get("ends"), cap=_cap(opts),
                              conformation=opts.get("conformation") or "coil", seed=opts.get("seed", 0), library=lib)
    return structure, _stem_with_initiator(opts, stem, lib)


def _cap(opts) -> str:
    return "CH3" if str(opts.get("cap") or "ch3").lower() == "ch3" else "H"


def _times(opts) -> dict | None:
    items = opts.get("times")
    if items is None:
        return None
    if len(items) % 2:
        raise InputError(f"TIMES needs name/value pairs ({KEYS['TIMES'].value}), got {' '.join(items)!r}",
                         key="TIMES")
    out = {}
    for name, value in zip(items[::2], items[1::2], strict=True):
        try:
            out[name.lower()] = float(value)
        except ValueError:
            raise InputError(f"TIMES: {value!r} for {name} is not a number (ns)", key="TIMES") from None
    return out


def _salt(opts):
    items = opts["salt"]
    if isinstance(items[0], str) and items[0].lower() == "none":
        if len(items) != 1:
            raise InputError("SALT = none takes no molarity", key="SALT")
        return None, 0.15
    if len(items) > 2:
        raise InputError(f"SALT must be {KEYS['SALT'].value}, got {' '.join(map(str, items))!r}", key="SALT")
    name = {"nacl": "NaCl", "kcl": "KCl"}.get(str(items[0]).lower())
    if name is None:
        raise InputError(f"SALT must be NaCl, KCl or none, got {items[0]!r}", key="SALT")
    try:
        molar = float(items[1]) if len(items) == 2 else 0.15
    except ValueError:
        raise InputError(f"SALT molarity must be a number, got {items[1]!r}", key="SALT") from None
    return name, molar


def _box(opts):
    box = opts.get("box")
    if box is None:
        return None
    return box[0] if len(box) == 1 else box


def _configure_ff(opts) -> None:
    if opts.get("oplsaa_ff") or opts.get("opls_rules"):
        register_backend(OplsaaBackend(oplsaa_ff=opts.get("oplsaa_ff"), rules=opts.get("opls_rules")))


def _job_build(opts, lib):
    from .simple import BuildResult

    structure, stem = _structure(opts, lib)
    out = Path(opts.get("output") or stem)
    stem, folder = (out.name[:-4] if out.name.endswith(".pdb") else out.name), out.parent
    chain = structure.chain
    lines = []
    if opts["ff"].lower() == "none":
        pdb, top, ff_dir, topology, dft_blocks = write_pdb(structure, folder / f"{stem}.pdb"), None, None, None, ()
    else:
        topology = parametrize(structure, opts["ff"], name=stem, charges=opts.get("charges"),
                               params=opts.get("params"), dft_charges=opts.get("dft"),
                               dft_parts=opts["dft_parts"], charge_scale=opts["charge_scale"],
                               builtin_charges=opts["builtin_charges"] == "yes")
        if topology.dft_missing:
            print(f"note: no DFT charges for {', '.join(topology.dft_missing)}; using {topology.force_field} "
                  "charges (write the DFT molecules with JOB = dft)", file=sys.stderr)
        top = write_gromacs(topology, folder, stem)
        pdb, ff_dir, dft_blocks = folder / f"{stem}.pdb", folder / f"{stem}.ff", topology.dft_blocks
        lines.append(f"force field {topology.force_field}: {top}, {ff_dir}/")
        if dft_blocks:
            lines.append(f"DFT charges: {len(dft_blocks)} residue blocks from {opts['dft']}")
        if topology.builtin_blocks:
            lines.append(f"built-in RESP charges: {', '.join(topology.builtin_blocks)} "
                         "(BUILTIN_CHARGES = no to switch off)")
        if topology.explicit:
            lines.append(f"built-in bonded terms (sulfonate, GAFF2 for SO4-): {len(topology.explicit)}, "
                         "written on their lines in the .top and .rtp")
        if topology.ion_blocks:
            lines.append(f"ion fragments: {len(topology.ion_blocks)} residue blocks")
        if topology.charge_scale != 1:
            lines.append(f"ionic groups scaled by {topology.charge_scale:g}: scale the counterions/salt too "
                         "(JOB = bulk / slab does this)")
        if topology.user_parameters:
            lines.append(f"user parameters: {len(topology.user_parameters)}")
        for rep in topology.charge_report:
            lines.append(f"user charges: {rep.residue} {rep.user_atoms} atoms, correction "
                         f"{rep.correction_per_atom:+.4f} e/atom" + (" (large)" if rep.large else ""))
    counter = f", counterion {chain.counterion}" if chain.counterion else ""
    dmin = structure.min_nonbonded_distance
    contact = f"{dmin:.2f} Å" if dmin != float("inf") else "> 3.00 Å"
    print(f"{pdb}: {len(structure.sequence)} units, {structure.mol.GetNumAtoms()} atoms, "
          f"net charge {chain.net_charge:+d}{counter}, Rg {structure.radius_of_gyration:.1f} Å, "
          f"closest non-bonded contact {contact}")
    for line in lines:
        print(line)
    return BuildResult(pdb, top, ff_dir, structure, topology, chain.net_charge, chain.counterion, dft_blocks)


def _job_dft(opts, lib):
    folder = Path(opts.get("dft") or "dft")
    parts = opts["dft_parts"]
    if opts.get("molecule") is not None:
        manifest = export_dft_molecule(opts["molecule"], folder, parts=parts, library=lib)
    else:
        seq, _ = sequence_and_stem(opts)
        units, _ = resolve_units(seq, None, opts.get("spacer"), opts.get("alkyl"), lib)
        initiator = opts.get("initiator") if opts.get("ends") != "none" else None
        manifest = None
        for k, unit in enumerate(dict.fromkeys(units)):
            manifest = export_dft_inputs(unit, folder, initiator=initiator if k == 0 else None, library=lib,
                                         parts=parts)
    for mol in yaml.safe_load(manifest.read_text(encoding="utf-8"))["molecules"]:
        print(f"{folder / mol['files'][0]}: {len(mol['atoms'])} atoms, charge {mol['net_charge']}, "
              f"multiplicity {mol['multiplicity']}")
    use = f"DFT = {folder}" + ("  and  DFT_PARTS = ions" if parts == "ions" else "")
    print(f"Next: run DFT (e.g. B3LYP or HF, ESP/RESP) on each .xyz, keep the atom order, and save the charges\n"
          f"in {folder} as <name>.charges, .mol2, .log (Gaussian) or .out (ORCA). Then add  {use}\n"
          f"to your build, bulk or slab input.")
    return manifest


def _job_system(opts, lib):
    seq, stem = sequence_and_stem(opts)
    stem = _stem_with_initiator(opts, stem, lib)
    job = opts["job"]
    common = dict(chains=opts["chains"], box_nm=_box(opts), initiator=opts.get("initiator"), ends=opts.get("ends"),
                  spacer=opts.get("spacer"), alkyl=opts.get("alkyl"), dft=opts.get("dft"), ff=opts["ff"],
                  dft_parts=opts["dft_parts"], charge_scale=opts["charge_scale"], cap=_cap(opts),
                  solvent=opts["solvent"],
                  builtin_charges=opts["builtin_charges"] == "yes",
                  packer=opts["packer"], times_ns=_times(opts), seed=opts["seed"], solvent_ff=opts.get("solvent_ff"),
                  out=opts.get("output") or f"{stem}_{job}", library=lib)
    if job == "bulk":
        t_high, t_low, t_step = opts["tg"]
        if opts.get("density") is not None:
            common.update(box_nm=None, density=opts["density"])
        folder = bulk_system(seq, None, t_high=t_high, t_low=t_low, t_step=t_step, **common)
    else:
        salt, molar = _salt(opts)
        if opts.get("density") is not None:
            common.update(box_nm=None, density=opts["density"], aspect=opts["aspect"])
        folder = slab_system(seq, None, slab_nm=opts.get("slab"), salt=salt, salt_molar=molar, t_equil=opts["tequil"],
                             temperatures=opts["temperatures"], **common)
    print(f"{job} system written to {folder}/\nNext: bash {folder}/run.sh   (or sbatch {folder}/run_slurm.sh), "
          f"then run an input with  JOB = analyze  and  OUTPUT = {folder}")
    return folder


def _job_analyze(opts, lib):
    folder = Path(opts["output"])
    try:
        kind = yaml.safe_load((folder / "protocol.yaml").read_text(encoding="utf-8")).get("kind")
    except (OSError, AttributeError, yaml.YAMLError):
        raise InputError(f"{folder} has no readable protocol.yaml (was it written by JOB = bulk or slab?)",
                         key="OUTPUT") from None
    if kind == "bulk":
        result = analyze_tg(folder)
        density = f"{result.density_300:.3f} g/cm3" if result.density_300 is not None else "n/a"
        print(f"Tg = {result.tg:.1f} K, density at 300 K = {density} -> {folder / 'tg_result.yaml'}")
    elif kind == "slab":
        result = analyze_slab(folder)
        for r in result:
            dilute = f"{r.dilute:.1f}" if r.dilute is not None else "n/a"
            print(f"{r.temperature:g} K: dense {r.dense:.1f} g/L, dilute {dilute} g/L, width {r.width_nm:.2f} nm")
        print(f"-> {folder / 'phase_diagram.csv'}")
    else:
        raise InputError(f"{folder}/protocol.yaml has unknown kind {kind!r}", key="OUTPUT")
    return result


def _job_template(opts, lib):
    structure, stem = _structure(opts, lib)
    stem = opts.get("output") or stem
    ff = opts["ff"].lower()
    if opts["job"] == "charges":
        text = charge_template(structure, parametrize(structure, ff) if ff != "none" else None)
    else:
        text = ionic_parameter_template(parametrize(structure, ff, strict=False))
    path = Path(f"{stem}_{opts['job']}.yaml")
    path.write_text(text, encoding="utf-8")
    use = "CHARGES" if opts["job"] == "charges" else "PARAMS"
    print(f"wrote {path}; edit it, then add  {use} = {path}  to your build input")
    return path


_RUNNERS = {"build": _job_build, "dft": _job_dft, "bulk": _job_system, "slab": _job_system,
            "analyze": _job_analyze, "charges": _job_template, "params": _job_template}


# rq-130422bb rq-df9384f8
def run_input(path: str | os.PathLike):
    """Parse an input file and run its job; returns the job's result"""
    source = Path(path)
    try:
        text = source.read_text(encoding="utf-8")
    except OSError as exc:
        raise InputError(f"cannot read input file {source}: {exc.strerror or exc}") from None
    opts = parse_input(text, str(source))
    lib = _library(opts)
    if opts["job"] in WITH_FF:
        _configure_ff(opts)
    return _RUNNERS[opts["job"]](opts, lib)
