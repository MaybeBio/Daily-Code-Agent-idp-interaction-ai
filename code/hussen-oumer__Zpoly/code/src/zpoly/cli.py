"""`zpoly`"""

# rq-80888a1d rq-5029ce48

from __future__ import annotations

import sys
from collections.abc import Sequence
from pathlib import Path

from . import __version__
from .forcefield import ForceFieldError, available_backends, get_backend
from .inputfile import JOBS, KEYS, InputError, run_input, suggest
from .structure import StructureError
from .zwitterions import ZwitterionError, default_library

USAGE = """\
Zpoly {version}: zwitterionic polymers -> PDB, GROMACS force field, MD

  zpoly list [monomers|initiators|forcefields]   names you can use
  zpoly help [KEY]                                keywords, or one keyword
  zpoly template [{jobs}]   a start file
  zpoly INPUT                                     run an input file

Input file: one KEY = value per line, # starts a comment."""

TEMPLATES = {
    "build": """\
# One chain + its GROMACS force field.  Run: zpoly build.inp
JOB          = build
MONOMER      = DMAPS        # see: zpoly list
N            = 12           # units per chain
INITIATOR    = KPS          # SO4- end group; see: zpoly list initiators
ENDS         = head         # head | tail | both | none (ends without KPS are CH3)
FF           = oplsaa       # oplsaa | gaff2 | none (PDB only)
CHARGE_SCALE = 0.85         # ionic groups x 0.85 (1 = unscaled)
# SPACER     = 4            # CH2 groups between cation and anion
# OUTPUT     = dmaps12      # default: DMAPS_n12_KPS
""",
    "dft": """\
# Molecules for your own DFT charges: whole monomer + cation, anion and ester fragments.
JOB          = dft
MONOMER      = DMAPS
N            = 12
INITIATOR    = KPS
DFT          = dft          # folder to write the .xyz/.pdb files to
DFT_PARTS    = ions         # whole | ions
""",
    "bulk": """\
# Dry bulk system + GROMACS protocol for density and Tg.
JOB          = bulk
MONOMER      = DMAPS
N            = 12
INITIATOR    = KPS
ENDS         = head
FF           = oplsaa
CHARGE_SCALE = 0.85
SOLVENT      = tip4p2005    # Madrid-2019 K+ counterions (tip3p: those of the force field)
CHAINS       = 50           # you choose
DENSITY      = 1.0          # g/cm3 start density: the cube edge follows (or BOX = edge | Lx Ly Lz, nm)
TG           = 600 200 20   # cooling T_high T_low step (K)
PACKER       = emc          # emc (Monte Carlo) | builtin | packmol
""",
    "slab": """\
# LLPS slab (Lx = Ly, Lz >= 3 Lx) in water + salt, run at several temperatures.
JOB          = slab
MONOMER      = DMAPS
N            = 12
INITIATOR    = KPS
ENDS         = head
FF           = oplsaa
CHARGE_SCALE = 0.85
SOLVENT      = tip4p2005    # TIP4P/2005 + Madrid-2019 ions | tip3p
CHAINS       = 50           # you choose
DENSITY      = 1.0          # g/cm3 of the starting dense phase: cube Lx = Ly = slab (or BOX + SLAB)
ASPECT       = 3            # Lz = ASPECT x Lx (at least 3)
SALT         = NaCl 0.15    # NaCl | KCl [molar] | none
TEMPERATURES = 280 300 320 340 360
PACKER       = emc          # emc (Monte Carlo) | builtin | packmol
""",
    "analyze": """\
# Analyse a finished bulk (Tg, density) or slab (phase diagram) folder.
JOB          = analyze
OUTPUT       = DMAPS_n12_KPS_slab
""",
}

# help: keyword sections
SECTIONS = (
    ("Job", ("JOB", "OUTPUT")),
    ("Polymer", ("MONOMER", "N", "SEQUENCE", "COPOLYMER", "INITIATOR", "ENDS", "CAP", "SPACER", "ALKYL")),
    ("Force field and charges", ("FF", "CHARGE_SCALE", "BUILTIN_CHARGES", "DFT", "DFT_PARTS", "CHARGES", "PARAMS")),
    ("System (bulk, slab)", ("CHAINS", "DENSITY", "ASPECT", "BOX", "SLAB", "SALT", "SOLVENT", "TEMPERATURES",
                             "TG", "PACKER")),
    ("Advanced", ("SEED", "CONFORMATION", "TIMES", "TEQUIL", "SOLVENT_FF", "OPLSAA_FF", "OPLS_RULES", "MOLECULE",
                  "LIBRARY")),
)


def _help() -> str:
    lines = [USAGE.format(version=__version__, jobs="|".join(TEMPLATES))]
    for title, names in SECTIONS:
        lines.append(f"\n {title}")
        lines += [f"  {n:<16s}{KEYS[n].short}" for n in names]
    lines.append("\nJOB = " + " | ".join(JOBS) + "\nMore about one keyword: zpoly help KEY")
    return "\n".join(lines)


def _key_help(name: str) -> str:
    key = KEYS[name.upper()]
    default = key.default
    if isinstance(default, tuple):
        default = " ".join(f"{v:g}" if isinstance(v, float) else str(v) for v in default)
    return (f"{key.name}: {key.short}\n"
            f"  value:   {key.value}\n"
            f"  default: {'none' if default is None else default}\n"
            f"  jobs:    {', '.join(key.jobs)}\n"
            f"  example: {key.example}")


def _list(what: str) -> None:
    lib = default_library()
    if what == "monomers":
        for z in lib.entries(polymerizable=True):
            print(f"{z.name:<14}{', '.join(z.aliases) or '-':<28}{z.family}")
    elif what == "initiators":
        for i in lib.initiators():
            print(f"{i.name:<8}{', '.join(i.aliases) or '-':<36}counterion {i.counterion or '-'}")
    else:
        for name in available_backends():
            try:
                get_backend(name).check_available()
                status = "available"
            except ForceFieldError as exc:
                status = f"not installed: {getattr(exc, 'hint', exc)}"
            print(f"{name:<10}{status}")


def _error(message: str) -> int:
    print(f"error: {message}", file=sys.stderr)
    return 1


def main(argv: Sequence[str] | None = None) -> int:
    args = list(sys.argv[1:] if argv is None else argv)
    command = args[0].lower() if args else "help"
    if command in ("help", "-h", "--help"):
        if len(args) > 2:
            return _usage("help takes at most one keyword")
        if len(args) == 2:
            if args[1].upper() not in KEYS:
                return _error(f"unknown keyword {args[1].upper()}{suggest(args[1])}")
            print(_key_help(args[1]))
        else:
            print(_help())
        return 0
    if command in ("-v", "--version"):
        print(__version__)
        return 0
    if command == "list":
        what = args[1].lower() if len(args) > 1 else "monomers"
        if len(args) > 2 or what not in ("monomers", "initiators", "forcefields"):
            return _usage("list monomers | initiators | forcefields")
        _list(what)
        return 0
    if command == "template":
        job = args[1].lower() if len(args) > 1 else "build"
        if len(args) > 2 or job not in TEMPLATES:
            return _usage(f"template {' | '.join(TEMPLATES)}")
        print(TEMPLATES[job], end="")
        return 0
    if len(args) != 1:
        return _usage("give one input file, e.g. zpoly my_polymer.inp")
    if not Path(args[0]).is_file():
        return _error(f"no input file {args[0]!r} (commands: list, help, template; or an input file)")
    try:
        run_input(args[0])
    except (InputError, ZwitterionError, StructureError, ForceFieldError, ValueError, OSError) as exc:
        return _error(str(exc))
    return 0


def _usage(message: str) -> int:
    print(f"usage: zpoly {message}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(main())
