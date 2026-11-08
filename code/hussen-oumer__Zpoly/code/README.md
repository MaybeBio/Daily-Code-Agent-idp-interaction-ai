<p align="center">
  <img src="https://github.com/hussen-oumer/Zpoly/raw/HEAD/docs/assets/Zpoly_horizontal_white_bg.png" alt="Zpoly — Zwitterionic Polymer Packer" width="620">
</p>

<p align="center">
  <b>Zwitterionic polymers, from a monomer name to a ready-to-run GROMACS simulation.</b>
</p>

<p align="center">
  <img src="https://img.shields.io/badge/python-3.10%2B-blue" alt="Python 3.10+">
  <img src="https://img.shields.io/badge/GROMACS-any%20version-green" alt="GROMACS">
  <img src="https://img.shields.io/badge/force%20fields-OPLS--AA%20%7C%20GAFF2-orange" alt="OPLS-AA | GAFF2">
  <img src="https://img.shields.io/badge/license-MIT-lightgrey" alt="MIT">
</p>

---

Zpoly builds all-atom models of **zwitterionic polymers** for molecular-dynamics studies of
**liquid–liquid phase separation (LLPS)**. You write one short input file 
and Zpoly gives you the polymer, its GROMACS force field, and a complete simulation box with its
protocol: a dry bulk system for density and Tg, or a direct-coexistence slab in water and salt for
the phase diagram.

## Features

- **Monomer library** — 61 polymerizable zwitterions and comonomers (sulfobetaines,
  carboxybetaines, phosphorylcholines, imidazolium and pyridinium betaines, DABCO, NIPAM, …) and
  9 initiators (KPS, APS, AIBN, …), by name; your own monomers in a small YAML file.
- **Chain builder** — homopolymers, alternating, block, random and pattern copolymers; KPS (SO4⁻)
  ends on one or both sides, CH3 elsewhere; spacer and N-alkyl variants; fast unit-by-unit growth
  for long chains.
- **Ready force field** — OPLS-AA (stock atom types) or GAFF2, written as a `.top` that works with
  any GROMACS version: residues named F/R/E, every special term with its values on its line.
- **Trusted charges** — built-in RESP charges for DMAPS, A3361, A3367 (and their spacer variants)
  and KPS; or your own DFT charges (whole monomer, or cation / anion / ester fragments);
  `CHARGE_SCALE` (e.g. 0.85) for the ionic groups, counterions and salt.
- **Water and ions** — TIP4P/2005 with the Madrid-2019 ions (their ion–ion and ion–water
  parameters built in), or TIP3P with the force field's ions.
- **Systems** — dry bulk cube (density and Tg by stepwise cooling) and LLPS slab
  (Lx = Ly, Lz ≥ 3 Lx, temperature series), from the number of chains and a density. An
  EMC-style Monte Carlo packer, then GROMACS compresses the box to the exact density.
- **Analysis** — Tg (bilinear fit) and density at 300 K; dense and dilute concentrations per
  temperature (phase diagram).
- **HPC ready** — prepare on your computer, run `run_slurm.sh` on the cluster, analyse anywhere.

## How it works

```
  monomer name ──► chain (PDB) ──► force field (.top, .ff) ──► box + protocol ──► GROMACS ──► analysis
  JOB = build                                    JOB = bulk | slab    run.sh       JOB = analyze
```

| `JOB =` | You give | You get |
|---|---|---|
| `build` | the polymer | `.pdb` and a ready GROMACS `.top` (+ `.ff/`) |
| `dft` | the polymer | molecules for your own DFT charges (`.xyz`) |
| `bulk` | `CHAINS`, `DENSITY` | dry bulk cube + density/Tg protocol (`run.sh`) |
| `slab` | `CHAINS`, `DENSITY`, `ASPECT` | LLPS slab in water + salt + temperature series |
| `analyze` | `OUTPUT = <folder>` | Tg and density, or the phase diagram |

## Installation

```
git clone <this repository> zpoly && cd zpoly
python -m venv ~/envs/zpoly && source ~/envs/zpoly/bin/activate
pip install -e .
zpoly help
```

You also need GROMACS (`gmx`). Installing on an HPC cluster in your own space, and moving
simulations between your computer and the cluster: **[INSTALLATION.md](https://github.com/hussen-oumer/Zpoly/blob/HEAD/INSTALLATION.md)**.

## Quick start

```
zpoly list                          # monomer names   (zpoly list initiators)
zpoly help                          # every keyword   (zpoly help DENSITY: one in detail)
zpoly template slab > slab.inp      # a ready input file: edit it
zpoly slab.inp                      # -> DMAPS_n12_KPS_slab/
bash DMAPS_n12_KPS_slab/run.sh      # or: sbatch DMAPS_n12_KPS_slab/run_slurm.sh
```

`slab.inp`:

```
JOB          = slab
MONOMER      = DMAPS        # zpoly list
N            = 12           # units per chain
INITIATOR    = KPS
ENDS         = head         # KPS on one end, CH3 on the other
FF           = oplsaa
CHARGE_SCALE = 0.85         # built-in RESP charges x 0.85
SOLVENT      = tip4p2005    # TIP4P/2005 + Madrid-2019 ions
CHAINS       = 50
DENSITY      = 1.0          # g/cm3 of the starting dense phase: the box follows
ASPECT       = 3            # Lz = 3 Lx
SALT         = NaCl 0.15
TEMPERATURES = 280 300 320 340 360
PACKER       = emc
```

Then the analysis:

```
JOB    = analyze
OUTPUT = DMAPS_n12_KPS_slab        # -> phase_diagram.csv, density_profile_<T>K.csv
```

Ready input files for every job: [`examples/inputs/`](https://github.com/hussen-oumer/Zpoly/blob/HEAD/examples/inputs).

### Common changes

| You want | Write |
|---|---|
| another monomer | `MONOMER = A3361` (see `zpoly list`) |
| one more CH2 between cation and anion | `SPACER = 4` |
| a copolymer | `COPOLYMER = random DMAPS:0.5 A3367:0.5` (or `alternating`, `block`, `pattern`) |
| KPS on both ends / none | `ENDS = both` / `ENDS = none` |
| unscaled charges | `CHARGE_SCALE = 1` |
| TIP3P water | `SOLVENT = tip3p` |
| GAFF2 | `FF = gaff2` (needs AmberTools) |
| your own box | `BOX = 8 8 24` and `SLAB = 6` instead of `DENSITY` |

### Python

```python
import zpoly as zp

zp.build("DMAPS", n=12, initiator="KPS", ends="head", ff="oplsaa", charge_scale=0.85, out="dmaps12")
zp.slab_system("DMAPS", n=12, initiator="KPS", ends="head", ff="oplsaa", charge_scale=0.85,
               solvent="tip4p2005", chains=50, density=1.0, packer="emc", out="dmaps_slab")
zp.run_input("slab.inp")
```

## Adding a New Monomer

**Without changing the software:** write the monomer in a small file and point the input to it.
Example: [`examples/inputs/my_monomers.yaml`](https://github.com/hussen-oumer/Zpoly/blob/HEAD/examples/inputs/my_monomers.yaml) (used by
`examples/inputs/new_monomer.inp`).

```yaml
version: 1
monomers:
  - name: SBMA-C3L                        # MONOMER = SBMA-C3L
    full_name: '[3-(Methacryloyloxy)propyl]dimethyl(3-sulfopropyl)ammonium'
    family: Sulfobetaine
    status: specialty                     # common | specialty | combinatorial | computed
    residue: SL                           # optional 2-character residue prefix
    blocks: [MA, C3, NMe2, C3, SO3]       # or: smiles: C=C(C)C(=O)OCCC[N+](C)(C)CCCS(=O)(=O)[O-]
```

```
LIBRARY = my_monomers.yaml       # in the input file
MONOMER = SBMA-C3L
```

- `smiles`: the monomer with its polymerizable C=C (vinyl, acrylate, methacrylate, acrylamide,
  styrene, …), net charge 0 with a + and a − centre. Non-zwitterionic comonomers go under
  `comonomers:` (with `counterion: Cl-` etc. when charged).
- `blocks`: the same monomer from building blocks (polymerizable group, linker, cation, spacer,
  anion; listed in [`docs/zwitterion_catalog.md`](https://github.com/hussen-oumer/Zpoly/blob/HEAD/docs/zwitterion_catalog.md)). Only this form
  allows `SPACER` / `ALKYL` variants.
- Charges: compute them with `JOB = dft` (`DFT_PARTS = ions` adds the cation, anion and ester
  fragments) and use `DFT = dft`; or give them in a `CHARGES` file.

**In the software (for everyone):**

| What | File |
|---|---|
| the monomer | `src/zpoly/zwitterions/data/zwitterion_library.yaml`, section `monomers:` |
| its RESP charges (optional) | `src/zpoly/forcefield/data/builtin_charges.yaml`, section `monomers:` |
| missing OPLS-AA atom types (optional) | `src/zpoly/forcefield/data/oplsaa_rules.yaml` (or `OPLS_RULES = file`) |

Then regenerate the catalog (`python devtools/build_zwitterion_catalog.py > docs/zwitterion_catalog.md`)
and run the tests (`pytest`), which check every built-in monomer.

## Project layout

```
src/zpoly/
  zwitterions/   monomer library (data/zwitterion_library.yaml), variants, chains
  structure/     3D chain builder, PDB output
  forcefield/    OPLS-AA and GAFF2 backends, DFT and built-in charges, GROMACS writers
  system/        bulk and slab systems, EMC packer, water/ions (Madrid-2019), protocols, analysis
  inputfile.py   the input file;  cli.py: the zpoly command;  simple.py: short Python API
examples/inputs/ ready input files          rqm/  specifications of every feature
tests/           test suite (pytest)        devtools/  catalog, GAFF2 extraction, topology comparison
```

## References

The data and methods Zpoly uses (monomer set, Madrid-2019, TIP4P/2005, EMC and packing methods) are
listed in [REFERENCES.md](https://github.com/hussen-oumer/Zpoly/blob/HEAD/REFERENCES.md). Please cite them together with Zpoly.

## Acknowledgements

🌊 **deSalSea** · 🇪🇺 **M-ERA.NET 2024**

This project is developed for deSalSea with funding from the M-ERA.NET Joint Call 2024.

🤖 Zpoly was developed with Claude (Anthropic) through Claude Code, using the
[Riprap](https://github.com/MolSSI/riprap) template (MolSSI) in a rootless Podman container. Its
requirements-first workflow keeps a specification of every feature in [`rqm/`](https://github.com/hussen-oumer/Zpoly/blob/HEAD/rqm), and every
specification scenario is checked by the test suite.

## License

MIT — see [LICENSE](https://github.com/hussen-oumer/Zpoly/blob/HEAD/LICENSE). Developing Zpoly itself: [DEVELOPMENT.md](https://github.com/hussen-oumer/Zpoly/blob/HEAD/DEVELOPMENT.md).

Hussen O. Mohammed
