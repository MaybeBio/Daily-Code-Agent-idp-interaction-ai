# Feature: Simple Python API (Polymer + GROMACS Force Field in a Few Lines) <!-- rq-b39965b3 -->

The everyday workflow in as few lines as possible. Users give the monomer, the chain length, the
end groups and their DFT charges, and get the PDB plus a ready GROMACS force field:

```python
import zpoly as pp

pp.dft_inputs("DMAPS", initiator="KPS", folder="dft")                 # 1. molecules for your DFT
# run your DFT (e.g. B3LYP or HF, ESP/RESP) on dft/DMAPS.xyz and dft/KPS_end.xyz, save the charges there
pp.build("DMAPS", n=12, initiator="KPS", dft="dft", out="dmaps12")   # 2. dmaps12.pdb + .top + .ff/
```

These two functions wrap the existing layers (library, chain builder, force-field core, DFT
charges). Everything they do is also available through the detailed APIs and the CLI.

## Feature API <!-- rq-735d2b77 -->

Importable directly from `zpoly`.

### `dft_inputs(monomer, *, initiator=None, folder="dft", spacer=None, alkyl=None, library=None) -> list[Path]` <!-- rq-9cc74219 -->

- Writes the DFT molecules exactly like `export_dft_inputs` (`rqm/forcefield/dft_charges.md`):
  the methyl-capped monomer, and the initiator end group(s) when `initiator` is given.
- `spacer` / `alkyl` apply `make_variant` first. <!-- rq-77f7a7bc -->
- `parts="ions"` also writes the ion fragments (`rqm/forcefield/dft_charges.md`). <!-- rq-0a935855 -->
- Returns the `.xyz` paths written, and prints one line per molecule (file, atoms, net charge, multiplicity).

### `build(monomer, n=None, *, initiator=None, ends=None, spacer=None, alkyl=None, dft=None, ff="gaff2", out="polymer", seed=0, library=None) -> BuildResult` <!-- rq-82b577fa -->

- `monomer` is a monomer name (homopolymer of `n` units), or a list of monomer names / `Zwitterion` <!-- rq-9876a6b6 -->
  objects (an explicit sequence; `n` must then be `None` or equal to its length). For alternating,
  block and random copolymers, pass the list from the sequence helpers (`alternating`, `block`,
  `random_copolymer`, `pattern`; also importable from `zpoly`).
- `spacer` / `alkyl` apply `make_variant` to every monomer, as the CLI options do. <!-- rq-4568c7e5 -->
- `dft` is one of: <!-- rq-4d7e893b -->
  - `None`: the force field's own charges (AM1-BCC for GAFF2), with the usual notice;
  - a folder written by `dft_inputs` / `dft-export` that contains the charge files;
  - a mapping from names to plain charge lists in the exported atom order. Keys are monomer names
    and the initiator name (`"KPS"`), or the exported molecule names (`"KPS_end"`, `"KPS_head"`,
    `"KPS_tail"`). Each list is checked like a charge file: the count must equal the exported
    molecule's atom count, and the sum must equal its net charge within 0.01 e (`DftChargeError`).
    A name that does not belong to the chain raises `DftChargeError`.
- `ff` is `"gaff2"` (default), `"oplsaa"` or `"none"` (PDB only). <!-- rq-facb859d -->
- `dft_parts="whole" | "ions"` and `charge_scale=1.0` as in `rqm/forcefield/dft_charges.md`; <!-- rq-f2870deb -->
  `bulk_system` / `slab_system` take the same two arguments.
- Writes `<out>.pdb`, and with a force field also `<out>.top` and `<out>.ff/`, into the folder of
  `out` (created if needed). Nothing is written if building or parametrizing fails.
- Returns `BuildResult` (frozen dataclass): `pdb: Path`, `top: Path | None`, `ff_dir: Path | None`,
  `structure: PolymerStructure`, `topology: Topology | None`, `net_charge: int`,
  `counterion: str | None`, `dft_blocks: tuple[str, ...]`.
- Prints one summary line, like the CLI.

## Engine Extensions (Planned) <!-- rq-567ac589 -->

GROMACS is the engine of this feature. **OpenMM** and **LAMMPS** writers are planned as extensions
(engine-writer plugins on the `Topology`, see `rqm/ARCHITECTURE.md`) and are out of scope here.

## Gherkin Scenarios <!-- rq-907bff89 -->

```gherkin
Feature: Simple Python API

  @rq-a3b3b51e
  Scenario: Top-level imports
    When "import zpoly as pp" is executed
    Then pp.build, pp.dft_inputs, pp.alternating, pp.block, pp.random_copolymer and pp.pattern exist

  @rq-621a8723
  Scenario: Export DFT molecules
    When pp.dft_inputs("DMAPS", initiator="KPS", folder=tmp/"dft") is called
    Then it returns the paths of DMAPS.xyz and KPS_end.xyz, which exist

  @rq-bd003416
  Scenario: Build with a DFT folder
    Given the fake backend and a dft folder with charge files for DMAPS and KPS_end
    When pp.build("DMAPS", n=4, initiator="KPS", dft=folder, ff="fake", out=tmp/"o/d4") is called
    Then tmp/o contains d4.pdb, d4.top and d4.ff/
    And the result's dft_blocks are KPH, DMF, DMR, DME, KPT and its net_charge is -2

  @rq-a7d55590
  Scenario: Build with pasted charge values
    Given the fake backend and charge lists for "DMAPS" and "KPS" in the exported atom order
    When pp.build("DMAPS", n=3, initiator="KPS", dft={"DMAPS": [...], "KPS": [...]}, ff="fake", out=tmp/"d3") is called
    Then the result's dft_blocks include DMR and KPH

  @rq-1fe6b25c
  Scenario: Pasted charges with the wrong count
    When pp.build("DMAPS", n=3, dft={"DMAPS": [0.0, 0.0]}, ff="fake", out=tmp/"x") is called
    Then DftChargeError is raised and no file is written

  @rq-da585695
  Scenario: Pasted charges for a monomer not in the chain
    When pp.build("DMAPS", n=3, dft={"CBMA": [...]}, ff="fake", out=tmp/"x") is called
    Then DftChargeError is raised naming "CBMA"

  @rq-78cb44ef
  Scenario: Copolymer from a sequence helper
    Given the fake backend
    When pp.build(pp.alternating("DMAPS", "CBMA", 4), ff="fake", out=tmp/"alt") is called
    Then the PDB residues are DMF, CBR, DMR, CBE

  @rq-7896b345
  Scenario: PDB only
    When pp.build("DMAPS", n=3, ff="none", out=tmp/"p") is called
    Then only p.pdb is written and the result's top is None

  @rq-ca771be1
  Scenario: Length mismatch
    When pp.build(["DMAPS", "CBMA"], n=5, ff="none", out=tmp/"m") is called
    Then ValueError is raised
```
