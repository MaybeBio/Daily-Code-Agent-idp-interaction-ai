# Feature: One Input File (VASP/Gaussian-Style Interface) <!-- rq-8ecf5b05 -->

The command-line program has one way of working, like VASP, Gaussian or Schrödinger: the user writes
a short input file of `KEY = value` lines and runs `zpoly file.inp`. Three helper commands
teach it:

```
zpoly list [monomers|initiators|forcefields]   # names to use in the input file
zpoly help [KEY]                               # keyword table, or one keyword in detail
zpoly template [build|dft|bulk|slab|analyze]   # a commented input file to start from
zpoly INPUT                                    # run an input file
```

Everything the input file does is also available from Python (`rqm/simple_api.md`).

## Input File Format <!-- rq-b18197c0 -->

- One `KEY = value` per line. Keys are case-insensitive. `#` and `!` start comments. Blank lines are
  ignored. Values are separated by spaces (or commas), e.g. `BOX = 8 8 24`.
- Unknown key → error naming the line and the closest known keys. A key given twice → error.
- Values are checked against their type (integer, number, list, choice); errors name the key, the line
  and the expected form.
- Keys that do not apply to the chosen `JOB` → error (for example `SLAB` with `JOB = bulk`).

## Keywords <!-- rq-f2d8deab -->

| Key | Jobs | Value | Default |
|---|---|---|---|
| `JOB` | all | `build`, `dft`, `bulk`, `slab`, `analyze`, `charges`, `params` | `build` |
| `MONOMER` | polymer jobs | a monomer name (homopolymer of `N` units) | — |
| `SEQUENCE` | polymer jobs | monomer names in chain order | — |
| `COPOLYMER` | polymer jobs | `alternating A B`, `block A:10 B:20`, `random A:0.5 B:0.5`, `pattern AAB A=X B=Y` | — |
| `N` | polymer jobs | units per chain (homopolymer, alternating, random, pattern) | — |
| `INITIATOR` | polymer jobs | initiator name (e.g. `KPS`) | none |
| `ENDS` | polymer jobs | `both`, `head`, `tail`, `none` | `both` with an initiator |
| `CAP` | polymer jobs | `CH3` or `H`: the chain ends without an initiator end group | `CH3` |
| `BUILTIN_CHARGES` | `build`, `bulk`, `slab` | `yes` / `no`: use the built-in RESP charges (`rqm/forcefield/builtin_charges.md`) | `yes` |
| `SPACER` | polymer jobs | CH2 groups between cation and anion | monomer's own |
| `ALKYL` | polymer jobs | N-alkyl length (1 methyl, 2 ethyl, …) | monomer's own |
| `CONFORMATION` | `build` | `coil`, `extended` | `coil` |
| `SEED` | polymer jobs | integer | 0 |
| `FF` | `build`, `bulk`, `slab`, `charges`, `params` | `gaff2`, `oplsaa`, `none` | `gaff2` |
| `DFT` | `build`, `bulk`, `slab`, `dft` | folder with DFT charge files (for `dft`: where to write the molecules) | `dft` for `JOB = dft` |
| `DFT_PARTS` | `build`, `dft`, `bulk`, `slab` | `whole`, or `ions` (+ cation/anion fragments, `rqm/forcefield/dft_charges.md`) | `whole` |
| `CHARGE_SCALE` | `build`, `bulk`, `slab` | number in (0, 1]: scales ionic groups, counterions and salt | 1 |
| `MOLECULE` | `dft` | a zwitterionic liquid (small molecule) instead of a polymer | — |
| `CHARGES` | `build` | charge file (per-atom overrides, from `JOB = charges`) | none |
| `PARAMS` | `build` | parameter file (bonded overrides, from `JOB = params`) | none |
| `OPLSAA_FF` | with `FF = oplsaa` | your `oplsaa.ff` folder | found via GROMACS |
| `OPLS_RULES` | with `FF = oplsaa` | extra typing rules | none |
| `CHAINS` | `bulk`, `slab` | number of chains (required) | — |
| `BOX` | `bulk`, `slab` | nm: one edge (cube) or `Lx Ly Lz` (required) | — |
| `SLAB` | `slab` | starting slab thickness in nm (required) | — |
| `DENSITY` | `slab` | starting dense-phase density, g/cm³: the box follows (cube Lx = Ly = slab, Lz = ASPECT × Lx); instead of `BOX` and `SLAB` | — |
| `ASPECT` | `slab` | Lz / Lx with `DENSITY` (at least 3) | 3 |
| `SALT` | `slab` | `NaCl 0.15`, `KCl 0.15` (molar) or `none` | `NaCl 0.15` |
| `TEMPERATURES` | `slab` | production temperatures in K | `280 300 320 340 360` |
| `TEQUIL` | `slab` | equilibration temperature in K | 300 |
| `TG` | `bulk` | cooling `T_high T_low step` in K | `600 200 20` |
| `TIMES` | `bulk`, `slab` | step times in ns, e.g. `cool 2 anneal 5` or `production 200` | protocol defaults |
| `SOLVENT` | `bulk`, `slab` | `tip3p` (water/ions of the GROMACS force field) or `tip4p2005` (TIP4P/2005 + Madrid-2019 ions, needs `CHARGE_SCALE = 0.85`) | `tip3p` |
| `PACKER` | `bulk`, `slab` | `packmol`, `builtin` | `packmol` |
| `SOLVENT_FF` | `bulk`, `slab` | GROMACS `.ff` folder for water and ions | matching one |
| `OUTPUT` | all but `dft` | output name (`build`, `charges`, `params`) or folder (`bulk`, `slab`); for `analyze`, the folder to analyse (required) | from the polymer |
| `LIBRARY` | all | extra monomer/initiator YAML files | none |

- Polymer jobs (`build`, `dft`, `bulk`, `slab`, `charges`, `params`) need exactly one of `MONOMER`,
  `SEQUENCE`, `COPOLYMER`.
  `JOB = dft` takes `MOLECULE` instead (and then no polymer keys).
- `N` is required with `MONOMER` and with the `alternating`, `random` and `pattern` forms, and is an
  error with `SEQUENCE` and `block` (their length is given).
- Defaults that describe the physical system (`CHAINS`, `BOX`, `SLAB`) do not exist: the user gives
  them.
- `FF = none` is allowed for `build` (PDB only, then `DFT`, `CHARGES`, `PARAMS` are errors) and
  `charges` (template without force-field charges); `bulk`, `slab` and `params` need a force field.
  `OPLSAA_FF` / `OPLS_RULES` need `FF = oplsaa`.

## Jobs <!-- rq-df9384f8 -->

| `JOB` | Does | Writes |
|---|---|---|
| `build` | one chain + force field | `<OUTPUT>.pdb`, `<OUTPUT>.top`, `<OUTPUT>.ff/` (only the PDB with `FF = none`) |
| `dft` | DFT molecules of every monomer of the polymer (+ initiator end groups) | `.xyz`/`.pdb` + `dft_manifest.yaml` in the `DFT` folder |
| `bulk` | dry bulk system + density/Tg protocol | the folder `OUTPUT` (`rqm/system/bulk.md`) |
| `slab` | LLPS slab + temperature-series protocol | the folder `OUTPUT` (`rqm/system/slab.md`) |
| `analyze` | analyses a finished bulk or slab folder (`OUTPUT`), detected from its `protocol.yaml` | `tg_result.yaml` / `phase_diagram.csv` |
| `charges` | editable charge-file template for the polymer's residues | `<OUTPUT>_charges.yaml` |
| `params` | editable parameter-file template for the ionic parts | `<OUTPUT>_params.yaml` |

- The default `OUTPUT` is derived from the polymer: `DMAPS_n12_KPS` for `MONOMER = DMAPS`, `N = 12`,
  `INITIATOR = KPS`; `DMAPS-alt-CBMA_n6`; `DMAPS10-CBMA20`; `DMAPS-ran-CBMA_n30_s1`; `AAB_n7`;
  `seq_n5`. `SPACER`/`ALKYL` add `_sp<n>`/`_alk<n>`. For `bulk`/`slab` the folder gets a `_bulk`/
  `_slab` suffix.
- Each job prints a short summary and exits 0. Any error prints `error: <message>` and exits 1 (2 for a
  wrong command line). Nothing is written when building or parametrizing fails.

## Help, List and Template <!-- rq-5029ce48 -->

- `zpoly` with no arguments, or `help`, prints the program line, the four commands, and the
  keywords in five sections (Job; Polymer; Force field and charges; System (bulk, slab); Advanced),
  one line each with its values or a short description, every keyword exactly once, within 80 columns.
- `help KEY` prints the key's value form, default, the jobs it applies to, and one example line.
  An unknown key → error with suggestions.
- `list` / `list monomers`: one line per polymerizable monomer, zwitterions and comonomers (name,
  aliases, family). `list initiators`:
  name, aliases, counterion. `list forcefields`: name and whether it is available here.
- `template JOB` prints a complete, commented input file for that job that runs as printed: DMAPS,
  KPS on one end (CH3 on the other), OPLS-AA with the built-in RESP charges and CHARGE_SCALE = 0.85,
  and for `bulk`/`slab` TIP4P/2005 with the Madrid-2019 ions and the built-in packer (values the
  built-in packer can place; the user sets the real ones). `template` alone prints the `build` template.
  `examples/inputs/` holds these files and a copolymer example.

## Feature API <!-- rq-44600338 -->

- `zpoly.inputfile.parse_input(text: str, source: str = "<input>") -> dict`: parses and checks <!-- rq-cfe286dc -->
  an input file; returns normalized values (lower-case key names, typed values, defaults filled in).
  Raises `InputError` (carries `line` and `key` when known).
- `zpoly.inputfile.run_input(path) -> object`: parses and runs; returns the job's result <!-- rq-130422bb -->
  (`BuildResult`; the manifest `Path` for `dft`; the folder `Path` for `bulk`/`slab`; `TgResult` or a
  list of `SlabResult` for `analyze`; the written file `Path` for `charges`/`params`). Both are also
  importable from `zpoly`.
- `zpoly.cli.main(argv=None) -> int`: the four commands above. <!-- rq-80888a1d -->

## Out of Scope <!-- rq-f0229af5 -->

- Several jobs in one file; variables or loops in input files (use Python for that).

## Gherkin Scenarios <!-- rq-7128a5a0 -->

```gherkin
Feature: One input file

  @rq-b4454a7b
  Scenario: Help lists every keyword
    When "zpoly help" is run
    Then the output lists JOB, MONOMER, N, INITIATOR, FF, DFT, CHAINS, BOX, SLAB and OUTPUT
    And every line is at most 80 characters wide

  @rq-0af106e4
  Scenario: Help for one keyword
    When "zpoly help box" is run
    Then the output shows the value form, that it applies to bulk and slab, and an example line "BOX = "

  @rq-73efd1ac
  Scenario: Help for an unknown keyword
    When "zpoly help boxx" is run
    Then it exits 1 and stderr suggests "BOX"

  @rq-1c00331e
  Scenario: List monomers and initiators
    When "zpoly list" and "zpoly list initiators" are run
    Then the first output has a line starting with "DMAPS" and the second a line starting with "KPS"

  @rq-4ae24f96
  Scenario: Templates run as printed
    Given the fake backend is the force field
    When "zpoly template build" is written to build.inp, FF is replaced by fake, and the file is run
    Then it exits 0 and writes DMAPS_n12_KPS.pdb, .top and .ff/

  @rq-ff28afbc
  Scenario: Build job
    Given the fake backend
    When an input "MONOMER = DMAPS / N = 3 / INITIATOR = KPS / FF = fake / OUTPUT = d3" is run
    Then d3.pdb, d3.top and d3.ff/ exist

  @rq-c68591e1
  Scenario: PDB only
    When an input "MONOMER = DMAPS / N = 3 / FF = none" is run
    Then only DMAPS_n3.pdb is written

  @rq-244b2f07
  Scenario: Copolymer forms
    Given the fake backend
    When inputs with COPOLYMER = "alternating DMAPS CBMA" (N = 4), "block DMAPS:2 CBMA:2", "random DMAPS:0.5 CBMA:0.5" (N = 4) and "pattern AAB A=DMAPS B=CBMA" (N = 3) are run with FF = none
    Then they write DMAPS-alt-CBMA_n4.pdb, DMAPS2-CBMA2.pdb, DMAPS-ran-CBMA_n4_s0.pdb and AAB_n3.pdb

  @rq-ff7e6eb6
  Scenario: Spacer and alkyl
    When an input "MONOMER = DMAPS / N = 3 / SPACER = 5 / ALKYL = 2 / FF = none" is run
    Then DMAPS_n3_sp5_alk2.pdb is written

  @rq-eda6a6c5
  Scenario: DFT job
    When an input "JOB = dft / MONOMER = DMAPS / N = 3 / INITIATOR = KPS / DFT = q" is run
    Then q/DMAPS.xyz, q/KPS_end.xyz and q/dft_manifest.yaml exist

  @rq-6eca734c
  Scenario: Bulk job
    Given the fake backend and the fixture solvent force field
    When an input "JOB = bulk / MONOMER = DMAPS / N = 3 / INITIATOR = KPS / FF = fake / CHAINS = 2 / BOX = 3 / PACKER = builtin / SOLVENT_FF = <fixture> / TG = 500 300 50" is run
    Then the folder DMAPS_n3_KPS_bulk has system.top, run.sh and five cool_<T>K.mdp files

  @rq-33ff1048
  Scenario: Slab job
    Given the fake backend and the fixture solvent force field
    When an input "JOB = slab / ... / CHAINS = 4 / BOX = 3.4 3.4 10.2 / SLAB = 3.4 / SALT = KCl 0.1 / TEMPERATURES = 300 320 / PACKER = builtin" is run
    Then the folder has prod_300K.mdp and prod_320K.mdp, and run.sh uses "-pname K -nname CL -conc 0.1"

  @rq-1c3e7626
  Scenario: Analyze job detects the folder kind
    Given a bulk folder without energy files
    When an input "JOB = analyze / OUTPUT = <folder>" is run
    Then it exits 1 with an error naming the first missing cooling step

  @rq-390e287d
  Scenario: Charges and params templates
    Given the fake backend
    When inputs "JOB = charges / MONOMER = DMAPS / N = 3 / FF = fake" and "JOB = params / MONOMER = DMAPS / N = 3 / FF = fake" are run
    Then DMAPS_n3_charges.yaml and DMAPS_n3_params.yaml are written

  @rq-d8db093e
  Scenario: Unknown keyword
    When an input with the line "MONMER = DMAPS" is run
    Then it exits 1 and the error names line 1 and suggests "MONOMER"

  @rq-1438d114
  Scenario: Repeated keyword
    When an input with "N = 3" on two lines is run
    Then it exits 1 and the error names the second line

  @rq-20746e76
  Scenario: Bad value
    When an input "MONOMER = DMAPS / N = three" is run
    Then it exits 1 and the error says N needs an integer

  @rq-f5c6245f
  Scenario: Keyword for another job
    When an input "JOB = bulk / ... / SLAB = 3" is run
    Then it exits 1 and the error says "SLAB is for JOB = slab"

  @rq-3d1ea632
  Scenario: Missing polymer or required system keys
    When an input "JOB = bulk / MONOMER = DMAPS / N = 3" (no CHAINS, BOX) is run
    Then it exits 1 and the error names CHAINS and BOX
    And an input with both MONOMER and SEQUENCE exits 1

  @rq-7474fcd4
  Scenario: Errors are reported, not raised
    When an input "MONOMER = NOPE / N = 3 / FF = none" is run
    Then it exits 1, stderr starts with "error:" and mentions "NOPE", and no file is written
```
