# Feature: Polymer Chain Builder (3D Structure and PDB) <!-- rq-5313ed51 -->

Builds a 3D structure of one polymer chain from a user's description, for example
"50 DMAPS with KPS at both ends" or "30 units of a random DMAPS/CBMA copolymer", and writes it as a
PDB file. It is the construction layer of `rqm/ARCHITECTURE.md`. It uses the monomers, initiators
and `build_chain` connectivity of `rqm/zwitterion_library.md` and needs no force field.

The feature has three parts:

1. **Sequence specification**: helpers that turn a description (homopolymer, alternating, block,
   random, explicit pattern or list) into an ordered list of monomer names.
2. **3D construction by unit-by-unit growth**: each distinct repeat unit and end group is
   embedded once with RDKit, and copies are placed along the backbone with clash checks. The chain
   is either a self-avoiding random coil (default) or extended. An optional UFF cleanup follows.
3. **Output**: a GROMACS-friendly PDB file, through a Python API and a `zpoly` command.

## Dependencies <!-- rq-93aa97d4 -->

No new runtime dependencies: RDKit, NumPy (installed with RDKit) and PyYAML. The command-line
interface uses `argparse`.

## Sequence Specification <!-- rq-b08ea041 -->

All counts are **total units in the chain**.

- **Homopolymer**: `n` copies of one monomer.
- **Alternating**: `n` units alternating A, B, A, B, …, starting with A. An odd `n` has one more A.
- **Block**: a list of `(monomer, count)` blocks in order, for example `[("DMAPS", 10), ("CBMA", 20)]`.
  Repeating a multiblock is written by listing the blocks again.
- **Random (statistical)**: a composition `{monomer: fraction}` and a total `n`. Fractions must be
  positive and sum to 1 within 1e-6. Unit counts are the fractions × `n` rounded by the
  largest-remainder method, so they sum exactly to `n`. The units are then shuffled with a seeded
  generator, so the same seed gives the same sequence.
- **Pattern**: a pattern string of single-letter symbols (for example `"AAB"`), a mapping from each
  symbol to a monomer, and a total `n`. The pattern repeats and is cut at `n` units.
- **Explicit list**: any list of monomer names or `Zwitterion` objects, used as given.

## 3D Construction <!-- rq-ee0a0744 -->

### Repeat-unit templates <!-- rq-57e067c6 -->

- The template for a monomer is its `repeat_unit_smiles()`. Each attachment point is replaced by a
  temporary carbon cap, so the attachment atoms get correct hydrogen counts and geometry. The
  template is embedded with hydrogens using RDKit ETKDGv3, giving up to 20 conformers, and each
  conformer is optimized with MMFF94 (UFF if MMFF94 has no parameters).
- **Second backbone torsion**: when the head atom has exactly one heavy neighbour besides its cap,
  and that bond is not in a ring, each conformer is also copied with everything beyond that bond
  rotated by 120° and 240°. This samples the head–Cα torsion, which conformer embedding alone
  covers poorly in crowded methacrylate backbones.
- Initiator head and tail fragments are templates the same way: a head fragment has an outgoing
  attachment, and a tail fragment has an incoming one.
- A chain end without an initiator end group is capped by `cap` (input key `CAP`):
  - `"CH3"` (default): every such end becomes a methyl group. The head CH2 gains one H (it becomes
    CH3), and the tail's substituted carbon gets a CH3 (the template's own methyl cap, with its
    coordinates). Example: `CH3–C(CH3)(R)–CH2–…–CH2–C(CH3)(R)–CH3`.
  - `"H"`: the head CH2 becomes CH3 and the tail carbon gains an H.
  - The cap atoms belong to the end unit's residue and come **after** its other atoms, so the
    heavy-atom names of the F, R and E residues are the same (`C4`, `N1`, `S1`, …).
  - With `INITIATOR` on one end only (`ENDS = head` / `tail`), the other end gets the cap.
- A template is embedded once per build and reused for every copy.

### Growth <!-- rq-7414a764 -->

- Units are placed in sequence order. The new unit's head atom goes where the previous unit's
  outgoing cap was, and its incoming cap points at the previous unit's tail atom. This makes the
  new backbone bond lie along the previous cap bond.
- The new unit is then rotated about that bond: 24 angles for each conformer (including the
  second-torsion copies). The **scored points** are the new unit's heavy atoms **plus its outgoing
  cap position**. Including the cap position matters, because it fixes where the next unit's head
  atom will be.
- Each placement is scored against every heavy atom already placed. A pair is **not** scored when
  the two atoms are 3 or fewer bonds apart across the junction, since bond geometry fixes their
  distance. This is the same rule `min_nonbonded_distance` uses.
  - The **clash score** is the sum of (3.2 Å − d)² over all scored distances d < 3.2 Å.
  - A placement is **clash-free** when no scored distance is below 2.6 Å.
- **Coil** conformation: conformers are scanned in a random (seeded) order, and the scan stops once
  48 clash-free placements have been found. One placement is chosen at random from the roomier half
  of the clash-free placements (clash score at or below their median). If none is clash-free, the
  placement with the lowest clash score is chosen.
- **Extended** conformation: all conformers are scanned. The clash-free placement whose outgoing cap
  (or, for the last unit, centroid) is farthest from the chain's first heavy atom is chosen. If none
  is clash-free, the lowest clash score wins.
- **Backtracking**: if the chosen placement has a scored distance below the attempt's limit, the
  builder removes placed units and grows them again with fresh random choices. The first failure at
  a unit removes 3 units, and each repeated failure at that unit doubles the depth (3, 6, 12, …).
  One growth attempt allows 30 backtracks.
- **Restarts**: each attempt uses a seed derived from `seed` and the attempt number. Up to 3
  **strict** attempts use a limit of 1.8 Å. If all fail, up to 2 **loose** attempts use a limit of
  1.0 Å. A loose result always gets a UFF cleanup of at least 500 steps, which pushes the close
  contacts apart; if UFF has no parameters, `BuildError` is raised. If every attempt fails,
  `BuildError` is raised with advice to try another seed or the extended conformation.
- After strict growth, the whole chain can be minimized with UFF for `relax_iterations` steps
  (default 0, meaning no cleanup: strict growth already places units clash-free, and whole-chain UFF
  costs about 0.1 s per step for 2,000 atoms). If UFF has no parameters for some atom, a requested
  relaxation is skipped with a logged warning.
- `PolymerStructure.relaxed` records whether a UFF cleanup ran.
- The assembled structure's heavy-atom connectivity is checked against `build_chain`'s SMILES; a
  mismatch is an internal error (`BuildError`).
- The final structure's `min_nonbonded_distance` is the smallest distance between heavy atoms more
  than 3 bonds apart. If it is below 1.5 Å after relaxation, the builder raises `BuildError`.
- The same inputs and seed give identical coordinates.

### Residues and atom names <!-- rq-65ea2a82 -->

- Each repeat unit is one residue; residues are numbered 1… in chain order, with the head end group
  first and the tail end group last when present.
- Monomer residue names are three characters: a **two-character prefix** plus a **chain-position
  letter**: `F` (Front) for the first monomer unit, `E` (End) for the last, and `R` (Repeat) for
  every unit in between. A one-unit chain uses `F`. For example, 12 DMAPS units give
  `DMF DMR … DMR DME`. Topology files (`rqm/forcefield/core.md`) use the same names, so
  `gmx pdb2gmx` matches them.
- A monomer's prefix is its library `residue` field when set (see `rqm/zwitterion_library.md`).
  Otherwise it is the first two alphanumeric characters of the monomer name, uppercased. If two
  different monomers in one chain would get the same prefix, the later one gets its second character
  replaced by `1`–`9`, then `A`–`Z`, choosing the first free one.
- End-group residues are the initiator's first two alphanumeric characters (uppercased) plus `H`
  (head) or `T` (tail); for example KPS gives `KPH` and `KPT`.
- Atom names are the element symbol plus a running number per element within the residue
  (`C1`, `C2`, …, `O1`, …, `H1`, …), at most 4 characters.

## PDB Output <!-- rq-18964b06 -->

- `ATOM` records, chain ID `A`, residue by residue in chain order. Within a residue, each heavy
  atom (in template order) is followed directly by its hydrogens. This order is what
  `gmx pdb2gmx` expects.
- One `CONECT` record set covering every bond, with up to 4 partners per line.
- A `CRYST1` record with an orthorhombic box equal to the coordinate bounding box plus
  `padding` (default 10 Å) on every side. The coordinates are translated to the box centre.
- The file ends with `TER` and `END`.
- PDB field limits: more than 99,999 atoms or more than 9,999 residues raises `PdbLimitError`.

## Feature API <!-- rq-6cd2ece7 -->

All public names are importable from `zpoly.structure`.

### Functions <!-- rq-5c8da3b8 -->

- `homopolymer(monomer: str, n: int) -> list[str]` <!-- rq-d34a145b -->
- `alternating(a: str, b: str, n: int) -> list[str]` <!-- rq-9c67798c -->
- `block(blocks: Sequence[tuple[str, int]]) -> list[str]` <!-- rq-88a61c9d -->
- `random_copolymer(composition: Mapping[str, float], n: int, seed: int = 0) -> list[str]` <!-- rq-ff223b41 -->
- `pattern(pattern: str, mapping: Mapping[str, str], n: int) -> list[str]` <!-- rq-6652f857 -->
  - All sequence helpers raise `SequenceError` for `n < 1`, a block count below 1, an empty
    block list, fractions that are not positive or do not sum to 1, an empty pattern, or a pattern
    symbol missing from the mapping. They do not check that monomer names exist; `build_polymer`
    does.
- `build_polymer(sequence: Sequence[str | Zwitterion], *, initiator: Initiator | str | None = None, <!-- rq-b9773ded -->
  ends: Ends | str | None = None, conformation: str = "coil", seed: int = 0,
  relax_iterations: int = 0, library: Library | None = None, cap: str = "CH3") -> PolymerStructure`
  (`build_chain` takes the same `cap`; `Chain.cap` records it; anything but `"CH3"`/`"H"` is a
  `ValueError`)
  - Resolves names and end groups exactly as `build_chain`, and raises the same errors (for
    example `UnknownEntryError`, `EndGroupError`, `ChainError`).
  - `conformation` is `"coil"` or `"extended"`; anything else raises `ValueError`.
  - The structure's connectivity is identical to `build_chain(...)` for the same arguments.
- `write_pdb(structure: PolymerStructure, path: str | os.PathLike, padding: float = 10.0) -> Path`: <!-- rq-671d8f4c -->
  writes the file, creating parent directories, and returns its path. It raises `PdbLimitError`
  for oversized structures and `OSError` for filesystem failures.

### Types <!-- rq-8e5557de -->

- `PolymerStructure` (frozen dataclass): <!-- rq-12389ffb -->
  - `mol: rdkit.Chem.Mol`: all atoms including hydrogens, with one conformer (Å). Each atom carries
    the properties `res` (residue name), `resnum` (int) and `name` (atom name).
  - `chain: Chain`: the `build_chain` result (SMILES, net charge, counterion, formula, …).
  - `sequence: tuple[str, ...]`: the monomer names in order.
  - `units: tuple[Zwitterion, ...]`: the monomer objects in order (variants and user-library
    monomers included), so later steps can rebuild fragments of the same chain.
  - `initiator: Initiator | None`: the resolved initiator, if any.
  - `templates: tuple[str, ...]`: the template SMILES used (repeat units and end-group fragments,
    with `[*:1]`/`[*:2]` attachment points), in order of first use. Every atom of `mol` also carries
    the integer properties `template` (index into `templates`) and `template_atom` (the atom's index
    in that template's embedded molecule, whose attachment points are methyl caps; −1 for the
    hydrogen added to an uncapped chain end). DFT charges are mapped through these
    (`rqm/forcefield/dft_charges.md`).
  - `residue_names: Mapping[str, str]`: monomer name → its two-character residue prefix, and
    end-group key (`"<INITIATOR>:head"` / `"<INITIATOR>:tail"`) → its residue name (`KPH`/`KPT`).
  - `conformation: str`, `seed: int`.
  - `min_nonbonded_distance: float` (Å).
  - `radius_of_gyration: float` (Å, heavy atoms).
  - `end_to_end_distance: float` (Å): from the first to the last backbone atom.
  - `relaxed: bool`: whether a UFF cleanup ran (requested, or needed after loose growth).
- `StructureError(Exception)`: base class. <!-- rq-280c4e6d -->
  - `SequenceError(StructureError, ValueError)`: invalid sequence description.
  - `BuildError(StructureError)`: no clash-free structure could be built. Carries
    `min_distance`.
  - `PdbLimitError(StructureError)`: the structure exceeds PDB field limits.

## Input File <!-- rq-1eaf96c1 -->

The console script `zpoly` (`zpoly.cli:main`) runs one input file
(`rqm/input_file.md`). The chain is given by `MONOMER` + `N`, `SEQUENCE`, or `COPOLYMER`
(`alternating A B`, `block A:10 B:20`, `random A:0.5 B:0.5`, `pattern AAB A=X B=Y`), with
`INITIATOR`, `ENDS`, `CONFORMATION` (`coil` | `extended`), `SEED`, `SPACER`, `ALKYL`, `LIBRARY` and
`OUTPUT`. `FF = none` writes only the PDB. `zpoly list [monomers|initiators]` lists the names.

- `SPACER` sets the number of CH2 groups between the cation and the anion, and `ALKYL` sets the
  N-alkyl length, for every monomer of the chain (`make_variant`); the default name gets
  `_sp<n>` / `_alk<n>`.
- The default output name is derived from the polymer (`rqm/input_file.md`, *Jobs*), plus
  `_<INITIATOR>` when end groups are attached.
- The summary line gives the file, units, atoms, net charge, counterion, Rg and the closest
  non-bonded contact. Library and build errors print `error: …`, exit 1, and write nothing.

## Out of Scope <!-- rq-2a066223 -->

- Force-field atom types, charges and topology files (the parametrization feature).
- Multiple chains, solvent, ions and boxes for simulation (the system builder).
- Branched, star, brush and cyclic architectures; tacticity control (the side groups get the
  orientation of the chosen template conformer and rotation).
- Recipe files and the web app (the front-end feature).

## Gherkin Scenarios <!-- rq-46e95474 -->

```gherkin
Feature: Polymer chain builder

  @rq-35976189
  Scenario: Homopolymer sequence
    When homopolymer("DMAPS", 50) is called
    Then it returns 50 times "DMAPS"

  @rq-95c60185
  Scenario: Alternating sequence with odd length
    When alternating("DMAPS", "CBMA", 5) is called
    Then it returns ["DMAPS", "CBMA", "DMAPS", "CBMA", "DMAPS"]

  @rq-2aa4244e
  Scenario: Block sequence
    When block([("DMAPS", 2), ("CBMA", 3)]) is called
    Then it returns ["DMAPS", "DMAPS", "CBMA", "CBMA", "CBMA"]

  @rq-c073d054
  Scenario: Random copolymer composition and reproducibility
    When random_copolymer({"ZB2": 0.56, "SB": 0.44}, 30, seed=1) is called twice
    Then both results are identical
    And they contain exactly 17 "ZB2" and 13 "SB"
    And random_copolymer with seed=2 returns the same counts

  @rq-d87d9e74
  Scenario: Pattern sequence
    When pattern("AAB", {"A": "DMAPS", "B": "CBMA"}, 7) is called
    Then it returns ["DMAPS", "DMAPS", "CBMA", "DMAPS", "DMAPS", "CBMA", "DMAPS"]

  @rq-826f003e
  Scenario: Invalid sequence descriptions
    When each of these is called:
      homopolymer("DMAPS", 0)
      block([])
      block([("DMAPS", 0)])
      random_copolymer({"DMAPS": 0.5, "CBMA": 0.4}, 10)
      random_copolymer({"DMAPS": 1.2, "CBMA": -0.2}, 10)
      pattern("", {"A": "DMAPS"}, 3)
      pattern("AB", {"A": "DMAPS"}, 3)
    Then each raises SequenceError

  @rq-06424eaf
  Scenario: Build a 50-unit homopolymer with KPS ends
    When build_polymer(homopolymer("DMAPS", 50), initiator="KPS") is called
    Then the structure has 52 residues: "KPH", "DMF", 48 × "DMR", "DME", "KPT"
    And its connectivity equals build_chain(["DMAPS"] * 50, initiator="KPS")
    And its net charge is -2
    And its min_nonbonded_distance is at least 1.5 Å

  @rq-6d35271c
  Scenario: Build a 30-unit random copolymer
    When build_polymer(random_copolymer({"DMAPS": 0.5, "CBMA": 0.5}, 30, seed=3)) is called
    Then the structure has 30 monomer residues using two distinct prefixes ("DM" and "CB")
    And the first residue ends in "F", the last in "E", and all others in "R"
    And its connectivity equals build_chain of the same sequence
    And its net charge is 0

  @rq-0660a451
  Scenario: Every atom has 3D coordinates and correct hydrogens
    When build_polymer(homopolymer("DMAPS", 5)) is called
    Then the molecule has one conformer
    And every heavy atom's total hydrogen count equals that in the chain SMILES
    And every bonded heavy-atom pair is between 1.1 Å and 2.0 Å apart

  @rq-d497f481
  Scenario: Same seed gives identical coordinates
    When build_polymer(homopolymer("DMAPS", 10), seed=5) is called twice
    Then both structures have identical coordinates

  @rq-3e9f3397
  Scenario: Extended chain is longer than a coil
    When build_polymer(homopolymer("DMAPS", 20), conformation="extended") is called
    And build_polymer(homopolymer("DMAPS", 20), conformation="coil") is called
    Then the extended end-to-end distance is greater than the coil's

  @rq-7f3fc73b
  Scenario: Unknown conformation
    When build_polymer(homopolymer("DMAPS", 3), conformation="helix") is called
    Then ValueError is raised

  @rq-61c06a74
  Scenario: Unknown monomer in a sequence
    When build_polymer(["DMAPS", "NOPE"]) is called
    Then UnknownEntryError is raised

  @rq-050546f1
  Scenario: Residue name collision is resolved
    When build_polymer(["CBMA", "CBMA-1"]) is called
    Then the residue names are "CBF" and "C1E"

  @rq-c1a4dfaa
  Scenario: Library residue field is used
    Given a user library where monomer "MYSB" has residue "MS"
    When build_polymer(["MYSB"] * 3, library=that library) is called
    Then the residues are named "MSF", "MSR", "MSE"

  @rq-8994396e
  Scenario: Uncapped ends are hydrogen-terminated
    When build_polymer(homopolymer("DMAPS", 3), cap="H") is called
    Then the first residue contains a CH3 group bonded to the first backbone carbon
    And the structure has net charge 0

  @rq-b1961917
  Scenario: Methyl caps
    When build_polymer(homopolymer("DMAPS", 3)) is called (cap="CH3", the default)
    Then the last residue's substituted backbone carbon carries a CH3 that is the residue's last heavy atom
    And DMF, DMR and DME list the same heavy-atom names in the same order before the cap
    And build_polymer(homopolymer("DMAPS", 3), initiator="KPS", ends="head") has KPH at the head and a CH3 at the tail

  @rq-cfb82d2e
  Scenario: Write a PDB file
    Given a structure from build_polymer(homopolymer("DMAPS", 3), initiator="KPS")
    When write_pdb(structure, tmp/"out/dmaps.pdb") is called
    Then the file exists, creating the "out" directory
    And it has one ATOM record per atom, a CRYST1 record, CONECT records, TER and END
    And residue numbers in ATOM records never decrease
    And every hydrogen record directly follows its parent heavy atom or another hydrogen of that parent
    And RDKit reads the file back with the same number of atoms and bonds

  @rq-92413a81
  Scenario: PDB box and centring
    When write_pdb(structure, path, padding=10.0) is called
    Then the CRYST1 box edges equal the bounding-box extent plus 20 Å
    And every coordinate lies at least 10 Å inside the box

  @rq-8b61e864
  Scenario: PDB atom limit
    Given a structure with more than 99,999 atoms
    When write_pdb is called
    Then PdbLimitError is raised and no file is written

  @rq-08692f3c
  Scenario: CLI builds a homopolymer
    When an input "MONOMER = DMAPS / N = 50 / INITIATOR = KPS / FF = none / OUTPUT = dmaps" is run
    Then it exits 0
    And tmp/dmaps.pdb exists
    And the summary line reports 50 units and net charge -2

  @rq-b5e5f3b3
  Scenario: CLI default output name
    When an input "COPOLYMER = alternating DMAPS CBMA / N = 6 / FF = none" is run in an empty directory
    Then the file "DMAPS-alt-CBMA_n6.pdb" is created

  @rq-e83ffc40
  Scenario: CLI random copolymer and block
    When an input "COPOLYMER = random DMAPS:0.5 CBMA:0.5 / N = 10 / SEED = 2 / FF = none" is run
    Then it exits 0 and writes "DMAPS-ran-CBMA_n10_s2.pdb"
    And an input "COPOLYMER = block DMAPS:3 CBMA:2 / FF = none" writes "DMAPS3-CBMA2.pdb"

  @rq-a405b58f
  Scenario: CLI reports a library error
    When an input "MONOMER = NOPE / N = 5 / FF = none" is run
    Then it exits 1 and stderr starts with "error:" and mentions "NOPE"

  @rq-94f765ba
  Scenario: CLI lists initiators
    When "zpoly list initiators" is run
    Then the output contains a line starting with "KPS"

  @rq-35c1237e
  Scenario: CLI spacer and alkyl options
    When an input "MONOMER = DMAPS / N = 5 / SPACER = 5 / ALKYL = 2 / FF = none" is run in an empty directory
    Then it exits 0 and writes "DMAPS_n5_sp5_alk2.pdb"
    And every unit has five CH2 groups between N and S and two ethyl groups on N

  @rq-32089525
  Scenario: CLI variant option a monomer cannot take
    When an input "COPOLYMER = alternating DMAPS SPVI / N = 4 / ALKYL = 2 / FF = none" is run in an empty directory
    Then it exits 1, stderr starts with "error:", and no PDB file is written
```
