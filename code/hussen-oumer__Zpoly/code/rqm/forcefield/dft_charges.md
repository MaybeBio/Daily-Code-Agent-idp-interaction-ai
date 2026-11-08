# Feature: DFT Charges (Export Monomers, Import QM Charges) <!-- rq-84dc267b -->

Charges for zwitterionic polymers should come from quantum chemistry (B3LYP, HF, RESP, …). This
feature splits the charge step out of force-field generation:

1. **Export**: the program writes each monomer, and each initiator end group, as a small
   methyl-capped molecule (`.xyz` and `.pdb`) for the user's DFT program, plus a manifest.
2. **DFT**: the user runs their own calculation and saves the charges next to each exported file.
3. **Import**: `DFT = folder` in the input file reads those charges, maps them onto every
   residue built from that monomer (F, R, E, and the initiator residues), and balances each
   residue to its integer charge.

When every residue of the chain has DFT charges, the backend only assigns atom types (for GAFF2,
`antechamber` without AM1-BCC), which makes force-field generation fast. Residues without DFT
charges fall back to the backend's charges (AM1-BCC for GAFF2) with a notice.

GAFF2 is the default force field (`FF = gaff2`, see *Input File*).

## Exported Molecules <!-- rq-9a042027 -->

- A monomer's DFT molecule is its **repeat unit with methyl caps**: each attachment point of
  `repeat_unit_smiles()` becomes a CH3 group, so the backbone carbons see a carbon neighbour as in
  the chain. For DMAPS: `CH3–CH2–C(CH3)(COO–CH2CH2–N⁺(CH3)2–(CH2)3–SO3⁻)–CH3`. This is exactly the
  molecule the chain builder uses as its template (`rqm/polymer_structure.md`), with the same atom
  order and hydrogens.
- An initiator's DFT molecule is its end-group fragment with a methyl cap on the attachment point
  (for KPS: `CH3–O–SO3⁻`, net charge −1). When the head and tail fragments are the same molecule,
  one export serves both ends.
- Geometry: the builder's first template conformer (ETKDG + MMFF94), in Å. This is a starting
  geometry; the user may optimize it, but **the atom order must stay the same**.
- Files per molecule `<name>`: `<name>.xyz` (comment line: name, net charge, multiplicity 1) and
  `<name>.pdb` (one residue `DFT`, atom names element + number, CONECT records). Names are
  `<MONOMER>` (e.g. `DMAPS`, `DMAPS-sp4`) and `<INITIATOR>_head` / `<INITIATOR>_tail`, or
  `<INITIATOR>_end` when both ends share one molecule.
- `dft_manifest.yaml` in the same folder lists each molecule: `name`, `files`, `templates` (the
  template SMILES with `[*:1]`/`[*:2]` attachment points that this molecule covers; two for a
  shared initiator molecule), `net_charge`, `multiplicity`, `atoms` (element per atom, in file
  order), and `caps` (1-based indices of the cap atoms: cap carbons and their hydrogens). It also
  records the zpoly version.

## Charge Files <!-- rq-17608b5f -->

For molecule `<name>`, the importer looks next to `<name>.xyz` for the first of these files:

| File | Format |
|---|---|
| `<name>.charges` | Plain list: one charge per atom in file order; each line is `q`, `element q`, or `index element q`. Blank lines and `#` comments are ignored. |
| `<name>.mol2` | Tripos mol2 (e.g. from RESP / antechamber): column 9 of `@<TRIPOS>ATOM`. |
| `<name>.log` | Gaussian output: the **last** `ESP charges:` block (RESP/ChelpG/MK); if none, the last `Mulliken charges:` block. |
| `<name>.out` | ORCA output: the **last** `CHELPG Charges` block; if none, the last `MULLIKEN ATOMIC CHARGES` block. |

- A manifest entry may set `charges: <file>` to use a specific file instead (relative to the
  folder), and `charge_type: esp | chelpg | mulliken` to pick the block from a Gaussian/ORCA output.
- The number of charges must equal the molecule's atom count, and the elements (where the format
  gives them) must match the exported atoms in order; otherwise `DftChargeError` names the file and
  the first mismatching atom.
- The sum of the DFT charges must equal the molecule's net charge within 0.01 e; otherwise
  `DftChargeError` (this usually means a different charge state was computed).
- A molecule with no charge file is reported as "missing" and its residues fall back to the
  backend's charges.

## Mapping Charges onto the Chain <!-- rq-31ab6f96 -->

- Every chain atom remembers its template and template atom index (see `rqm/polymer_structure.md`).
  A chain residue built from a template that has DFT charges takes, for each of its atoms, the DFT
  charge of the same template atom.
- Cap atoms are not part of any residue, so their charge is dropped. The one hydrogen that the
  builder adds to an uncapped chain end has no DFT counterpart and starts at 0.
- Each residue is then shifted equally on all its atoms to its integer target charge and rounded
  as usual (`rqm/forcefield/core.md`). The shift is recorded in the charge report as a DFT entry.
  It is normally small (about the cap charge spread over the residue's atoms).
- Identical residue names still have identical charges, because they come from the same template
  and the same shift.
- A charge file (`CHARGES`) is applied **after** the DFT charges, so single atoms can still be
  overridden by hand.
- If every residue of the chain has DFT charges, the backend is called with `charges=False` (types
  and parameters only). Otherwise the backend computes charges, and the CLI prints one notice:
  `note: no DFT charges for <names>; using <backend> charges (write the DFT molecules with JOB = dft)`.

## Backend Change <!-- rq-8a4726a9 -->

- `ForceFieldBackend.parametrize(mol, net_charge, work_dir, *, charges=True)`: with
  `charges=False` the backend may skip its charge model and return zero charges.
- GAFF2 with `charges=False` runs `antechamber … -at gaff2 -pf y -s 0` **without** `-c bcc` (no `sqm`),
  then `parmchk2` as usual.

## Ion Fragments (`DFT_PARTS = ions`) <!-- rq-d4bd3142 -->

A starting point for charge scaling and for ion-pairing studies, as for ionic liquids: each charged
group of a monomer is also computed on its own, as a small, exactly charged ion.

- **Ionic groups** (`ionic_groups(mol)`), found on the molecule's formal charges:
  - the *core*: each charged heavy atom; for an anion on O/S bonded to a central atom (S, P, C), that
    central atom and all its terminal O/S atoms (SO3⁻, COO⁻, PO4⁻); for a charged atom in an
    aromatic ring, the whole aromatic ring system (imidazolium, pyridinium);
  - the *shell*: the heavy atoms bonded to the core; a shell atom that is not carbon (an ester O of a
    phosphate) brings its carbon neighbours too;
  - the hydrogens of core and shell atoms;
  - groups that share atoms are merged (their charges add up, e.g. a DABCO dication, or a spacer of
    one CH2 between cation and anion).
- **Fragment molecules**: the group's atoms cut out of the methyl-capped monomer, keeping its
  coordinates; every cut bond gets a hydrogen cap at 1.09 Å along the bond. Examples: DMAPS →
  `DMAPS_cation` (an N(CH2R)(CH3)3⁺-like ion with the four α carbons as methyls) and
  `DMAPS_anion` (CH3SO3⁻). Names are `<monomer>_cation`, `<monomer>_anion` or `<monomer>_ion`
  (net charge 0, merged groups), numbered `_cation2`, … when a monomer has several.
- **Linker fragment** (the methacrylate / acrylate ester or acrylamide): the carbonyl C, its =O, the
  ester O or amide N (with its H), and the carbon bonded to that O/N (with its H). The DFT molecule
  also keeps the backbone carbon bonded to the carbonyl and that carbon's carbon neighbours, all
  capped with H: DMAPS (methacrylate) → methyl pivalate (CH3)3C–C(=O)O–CH3; an acrylate →
  methyl isobutyrate; an acrylamide → N-methyl isobutyramide. Only the linker atoms take its
  charges (the backbone atoms keep the force field's), listed in the manifest as `apply`. Its name is
  `<monomer>_ester` or `<monomer>_amide`. Monomers without such a linker (styrenics) have none.
- `JOB = dft` with `DFT_PARTS = ions` writes the whole methyl-capped monomer **and** its fragments
  (`.xyz`/`.pdb`); the end groups are unchanged (they are already ions). The manifest records, per
  fragment, `parent` (the monomer molecule), `parent_atoms` (1-based atom of the parent for each
  fragment atom, 0 for a cap), `caps` (`[cap, atom it is bonded to]`) and `apply` (the fragment
  atoms whose charges are used; all mapped atoms for the ions).
- Loading checks a fragment file like any charge file (count, elements, total = its net charge
  within 0.01 e). A cap's charge is added to the atom it is bonded to, so the group sums to its
  net charge.
- Building with `DFT_PARTS = ions`: the ionic and linker atoms take the fragment charges and stay fixed; the
  other atoms of the residue take the whole-monomer DFT charges (or the force field's, with the
  usual notice, when the whole monomer has no charge file) and absorb the balance, so the residue
  keeps its integer total. A user `CHARGES` file still overrides single atoms. A monomer with
  ionic groups but no fragment charges is reported as missing, like a missing whole molecule.
- `DFT_PARTS = whole` (default) ignores fragment files.
- Zwitterionic liquids (the library's small molecules): `JOB = dft` with `MOLECULE = <name>` writes
  the molecule (no caps) and, with `DFT_PARTS = ions`, its fragments, into the same folder.
  (Building their topology is a separate feature.)

## Charge Scaling (`CHARGE_SCALE`) <!-- rq-12060285 -->

`CHARGE_SCALE = s` (0 < s ≤ 1, default 1; e.g. 0.8 for the ECC approach) in `build`, `bulk` and
`slab`, and `parametrize(..., charge_scale=s)`:

- Every atom of an ionic group of the chain (the zwitterion's cation and anion groups, charged
  comonomers, the initiator's SO4⁻ end groups) has its charge multiplied by s. For scaling, a group
  stays within the residue of its charged atom: the SO4⁻ end group scales its own O–SO3 atoms, not
  the backbone CH2 it is bonded to. The linker (ester) is not ionic and is not scaled. Each residue's
  target becomes s × its formal charge; the residue's other atoms absorb the small rest (or the group
  atoms, if the residue has no others), with the same 4-decimal rounding to the exact target.
- Free ions are scaled by the same factor: the counterions and the salt (`NA`, `K`, `CL`) in the
  `ions.itp` written for `bulk` and `slab`, so the system stays exactly neutral. Water and the
  neutral parts of the polymer are not scaled.
- The topology records `charge_scale`; `.top` and `system.top` say `charges scaled by s`. A
  single-chain `build` notes that its counterions must be scaled by s too.
- `protocol.yaml` of `bulk`/`slab` records `charge_scale`.

## Feature API <!-- rq-905bf2b0 -->

Importable from `zpoly.forcefield`:

- `export_dft_inputs(monomer: Zwitterion | str, out_dir, *, initiator: Initiator | str | None = None, <!-- rq-396afb23 -->
  library: Library | None = None) -> Path`: writes the molecules and `dft_manifest.yaml`, returns
  the manifest path. Existing files with the same names are overwritten; other files in the folder
  (for example charge files from an earlier run) are left alone. If the manifest already exists,
  its entries for other molecules are kept, so several monomers can share one folder.
- `load_dft_charges(folder) -> DftCharges`: reads the manifest and every available charge file. <!-- rq-8c171443 -->
  `DftCharges` has `by_template: Mapping[str, tuple[float, ...]]` (charges per template, in file
  order), `sources: Mapping[str, str]` (molecule name → file used), and `missing: tuple[str, ...]`.
- `export_dft_inputs(..., parts="whole" | "ions")`; `export_dft_molecule(name_or_entry, out_dir, *, <!-- rq-d0c21c6a -->
  parts="whole", library=None) -> Path` for small molecules. `ionic_groups(mol) -> list[IonGroup]`
  (`atoms: tuple[int, ...]`, `charge: int`). `DftCharges.ions_by_template: Mapping[str, Mapping[int,
  float]]` (template atom → fragment charge).
- `parametrize(..., dft_parts="whole", charge_scale=1.0)`; `Topology.charge_scale`, <!-- rq-aa6ed061 -->
  `Topology.ion_blocks` (residues whose ionic atoms came from fragments).
- `parametrize(..., dft_charges: DftCharges | str | os.PathLike | None = None)`: applies them as <!-- rq-c0666681 -->
  described above.
- `DftChargeError(ForceFieldError)`: carries `path` and `atom` (1-based) when known. <!-- rq-9e291eae -->

## Input File <!-- rq-d125cfec -->

- `JOB = dft` writes the molecules of every monomer of the polymer (and the initiator end groups)
  and the manifest into the `DFT` folder (default `dft`), and prints, per molecule, its file, atom
  count, net charge and multiplicity, plus the next steps.
- `DFT = folder` in a `build`, `bulk` or `slab` input applies them (needs a force field). The
  summary reports `DFT charges: <n> residue blocks from <folder>` (and `ion fragments: <n>` with
  `DFT_PARTS = ions`).
- `DFT_PARTS = whole | ions` (jobs `dft`, `build`, `bulk`, `slab`) and `CHARGE_SCALE` (`build`,
  `bulk`, `slab`) as described above; `MOLECULE` (job `dft` only) exports a zwitterionic liquid.
- `FF` defaults to **`gaff2`**; `FF = none` writes only the PDB.

## Out of Scope <!-- rq-d01cc120 -->

- Running the DFT calculation itself; choosing the functional, basis set or ESP fitting scheme.
- Multi-conformer charge averaging; charge derivation for the whole chain at once.

## Gherkin Scenarios <!-- rq-ca38b1b0 -->

```gherkin
Feature: DFT charges

  @rq-0a9d9a0e
  Scenario: Export a monomer with methyl caps
    When export_dft_inputs("DMAPS", tmp) is called
    Then tmp contains DMAPS.xyz, DMAPS.pdb and dft_manifest.yaml
    And DMAPS.xyz has as many atoms as the builder's DMAPS template, with net charge 0 on its comment line
    And the manifest's caps for DMAPS are the two cap carbons and their six hydrogens

  @rq-ca2b6f77
  Scenario: Export an initiator end group
    When export_dft_inputs("DMAPS", tmp, initiator="KPS") is called
    Then tmp also contains KPS_end.xyz with 9 atoms (CH3-O-SO3-) and net charge -1
    And the manifest has one KPS entry used for both chain ends

  @rq-dea7e8ed
  Scenario: Export keeps other manifest entries
    Given export_dft_inputs("DMAPS", tmp) was called
    When export_dft_inputs("CBMA", tmp) is called
    Then the manifest lists both DMAPS and CBMA

  @rq-5abcd130
  Scenario: Plain charge list
    Given an exported DMAPS folder and a DMAPS.charges file with "element q" lines summing to 0
    When load_dft_charges(folder) is called
    Then by_template has those charges for the DMAPS template and sources names DMAPS.charges

  @rq-1132056a
  Scenario: Charge list with the wrong count
    Given a DMAPS.charges file with one charge too few
    When load_dft_charges(folder) is called
    Then DftChargeError is raised naming the file

  @rq-0bf14ca4
  Scenario: Element mismatch
    Given a DMAPS.charges file whose first line says "O" while atom 1 is C
    When load_dft_charges(folder) is called
    Then DftChargeError is raised with atom 1

  @rq-e24ad818
  Scenario: Wrong total charge
    Given a DMAPS.charges file whose charges sum to 1.0
    When load_dft_charges(folder) is called
    Then DftChargeError is raised mentioning the net charge 0

  @rq-d2d065e8
  Scenario: mol2 charges
    Given a DMAPS.mol2 with the exported atoms in order and charges summing to 0
    When load_dft_charges(folder) is called
    Then the charges are read from column 9

  @rq-108f0cc1
  Scenario: Gaussian ESP charges
    Given a DMAPS.log with a "Mulliken charges:" block followed by an "ESP charges:" block
    When load_dft_charges(folder) is called
    Then the ESP charges are used

  @rq-b599705d
  Scenario: ORCA CHELPG charges
    Given a DMAPS.out with "MULLIKEN ATOMIC CHARGES" and "CHELPG Charges" blocks
    When load_dft_charges(folder) is called
    Then the CHELPG charges are used

  @rq-b90af318
  Scenario: Missing charge file
    Given an exported folder with no charge file for DMAPS
    When load_dft_charges(folder) is called
    Then missing contains "DMAPS"

  @rq-90737175
  Scenario: DFT charges mapped onto F, R and E
    Given DFT charges for DMAPS and KPS_end in an exported folder
    When a 5-unit KPS-capped DMAPS chain is parametrized with the fake backend and those DFT charges
    Then every DMR atom's charge equals its template atom's DFT charge plus the DMR shift
    And every block sums exactly to its integer charge
    And the backend was called with charges=False every time

  @rq-fa62c105
  Scenario: Partial DFT charges fall back
    Given DFT charges for DMAPS only
    When a 3-unit KPS-capped DMAPS chain is parametrized with them
    Then the DM* blocks use DFT charges, the KP* blocks use the backend's charges
    And the backend was called with charges=True

  @rq-a8df24ac
  Scenario: Charge file after DFT charges
    Given DFT charges for DMAPS and a charge file {"DMR": {"N1": 0.2}}
    When a 3-unit DMAPS chain is parametrized with both
    Then N1 in DMR is exactly 0.2

  @rq-98ce435e
  Scenario: GAFF2 without charges skips AM1-BCC
    Given fake antechamber and parmchk2 programs that record their arguments
    When the GAFF2 backend parametrizes a fragment with charges=False
    Then antechamber was called without "-c"

  @rq-684b318c
  Scenario: CLI export, check and build
    Given the fake backend is registered
    When an input "JOB = dft / MONOMER = DMAPS / N = 3 / INITIATOR = KPS / DFT = dft" is run
    And charge files are written for DMAPS and KPS_end
    And an input "MONOMER = DMAPS / N = 3 / INITIATOR = KPS / FF = fake / DFT = dft / OUTPUT = out/d" is run
    Then both exit 0 and the build reports "DFT charges: 5 residue blocks"

  @rq-4d98df6a
  Scenario: CLI notice when DFT charges are incomplete
    Given the fake backend and an exported folder with charges for DMAPS only
    When an input "MONOMER = DMAPS / N = 3 / INITIATOR = KPS / FF = fake / DFT = dft / OUTPUT = out/d" is run
    Then it exits 0 and stderr contains "note: no DFT charges for KPS_end"

  @rq-10a45821
  Scenario: GAFF2 is the default force field
    When an input "MONOMER = DMAPS / N = 3" is run on a machine without AmberTools
    Then it exits 1 with an error about AmberTools
    And the same input with "FF = none" exits 0 and writes only the PDB

  @rq-013e1066
  Scenario: Ionic groups of DMAPS
    When ionic_groups is called on the methyl-capped DMAPS molecule
    Then there are two groups: charge +1 with N and its 4 carbon neighbours, and charge -1 with S, its 3 O and the CH2 bonded to S

  @rq-39505fe6
  Scenario: Merged groups
    When ionic_groups is called on the VBD comonomer and on a sulfobetaine with a one-CH2 spacer
    Then VBD has one group of charge +2 and the sulfobetaine one group of charge 0

  @rq-cefbfada
  Scenario: Export ion fragments
    When an input "JOB = dft / MONOMER = DMAPS / N = 3 / INITIATOR = KPS / DFT = q / DFT_PARTS = ions" is run
    Then q has DMAPS.xyz, DMAPS_cation.xyz (charge +1, 17 atoms: N, 4 C, 12 H), DMAPS_anion.xyz (CH3SO3-, charge -1, 8 atoms), DMAPS_ester.xyz (methyl pivalate, 20 atoms, 6 applied: C, 2 O, C and its 2 H) and KPS_end.xyz
    And every fragment atom of the manifest maps to a parent atom of the same element, or is a cap

  @rq-f8346d02
  Scenario: Fragment charges fixed in the build
    Given the fake backend, whole-monomer charges and fragment charges for DMAPS
    When the chain is parametrized with dft_parts="ions"
    Then the ionic atoms of each DM residue have the fragment charges (cap charge added to its atom, within 0.0001)
    And every residue still sums to its integer target, and ion_blocks lists the DM residues

  @rq-7670d71d
  Scenario: Missing fragment charges
    Given fragment files are missing
    When the chain is parametrized with dft_parts="ions"
    Then dft_missing names DMAPS_cation and DMAPS_anion

  @rq-e14e60c1
  Scenario: Charge scaling of the chain
    Given the fake backend
    When a DMAPS 3-mer with KPS ends is parametrized with charge_scale=0.8
    Then the KPS end residues sum to -0.8, the DM residues to 0, and the ionic atoms are 0.8 times the unscaled ones (within the balancing rest)

  @rq-db0975d0
  Scenario: Charge scaling of the system
    Given the fake backend and the fixture solvent force field
    When a bulk system with CHARGE_SCALE = 0.8 is written
    Then ions.itp has K with charge 0.8, system.top says "charges scaled by 0.8", and the total charge is 0

  @rq-c7326013
  Scenario: Zwitterionic liquid export
    When an input "JOB = dft / MOLECULE = <a library small molecule> / DFT_PARTS = ions / DFT = z" is run
    Then the molecule and its fragments are written to z

  @rq-a3491212
  Scenario: Bad scale
    When an input with CHARGE_SCALE = 0 or 1.5 is run
    Then it exits 1 with an error that CHARGE_SCALE must be in (0, 1]
```
