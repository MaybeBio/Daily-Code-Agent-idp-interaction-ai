# Feature: Force-Field Topology Core (F/R/E Residues and GROMACS Output) <!-- rq-2d53dd42 -->

Turns a built chain (`PolymerStructure`, `rqm/polymer_structure.md`) into a GROMACS-ready
force-field description. The design aims to be simple for beginners and complete for advanced
users:

- every monomer appears as **Front / Repeat / End** residues (`DMF`, `DMR`, `DME`), and initiator
  end groups as their own residues (`KPH`, `KPT`);
- the output is a self-contained **`<name>.ff/` folder** with a residue topology (`.rtp`), so
  `gmx pdb2gmx` works, plus a **ready `<name>.top`** for the built chain, so beginners can skip
  `pdb2gmx`;
- the chemistry-specific work (atom types, charges, parameters) is done by a **backend**: GAFF2
  (`rqm/forcefield/gaff2.md`) or OPLS-AA (`rqm/forcefield/oplsaa.md`, a later part). This file
  defines everything the backends share.

## Residue Blocks <!-- rq-f928598f -->

- A **residue block** is one distinct residue name in the chain with its atoms (names, types,
  charges), its bonds, and its links to the previous and next residue. All residues with the same
  name are identical: every `DMR` has the same atoms, types and charges.
- Residue names and atom names come from the built structure (`rqm/polymer_structure.md`). Blocks
  are listed in chain order of first appearance, for example `KPH, DMF, DMR, DME, KPT`.
- Because the `.ff` folder is generated for one built chain, `F` and `E` describe exactly that
  chain's ends: with a KPS head, `DMF`'s head CH2 bonds to `KPH`; with no head group, `DMF`'s head is
  a CH3.

## Charges from Context Oligomers <!-- rq-2dccba6a -->

Charges depend on chemical surroundings, so each block is parametrized inside a short **context
oligomer** built with `build_polymer` (same seed, coil conformation, no relaxation):

- **R** block of monomer M: three units `M M M`, no end groups. The middle residue is used.
- **F** block: the chain's first monomer M followed by two more units of M, with the chain's head
  group (if any) and no tail group. The first monomer residue is used, plus the head group's
  residue (e.g. `KPH`).
- **E** block: two units of the chain's last monomer M followed by M, with the chain's tail group
  (if any) and no head group. The last monomer residue is used, plus the tail group's residue
  (e.g. `KPT`).
- A one-unit chain is its own context.
- A context is parametrized once per build even if several blocks need it. Context residues are
  matched to chain residues by atom name; a mismatch is an internal error (`ForceFieldError`).
- **Charge normalization**: each residue's target charge is the sum of its atoms' formal charges
  (0 for a zwitterionic unit, −1 for `KPH`). The difference between the target and the backend's
  charges is spread equally over the residue's atoms. Charges are then rounded to 4 decimals with
  the largest-remainder method, so each written charge is within 0.0001 of its exact value and
  each block's written charges sum **exactly** to its target. The chain's total charge equals its net
  charge.

## User Charges (QM / RESP) for Selected Atoms <!-- rq-cb6ca230 -->

By default every charge comes from the backend (AM1-BCC for GAFF2, the OPLS-AA rules for
OPLS-AA). Users who have better charges for part of a residue, typically the ionic groups from a
QM calculation (HF/DFT ESP or RESP), can supply them in a **charge file**. This works the same way
for every backend.

### Charge file <!-- rq-31b6ec98 -->

A YAML mapping from residue name to a mapping of atom name → charge (e). A heavy atom may also
carry the charge of **its hydrogens** in one entry, `{charge: q, H: qH}`, which sets every
hydrogen bonded to that atom to `qH`. This matches how DFT/RESP charges for CH2/CH3 groups are
usually reported:

```yaml
# my_dft.yaml: DFT/RESP charges for the ionic groups of DMAPS
DMR:
  C6: {charge: -0.1976, H: 0.2042}   # N+-CH2 and its two H
  N1: 0.1573
  C7: {charge: -0.4018, H: 0.2042}   # N+-CH3 and its three H
  C8: {charge: -0.4018, H: 0.2042}
  S1: 1.4022
  O3: -0.7527
  O4: -0.7527
  O5: -0.7527
DMF:
  N1: 0.1573
```

- An entry's value is a number (the atom's charge), or a mapping with `charge` and/or `H`
  (either may be left out: `{H: 0.06}` changes only the hydrogens).
- Hydrogens may also be listed one by one by name (`H6: 0.2042`). Giving a hydrogen both ways with
  different values is an error.

- Residue names are the chain's residue names (`DMF`, `DMR`, `DME`, `KPH`, …). Atom names are the
  builder's names, which `JOB = charges` writes (see below).
- Only the listed atoms are changed; any subset of atoms, hydrogens included, may be given. A
  residue with no atoms listed (for example a template block left fully commented out) and an empty
  file mean "no user charges".
- A residue name that does not occur in the chain is ignored with a logged warning (the CLI prints
  it), so one file can serve several chains.
- These are errors (`ChargeFileError`, naming the file, residue and atom): an atom name that does
  not exist in that residue, a value that is not a finite number, a mapping value with keys other
  than `charge`/`H` or with neither, `H` on an atom that has no hydrogens, a hydrogen given two
  different values, a file that is not a YAML mapping of mappings, or a missing or unreadable file.

### Balancing <!-- rq-3541c593 -->

- Listed charges are kept **exactly**. The difference between the residue's integer target (see
  *Charges from Context Oligomers*) and the sum of all its charges is spread equally over the
  residue's **other** atoms. Then charges are rounded to 4 decimals; the rounding is absorbed by
  the unlisted atoms with the largest-remainder method, so each stays within 0.0001 of its exact
  value and the residue sums exactly to its target.
- If every atom of a residue is listed, their sum must equal the target within 0.0001, otherwise
  `ChargeFileError`.
- A **charge report** records, per residue with user charges: the number of user atoms, and the
  correction added to each other atom. If that correction is larger than 0.1 e in magnitude, a
  warning is logged and the report marks it, because it usually means the user charges were fitted
  for a different charge state or fragment.

### Charge template <!-- rq-194370d8 -->

`charge_template(structure, topology=None) -> str` returns a ready-to-edit charge file for the
structure. Each residue block of the chain is listed once, in chain order. Each **heavy atom** is
a commented line, so the user uncomments only the atoms they have charges for. A heavy atom with
hydrogens is written in the `{charge, H}` form, and its comment names those hydrogens:

```yaml
# DMR: repeat unit of DMAPS. Uncomment the atoms you have QM/RESP charges for.
DMR:
  # C6: {charge: -0.1203, H: 0.0600}   # C, 2 H (H6 H7), bonded to N1 C5
  # N1: 0.0821   # N, bonded to C6 C7 C8 C9
```

- The values shown are the backend's charges when a topology is given (for `H`, the mean of that
  atom's hydrogens), otherwise `0.0`.
- The file header explains that hydrogens can also be given one by one by name.
- Each comment gives the element and the bonded atoms (atoms in a neighbouring residue are written
  `-NAME` / `+NAME`), so the ionic atoms are easy to recognise.

## User Bonded Parameters <!-- rq-5a64d41e -->

The bonded parameters of the **ionic parts** (the cationic and anionic groups) are the ones that
cause problems: the force field may lack them (for example, the stock GROMACS `oplsaa.ff` has no
bonds for its sulfone S bond type `SY2`), or the user has better values from DFT. Everything else
comes from the force-field files as usual. A **parameter file** adds or **overrides** bonded
parameters for every backend, the same way a charge file overrides charges, and
`JOB = params` writes one for the ionic parts with the current values filled in
(see *Ionic parameter template*).

```yaml
# my_params.yaml (GROMACS units: nm, degrees, kJ/mol)
version: 1
bonds:
  - {types: [CT, SY2], same_as: [CT, SY]}          # reuse an existing parameter
  - {types: [SY2, OY], b0: 0.1440, kb: 585760.0}
angles:
  - {types: [OY, SY2, OY], theta0: 119.0, k: 870.272}
dihedrals:
  - {types: [X, CT, SY2, X], rb: [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]}   # OPLS-AA (RB, funct 3)
  - {types: [X, c3, s6, X], terms: [[0.0, 0.6508, 3]]}             # GAFF2 (Fourier, funct 9)
impropers:
  - {types: [X, X, C_2, X], phase: 180.0, k: 43.932, mult: 2}
```

- `types` are **bond types** (what `MissingParameterError` prints): GAFF2 type names, or OPLS-AA
  classes such as `CT`. Dihedral and improper keys may use `X`.
- Each entry either gives its values or `same_as: [types…]` (same length), which copies the
  parameter for that other key (wildcards and reverse order allowed, as in matching). The key is
  looked up in the **whole** force-field file (a backend's optional `full_parameters()`: all of
  `ffbonded.itp` for OPLS-AA, all of `gaff2.dat` for GAFF2), not only in the types the chain uses,
  so `same_as: [CT, SY]` works even when no `SY` atom is present.
- Values: bonds `b0` (nm), `kb` (kJ/mol/nm²); angles `theta0` (degrees), `k` (kJ/mol/rad²);
  impropers `phase` (degrees), `k` (kJ/mol), `mult`; dihedrals `rb` (six coefficients, for force
  fields with RB dihedrals) or `terms` (a list of `[phase, k, mult]`, for Fourier dihedrals).
- A user entry replaces a backend parameter with the same key (or its reverse), and otherwise adds
  a new one. User entries are applied after the contexts' parameters are merged and before the
  missing-parameter check, so they also fill gaps.
- Errors (`ParameterFileError`, naming the file and entry): unreadable file or invalid YAML, a
  missing or unsupported `version`, an unknown section or key, a wrong `types` length, missing or
  non-numeric values, both values and `same_as`, a `same_as` key with no parameter, or a dihedral
  in the other style than the force field's (`rb` for a Fourier force field or `terms` for an RB
  one).
- The topology records the keys that came from the user (`Topology.user_parameters`), and the
  generated `ffbonded.itp` marks those lines with `; user`.

### Ionic parameter template <!-- rq-702d4a58 -->

`ionic_parameter_template(topology) -> str` returns a parameter file for the ionic parts of the
chain:

- **Ionic atoms** are atoms with a non-zero formal charge plus the heavy atoms within two bonds of
  them. For SO₃⁻ that is the O⁻, the S, the other O and the C on S; for N⁺ its carbons and the next
  carbons along each chain; for a KPS end group the ester-like O as well.
- Listed terms: every bond, angle and proper dihedral that contains at least one ionic atom, each
  bond-type key once.
- Each line has the force field's current values, or zeros with the comment `# missing in the
  force field` when there is none. Each line's comment names one example of the atoms involved
  (residue and atom names), so the user can match it to their DFT data.
- Loading the template unchanged and applying it changes no value (only the missing entries need
  editing first).
- To build the template even when parameters are missing, `parametrize(..., strict=False)` returns
  the topology without raising `MissingParameterError`; the missing combinations are recorded in
  `Topology.missing_parameters`. Such a topology cannot be written (`write_gromacs` raises
  `MissingParameterError`).

## Backends <!-- rq-17848072 -->

- A backend has a `name` (`"gaff2"`, `"oplsaa"`), `check_available()` (raises
  `ForceFieldUnavailableError` with an installation hint), and
  `parametrize(mol, net_charge, work_dir) -> FragmentParameters`. `mol` is an RDKit molecule with
  explicit hydrogens and one conformer (Å).
- Backends are found in a registry: built-ins register themselves, `register_backend()` adds one
  at runtime, and installed packages can provide backends through the entry-point group
  `zpoly.forcefields` (each entry point loads a backend object). Names are
  case-insensitive (`ForceField` from `rqm/zwitterion_library.md` resolves them).
- **Merging** across contexts: atom types with identical definitions merge. An atom type defined
  twice with a different mass, σ or ε raises `ForceFieldError`. For bonded parameters with the same
  type key, the first context's value wins, and a differing later value is logged as a warning.

## Bonded Terms of the Chain <!-- rq-79e07d8e -->

- **Bonds**: every bond. **Angles**: every i–j–k path. **Proper dihedrals**: every i–j–k–l path
  (each central bond once). **Pairs**: the 1–4 atom pairs of the proper dihedrals, each once, and
  never a pair that is also 1–2 or 1–3 (rings).
- **Impropers**: for every atom with exactly three neighbours and SP2 hybridization (RDKit), one
  improper with the central atom third: (n1, n2, centre, n3), the neighbours ordered by atom type
  name then index. It is kept only if a matching improper parameter exists.
- **Parameter matching** uses each atom's **bond type**: the atom type's `bond_type`, or the type
  name itself when the backend sets none (GAFF2). OPLS-AA uses classes such as `CT`; dihedral and improper keys may contain the wildcard
  `X`, and a key also matches in reverse order. Every bond, angle and proper dihedral of the chain
  must have a matching parameter, otherwise `MissingParameterError` lists **all** missing type
  combinations (bonds, angles and dihedrals together). Its message ends with a ready-to-edit
  parameter-file block for exactly those combinations (see *User Bonded Parameters*). This
  guarantees that `gmx grompp` finds all parameters.

## GROMACS Output <!-- rq-39315c87 -->

`write_gromacs(topology, out_dir, name)` writes, using `name` for every file:

- `<name>.pdb`: the structure (`write_pdb`), with the F/R/E residue names.
- `<name>.ff/forcefield.itp`: a `[ defaults ]` section from the backend (nbfunc, comb-rule,
  gen-pairs, fudgeLJ, fudgeQQ), then `#include "ffnonbonded.itp"` and `#include "ffbonded.itp"`.
- `<name>.ff/forcefield.doc`: one line naming the force field and the generator.
- `<name>.ff/ffnonbonded.itp`: `[ atomtypes ]` with name, atomic number, mass, charge 0.0, ptype
  `A`, σ (nm) and ε (kJ/mol). When any atom type has a bond type different from its name, a
  bond-type column follows the name (GROMACS 8-column format).
- `<name>.ff/ffbonded.itp`: `[ bondtypes ]` (funct 1, b0 nm, kb kJ/mol/nm²), `[ angletypes ]`
  (funct 1, θ0 degrees, kθ kJ/mol/rad²), `[ dihedraltypes ]` for propers and impropers (funct 4). Propers use the
  parameter set's `dihedral_funct`: funct 9 writes one line per Fourier term (φs degrees, kφ
  kJ/mol, multiplicity); funct 3 writes the six RB coefficients on one line.
- `<name>.ff/<name>.rtp`: a `[ bondedtypes ]` header from the backend, then one block per residue
  with `[ atoms ]` (name, type, charge, charge-group number), `[ bonds ]` (intra-residue bonds plus
  the link to the previous residue written as `-<atom>`), and `[ impropers ]` when the residue has
  any.
- `<name>.top`: `#include "<name>.ff/forcefield.itp"`, a `[ moleculetype ]` (nrexcl 3) with
  `[ atoms ]` (nr, type, resnr, residue, atom, cgnr = residue number, charge, mass), `[ bonds ]`,
  `[ pairs ]`, `[ angles ]`, `[ dihedrals ]` (propers with the parameter set's `dihedral_funct`, then funct 4 impropers), all without
  explicit parameters (grompp looks them up), then `[ system ]` and `[ molecules ]` with one
  molecule. A comment states the net charge and the counterions needed to neutralize it (for
  example `; net charge -2: add 2 K+`).
- Output is deterministic: the same inputs give byte-identical files (no dates).
- Parent directories are created. If any step before writing fails, no file is written.

## Feature API <!-- rq-b45bedcc -->

All public names are importable from `zpoly.forcefield`.

### Types <!-- rq-2691f166 -->

- `AtomType` (frozen dataclass): `name`, `atomic_number: int`, `mass: float` (u), <!-- rq-a68ca1c4 -->
  `sigma: float` (nm), `epsilon: float` (kJ/mol), `bond_type: str = ""` (empty means the name).
- `ParameterSet` (frozen dataclass): `defaults: tuple[int, int, str, float, float]` <!-- rq-840650f8 -->
  (nbfunc, comb-rule, gen-pairs `"yes"`/`"no"`, fudgeLJ, fudgeQQ), `rtp_bondedtypes:
  tuple[int, ...]` (the 8 `[ bondedtypes ]` integers), `atomtypes: Mapping[str, AtomType]`,
  `bondtypes: Mapping[tuple[str, str], tuple[float, float]]`,
  `angletypes: Mapping[tuple[str, str, str], tuple[float, float]]`,
  `dihedraltypes: Mapping[tuple[str, str, str, str], tuple]`, `impropertypes:
  Mapping[tuple[str, str, str, str], tuple[float, float, int]]`, `dihedral_funct: int = 9`. With
  `dihedral_funct` 9 a dihedral value is a tuple of Fourier terms (φs, kφ, multiplicity); with 3 it
  is the six Ryckaert–Bellemans coefficients C0–C5 (kJ/mol). Bonded keys use bond types.
- `FragmentParameters` (frozen dataclass): `atom_types: tuple[str, ...]` and `charges: tuple[float, ...]` <!-- rq-81766163 -->
  (one per atom of the parametrized molecule, in its atom order), and `parameters: ParameterSet`.
- `ResidueBlock` (frozen dataclass): `name`, `atom_names`, `atom_types`, `charges` (tuples, block <!-- rq-9616f2d1 -->
  atom order), `bonds: tuple[tuple[str, str], ...]` (a link atom is written `-NAME`),
  `impropers: tuple[tuple[str, str, str, str], ...]`, `charge: int`.
- `Topology` (frozen dataclass): `name`, `force_field: str`, `structure: PolymerStructure`, <!-- rq-853c44c5 -->
  `blocks: tuple[ResidueBlock, ...]`, `atom_types: tuple[str, ...]` and `charges: tuple[float, ...]`
  (per atom of `structure.mol`), `parameters: ParameterSet`, `bonds`, `pairs`, `angles`,
  `propers`, `impropers` (tuples of 0-based atom-index tuples), `net_charge: int`,
  `charge_report: tuple[ChargeReport, ...]` (empty without user charges),
  `user_parameters: tuple[tuple[str, tuple[str, ...]], ...]` ((kind, key) of every user-supplied
  parameter; empty without a parameter file), `missing_parameters: Mapping[str, tuple]` (empty
  unless built with `strict=False`).
- `ChargeReport` (frozen dataclass): `residue: str`, `user_atoms: int`, `correction_per_atom: float`, <!-- rq-ea59ede3 -->
  `large: bool` (correction above 0.1 e).
- `ForceFieldBackend` (Protocol): `name: str`, `check_available() -> None`, <!-- rq-54f83eba -->
  `parametrize(mol, net_charge: int, work_dir: Path, *, charges: bool = True) -> FragmentParameters`.
  With `charges=False` (all charges come from DFT, `rqm/forcefield/dft_charges.md`) a backend may
  skip its charge model and return zero charges.

### Functions <!-- rq-d843e829 -->

- `parametrize(structure: PolymerStructure, force_field: str = "gaff2", *, name: str | None = None, <!-- rq-1229e2e7 -->
  work_dir: str | os.PathLike | None = None, charges: str | os.PathLike | Mapping | None = None,
  params: str | os.PathLike | Mapping | None = None, strict: bool = True) -> Topology`: runs the
  backend on the context
  oligomers and assembles the topology. `name` defaults to `"polymer"`. Intermediate files go to
  `work_dir` (default: a temporary directory that is removed afterwards). `charges` is a
  charge-file path or an already-loaded mapping (see *User Charges*); `params` is a parameter-file
  path or mapping (see *User Bonded Parameters*).
- `load_parameter_file(path: str | os.PathLike) -> dict`: reads and checks a parameter file <!-- rq-7245ccbb -->
  (`same_as` keys are resolved later, against the backend's parameters).
- `ionic_parameter_template(topology: Topology) -> str` (see *Ionic parameter template*). <!-- rq-2e3f1ce0 -->
- `load_charge_file(path: str | os.PathLike) -> dict[str, dict[str, float]]`: reads and checks the <!-- rq-2b7983f7 -->
  file format (not the atom names, which depend on the chain).
- `charge_template(structure: PolymerStructure, topology: Topology | None = None) -> str`. <!-- rq-e560bace -->
- `write_gromacs(topology: Topology, out_dir: str | os.PathLike, name: str | None = None, <!-- rq-58b0a34d -->
  padding: float = 10.0) -> Path`: writes the files above and returns the `.top` path. `name`
  defaults to `topology.name`; `padding` is passed to `write_pdb`.
- `get_backend(name: str) -> ForceFieldBackend`, `available_backends() -> list[str]`, <!-- rq-adb49473 -->
  `register_backend(backend: ForceFieldBackend) -> None`.

### Errors <!-- rq-ad30757b -->

- `ForceFieldError(Exception)`: base class. <!-- rq-2802a0d0 -->
  - `UnknownForceFieldError`: no backend has this name. Carries `name` and `available`.
  - `ForceFieldUnavailableError`: the backend's external tools are missing. Carries `name` and
    `hint` (how to install them).
  - `ParametrizationError`: an external tool failed. Carries `tool`, `returncode` and `log`
    (the tool's output, last 40 lines).
  - `ChargeFileError`: the charge file is unreadable or invalid, names an atom that does not
    exist, or fully specifies a residue with the wrong total. Carries `path` (or `None` for a
    mapping), `residue` and `atom` when known.
  - `MissingParameterError`: the chain has bonded terms without parameters. Carries `missing`
    (a mapping from `"bond"`, `"angle"`, `"dihedral"` to the missing type combinations), plus
    `kind` and `types` for the first kind that is missing. `template` is the ready-to-edit
    parameter-file text, also included at the end of the message.
  - `ParameterFileError`: the parameter file is invalid (see *User Bonded Parameters*). Carries
    `path` and `entry` (e.g. `"bonds[2]"`) when known.

## Input File <!-- rq-0e26d39f -->

In the input file (`rqm/input_file.md`), `FF = gaff2 | oplsaa | none` (default `gaff2`; `none`
writes only the PDB):

- With a force field, `JOB = build` parametrizes the chain before writing anything, then writes
  `<OUTPUT>.pdb`, `<OUTPUT>.top` and `<OUTPUT>.ff/` next to each other, and prints the summary line
  plus `force field <name>: <OUTPUT>.top, <OUTPUT>.ff/`.
- `PARAMS = FILE` applies a parameter file. The summary reports `user parameters: N`.
  `JOB = params` writes the ionic parameter template `<OUTPUT>_params.yaml` for the chain.
- `CHARGES = FILE` applies a charge file. Residues with user charges are reported one line each
  (`user charges: DMR 5 atoms, correction +0.0123 e/atom`, `(large)` above 0.1 e).
  `JOB = charges` writes the charge template `<OUTPUT>_charges.yaml` (use `N = 3` so the F, R and
  E blocks all appear); with a force field its default charges are filled in.
- Force-field errors, including charge-file errors, print `error: …` and exit 1, and no file is
  written.
- `zpoly list forcefields` prints each registered backend's name and whether it is
  available on this machine (`available` / `not installed: <hint>`).

## Out of Scope <!-- rq-07bf5734 -->

- The GAFF2 and OPLS-AA chemistry itself (their own requirement files).
- Solvent, ions and multi-chain systems; `.mdp` run files (the system builder).
- Validation by running GROMACS in this repository's CI (no GROMACS here; users run `gmx grompp`).

## Gherkin Scenarios <!-- rq-6d13caec -->

These scenarios use a deterministic test backend registered as `"fake"`: it types every atom by
element (`C`, `H`, …), gives every atom charge 0.1, and supplies parameters for every type
combination unless told otherwise.

```gherkin
Feature: Force-field topology core

  @rq-be5f3279
  Scenario: Residue blocks of a KPS-capped homopolymer
    Given a structure from build_polymer(homopolymer("DMAPS", 5), initiator="KPS")
    When parametrize(structure, "fake") is called
    Then the blocks are named "KPH", "DMF", "DMR", "DME", "KPT" in that order

  @rq-e2efc78d
  Scenario: Repeat residues are identical
    Given the topology of a 5-unit KPS-capped DMAPS chain
    Then every residue named "DMR" in the chain has the same atom names, types and charges as the "DMR" block

  @rq-998720ac
  Scenario: Block charges are normalized to integer targets
    Given the topology of a 5-unit KPS-capped DMAPS chain from the fake backend
    Then the charges of each DM* block sum exactly to 0 at 4 decimals
    And the charges of "KPH" and "KPT" each sum exactly to -1
    And the sum of all atom charges equals -2

  @rq-882e553a
  Scenario: Context oligomers
    Given a structure from build_polymer(homopolymer("DMAPS", 5), initiator="KPS")
    When parametrize(structure, "fake") is called
    Then the backend is called exactly 3 times, with net charges 0, -1 and -1

  @rq-71d21eb1
  Scenario: One-unit chain
    Given a structure from build_polymer(["DMAPS"])
    When parametrize(structure, "fake") is called
    Then the backend is called once and the only block is "DMF"

  @rq-c3c66893
  Scenario: Copolymer repeat blocks
    Given a structure from build_polymer(["DMAPS", "CBMA", "DMAPS", "CBMA"])
    When parametrize(structure, "fake") is called
    Then the blocks are "DMF", "CBR", "DMR", "CBE"

  @rq-5ddf514c
  Scenario: Conflicting atom type definitions
    Given a fake backend that defines type "C" with a different mass in the second context
    When parametrize is called on a 5-unit KPS-capped DMAPS chain
    Then ForceFieldError is raised naming type "C"

  @rq-26506f0a
  Scenario: Missing bonded parameter
    Given a fake backend that supplies no angle parameter for "C-S-O"
    When parametrize is called on a 3-unit DMAPS chain
    Then MissingParameterError is raised with kind "angle" listing ("C", "S", "O")

  @rq-37d656d6
  Scenario: Wildcard dihedral parameters match
    Given a fake backend whose only proper dihedral parameter is ("X", "C", "C", "X")
    And parameters for every other dihedral central pair
    When parametrize is called on a 3-unit DMAPS chain
    Then no MissingParameterError is raised

  @rq-ba87a649
  Scenario: Write GROMACS files
    Given the topology of a 3-unit KPS-capped DMAPS chain named "dmaps3"
    When write_gromacs(topology, tmp) is called
    Then tmp contains "dmaps3.pdb", "dmaps3.top" and "dmaps3.ff/" with forcefield.itp, forcefield.doc, ffnonbonded.itp, ffbonded.itp and dmaps3.rtp
    And forcefield.itp includes ffnonbonded.itp and ffbonded.itp after [ defaults ]
    And the .top includes "dmaps3.ff/forcefield.itp"
    And the .top [ atoms ] count equals the PDB ATOM count, with the same residue and atom names in the same order

  @rq-fecefe8b
  Scenario: RTP blocks and links
    Given the GROMACS files of a 3-unit KPS-capped DMAPS chain
    Then the .rtp has blocks [ KPH ], [ DMF ], [ DMR ], [ DME ], [ KPT ]
    And [ KPH ] has no bond to a previous residue
    And [ DMF ] has exactly one bond containing a "-" atom, to the KPH oxygen
    And every [ atoms ] line has a name, a type, a charge and a charge-group number

  @rq-e3e93468
  Scenario: Topology bonded sections
    Given the GROMACS files of a 3-unit DMAPS chain
    Then the .top [ bonds ] count equals the number of bonds of the molecule
    And [ pairs ] contains no pair that is also in [ bonds ] or 1-3 related
    And proper dihedrals use funct 9 and impropers use funct 4

  @rq-0e0b7ea5
  Scenario: Net charge comment
    Given the GROMACS files of a 3-unit KPS-capped DMAPS chain
    Then the .top contains "; net charge -2: add 2 K+"

  @rq-7348427c
  Scenario: Deterministic output
    When the same topology is written twice to different folders
    Then the corresponding files are byte-identical

  @rq-23d1bfb8
  Scenario: Unknown force field
    When parametrize(structure, "charmm") is called
    Then UnknownForceFieldError is raised listing the available backends

  @rq-d876e1b0
  Scenario: Backend registered at runtime
    When register_backend(a backend named "mine") is called
    Then "mine" is in available_backends() and get_backend("MINE") returns it

  @rq-ba2adece
  Scenario: CLI writes force-field files
    Given the fake backend is registered
    When an input "MONOMER = DMAPS / N = 3 / INITIATOR = KPS / FF = fake / OUTPUT = out/dm" is run
    Then it exits 0 and out/ contains dm.pdb, dm.top and dm.ff/
    And the output mentions "force field fake"

  @rq-07b1387b
  Scenario: CLI with an unavailable force field writes nothing
    Given a registered backend "absent" whose check_available raises ForceFieldUnavailableError
    When an input "MONOMER = DMAPS / N = 3 / FF = absent / OUTPUT = out/dm" is run
    Then it exits 1, stderr starts with "error:", and out/ contains no files

  @rq-dca9b6e9
  Scenario: CLI lists force fields
    When "zpoly list forcefields" is run
    Then the output has a line starting with "gaff2"

  @rq-56f6f132
  Scenario: User charges are kept exactly and the rest is balanced
    Given the fake backend and a charge file {"DMR": {"N1": 0.5, "S1": 1.2}}
    When a 5-unit DMAPS chain is parametrized with that charge file
    Then in the DMR block N1 is 0.5 and S1 is 1.2 exactly
    And the DMR charges sum exactly to 0
    And every other DMR atom changed by the same amount (within 0.0001)
    And the charge report for DMR has 2 user atoms

  @rq-d3724017
  Scenario: Same user charges in every repeat residue
    Given the fake backend and a charge file {"DMR": {"N1": 0.5}}
    When a 5-unit DMAPS chain is parametrized with that charge file
    Then atom N1 of every DMR residue in the topology has charge 0.5

  @rq-e3a90ed3
  Scenario: Residue not in the chain is ignored
    Given the fake backend and a charge file {"XXR": {"C1": 0.1}, "DMR": {"N1": 0.5}}
    When a 3-unit DMAPS chain is parametrized with that charge file
    Then no error is raised and N1 in DMR is 0.5

  @rq-ab5245a5
  Scenario: Unknown atom name
    Given a charge file {"DMR": {"N99": 0.5}}
    When a 3-unit DMAPS chain is parametrized with it
    Then ChargeFileError is raised with residue "DMR" and atom "N99"

  @rq-1ea0268f
  Scenario: Fully specified residue with the wrong total
    Given a charge file that gives every atom of KPH a charge, summing to -0.5
    When a 3-unit KPS-capped DMAPS chain is parametrized with it
    Then ChargeFileError is raised with residue "KPH"

  @rq-e4d525a8
  Scenario: Large correction is flagged
    Given the fake backend and a charge file {"DMR": {"N1": 1.0, "S1": 2.0}}
    When a 3-unit DMAPS chain is parametrized with that charge file
    Then the charge report for DMR is marked large

  @rq-58d9d01e
  Scenario: Invalid charge file
    Given a charge file containing "DMR: [1, 2]"
    When load_charge_file is called
    Then ChargeFileError is raised naming the file

  @rq-4700787b
  Scenario: Charge template lists every block and atom
    Given a structure from build_polymer(homopolymer("DMAPS", 3), initiator="KPS")
    When charge_template(structure) is called
    Then the text has the blocks KPH, DMF, DMR, DME, KPT in that order
    And every heavy atom of each block appears on a commented line with its element and bonded atoms
    And every hydrogen is named in its parent atom's line
    And after uncommenting every line it loads with load_charge_file

  @rq-83e95ea1
  Scenario: CLI charges option
    Given the fake backend and a charge file {"DMR": {"N1": 0.5}}
    When an input "MONOMER = DMAPS / N = 3 / FF = fake / CHARGES = FILE / OUTPUT = out/d" is run
    Then it exits 0, the output has a line "user charges: DMR 1 atoms", and in out/d.top every DMR atom N1 has charge 0.5000

  @rq-3618722f
  Scenario: CLI charges-template
    When an input "JOB = charges / MONOMER = DMAPS / N = 3 / INITIATOR = KPS / FF = none" is run
    Then it exits 0 and DMAPS_n3_KPS_charges.yaml has blocks for KPH, DMF, DMR, DME and KPT

  @rq-0bd185f8
  Scenario: Heavy atom and its hydrogens in one entry
    Given the fake backend and a charge file {"DMR": {"C7": {"charge": -0.4018, "H": 0.2042}, "N1": 0.1573}}
    When a 3-unit DMAPS chain is parametrized with that charge file
    Then in the DMR block C7 is -0.4018, N1 is 0.1573, and each of C7's three hydrogens is 0.2042
    And the charge report for DMR has 5 user atoms
    And the DMR charges sum exactly to 0

  @rq-27bf349e
  Scenario: Only the hydrogens of an atom
    Given the fake backend and a charge file {"DMR": {"C7": {"H": 0.1}}}
    When a 3-unit DMAPS chain is parametrized with that charge file
    Then C7's hydrogens are 0.1 and C7 itself was balanced like the other unlisted atoms

  @rq-c5f83c02
  Scenario: Hydrogen given two different values
    Given a charge file that sets C7's hydrogens with {"H": 0.2} and one of them by name to 0.3
    When a 3-unit DMAPS chain is parametrized with it
    Then ChargeFileError is raised naming that hydrogen

  @rq-d93f2f21
  Scenario: H on an atom without hydrogens
    Given a charge file {"DMR": {"N1": {"charge": 0.15, "H": 0.2}}}
    When a 3-unit DMAPS chain is parametrized with it
    Then ChargeFileError is raised with atom "N1"

  @rq-a2452483
  Scenario: Invalid compact entry
    Given a charge file {"DMR": {"C7": {"q": 0.1}}}
    When load_charge_file is called
    Then ChargeFileError is raised naming atom "C7"

  @rq-38038e4b
  Scenario: Missing parameters produce a ready-to-edit block
    Given a fake backend that supplies no angle parameter for "C-S-O" and no bond for "C-S"
    When parametrize is called on a 3-unit DMAPS chain
    Then MissingParameterError is raised whose missing has kinds "bond" and "angle"
    And its template loads with load_parameter_file after values are filled in

  @rq-3072b091
  Scenario: User parameters fill a gap
    Given a fake backend that supplies no angle parameter for "C-S-O"
    And a parameter file with angle [C, S, O] theta0 108.9 k 619.232
    When parametrize is called with that parameter file
    Then no error is raised and the angle (C, S, O) is (108.9, 619.232)
    And user_parameters contains ("angle", ("C", "S", "O"))

  @rq-96dc51b6
  Scenario: same_as copies an existing parameter
    Given a fake backend that supplies no bond for "C-S" but has bond "C-C" (0.15, 250000.0)
    And a parameter file with bond [C, S] same_as [C, C]
    When parametrize is called with that parameter file
    Then the bond (C, S) is (0.15, 250000.0)

  @rq-74842904
  Scenario: User parameters override the backend
    Given the fake backend and a parameter file with bond [C, C] b0 0.2 kb 1000.0
    When parametrize is called with that parameter file
    Then the bond (C, C) is (0.2, 1000.0)
    And the ffbonded.itp line for C C ends with "; user"

  @rq-ba3008a1
  Scenario: Wrong dihedral style
    Given the fake backend (Fourier dihedrals) and a parameter file with dihedral [X, C, C, X] rb [0, 0, 0, 0, 0, 0]
    When parametrize is called with that parameter file
    Then ParameterFileError is raised naming "dihedrals[1]"

  @rq-db609516
  Scenario: same_as to a missing key
    Given a parameter file with bond [C, S] same_as [Q, Q]
    When parametrize is called with that parameter file
    Then ParameterFileError is raised naming "bonds[1]"

  @rq-ff8c4f34
  Scenario: Invalid parameter file
    Given a parameter file with a bond entry that has types [C] and b0 0.1
    When load_parameter_file is called
    Then ParameterFileError is raised naming the file and "bonds[1]"

  @rq-0d2848e7
  Scenario: CLI parameter file
    Given the fake backend and a parameter file with bond [C, C] b0 0.2 kb 1000.0
    When an input "MONOMER = DMAPS / N = 3 / FF = fake / PARAMS = FILE / OUTPUT = out/d" is run
    Then it exits 0 and the output has a line "user parameters: 1"

  @rq-38b4e442
  Scenario: Ionic parameter template lists only the ionic parts
    Given the fake backend and the topology of a 3-unit DMAPS chain
    When ionic_parameter_template(topology) is called
    Then it lists the bonds (C, N) and (O, S) and the angle (O, S, O)
    And it does not list the angle (O, C, O) of the ester
    And loading it and parametrizing again with it changes no parameter value

  @rq-0df95ac2
  Scenario: Template marks missing ionic parameters
    Given a fake backend that supplies no bond for "C-S"
    When parametrize(structure, "nobondcs", strict=False) and ionic_parameter_template are called
    Then the template lists bond [C, S] with zeros and "# missing in the force field"
    And write_gromacs on that topology raises MissingParameterError

  @rq-d76b560f
  Scenario: CLI params-template
    Given the fake backend is registered
    When an input "JOB = params / MONOMER = DMAPS / N = 3 / FF = fake" is run (written to DMAPS_n3_params.yaml)
    Then it exits 0 and prints a parameter file with "bonds:" and "angles:"
```
