# Feature: OPLS-AA Force-Field Backend (SMARTS Typing on Your oplsaa.ff) <!-- rq-f6d54673 -->

The OPLS-AA backend for the force-field core (`rqm/forcefield/core.md`). It works like the GAFF2
backend from the user's point of view: `FF = oplsaa` or `parametrize(structure, "oplsaa")` produces
the same F/R/E `.rtp`, `.ff` folder and ready `.top`, and the same charge files (QM/RESP) apply.

Two inputs drive it:

1. **An `oplsaa.ff` folder** with the OPLS-AA data: the user's own folder (for example a
   group-extended one that defines `opls_103`, `opls_493`, `opls_494`), or the stock GROMACS one.
   Nothing from it is copied into this package.
2. **Typing rules**: a YAML table mapping chemical groups (SMARTS) to `opls_XXX` types. A starter
   table ships with the package; users extend or override it with their own file.

## Locating oplsaa.ff <!-- rq-040346db -->

In order: the backend's `oplsaa_ff` argument (input key `OPLSAA_FF`); the environment variable
`ZPOLY_OPLSAA_FF`; `$GMXDATA/top/oplsaa.ff`; and `<prefix>/share/gromacs/top/oplsaa.ff` next to
a `gmx` (or `gmx_mpi`) found on `PATH`. The folder must contain `forcefield.itp`,
`ffnonbonded.itp` and `ffbonded.itp`. If none is found, `check_available()` raises
`ForceFieldUnavailableError` with the hint "set OPLSAA_FF = PATH to your oplsaa.ff folder, or
install GROMACS".

## Reading the GROMACS Files <!-- rq-fb24e397 -->

- `forcefield.itp`: the `[ defaults ]` line (nbfunc, comb-rule, gen-pairs, fudgeLJ, fudgeQQ).
- `ffnonbonded.itp`: `[ atomtypes ]` with 8 columns: name, bond type, atomic number, mass, charge,
  ptype, σ, ε.
- `ffbonded.itp`: `[ bondtypes ]` (funct 1: b0, kb), `[ angletypes ]` (funct 1: θ0, kθ),
  `[ dihedraltypes ]` with funct 3 (Ryckaert–Bellemans, 6 coefficients C0–C5; keys may use `X`),
  and `#define` macros (`#define improper_O_C_X_Y 180.0 43.93200 2`: φs, kφ, multiplicity).
  Other functions and sections are ignored.
- `atomtypes.atp` (optional): `name mass ; description` lines. The descriptions are used to check
  rules that declare `expect` keywords (see *Typing Rules*).
- Comments after `;` and blank lines are skipped, and `#include` lines are ignored.
  `#ifdef NAME`, `#ifndef NAME`, `#else` and `#endif` are evaluated like the GROMACS preprocessor
  with no `-D` options: a name counts as defined only if an earlier `#define NAME` in the files
  read so far defines it. So optional variants such as `#ifdef HEAVY_H` are skipped, and their
  `#else` branches are read.
- Bonded parameters are keyed by **bond type** (the second `[ atomtypes ]` column), as in GROMACS.
- A file that cannot be read or parsed raises `ForceFieldError` naming the file and line.

## Typing Rules <!-- rq-87cd66f0 -->

```yaml
version: 1
rules:                          # tried in order; the first match types an atom
  - name: ester alkoxy CH2
    smarts: "[CX4;H2,H1][OX2][CX3]=O"   # the FIRST atom of the pattern gets the type
    type: opls_490
  - name: carboxylate C
    smarts: "[CX3](=O)[O-]"
    type: opls_271
    expect: [carboxyl, COO]     # optional: checked against atomtypes.atp
impropers:
  - name: carbonyl centre
    smarts: "[CX3](=O)"         # the first atom is the improper's centre
    define: improper_O_C_X_Y    # a #define in ffbonded.itp
```

- Matching is done on the fragment with explicit hydrogens. Each atom gets the type of the
  **first** rule whose pattern matches with that atom as the pattern's first atom.
- User rules (`rules=` / `OPLS_RULES = FILE`) are tried **before** the starter rules, so they
  override them. User impropers are likewise checked before the starter ones.
- **Charges**: by default each atom gets its type's charge from `ffnonbonded.itp`, then the core
  balances each residue to its integer charge as usual. A charge file overrides selected atoms, as
  for every backend.
- Every rule's type must exist in the loaded `ffnonbonded.itp`, and every improper's `define` in
  `ffbonded.itp`, otherwise `TypingError` names the rule.
- An atom that no rule types raises `TypingError`, listing each untyped atom with its element, its
  neighbours and the fragment SMILES, so the user knows which rule to add.
- `expect` (optional) is a list of keywords. When the `oplsaa.ff` folder has `atomtypes.atp`,
  the rule's type description must contain at least one keyword (case-insensitive), otherwise
  `TypingError` names the rule and shows the folder's description, so a wrong type number is caught
  before anything is built. Without `atomtypes.atp`, the check is skipped with one logged warning.
- An unknown key, a missing `smarts`/`type`, an unparsable SMARTS, an `expect` that is not a list of
  strings, or an unsupported `version` raises `TypingError` naming the file and rule.

### Starter rules <!-- rq-11b0641d -->

Two kinds of starter rules ship. **Confirmed** rules use the types in the group's polybuilder
OPLS-AA residues (DMAPS/A3316/KPS and the catalog). **Standard** rules use standard OPLS-AA type
numbers and declare `expect` keywords, so they are checked against the user's `atomtypes.atp`.
Order matters (specific first):

| Group | SMARTS (first atom typed) | Type | Kind (`expect`) |
|---|---|---|---|
| Sulfonate / sulfate S | `[SX4](https://github.com/hussen-oumer/Zpoly/blob/HEAD/rqm/forcefield/%3DO)(=O)[O-]` | opls_493 | confirmed |
| Sulfonate / sulfate O (incl. O⁻) | `[OX1]~[SX4]` | opls_494 | confirmed |
| Sulfate ester O (KPS attachment) | `[OX2](https://github.com/hussen-oumer/Zpoly/blob/HEAD/rqm/forcefield/%5B#6])[SX4]` | opls_467 | confirmed |
| Phosphate diester P | `[PX4](https://github.com/hussen-oumer/Zpoly/blob/HEAD/rqm/forcefield/%3DO)([O-])([OX2])[OX2]` | opls_440 | standard (`PO4`, `phosph`) |
| Phosphate terminal O (O=, O⁻) | `[OX1]~[PX4]` | opls_441 | standard (`PO4`, `phosph`) |
| Phosphate ester O | `[OX2](https://github.com/hussen-oumer/Zpoly/blob/HEAD/rqm/forcefield/%5B#6])[PX4]` | opls_442 | standard (`PO4`, `phosph`) |
| Quaternary ammonium N⁺ | `[NX4+;H0]` | opls_103 | confirmed |
| Carboxylate C | `[CX3](https://github.com/hussen-oumer/Zpoly/blob/HEAD/rqm/forcefield/%3DO)[O-]` | opls_271 | standard (`carboxyl`, `COO`) |
| Carboxylate O | `[OX1]~[CX3]~[OX1]` | opls_272 | standard (`carboxyl`, `COO`) |
| Ester carbonyl C | `[CX3](https://github.com/hussen-oumer/Zpoly/blob/HEAD/rqm/forcefield/%3DO)[OX2][#6]` | opls_465 | confirmed |
| Ester C=O oxygen | `[OX1]=[CX3][OX2][#6]` | opls_466 | confirmed |
| Ester –O– | `[OX2](https://github.com/hussen-oumer/Zpoly/blob/HEAD/rqm/forcefield/%5BCX3%5D%3DO)[#6]` | opls_467 | confirmed |
| Ester alkoxy C (O–CH2 / O–CH) | `[CX4;H2,H1][OX2][CX3]=O` | opls_490 | confirmed |
| Amide C | `[CX3](https://github.com/hussen-oumer/Zpoly/blob/HEAD/rqm/forcefield/%3DO)[NX3]` | opls_235 | confirmed |
| Amide O | `[OX1]=[CX3][NX3]` | opls_236 | confirmed |
| Amide N–H | `[NX3;H1][CX3]=O` | opls_238 | confirmed |
| H on amide N | `[#1][NX3][CX3]=O` | opls_241 | confirmed |
| Alcohol O | `[OX2;H1][CX4]` | opls_154 | standard (`alcohol`) |
| H on alcohol O | `[#1][OX2;H1][CX4]` | opls_155 | standard (`alcohol`) |
| Dialkyl ether O | `[OX2](https://github.com/hussen-oumer/Zpoly/blob/HEAD/rqm/forcefield/%5BCX4%5D)[CX4]` | opls_180 | standard (`ether`) |
| Benzene-ring C (all-carbon 6-ring) | `[c;$(c1ccccc1)]` | opls_145 | standard (`benzene`, `aromatic`) |
| H on benzene-ring C | `[#1][c;$(c1ccccc1)]` | opls_146 | standard (`benzene`, `aromatic`) |
| Alkyl CH3 | `[CX4;H3]` | opls_135 | confirmed |
| Alkyl CH2 | `[CX4;H2]` | opls_136 | confirmed |
| Alkyl CH | `[CX4;H1]` | opls_137 | confirmed |
| Quaternary C | `[CX4;H0]` | opls_139 | confirmed |
| H on carbon | `[#1][#6]` | opls_140 | confirmed |

Starter impropers: `[CX3](https://github.com/hussen-oumer/Zpoly/blob/HEAD/rqm/forcefield/%3DO)` → `improper_O_C_X_Y` (ester, amide and carboxylate C);
`[NX3;H1][CX3]=O` → `improper_Z_N_X_Y`; `[c;$(c1ccccc1)]` → `improper_Z_CA_X_Y`.

Heteroaromatic cations (pyridinium, imidazolium, isoquinolinium), protonated amines (N–H⁺),
fluorinated groups, the amine N-oxide, phosphonium, sulfonium, phosphonates and sulfonylimides have
no starter rule; users add them in their rule file.

### Impropers <!-- rq-57822d06 -->

For each atom that an improper rule matches as centre, the backend supplies an improper parameter
from the rule's `define`, keyed `(X, X, <centre bond type>, X)`. The core then places the improper
on that atom (centre third) as for every backend. If two rules give different values for the same
centre bond type, the first one is kept.

## Output Conventions <!-- rq-32f43f89 -->

- `[ defaults ]` comes from the `oplsaa.ff` `forcefield.itp` (normally `1 3 yes 0.5 0.5`).
- The `.rtp` `[ bondedtypes ]` header is `1 1 3 4 1 3 1 0`: RB proper dihedrals (funct 3) and
  periodic impropers (funct 4, written into the generated `.ff`'s `[ dihedraltypes ]`).
- `ffnonbonded.itp` in the generated `.ff` keeps the bond-type column.

## Feature API <!-- rq-7e23c2b9 -->

Importable from `zpoly.forcefield`:

- `OplsaaBackend(oplsaa_ff: str | os.PathLike | None = None, rules: str | os.PathLike | None = None)`: <!-- rq-a47637a7 -->
  implements `ForceFieldBackend`, registered as `"oplsaa"`. `check_available()` locates and parses
  the folder and checks the rules against it.
- `read_gromacs_forcefield(path: str | os.PathLike) -> GromacsForceField`: parses an `.ff` folder. <!-- rq-9a642c3a -->
  `GromacsForceField` holds `defaults`, `atomtypes` (name → (bond type, atomic number, mass,
  charge, σ, ε)), `bondtypes`, `angletypes`, `rb_dihedraltypes`, `defines` and `descriptions`
  (name → `atomtypes.atp` description; empty when the file is absent).
- `load_typing_rules(path: str | os.PathLike | None = None) -> TypingRules`: the starter rules <!-- rq-779e8915 -->
  (`None`) or a user file.
- `TypingError(ForceFieldError)`: carries `rule` (name or `None`) and `atoms` (untyped atom <!-- rq-7705a23e -->
  descriptions).

## Input File <!-- rq-58130a41 -->

`FF = oplsaa`, with optional `OPLSAA_FF = PATH` and `OPLS_RULES = FILE`, configures the `oplsaa`
backend for that run (`rqm/input_file.md`). `OPLSAA_FF` / `OPLS_RULES` without `FF = oplsaa` is an
input error.

## Out of Scope <!-- rq-1728b949 -->

- Downloading or shipping OPLS-AA data; LigParGen or CM1A charges.
- Starter rules for chemistry not confirmed by the group's files (see *Starter rules*).

## Gherkin Scenarios <!-- rq-e15a0d6a -->

These scenarios use a small fixture `oplsaa.ff` folder written by the tests. It defines every
starter type with a bond type, the bond, angle and RB dihedral parameters needed (some via `X`
wildcards), and the two improper macros. A test marked `gromacs` uses a real `oplsaa.ff` when one is
found and is skipped otherwise.

```gherkin
Feature: OPLS-AA backend

  @rq-c593d1ee
  Scenario: oplsaa.ff not found
    Given no OPLSAA_FF, no ZPOLY_OPLSAA_FF, no GMXDATA and no gmx on PATH
    When OplsaaBackend().check_available() is called
    Then ForceFieldUnavailableError is raised whose hint mentions "OPLSAA_FF"

  @rq-a641bad1
  Scenario: oplsaa.ff from the environment
    Given ZPOLY_OPLSAA_FF points to the fixture folder
    When OplsaaBackend().check_available() is called
    Then no error is raised

  @rq-de0e6dd9
  Scenario: Reading the GROMACS files
    When read_gromacs_forcefield(fixture folder) is called
    Then opls_135 has bond type "CT", atomic number 6, mass 12.011 and charge -0.18
    And the bond type ("CT", "CT") has b0 0.1529 and kb 224262.4
    And the RB dihedral ("X", "CT", "CT", "X") has 6 coefficients
    And the define "improper_O_C_X_Y" is (180.0, 43.932, 2)
    And the defaults are (1, 3, "yes", 0.5, 0.5)

  @rq-daf83e1c
  Scenario: Preprocessor lines and comments are skipped
    Given an ffbonded.itp with #ifdef/#endif lines and trailing ";" comments
    When it is read
    Then the parameters outside the skipped #ifdef block are read and no error is raised

  @rq-9dd1d715
  Scenario: Conditional atom type variants
    Given an ffnonbonded.itp defining opls_140 with mass 4.032 under "#ifdef HEAVY_H" and 1.008 under "#else"
    When it is read
    Then opls_140 has mass 1.008

  @rq-3b60921d
  Scenario: Typing DMAPS with the starter rules
    Given the fixture folder
    When the backend parametrizes a 3-unit DMAPS chain fragment
    Then the ammonium N is opls_103, the sulfonate S is opls_493 and its three O are opls_494
    And the ester C, =O, -O- and O-CH2 are opls_465, opls_466, opls_467 and opls_490
    And the backbone quaternary C is opls_139 and every H on carbon is opls_140

  @rq-7fac93b5
  Scenario: Typing an acrylamide sulfobetaine
    When the backend parametrizes an A3361 (SBAA) fragment
    Then the amide C, O, N and N-H are opls_235, opls_236, opls_238 and opls_241
    And the backbone CH is opls_137

  @rq-cffa5145
  Scenario: KPS end groups
    When the backend parametrizes a KPS-capped DMAPS fragment
    Then the sulfate S is opls_493, its O and O- are opls_494, and the attachment O is opls_467

  @rq-08faba23
  Scenario: Default charges are the type charges
    When the backend parametrizes a fragment
    Then each atom's charge equals its type's charge in ffnonbonded.itp

  @rq-e02bf89b
  Scenario: Untyped atoms are reported
    When the backend parametrizes an SPVI (vinylimidazolium) fragment with the starter rules
    Then TypingError is raised listing the aromatic atoms with their elements and neighbours

  @rq-ea516f0e
  Scenario: User rules override the starter rules
    Given a user rule file typing "[CX4;H3][NX4+]" as opls_136
    When the backend parametrizes a DMAPS fragment with those rules
    Then the N-methyl carbons are opls_136 and the backbone methyl is still opls_135

  @rq-c7ee0327
  Scenario: Rule with a type missing from oplsaa.ff
    Given a user rule file typing "[CX4;H3]" as opls_999
    When check_available() is called
    Then TypingError is raised naming that rule

  @rq-a2b7965f
  Scenario: Invalid rule file
    Given a user rule file with a rule that has no "type"
    When load_typing_rules is called
    Then TypingError is raised naming the file

  @rq-91413f6c
  Scenario: Impropers from defines
    When the backend parametrizes a DMAPS fragment
    Then its impropertypes contain ("X", "X", <ester C bond type>, "X") = (180.0, 43.932, 2)

  @rq-010b9f97
  Scenario: RB dihedrals in the output
    Given a 3-unit DMAPS chain parametrized with "oplsaa" and written with write_gromacs
    Then ffbonded.itp has [ dihedraltypes ] lines with funct 3 and six coefficients
    And the .top proper dihedrals use funct 3
    And the .rtp [ bondedtypes ] line is "1 1 3 4 1 3 1 0"
    And ffnonbonded.itp lines have a bond-type column

  @rq-34b90470
  Scenario: Bonded parameters are matched by bond type
    Given the fixture defines bond parameters only for bond-type pairs
    When a 3-unit DMAPS chain is parametrized with "oplsaa"
    Then no MissingParameterError is raised

  @rq-b2b4afa0
  Scenario: Charge file works with OPLS-AA
    Given a charge file {"DMR": {"N1": 0.1573}}
    When a 3-unit DMAPS chain is parametrized with "oplsaa" and that charge file
    Then N1 of every DMR residue has charge 0.1573 and each DM* block sums to 0

  @rq-0e9e4781
  Scenario: CLI with OPLS-AA
    When an input "MONOMER = DMAPS / N = 3 / INITIATOR = KPS / FF = oplsaa / OPLSAA_FF = FIXTURE / OUTPUT = out/d" is run
    Then it exits 0 and writes out/d.top and out/d.ff/

  @rq-db5f1df5
  Scenario: OPLSAA_FF without FF = oplsaa
    When an input "MONOMER = DMAPS / N = 3 / OPLSAA_FF = FIXTURE" is run
    Then it exits 1 with an error saying it needs FF = oplsaa

  @rq-f129ffeb
  Scenario: Real oplsaa.ff (runs only when found)
    Given a real oplsaa.ff folder is found
    When read_gromacs_forcefield is called on it
    Then it has at least 800 atom types and opls_135 has bond type "CT"

  @rq-e2e215ea
  Scenario: Typing a carboxybetaine
    When the backend parametrizes a CBMA fragment
    Then the carboxylate C is opls_271 and both carboxylate O are opls_272

  @rq-87d16db3
  Scenario: Typing alcohol and ether oxygens
    When the backend parametrizes an HSBMA fragment and a ZB4 (morpholinium) fragment
    Then the hydroxyl O and H are opls_154 and opls_155
    And the morpholine ring O is opls_180

  @rq-8f9bad27
  Scenario: Typing a phosphorylcholine
    When the backend parametrizes an MPC fragment
    Then the P is opls_440, its two terminal O are opls_441 and its two ester O are opls_442

  @rq-c823b8ce
  Scenario: Typing a benzene ring
    When the backend parametrizes an SBVB (styrenic) fragment
    Then the six ring C are opls_145 and the four ring H are opls_146

  @rq-3d69e50c
  Scenario: Expected description does not match atomtypes.atp
    Given the fixture atomtypes.atp describes opls_271 as "alkane CH3"
    When check_available() is called
    Then TypingError is raised naming the "carboxylate C" rule and showing "alkane CH3"

  @rq-2106839f
  Scenario: No atomtypes.atp
    Given the fixture folder has no atomtypes.atp
    When check_available() is called
    Then no error is raised and a warning says the expect checks were skipped
```
