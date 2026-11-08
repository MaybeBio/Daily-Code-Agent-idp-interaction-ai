# Feature: GAFF2 Force-Field Backend (AmberTools) <!-- rq-85e6069f -->

The GAFF2 backend for the force-field core (`rqm/forcefield/core.md`). It assigns GAFF2 atom types
and AM1-BCC charges with AmberTools, fills in missing parameters with `parmchk2`, and converts
everything to GROMACS units. Users select it with `FF = gaff2` or `parametrize(structure, "gaff2")`.

## Dependencies <!-- rq-c7ea04ff -->

- **AmberTools** (`antechamber`, `parmchk2`, and the `gaff2.dat` parameter file) is an **optional**
  external dependency, installed with `conda install -c conda-forge ambertools`. The package itself
  does not depend on it; the backend reports clearly when it is missing.

## Locating AmberTools <!-- rq-8c905ad1 -->

- `antechamber` and `parmchk2` are looked up on `PATH`.
- `gaff2.dat` is `$AMBERHOME/dat/leap/parm/gaff2.dat` when `AMBERHOME` is set. Otherwise it is
  derived from the `antechamber` location: `<prefix>/bin/antechamber` → `<prefix>/dat/leap/parm/gaff2.dat`.
- `check_available()` raises `ForceFieldUnavailableError` naming the first missing piece. Its hint is
  `conda install -c conda-forge ambertools` (and, for a missing `gaff2.dat`, to set `AMBERHOME`).

## Parametrizing One Context Oligomer <!-- rq-1b967d27 -->

`parametrize(mol, net_charge, work_dir)`:

1. Writes `frag.sdf` (RDKit MDL V2000 with formal charges and explicit hydrogens, in the molecule's
   atom order).
2. Runs `antechamber -i frag.sdf -fi mdl -o frag.mol2 -fo mol2 -c bcc -nc <net_charge> -at gaff2
   -pf y -s 0` in `work_dir`. AM1-BCC charges come from `sqm`, which antechamber runs.
3. Runs `parmchk2 -i frag.mol2 -f mol2 -o frag.frcmod -s gaff2`.
4. Reads atom types (column 6) and charges (column 9) of `frag.mol2`. The atom count must equal the
   input's, and the elements must match atom by atom, otherwise `ParametrizationError`.
5. Builds the `ParameterSet` from `gaff2.dat` plus `frag.frcmod`, keeping only parameters whose
   atom types all occur in the fragment (wildcard `X` always qualifies). A `frcmod` entry replaces
   a `gaff2.dat` entry with the same key.

- Each tool runs with a timeout (default 3600 s, configurable on the backend). A non-zero exit
  status, a timeout, or a missing output file raises `ParametrizationError` with the tool name,
  return code (−1 for a timeout) and the last 40 lines of its output.
- `[ defaults ]` is `1 2 yes 0.5 0.8333` (Amber combination rule and 1–4 scaling), and the `.rtp`
  `[ bondedtypes ]` header is `1 1 9 4 1 3 1 0`.

## Reading Amber Parameter Files <!-- rq-e2b2a86f -->

`gaff2.dat` sections, in order: a title line, `MASS` lines up to a blank line, one line of
hydrophilic types, then `BOND`, `ANGLE`, `DIHE`, `IMPROPER` and 10–12 H-bond sections (each up to a
blank line), equivalence lines up to a blank line, a `MOD4 RE` header, `NONBON` lines up to a blank
line, and `END`. `frcmod` files have keyword sections `MASS`, `BOND`, `ANGLE`, `DIHE`, `IMPROPER`,
`NONBON`, each up to a blank line. Type keys are the hyphen-separated, whitespace-stripped fields
(`c3-c3`, `X -c -c -X`).

Unit conversions (1 kcal = 4.184 kJ, 1 Å = 0.1 nm):

| Amber term | Amber fields | GROMACS |
|---|---|---|
| Bond | k (kcal/mol/Å²), r0 (Å) | b0 = r0 / 10 nm; kb = 2 · k · 418.4 kJ/mol/nm² |
| Angle | k (kcal/mol/rad²), θ0 (°) | θ0 °; kθ = 2 · k · 4.184 kJ/mol/rad² |
| Proper dihedral | IDIVF, PK (kcal/mol), phase (°), PN | funct 9 term: φs = phase, kφ = PK / IDIVF · 4.184, multiplicity = abs(PN). A negative PN means more terms for the same key follow on the next lines. |
| Improper | PK, phase, PN | funct 4: φs = phase, kφ = PK · 4.184, multiplicity = abs(PN) |
| Non-bonded | R* (Å), ε (kcal/mol) | σ = R* · 2^(5/6) / 10 nm; ε · 4.184 kJ/mol |
| Mass | mass (u) | mass (u) |

The atomic number of a type comes from the fragment atoms that carry it.

## Feature API <!-- rq-505753b3 -->

Importable from `zpoly.forcefield`:

- `Gaff2Backend` (implements `ForceFieldBackend`, registered as `"gaff2"`): constructor arguments <!-- rq-455dff93 -->
  `timeout: float = 3600.0`, `antechamber: str = "antechamber"`, `parmchk2: str = "parmchk2"`,
  `gaff2_dat: str | os.PathLike | None = None` (overrides the lookup).
- `read_amber_parameters(path, *, frcmod: bool = False) -> ParameterSet`: parses a `gaff2.dat` <!-- rq-2702f408 -->
  or `frcmod` file into GROMACS units (the `defaults` and `rtp_bondedtypes` are the Amber ones).

## Out of Scope <!-- rq-c613f8ae -->

- RESP or other QM charges; GAFF (version 1); Amber `prmtop` output.
- Running the real AmberTools in this repository's tests (it is not installed here). The tests use
  small fake `antechamber`/`parmchk2` programs and a fixture `gaff2.dat`; a test marked `amber`
  runs the real tools when they are on `PATH` and is skipped otherwise.

## Gherkin Scenarios <!-- rq-a54a6990 -->

```gherkin
Feature: GAFF2 backend

  @rq-8d5a0d31
  Scenario: AmberTools missing
    Given antechamber is not on PATH
    When Gaff2Backend().check_available() is called
    Then ForceFieldUnavailableError is raised whose hint contains "conda install -c conda-forge ambertools"

  @rq-5392e74f
  Scenario: gaff2.dat located from the antechamber prefix
    Given AMBERHOME is unset and <prefix>/bin/antechamber and <prefix>/dat/leap/parm/gaff2.dat exist
    When Gaff2Backend().check_available() is called
    Then no error is raised and the backend uses that gaff2.dat

  @rq-6c330b7c
  Scenario: gaff2.dat missing
    Given antechamber and parmchk2 are on PATH but no gaff2.dat can be found
    When Gaff2Backend().check_available() is called
    Then ForceFieldUnavailableError is raised whose hint mentions AMBERHOME

  @rq-bd8eb0c7
  Scenario: Tool command lines
    Given fake antechamber and parmchk2 programs that record their arguments
    When the backend parametrizes a fragment with net charge -1
    Then antechamber was called with "-c bcc", "-nc -1", "-at gaff2" and "-fi mdl"
    And parmchk2 was called with "-s gaff2"

  @rq-9b94cbfd
  Scenario: Types and charges from the mol2 file
    Given a fixture mol2 for a 5-atom fragment
    When the backend reads it for a matching 5-atom molecule
    Then the atom types and charges are returned in input atom order

  @rq-b018fb93
  Scenario: Mol2 atom mismatch
    Given a fixture mol2 whose second atom is an oxygen while the input's is a carbon
    When the backend reads it
    Then ParametrizationError is raised

  @rq-f5a523d0
  Scenario: Unit conversion
    Given a parameter file with bond "c3-c3 300.9 1.5375", angle "c3-c3-c3 62.9 111.51", nonbonded "c3 1.9080 0.1094" and mass "c3 12.01"
    When read_amber_parameters is called
    Then bond (c3, c3) is b0 0.15375 nm and kb 251793.12 kJ/mol/nm²
    And angle (c3, c3, c3) is θ0 111.51 and kθ 526.3472
    And atom type c3 has σ 0.339967 nm (to 6 decimals) and ε 0.4577296 kJ/mol and mass 12.01

  @rq-13487d5c
  Scenario: Multi-term dihedral
    Given DIHE lines "X -c -c -X 1 0.5 180.0 -2." and "X -c -c -X 1 0.3 0.0 3."
    When read_amber_parameters is called
    Then dihedral (X, c, c, X) has terms (180.0, 2.092, 2) and (0.0, 1.2552, 3)

  @rq-ed920f0f
  Scenario: Improper conversion
    Given IMPROPER line "c -o -o -o 1.1 180. 2."
    When read_amber_parameters is called
    Then improper (c, o, o, o) is (180.0, 4.6024, 2)

  @rq-2708fa8a
  Scenario: frcmod overrides gaff2.dat
    Given gaff2.dat and a frcmod that both define bond c3-n4 with different values
    When the backend builds the parameter set
    Then the frcmod value is used

  @rq-8aa9cd8f
  Scenario: Only the fragment's types are kept
    Given a gaff2.dat defining types c3, hc and os
    And a fragment typed only c3 and hc
    When the backend builds the parameter set
    Then atom type os is absent and every kept bonded key uses only c3, hc or X

  @rq-d50c2e10
  Scenario: Tool failure
    Given a fake antechamber that prints "sqm failed" and exits 1
    When the backend parametrizes a fragment
    Then ParametrizationError is raised with tool "antechamber", return code 1 and a log containing "sqm failed"

  @rq-d03cde6e
  Scenario: Tool timeout
    Given a fake antechamber that never finishes and a backend timeout of 1 s
    When the backend parametrizes a fragment
    Then ParametrizationError is raised with return code -1

  @rq-17e6c2e5
  Scenario: End-to-end with fake AmberTools
    Given fake antechamber and parmchk2 programs that type atoms by element and write a complete frcmod
    When parametrize(build_polymer(homopolymer("DMAPS", 3), initiator="KPS"), "gaff2") and write_gromacs are called
    Then the .ff folder, .rtp and .top are written and every DM* block charge sums to 0

  @rq-bf118372
  Scenario: Real AmberTools (runs only when installed)
    Given antechamber is on PATH
    When a 3-unit KPS-capped DMAPS chain is parametrized with "gaff2"
    Then the .top is written and its total charge is -2 within 1e-4
```
