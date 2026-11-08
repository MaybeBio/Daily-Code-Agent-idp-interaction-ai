# Feature: Zwitterion Monomer Library <!-- rq-494a5e56 -->

The zwitterion library is the chemistry layer of `zpoly` (see `rqm/ARCHITECTURE.md`). It
turns the zwitterion building-block catalog into data the software can use:

- a **built-in library** of named zwitterionic monomers and small-molecule zwitterions
  (zwitterionic liquids, osmolytes), looked up by name or alias;
- the **building blocks** (polymerizable groups, linkers, cations, spacers, anions) those monomers
  are made from;
- **composition** of new monomers from blocks, and **import** of arbitrary zwitterion SMILES;
- a **force-field option** (OPLS-AA or GAFF2) recorded per entry, naming where curated parameters
  exist;
- **initiators** (for example KPS, which leaves an SO4⁻ end group) with a choice of capping both
  chain ends, only the head, only the tail, or neither. Users can define their own initiators;
- **capped chain SMILES** for a sequence of monomers, so the end-group choice can be applied and
  checked (net charge, formula).

Every entry is validated on load: it parses, it is net-neutral, it has at least one positive and
one negative centre, and its polymerizable group is identified. The library describes structure
(2D connectivity) only. It does not assign atom types or partial charges. That is the job of the force-field
parametrization feature, which uses the `ForceField` vocabulary defined here.

The built-in data lives in one YAML file shipped inside the package. Users add their own entries in
YAML files with the same schema.

## Dependencies <!-- rq-561736a8 -->

- `rdkit` (structure parsing, validation, SMARTS matching) and `PyYAML` (library files) are runtime
  dependencies of `zpoly`.

## Data Files <!-- rq-14c6b03b -->

### Built-in library <!-- rq-b184f495 -->

`src/zpoly/zwitterions/data/zwitterion_library.yaml` is shipped as package data. It contains:

- every building block listed in `docs/zwitterion_catalog.md` (sections 1–5);
- every named monomer in the catalog (section 6);
- every monomer of the Schneider *et al.* set (section 9), merged with the section-6 entry when the
  structure is identical (for example `SB` is an alias of `DMAPS`, and `CB` is an alias of `CBMA-1`).
  `ZB4` is the as-synthesized structure (ethylene spacer); `ZB4-C3` is the structure as drawn in
  SI Table S1 (propylene spacer). Both are present and each one's notes mention the other.
  The paper's computed-only `M5` (N-dibenzyl) is deliberately **not** included; its `NBn2` cation
  block remains available for `compose`;
- every small-molecule zwitterion in the catalog (section 7);
- the initiators listed in *Built-in initiators*.

### Schema <!-- rq-32c147ad -->

A library file is a YAML mapping with these keys:

- `version` (required): integer schema version. The only supported value is `1`.
- `blocks` (optional): a mapping from block category (`polymerizable`, `linker`, `cation`,
  `spacer`, `anion`) to a mapping of block key → block record:
  - `name` (required), `smiles` (required), `smarts` (required for `polymerizable`, `cation`,
    `anion`; optional otherwise), `notes` (optional).
  - Attachment points in `smiles` are written `[*:1]` (toward the backbone) and `[*:2]` (away from
    the backbone). A block with only `[*:1]` is **terminal**. A polymerizable block has only
    `[*:2]`.
- `monomers` (optional): a list of entry records:
  - `name` (required): the primary identifier, for example `DMAPS`.
  - `aliases` (optional): list of alternative identifiers.
  - `full_name` (required), `family` (required, free text such as `Sulfobetaine`).
  - `status` (required): one of `common`, `specialty`, `combinatorial`, `computed`.
  - exactly one of `blocks` (an ordered list of block keys) or `smiles`. When an entry is written
    with `smiles` but a block decomposition is also known, `blocks` may be given as well, and
    the structure assembled from the blocks must be identical to the SMILES.
  - `residue` (optional): a two-character residue prefix (uppercase letters and digits) used in
    structure and topology files, for example `DM`. Builders append the chain-position letter
    `F`, `R` or `E` (see `rqm/polymer_structure.md`), giving `DMF`, `DMR`, `DME`. When it is
    absent, builders derive the prefix from the name.
  - `source` (optional): citation text for user entries. The built-in entries leave it empty;
    their references are listed once in `REFERENCES.md`.
  - `notes` (optional).
  - `properties` (optional): mapping of property name → number or string (for example
    `fmin_kbt: 2.77`, `fmin_kbt_error: 0.23`).
  - `parameters` (optional): mapping of force-field identifier (`oplsaa`, `gaff2`) → text naming a
    curated parameter source (for example `polybuilder:DMAPS`).
- `small_molecules` (optional): a list of entry records with the same fields as `monomers`, except
  that they need no polymerizable group.
- `comonomers` (optional): polymerizable monomers that are **not** zwitterions, for copolymers with
  the zwitterions. Same fields as `monomers`, plus `counterion`. They skip the zwitterion check
  (net charge 0 with both charge signs) and instead need: net charge 0 and no `counterion`
  (neutral, e.g. NIPAM), or a non-zero net charge and a `counterion` of the opposite sign
  (`Cl-`, `Br-`, `Na+`, `K+`), e.g. a DABCO dication with `Cl-`. Their entries have
  `zwitterionic = False`; `make_variant` rejects them (`VariantError`).
- `initiators` (optional): a list of initiator records:
  - `name` (required), `aliases` (optional), `full_name` (required).
  - `head` and `tail` (at least one required): end-group SMILES fragments, each with exactly one
    attachment point `[*:1]`. `head` caps the chain start (it bonds to the CH2 of the first unit,
    where the initiator radical adds). `tail` caps the chain end (it bonds to the substituted
    carbon of the last unit). Either may be omitted or `null` for a one-sided initiator.
  - `counterion` (optional): counterion that balances the charged end groups, for example `K+`.
  - `source`, `notes` (optional), `parameters` (optional, as for monomers).

Unknown top-level keys or unknown record fields are errors, so that typos are not ignored.

### Comonomers <!-- rq-ddcf8658 -->

Built in, taken over from polybuilder: `NIPAM` (N-isopropylacrylamide) and `IBOA` (isobornyl
acrylate), neutral; `VBD` (4-vinylbenzyl-DABCO with an N-propyl, +2) and `BVD`
(bis(4-vinylbenzyl)-DABCO, +2, one styrene polymerizes and the other stays as a pendant), both
with `Cl-`. A chain's counterion follows the sign of its net charge: for a negative chain the
initiator's counterion, for a positive one the counterion of its charged comonomers (the
neutralizing count is |net charge| per chain).

When several terminal C=C groups are symmetry-equivalent (as in BVD), the first one polymerizes
and no `vinyl_atoms` are needed.

### Curated parameters at the time of writing <!-- rq-9b2b9d78 -->

The built-in `parameters` entries record the OPLS-AA residues that polybuilder already carries:
monomers `DMAPS`, `A3367`, `M3295`, `A3361`, and initiators `KPS`, `APS`. No built-in entry has
curated GAFF2 parameters.

## Variants: N-Alkyl Length and Spacer Length <!-- rq-6516039f -->

Users often study a monomer family by changing two chain lengths. `make_variant` derives a new
monomer from any entry that has a block decomposition:

- **Spacer length** (`spacer=n`, n ≥ 1): the number of CH2 groups between the cation block and the
  anion block, for example N⁺–(CH2)ₙ–SO3⁻. The spacer block between the two charged blocks
  (in either order) is replaced by the linear `C<n>` spacer. If that spacer is not a plain
  alkylene chain (`C<k>`, for example the hydroxy spacer `C3OH`), or there is no spacer block
  between the charged blocks (for example the N-oxide), `VariantError` is raised.
- **N-alkyl length** (`alkyl=n`, n ≥ 1): the length of the identical linear alkyl groups on the
  ammonium nitrogen, so 1 = methyl, 2 = ethyl, 3 = propyl, and so on. This applies to the cation
  blocks `NMe2` and `NEt2` (two alkyl groups) and `NMe3` (three alkyl groups, as in MPC). For any
  other cation block (for example imidazolium or pyrrolidinium) `VariantError` is raised.
- A variant is validated like any composed structure (net-neutral, polymerizable group found).
- When an option equals the entry's existing value, it changes nothing; when no option is given,
  the entry is returned unchanged.
- The variant's name is the entry name plus `-sp<n>` and/or `-alk<n>` for each option that changes
  the structure, in that order, for example `DMAPS-sp5-alk2`. Its `status` is `user`, its
  `residue` is `None` (residue names are derived), its `blocks` are the new block list, and its
  `notes` say which entry it was derived from. Aliases, properties and curated `parameters` are
  not carried over, since they describe the original structure.

## Name Normalization <!-- rq-beab5e79 -->

Names and aliases are matched case-insensitively and ignoring `-`, `_`, spaces and `/`. For
example, `cbma1`, `CBMA-1` and `CBMA_1` all refer to the same entry. Two entries (built-in or user)
whose normalized names or aliases collide are a `DuplicateEntryError`.

## Composition Rules <!-- rq-458386a5 -->

`compose` joins blocks in the order given:

- The first block is a `polymerizable` block.
- Each middle block has both `[*:1]` and `[*:2]`. The last block is terminal (`[*:1]` only).
- A linker is optional: a cation or anion may follow the polymerizable group directly (for
  example vinylpyridinium).
- Block keys are unique across all categories. Linear `(CH2)n` chains are stored once, in the
  `spacer` category (`C1`–`C4`), and serve in the linker position as well. The `linker` category
  holds only linkers that are not plain alkylene chains.
- Spacer shorthand: `C<n>` with any integer `n ≥ 1` is a linear `(CH2)n` spacer or linker, even
  when `C<n>` is not a stored block key.
- The order *cation → spacer → anion* (betaines) and *anion → spacer → cation*
  (phosphorylcholine-type, "inverted") are both valid. Validity is decided by the charge rules
  below, not by the block order.
- The result must be net-neutral with at least one positive and one negative formal charge.

## Initiators and End Groups <!-- rq-c7899aba -->

### Built-in initiators <!-- rq-e283b2ae -->

| Name | Aliases | Head group | Tail group | Counterion |
|---|---|---|---|---|
| `KPS` (potassium persulfate) | `SO4`, `KSO4`, `K2S2O8` | sulfate `[*:1]OS(=O)(=O)[O-]` | sulfate | `K+` |
| `APS` (ammonium persulfate) | | sulfate | sulfate | `NH4+` |
| `NaPS` (sodium persulfate) | | sulfate | sulfate | `Na+` |
| `AIBA` (2,2'-azobis(2-methylpropionamidine)) | `V-50` | amidinium `[*:1]C(C)(C)C(N)=[NH2+]` | amidinium | `Cl-` |
| `AIBN` | | `[*:1]C(C)(C)C#N` | same | — |
| `ACVA` (4,4'-azobis(4-cyanovaleric acid)) | `V-501` | `[*:1]C(C)(CCC(=O)O)C#N` | same | — |
| `TBHP` (tert-butyl hydroperoxide) | | `[*:1]OC(C)(C)C` | same | — |
| `H2O2` (hydroxyl radical) | `hydroxyl` | `[*:1]O` | same | — |
| `CPA` (RAFT agent, 4-cyano-4-(phenylcarbonothioylthio)pentanoic acid) | | R group `[*:1]C(C)(CCC(=O)O)C#N` | dithiobenzoate `[*:1]SC(=S)c1ccccc1` | — |

Initiator names have their own namespace and normalization (see *Name Normalization*), separate
from monomer names.

### End selection <!-- rq-fbb9e221 -->

The `ends` choice is one of `both`, `head`, `tail`, `none`:

- `both`: cap both ends. The initiator needs both a head and a tail fragment.
- `head` / `tail`: cap only that end. The initiator needs that fragment.
- `none`: no caps.
- If an initiator is given and `ends` is not, `ends` defaults to `both` when the initiator has both
  fragments, otherwise to the side it has.
- Without an initiator, `ends` must be `none` or not given.
- An uncapped end is hydrogen-terminated (the head CH2 becomes CH3, the tail carbon gains an H).

The chain's net charge is the sum of the end-group charges, since every unit is net-neutral. For
example, a KPS chain capped at both ends has net charge −2.

## Feature API <!-- rq-a27c63b4 -->

All public names are importable from `zpoly.zwitterions`.

### Types <!-- rq-d6ec487d -->

- `ForceField(str, Enum)`: the force fields known to the project, `OPLSAA = "oplsaa"` and <!-- rq-a09ee60d -->
  `GAFF2 = "gaff2"`. `ForceField("OPLS-AA")` and `ForceField("gaff2")` both resolve, so lookup is
  case-insensitive and ignores `-`.
- `BlockCategory(str, Enum)`: `POLYMERIZABLE`, `LINKER`, `CATION`, `SPACER`, `ANION`. <!-- rq-158c925e -->
- `Block` (frozen dataclass): `key`, `category: BlockCategory`, `name`, `smiles`, <!-- rq-e896b7cf -->
  `smarts: str | None`, `formal_charge: int`, `terminal: bool`, `notes: str`.
- `Zwitterion` (frozen dataclass): one library entry or composed/imported structure. <!-- rq-5e9d599c -->
  - `name: str`, `aliases: tuple[str, ...]`, `full_name: str`, `family: str`, `status: str`
    (`common`, `specialty`, `combinatorial`, `computed`, or `user` for structures from
    `compose`/`from_smiles`).
  - `smiles: str`: RDKit canonical SMILES.
  - `blocks: tuple[str, ...] | None`: the block keys, if a decomposition is known.
  - `polymerizable_group: str | None`: key of the matched polymerizable block, or `None` for a
    small molecule.
  - `vinyl_atoms: tuple[int, int] | None`: heavy-atom indices (in `smiles` order) of the
    polymerizing C=C, terminal CH2 first.
  - `cation_classes: tuple[str, ...]`, `anion_classes: tuple[str, ...]`: keys of all cation and
    anion blocks whose SMARTS match.
  - `formula: str`, `molecular_weight: float`.
  - `source: str`, `notes: str`, `properties: Mapping[str, float | str]`.
  - `residue: str | None`: the two-character residue prefix from the library file, if any.
  - `parameters: Mapping[ForceField, str]`: curated parameter sources.
  - `polymerizable: bool` (property): `polymerizable_group is not None`.
  - `to_mol(add_hs: bool = True) -> rdkit.Chem.Mol`: a new RDKit molecule.
  - `repeat_unit_smiles() -> str`: the polymerized repeat unit. The C=C becomes a single bond, the
    terminal CH2 gets `[*:1]` and the other vinyl carbon gets `[*:2]`. Raises
    `NoPolymerizableGroupError` for small molecules.
- `Library`: an immutable, validated collection of blocks and entries. <!-- rq-226ae0c8 -->
  - `get(name: str) -> Zwitterion`: look up by name or alias after normalization. Raises
    `UnknownEntryError`.
  - `block(key: str) -> Block`: raises `UnknownBlockError`. Accepts `C<n>` spacer shorthand.
  - `blocks(category: BlockCategory | str | None = None) -> list[Block]`
  - `entries(*, family: str | None = None, status: str | None = None,
    polymerizable: bool | None = None, cation: str | None = None, anion: str | None = None,
    force_field: ForceField | str | None = None) -> list[Zwitterion]`: all filters are combined with AND.
    `family` matches case-insensitively as a substring. `cation`/`anion` match against
    `cation_classes`/`anion_classes`. `force_field` keeps only entries with curated parameters
    for that force field. Results are sorted by `name`.
  - `__contains__(name)`, `__len__()` (number of entries), `__iter__()` (entries sorted by name).
  - `initiator(name: str) -> Initiator`: look up by name or alias after normalization. Raises
    `UnknownInitiatorError`.
  - `initiators() -> list[Initiator]`: all initiators sorted by name.
- `Ends(str, Enum)`: `BOTH = "both"`, `HEAD = "head"`, `TAIL = "tail"`, `NONE = "none"`. <!-- rq-4711c4d5 -->
- `Initiator` (frozen dataclass): `name`, `aliases: tuple[str, ...]`, `full_name`, <!-- rq-6114615c -->
  `head_smiles: str | None`, `tail_smiles: str | None`, `head_charge: int`, `tail_charge: int`,
  `counterion: str | None`, `source`, `notes`, `parameters: Mapping[ForceField, str]`.
- `Chain` (frozen dataclass): the result of `build_chain`. <!-- rq-f3793f01 -->
  - `smiles: str` (RDKit canonical), `units: tuple[str, ...]` (unit names in order),
    `initiator: str | None`, `ends: Ends`, `head_group: str | None`, `tail_group: str | None`
    (the fragment SMILES used), `net_charge: int`, `counterion: str | None`, `formula: str`,
    `molecular_weight: float`.

### Functions <!-- rq-338a5f7b -->

- `default_library() -> Library`: the built-in library. It is loaded and validated once per <!-- rq-26d9d524 -->
  process and cached.
- `load_library(*paths: str | os.PathLike, include_builtin: bool = True) -> Library`: the built-in <!-- rq-54e8d4f2 -->
  library (unless `include_builtin=False`) merged with user YAML files, in the order given. User
  blocks and entries may not redefine built-in ones (`DuplicateEntryError`). User entries may
  use built-in blocks.
- `compose(*block_keys: str, name: str | None = None, library: Library | None = None) -> Zwitterion`: <!-- rq-3c5af765 -->
  assembles a structure per *Composition Rules* using blocks from `library` (default:
  `default_library()`). `name` defaults to the block keys joined by `-`. Status is `user`.
- `make_variant(entry: Zwitterion | str, *, spacer: int | None = None, alkyl: int | None = None, <!-- rq-c9621362 -->
  library: Library | None = None) -> Zwitterion`: derives a variant per *Variants*. A string
  `entry` is looked up in `library` (default: `default_library()`). Raises `VariantError`, or
  `ValueError` for an option below 1.
- `from_smiles(smiles: str, name: str, *, require_polymerizable: bool = True, <!-- rq-d2f618ff -->
  vinyl_atoms: tuple[int, int] | None = None, library: Library | None = None) -> Zwitterion`:
  validates an arbitrary structure. It identifies the polymerizable group with the polymerizable
  block SMARTS. If no group matches and `require_polymerizable` is false, the result is a small
  molecule. If more than one distinct C=C could polymerize, the caller must pass `vinyl_atoms`
  (except for the diallyl block, whose two allyl arms are one cyclopolymerizing group). Status is
  `user`.
- `custom_initiator(name: str, *, head: str | None = None, tail: str | None = None, <!-- rq-44802ad6 -->
  counterion: str | None = None, full_name: str | None = None) -> Initiator`: defines an initiator in
  code with the same validation as library files. At least one of `head`/`tail` is required, and each
  must parse and carry exactly one `[*:1]`.
- `build_chain(units: Sequence[Zwitterion | str], *, initiator: Initiator | str | None = None, <!-- rq-a2534bc5 -->
  ends: Ends | str | None = None, library: Library | None = None) -> Chain`: links the repeat units
  in order (unit *i*'s `[*:2]` to unit *i+1*'s `[*:1]`) and applies the end groups per *End
  selection*. Names are resolved in `library` (default: `default_library()`), and a string
  `initiator` is looked up with `Library.initiator`. Homopolymer: `build_chain(["DMAPS"] * 20,
  initiator="KPS")`.

### Errors <!-- rq-ef81db1c -->

- `ZwitterionError(Exception)`: base class. <!-- rq-687220b8 -->
  - `UnknownEntryError`: carries `name` and `suggestions: tuple[str, ...]` (up to 3 close
    names). The message includes the suggestions.
  - `UnknownBlockError`: carries `key`.
  - `InvalidStructureError`: SMILES does not parse or sanitize. Carries `smiles` and `reason`.
    - `ChargeError`: net charge is not 0, or there is no positive or no negative centre. Carries
      `net_charge`.
  - `IncompatibleBlocksError`: block attachment points do not chain (for example a terminal block
    in the middle, or a first block that is not polymerizable). Carries `block_keys` and
    `position`.
  - `NoPolymerizableGroupError`: no polymerizable group was found where one is required.
  - `AmbiguousPolymerizableGroupError`: several candidate C=C groups that are not
    symmetry-equivalent, and no `vinyl_atoms` given. Carries `candidates`.
  - `DuplicateEntryError`: carries `name` and the conflicting `sources`.
  - `VariantError(ZwitterionError, ValueError)`: the entry cannot take the requested variant
    option. Carries `name` and `option` (`"spacer"` or `"alkyl"`); the message says why.
  - `UnknownInitiatorError(UnknownEntryError)`: no initiator has this name or alias.
  - `ChainError(ZwitterionError, ValueError)`: `build_chain` was given an empty sequence or a unit
    that is not polymerizable.
    - `EndGroupError`: the `ends` choice needs a fragment the initiator does not have, or an end
      cap was requested without an initiator. Carries `ends` and `initiator`.
  - `LibraryFileError`: a library file is missing, unreadable, not valid YAML, violates the
    schema, or contains an entry that fails validation. Carries `path` and `reason`; the reason
    names the offending entry.

## Catalog Document <!-- rq-a3c277ca -->

`docs/zwitterion_catalog.md` is generated from the built-in library by
`devtools/build_zwitterion_catalog.py`, which imports `zpoly.zwitterions`. It includes
a table of the initiators. The YAML file
is the single source of truth; the Markdown is never edited by hand.

## Out of Scope <!-- rq-6945abea -->

- Atom typing, partial charges and bonded parameters for either force field (the force-field
  parametrization feature).
- 3D chain structures, mid-chain initiator bridges, initiator placement by percentage, and
  PDB/RTP/topology files (the polybuilder adapter feature).
- Command-line and web front ends.

## Gherkin Scenarios <!-- rq-14b51091 -->

```gherkin
Feature: Zwitterion monomer library

  @rq-ae25abea
  Scenario: Built-in library loads and every entry passes conformance
    When default_library() is called
    Then every entry parses and has net formal charge 0
    And every entry has at least one positive and one negative formal charge
    And every entry in "monomers" has a polymerizable group and vinyl atoms
    And every entry in "small_molecules" is not required to be polymerizable
    And every entry with blocks has a SMILES identical to the structure assembled from its blocks

  @rq-83c0765f
  Scenario: Built-in library contains the catalog and paper monomers
    When default_library() is called
    Then it contains the entries "DMAPS", "MPC", "CBAA", "SPV" by name or alias
    And it contains "ZB", "ZB2", "ZB3", "ZB4", "ZB4-C3", "ZB5" and "M1" through "M19" except "M5"
    And it does not contain "M5"
    And it contains the small molecules "TMAO" and "MIm-PS"

  @rq-471ad711
  Scenario: Look up an entry by alias with name normalization
    When get("sb") and get("SBMA") and get("DMAPS") are called on the default library
    Then all three return the same entry

  @rq-881ae9e4
  Scenario: Unknown name gives suggestions
    When get("DMAP") is called on the default library
    Then UnknownEntryError is raised with name "DMAP"
    And its suggestions include "DMAPS"

  @rq-ad1e99c1
  Scenario: Filter entries by family and status
    When entries(family="carboxybetaine", status="common") is called on the default library
    Then every returned entry has family containing "Carboxybetaine" and status "common"
    And the result includes the entry named "CBMA"

  @rq-c9865bcb
  Scenario: Filter entries by anion class
    When entries(anion="OSO3") is called on the default library
    Then the result includes "ZB" and "ZB5"
    And every returned entry has "OSO3" in its anion classes

  @rq-6171c7c9
  Scenario: Filter entries by force field with curated parameters
    When entries(force_field="oplsaa") is called on the default library
    Then the result includes "DMAPS"
    And every returned entry has an "oplsaa" parameter source
    And entries(force_field=ForceField.GAFF2) returns an empty list

  @rq-abf15540
  Scenario: Force field names are case-insensitive
    When ForceField("OPLS-AA") and ForceField("GAFF2") are evaluated
    Then they equal ForceField.OPLSAA and ForceField.GAFF2

  @rq-82d1c963
  Scenario: Paper properties are attached to entries
    When get("ZB5") is called on the default library
    Then its properties contain fmin_kbt 2.77
    And REFERENCES.md cites Schneider et al.

  @rq-e5af84ba
  Scenario: Diallyl monomer is recognised as cyclopolymerizing
    When get("DAMAPS") is called on the default library
    Then its polymerizable group is "Al"

  @rq-9769c65a
  Scenario: Compose DMAPS from blocks
    When compose("MA", "C2", "NMe2", "C3", "SO3") is called
    Then the result has the same SMILES as get("DMAPS")
    And its status is "user"
    And its name is "MA-C2-NMe2-C3-SO3"

  @rq-7141c236
  Scenario: Compose with spacer shorthand beyond stored keys
    When compose("MA", "C2", "NMe2", "C6", "SO3") is called
    Then the result contains a six-carbon chain between the ammonium nitrogen and the sulfur
    And it has net charge 0

  @rq-f26c0a78
  Scenario: Compose an inverted phosphorylcholine
    When compose("MA", "C2", "PO4i", "C2", "NMe3") is called
    Then the result has the same SMILES as get("MPC")

  @rq-e3d5290e
  Scenario: Compose without a linker
    When compose("V", "Py2", "C3", "SO3") is called
    Then the result has the same SMILES as get("SPV")

  @rq-30c6d511
  Scenario: Compose a single-group zwitterion
    When compose("MA", "C2", "NO") is called
    Then the result has net charge 0 and one positive and one negative centre

  @rq-fb61858f
  Scenario: Compose with an unknown block
    When compose("MA", "C2", "XYZ", "C3", "SO3") is called
    Then UnknownBlockError is raised with key "XYZ"

  @rq-bf757833
  Scenario: Compose a non-neutral structure
    When compose("MA", "C2", "NMe3") is called
    Then ChargeError is raised with net_charge 1

  @rq-3c28b422
  Scenario: Compose with a terminal block in the middle
    When compose("MA", "SO3", "C3", "NMe3") is called
    Then IncompatibleBlocksError is raised with position 1

  @rq-372fced8
  Scenario: Compose that does not start with a polymerizable block
    When compose("C2", "NMe2", "C3", "SO3") is called
    Then IncompatibleBlocksError is raised with position 0

  @rq-e70c4dde
  Scenario: Import a zwitterion from SMILES
    When from_smiles("C=C(C)C(=O)OCC[N+]1(CCCS(=O)(=O)[O-])CCCCC1", "MyPip") is called
    Then the result has polymerizable group "MA"
    And its cation classes contain "Pip" and its anion classes contain "SO3"
    And its status is "user"

  @rq-878e77a9
  Scenario: Import an unparsable SMILES
    When from_smiles("C=C(C", "Bad") is called
    Then InvalidStructureError is raised

  @rq-3e28fc65
  Scenario: Import a charged, non-zwitterionic structure
    When from_smiles("C=C(C)C(=O)OCC[N+](C)(C)C", "METAC") is called
    Then ChargeError is raised with net_charge 1

  @rq-bf11ad76
  Scenario: Import a small molecule
    When from_smiles("C[N+](C)(C)CC(=O)[O-]", "Betaine") is called
    Then NoPolymerizableGroupError is raised
    And from_smiles with require_polymerizable=False returns an entry that is not polymerizable

  @rq-9ea65b16
  Scenario: Import a structure with two candidate vinyl groups
    Given a zwitterion SMILES with a methacrylate and a pendant styrene
    When from_smiles is called without vinyl_atoms
    Then AmbiguousPolymerizableGroupError is raised listing both candidates
    And from_smiles with vinyl_atoms set to the methacrylate C=C succeeds

  @rq-456dad55
  Scenario: Repeat unit of DMAPS
    When repeat_unit_smiles() is called on get("DMAPS")
    Then the result is the canonical form of "[*:1]CC([*:2])(C)C(=O)OCC[N+](C)(C)CCCS(=O)(=O)[O-]"

  @rq-052fdbed
  Scenario: Repeat unit of a small molecule
    When repeat_unit_smiles() is called on get("TMAO")
    Then NoPolymerizableGroupError is raised

  @rq-3ae76166
  Scenario: Load a user library file
    Given a YAML file with version 1 and a monomer "MYSB" built from blocks ["AAm", "C3", "Pip", "C4", "SO3"]
    When load_library(path) is called
    Then the library contains "MYSB" and all built-in entries

  @rq-bd5de7cc
  Scenario: User file redefines a built-in name
    Given a YAML file with version 1 and a monomer whose alias is "cbma-1"
    When load_library(path) is called
    Then DuplicateEntryError is raised with name "cbma-1"

  @rq-6698acbe
  Scenario: User file is not valid YAML
    Given a file containing "monomers: [unclosed"
    When load_library(path) is called
    Then LibraryFileError is raised carrying that path

  @rq-a50b7bb2
  Scenario: User file does not exist
    When load_library("missing.yaml") is called
    Then LibraryFileError is raised carrying the path "missing.yaml"

  @rq-63630b82
  Scenario: User file has an unknown field
    Given a YAML file with version 1 and a monomer record containing the field "smile"
    When load_library(path) is called
    Then LibraryFileError is raised whose reason names "smile"

  @rq-134bdc27
  Scenario: User file has an unsupported version
    Given a YAML file with version 2
    When load_library(path) is called
    Then LibraryFileError is raised

  @rq-cf1496b1
  Scenario: User entry whose blocks and SMILES disagree
    Given a YAML file with a monomer giving blocks ["MA", "C2", "NMe2", "C3", "SO3"] and the SMILES of CBMA
    When load_library(path) is called
    Then LibraryFileError is raised whose reason names that monomer

  @rq-3008e140
  Scenario: User entry that is not neutral
    Given a YAML file with a monomer whose SMILES is "C=C(C)C(=O)OCC[N+](C)(C)C"
    When load_library(path) is called
    Then LibraryFileError is raised whose reason names that monomer

  @rq-33491788
  Scenario: Catalog document is up to date
    When devtools/build_zwitterion_catalog.py is run
    Then its output is identical to docs/zwitterion_catalog.md

  @rq-30331843
  Scenario: Built-in initiators include KPS with sulfate end groups
    When initiator("SO4") and initiator("kso4") are called on the default library
    Then both return the initiator named "KPS"
    And its head and tail fragments are sulfate with charge -1
    And its counterion is "K+"

  @rq-f5dfd922
  Scenario: Initiator lookup with an unknown name
    When initiator("KPSS") is called on the default library
    Then UnknownInitiatorError is raised whose suggestions include "KPS"

  @rq-9218ad81
  Scenario: KPS caps both chain ends by default
    When build_chain(["DMAPS"] * 3, initiator="KPS") is called
    Then the chain has ends "both" and net charge -2
    And the chain SMILES contains two sulfate end groups bonded through oxygen

  @rq-94a7ac1f
  Scenario: KPS caps only the head
    When build_chain(["DMAPS"] * 3, initiator="KPS", ends="head", cap="H") is called
    Then the chain has net charge -1
    And the head CH2 of the first unit is bonded to the sulfate oxygen
    And the last unit's tail carbon carries a hydrogen instead of an end group

  @rq-4fa24446
  Scenario: KPS caps only the tail
    When build_chain(["DMAPS"] * 3, initiator="KPS", ends="tail") is called
    Then the chain has net charge -1
    And the first unit starts with a CH3 group

  @rq-f9862ff3
  Scenario: Chain without end groups
    When build_chain(["DMAPS"] * 3, cap="H") is called
    Then the chain has ends "none" and net charge 0
    And its formula equals three DMAPS units plus two hydrogens
    And with the default cap="CH3" it equals three DMAPS units plus H and CH3 (the chain ends are methyls)

  @rq-0c9b91ab
  Scenario: End cap requested without an initiator
    When build_chain(["DMAPS"] * 3, ends="both") is called
    Then EndGroupError is raised

  @rq-3ed1d11c
  Scenario: Different head and tail groups (RAFT)
    When build_chain(["DMAPS"] * 2, initiator="CPA") is called
    Then the chain contains the cyanopentanoic acid head and the dithiobenzoate tail

  @rq-21f96ecb
  Scenario: Copolymer sequence
    When build_chain(["DMAPS", "CBMA", "DMAPS"], initiator="KPS") is called
    Then the chain units are ("DMAPS", "CBMA", "DMAPS") and its net charge is -2

  @rq-8d1fb6d8
  Scenario: Custom initiator defined in code
    When custom_initiator("MYINIT", head="[*:1]OC(=O)c1ccccc1") is created
    And build_chain(["DMAPS"] * 2, initiator=that initiator) is called
    Then the chain ends default to "head"
    And build_chain with the same initiator and ends="both" raises EndGroupError

  @rq-b042d5d7
  Scenario: Custom initiator with an invalid fragment
    When custom_initiator("BAD", head="[*:1]OC(=O)[*:2]") is called
    Then InvalidStructureError is raised

  @rq-82181c68
  Scenario: Custom initiator without fragments
    When custom_initiator("EMPTY") is called
    Then InvalidStructureError is raised

  @rq-c62fdb4d
  Scenario: Initiator defined in a user library file
    Given a YAML file with version 1 and an initiator "VA044" with head and tail "[*:1]C(C)(C)C1=[NH+]CCN1" and counterion "Cl-"
    When load_library(path) is called
    Then initiator("VA044") returns it with head charge 1 and tail charge 1
    And build_chain(["MPC"] * 2, initiator="VA044", library=that library) has net charge 2

  @rq-dec70936
  Scenario: User initiator redefines a built-in name
    Given a YAML file with version 1 and an initiator whose alias is "k2s2o8"
    When load_library(path) is called
    Then DuplicateEntryError is raised with name "k2s2o8"

  @rq-42202c0e
  Scenario: Building a chain from a small molecule
    When build_chain(["TMAO"]) is called
    Then ChainError is raised

  @rq-cffe1741
  Scenario: Building an empty chain
    When build_chain([]) is called
    Then ChainError is raised

  @rq-9fa59470
  Scenario: Longer spacer between cation and anion
    When make_variant("DMAPS", spacer=5) is called
    Then the result has exactly five CH2 groups between the ammonium N and the sulfonate S
    And its name is "DMAPS-sp5" and its status is "user"
    And it has net charge 0

  @rq-1c79d863
  Scenario: Spacer variant of an inverted (phosphorylcholine) monomer
    When make_variant("MPC", spacer=4) is called
    Then the result has four CH2 groups between the phosphate O and the ammonium N

  @rq-10b3fee1
  Scenario: Longer N-alkyl groups
    When make_variant("DMAPS", alkyl=2) is called
    Then the ammonium N carries two ethyl groups
    And the result has the same SMILES as get("SBMA-Et")

  @rq-f083b2fa
  Scenario: Spacer and alkyl together
    When make_variant("CBMA", spacer=3, alkyl=3) is called
    Then the result's name is "CBMA-sp3-alk3"
    And the ammonium N carries two propyl groups and the carboxylate is three CH2 groups away

  @rq-046183e6
  Scenario: Options equal to the current values change nothing
    When make_variant("DMAPS", spacer=3, alkyl=1) is called
    Then the result is the DMAPS entry itself

  @rq-942f1853
  Scenario: Variant of an entry without a plain spacer
    When make_variant("HSBMA", spacer=4) is called
    Then VariantError is raised with option "spacer"

  @rq-33cd927c
  Scenario: Alkyl variant of a non-ammonium cation
    When make_variant("SPVI", alkyl=2) is called
    Then VariantError is raised with option "alkyl"

  @rq-9ca58870
  Scenario: Variant of an entry without blocks
    When make_variant("M1", spacer=2) is called
    Then VariantError is raised

  @rq-5ec68583
  Scenario: Invalid variant length
    When make_variant("DMAPS", spacer=0) is called
    Then ValueError is raised

  @rq-5769e593
  Scenario: Neutral and charged comonomers
    When the built-in library is loaded
    Then NIPAM and IBOA have net charge 0 and no counterion, VBD and BVD have net charge +2 and counterion Cl-
    And all four are polymerizable with zwitterionic False

  @rq-914103fa
  Scenario: Comonomer charge rules
    Given a user file with a comonomer of net charge +1 and no counterion, or a neutral one with counterion Cl-, or charge +1 with counterion K+
    Then loading raises LibraryFileError each time

  @rq-b7a9f7b4
  Scenario: Copolymer with a cationic comonomer
    When build_chain(["DMAPS", "VBD", "DMAPS"], initiator="KPS") is called
    Then its net charge is 0 and its counterion is None
    And build_chain(["VBD", "VBD"], initiator="KPS") has net charge +2 and counterion Cl-

  @rq-21cde3c2
  Scenario: Symmetric divinyl comonomer
    When BVD is built into a chain of 3 units
    Then each unit keeps one C=C double bond
```
