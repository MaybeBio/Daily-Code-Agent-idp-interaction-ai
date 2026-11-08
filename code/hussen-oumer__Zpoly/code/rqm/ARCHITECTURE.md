# Architecture <!-- rq-cf532132 -->

## Purpose <!-- rq-0883fa2e -->

Zpoly is an open, general-purpose toolkit for building molecular dynamics models of
**zwitterionic polymers and zwitterionic liquids**, aimed at studies such as liquid–liquid phase
separation (LLPS). Its defining goal is ease of use: a user names a monomer (or supplies a SMILES
string), describes the polymer and the surrounding system, and receives a complete, simulation-ready
system for their MD engine of choice, with no hand-written atom types or charges.

The intended audience is the **broad public research community** working with polymer and soft-matter
simulations, not only the originating group. The project is therefore distributed as an installable
Python package (PyPI / conda-forge) with stable public interfaces, user documentation, and versioned
releases.

## Scope and Key Features <!-- rq-519df0f5 -->

- **General zwitterion chemistry.** The tool is not tied to one chemical family. It covers
  zwitterionic polymers (sulfobetaines, carboxybetaines, phosphorylcholines, and others),
  polyampholytes (chains mixing cationic and anionic monomers), and small-molecule zwitterionic
  liquids. The concrete list of chemistries and the command set will be specified by the user
  step by step through `/rr-plan`. The architecture must not assume a fixed set of families.
- **Monomer by name.** A built-in, curated library lets users refer to common monomers and
  zwitterionic liquids by name (for example DMAPS, SBMA, CBMA, MPC) with validated parameters.
- **Automatic parametrization.** Any new monomer or zwitterionic liquid given as SMILES is assigned
  atom types and partial charges automatically, under the active force field. This is the
  feature with the most architectural weight: chemistry is described once, engine-neutrally, and
  parameters are derived from that description, never typed by hand.
- **Polymer construction.** Homopolymers, copolymers (alternating, block, random, repeating
  patterns), end groups and initiators, with chain length and sequence chosen by the user.
- **Full LLPS-ready systems.** The tool goes beyond single chains. It builds a simulation box with
  multiple chains, solvent, salt/counter-ions and optional additional species. It writes a complete
  topology and template run-input files (minimization, equilibration, production).
- **MD engines.** GROMACS is the primary target and ships with the core. OpenMM and LAMMPS
  writers are planned as **extensions** (engine-writer plugins) built on the same `Topology`.
- **One input file, like VASP, Gaussian or Schrödinger.** The user writes a short `KEY = value`
  file and runs `zpoly file.inp`. `list`, `help [KEY]` and `template JOB` teach the
  keywords (`rqm/input_file.md`). The same jobs are a few lines of Python (`rqm/simple_api.md`). A
  local web app may later be added as another front end over the same layer.

### Non-Goals <!-- rq-ff620166 -->

- **Coarse-grained models** (for example Martini). Only all-atom resolution is in scope.
- **A hosted public web service.** The GUI runs locally. There are no servers, accounts or
  multi-user deployments.
- **Running production simulations.** The tool prepares systems and writes the GROMACS protocol
  and run scripts; the user runs them. Analysis is limited to the publication quantities it sets up
  (density/Tg, slab phase concentrations).
- **New force-field development.** The tool applies existing force fields. It does not run quantum
  chemistry; it takes the user's DFT charges and the user's bonded parameters for the ionic parts.

## Architectural Decisions <!-- rq-0e547541 -->

- **Language:** Python (≥ 3.10). Performance-critical pieces may later be moved
  to compiled code, but there are no compiled components at the outset.
- **Execution target:** CPU, single node. System building is not compute-bound, so no GPU or MPI
  support is planned.
- **Origin.** The project grew out of polybuilder. Everything it needed from polybuilder (its
  zwitterionic monomers, initiators, comonomers and OPLS-AA atom types) now lives in this package's
  data files, and polybuilder is no longer used. polybuilder embedded a whole chain at once with
  RDKit, which does not scale (a 50-unit DMAPS chain did not finish in 10 minutes). This project
  grows 3D chains **unit by unit**: each distinct repeat unit is embedded once, and copies are placed
  along the backbone with clash checks and backtracking, which takes seconds for 50 units.
- **Three layers.** Each layer uses only the layers below it:
  1. **User interface layer** (`cli`, `inputfile`, `simple`): the VASP/Gaussian/Schrödinger-style
     input file, its `list` / `help` / `template` commands, and the short Python API. It checks
     the user's input and reports errors in plain words. It contains no chemistry or file-format
     logic and makes no physical choices for the user (number of chains, box, slab are the
     user's).
  2. **Core geometry engine** (`zwitterions`, `structure`, the packing part of `system`): the
     monomer library (zwitterions, comonomers, initiators, variants) as RDKit molecules; chain
     sequences (homopolymer, alternating, block, random, pattern, explicit); 3D chains by
     unit-by-unit growth with F/R/E residue naming; and multi-chain boxes (bulk, slab) packed with
     Packmol or the built-in packer.
  3. **Topology and parameter layer, output ready** (`forcefield`, the writing part of `system`):
     force-field backends (**GAFF2 default**, OPLS-AA optional) behind one interface; charges from
     the user's DFT (B3LYP/HF, RESP) or the force field's own; user charge and bonded-parameter
     overrides for the ionic parts; integer-charge blocks. It writes ready-to-run GROMACS output:
     `.pdb`, `.gro`, `.top`, `.itp`, `.ff/`, `.mdp`, `run.sh`, plus the analysis of the finished
     runs. OpenMM and LAMMPS writers are planned as extension plugins over the same `Topology`.

  Engine-neutral data (the `Topology`) are stored as structures of arrays. There is no
  engine-to-engine conversion.
- **Packing:** multi-component boxes are packed with **Packmol**, an optional external dependency,
  with a simple built-in packer as a fallback when Packmol is not installed.
- **Dependencies:** core requirements are NumPy, RDKit and PyYAML. Heavier dependencies are
  optional: AmberTools (GAFF2), the user's GROMACS `oplsaa.ff` (OPLS-AA), Packmol, GROMACS itself.
  Without them the core still builds the PDB structures.
- **Charge integrity** is a cross-cutting invariant. Every built molecule and system has an
  integer net charge, zwitterionic units are net-neutral, and charged comonomers and initiators
  declare their counterions.
  Violations are errors, not warnings.

## Extensibility <!-- rq-fc36f6e1 -->

Users extend the tool at two levels:

1. **Data files (simple cases).** A YAML/TOML entry adds a monomer, zwitterionic liquid, ion or
   solvent by name and SMILES, plus optional metadata. It is parametrized automatically by the
   active force field. No Python is required.
2. **Python plugin packages (advanced cases).** New force fields, engine writers, packers,
   polymer-architecture generators, and hand-curated parameter sets are registered through standard
   Python **entry points**, so a third-party package can extend the tool just by being installed.

Every extension category has a clearly defined interface and a central registry. Built-in
components register themselves through the same mechanism as third-party plugins, so built-ins have
no special privileges.

## Testing Strategy <!-- rq-d1b83104 -->

- **Conformance harness for chemistry.** Every library entry, whether built-in or loaded from a
  plugin or data file, automatically runs the same standard checks: integer net charge, zwitterion
  neutrality, valence and connectivity, and complete parameter coverage under each supported force
  field. New entries are tested without writing new tests.
- **Conformance harness for engine writers.** Every writer runs against a shared set of reference
  systems. Generated inputs are validated by the real engines (`gmx grompp`, LAMMPS, OpenMM) with a
  short energy minimization in CI.
- **Physical sanity runs.** Short MD runs check properties such as density and hydration against
  literature values. These run nightly or on demand, not on every commit.
- Unit tests cover each layer in isolation. The front ends are tested through the public API they
  wrap.

## Open Questions <!-- rq-b18f0413 -->

- OPLS-AA typing rules for the monomers that do not have them yet.
- A local web app as a second front end.
- The package's public distribution name and the release/versioning policy.
