# Feature: Bulk System for Density and Tg (Dry Amorphous Polymer) <!-- rq-4a8e3ff0 -->

For publication, the polymer model is validated on two bulk properties of the **dry** amorphous
polymer: the **density** at 300 K and the **glass transition temperature Tg**. This feature builds
the bulk system and writes a complete, ready-to-run GROMACS protocol. After the user runs it (on a
workstation or HPC), one command analyses the results.

It also introduces the parts shared by all multi-chain systems (`rqm/system/slab.md` builds on them):
packing several chains (Packmol by default, a built-in packer as the alternative), counterions,
water/ion force-field files, and a system topology.

```python
import zpoly as pp
pp.bulk_system("DMAPS", n=12, initiator="KPS", dft="dft", chains=20, box_nm=6.0, out="dmaps_bulk")
# run dmaps_bulk/run.sh (or sbatch dmaps_bulk/run_slurm.sh), then:
pp.analyze_tg("dmaps_bulk")          # density at 300 K, Tg, CSV of density vs T
```

## The System <!-- rq-88cb8d61 -->

- **The user decides the system**: the number of chains (`chains`) and the box (`box_nm`, an edge
  length for a cubic box or `(Lx, Ly, Lz)`) are required; the code never chooses them. It reports the
  resulting starting density (polymer + counterion mass / box volume) so the user can judge it.
- `chains` copies of the chain, each built with its own seed (`seed + k`), so every chain has its
  own coil conformation but the same topology (one `moleculetype`).
- **Counterions** neutralize the system: for a chain with net charge q and counterion `K+`, |q| K⁺
  per chain. Supported counterions: `K+`, `Na+`, `Cl-` (GROMACS ion molecules `K`, `NA`, `CL`).
  Any other counterion raises `SystemBuildError`. Zwitterionic chains without charged end groups
  need none.
- Packed into the user's box; the protocol then compresses the system (NPT at `compress_bar`).

## Packing <!-- rq-f3321ea0 -->

- `packer="packmol"` (default): the chain and ion coordinates are written as PDB, Packmol packs them
  in the box with tolerance 2.0 Å (`seed` passed to Packmol), and its output is read back.
  If `packmol` is not on `PATH`, `SystemBuildError` says to install it
  (`conda install -c conda-forge packmol`) or use `packer="builtin"`.
- `packer="builtin"`: chains are placed one by one at a random position and random orientation
  (seeded); a placement is accepted when no heavy atom is closer than 2.0 Å to any atom already
  placed, with periodic images. Up to 1000 tries per molecule, otherwise `SystemBuildError`
  suggests a lower start density. Ions are placed the same way.
- Both packers give the same atom order: chain 1, chain 2, …, then the ions.

## Water and Ion Parameters <!-- rq-39bee73a -->

- The ions (and, for slab systems, the water) come from a GROMACS force-field folder that matches
  the polymer force field: `amber99sb-ildn.ff` for GAFF2 (same combination rule 2 and 1–4 scaling),
  `oplsaa.ff` for OPLS-AA. It is located like `oplsaa.ff` (GMXDATA, or next to `gmx`), or given as
  `solvent_ff`. Its `[ defaults ]` must match the polymer's, otherwise `SystemBuildError`.
- The needed atom types are appended to the generated `.ff`'s `ffnonbonded.itp` (respecting
  `#ifdef` variants), and the needed molecule blocks are written to `<name>.ff/ions.itp` (and
  `<name>.ff/tip3p.itp` for water). The `.ff` folder stays self-contained.

## System Topology <!-- rq-b2c69832 -->

- `<name>.itp`: the chain's `[ moleculetype ]` (the same content `write_gromacs` puts in a `.top`).
- `system.top`: `#include "<name>.ff/forcefield.itp"`, `#include "<name>.itp"`, the ion (and water)
  includes, `[ system ]`, and `[ molecules ]` listing the polymer and each ion species with counts.
  The system's total charge is 0.

## Protocol (GROMACS) <!-- rq-2d93d71c -->

Written as `.mdp` files, numbered in run order. All use the Verlet scheme, PME electrostatics, 1.2
nm cut-offs, h-bond constraints (LINCS), 2 fs steps, V-rescale thermostat (τ 0.1 ps) and C-rescale
barostat (τ 2 ps, compressibility 4.5 × 10⁻⁵ bar⁻¹) where NPT:

| Step | File | What |
|---|---|---|
| 1 | `01_em.mdp` | steepest-descent minimization |
| 2 | `02_nvt_hot.mdp` | NVT at T_high (default 600 K), 1 ns, velocities generated |
| 3 | `03_npt_compress.mdp` | NPT at T_high, 1000 bar, 2 ns |
| 4 | `04_npt_anneal.mdp` | NPT at T_high, 1 bar, 5 ns |
| 5… | `cool_<T>K.mdp` | NPT at each T from T_high down to T_low (default 200 K) in steps of 20 K, 2 ns each, continuing from the previous step |

- Every number is a parameter (`temperatures`, `times_ns`, `pressure_bar`, `compress_bar`).
- `run.sh` runs the steps in order with `gmx grompp`/`gmx mdrun` (GMX and thread count from the
  environment), stopping at the first error. `run_slurm.sh` is the same with an `#SBATCH` header
  (job name, 1 node, `--time` from the total simulated time at an assumed 50 ns/day, editable).
- `protocol.yaml` records the system, every parameter, and the cooling temperatures; the analysis
  reads it.

## Analysis <!-- rq-a5952119 -->

`analyze_tg(folder, average_last_ns=1.0)`:

- For each cooling step, the density is averaged over its last `average_last_ns` from
  `cool_<T>K.edr` (via `gmx energy`), with its standard deviation. If an `.edr` is missing,
  `SystemAnalysisError` names the step.
- **Tg** comes from a bilinear fit of density versus temperature: for every split between two
  neighbouring temperatures (at least 3 points on each side), two straight lines are fitted by least
  squares; the split with the smallest total squared error wins, and Tg is where the two lines
  cross.
- **Density at 300 K**: the cooling step at 300 K, if any, otherwise linear interpolation between
  the neighbouring steps.
- Writes `tg_density.csv` (T, density, std) and `tg_result.yaml` (Tg, both lines' slopes and
  intercepts, density at 300 K), prints a summary, and returns a `TgResult` (`tg`, `density_300`,
  `temperatures`, `densities`, `std`, `fit`).

## Density Mode (`DENSITY`) <!-- rq-f791d235 -->

Instead of `box_nm`, the user can give the start `density` (g/cm³, input key `DENSITY`): the cube edge
is (M / ρ)^(1/3) (M: chains and counterions). The molecules are packed in a looser cube at 0.35 g/cm³
(or at ρ if lower); after `01_em`, `01_heat` (10 ps NVT at T_high, velocities) and `01_deform`
(NVT continuing from it, `deform` shrinking the cube to the edge in `deform` ns, default 1) bring it
to the requested density; then the protocol continues as usual. `box_nm` and `density` are
exclusive (`ValueError`). `protocol.yaml` records the density and the packing box.

## Madrid-2019 Counterions (`SOLVENT = tip4p2005`) <!-- rq-0de75a80 -->

With `solvent="tip4p2005"` the counterions are the Madrid-2019 ions (K⁺, Na⁺, Cl⁻ at ±0.85 e, with
their LJ and `[ nonbond_params ]` pairs, see `rqm/system/slab.md`); `CHARGE_SCALE = 0.85` is required.
This keeps the dry bulk consistent with the TIP4P/2005 slab.

## Feature API <!-- rq-d8a94b51 -->

Importable from `zpoly` (and `zpoly.system`):

- `bulk_system(monomer, n=None, *, chains, box_nm=None, density=None, initiator=None, ends=None, spacer=None, alkyl=None, <!-- rq-b510a95a -->
  dft=None, ff="gaff2", packer="packmol", t_high=600.0, t_low=200.0, t_step=20.0,
  times_ns=None, pressure_bar=1.0, compress_bar=1000.0, seed=0, solvent_ff=None, out="bulk",
  library=None) -> Path`: builds and parametrizes the chain like `build`, packs the system, and
  writes the folder (returned). `chains` (≥ 1) and `box_nm` (a positive edge, or three positive
  edges) are required (`ValueError` otherwise). `times_ns` overrides the step times by name
  (`nvt_hot`, `compress`, `anneal`, `cool`).
- `analyze_tg(folder, *, average_last_ns=1.0, gmx="gmx") -> TgResult`. <!-- rq-da2f108a -->
- `fit_tg(temperatures, densities) -> dict` (the bilinear fit alone, for any data). <!-- rq-57646a02 -->
- `SystemBuildError`, `SystemAnalysisError` (subclasses of `zpoly.structure.StructureError`). <!-- rq-47e41beb -->

## Input File <!-- rq-b202ec73 -->

`JOB = bulk` with `CHAINS`, `BOX` (required), `TG = T_high T_low step`, `TIMES`, `PACKER`,
`SOLVENT_FF`, `DFT`, `FF` and `OUTPUT` (`rqm/input_file.md`). `JOB = analyze` with `OUTPUT = <folder>`
runs `analyze_tg` on a bulk folder.

## Out of Scope <!-- rq-436a610c -->

- Running the simulations from the package; hydrated Tg; other Tg methods (e.g. from mobility).
- Plots (the CSV is ready for any plotting tool).

## Gherkin Scenarios <!-- rq-766f7511 -->

```gherkin
Feature: Bulk system for density and Tg

  @rq-7f93393d
  Scenario: Bulk folder contents
    Given the fake backend
    When bulk_system("DMAPS", n=3, chains=4, box_nm=3.0, initiator="KPS", ff="fake", packer="builtin", out=tmp/"b") is called
    Then tmp/b contains system.gro, system.top, polymer.itp, polymer.ff/, 01_em.mdp, 02_nvt_hot.mdp, 03_npt_compress.mdp, 04_npt_anneal.mdp, run.sh, run_slurm.sh and protocol.yaml
    And there is one cool_<T>K.mdp for each of 600, 580, ..., 200 K (21 files)

  @rq-a09233c5
  Scenario: Counterions neutralize the system
    Given the same bulk system (4 chains with KPS at both ends)
    Then system.top [ molecules ] lists polymer 4 and K 8
    And the total charge of system.top is 0
    And system.gro has 4 chain copies followed by 8 K atoms

  @rq-f413df4d
  Scenario: Chains have different conformations
    Then no two chains in system.gro have identical internal coordinates (pairwise distance matrices differ)

  @rq-1dba016a
  Scenario: Built-in packer respects the minimum distance
    Then no two heavy atoms of different molecules in system.gro are closer than 2.0 Å (with periodic images)

  @rq-b96476fb
  Scenario: The box is the user's
    Then the box in system.gro is 3.0 x 3.0 x 3.0 nm
    And bulk_system with box_nm=(3.0, 3.5, 4.0) gives that box
    And bulk_system without chains or box_nm, or with a non-positive edge, raises an error

  @rq-db5c6f04
  Scenario: Packmol packer
    Given packmol is on PATH
    When bulk_system with packer="packmol" is called
    Then system.gro has the same atoms in the same order as with the built-in packer, inside the box

  @rq-af2c836d
  Scenario: Packmol missing
    Given packmol is not on PATH
    When bulk_system with packer="packmol" is called
    Then SystemBuildError is raised whose message mentions "conda install -c conda-forge packmol" and packer="builtin"

  @rq-72fe414e
  Scenario: Unsupported counterion
    Given a chain whose counterion is "NH4+"
    When bulk_system is called
    Then SystemBuildError is raised naming "NH4+"

  @rq-9ae6f9ab
  Scenario: Protocol parameters
    When bulk_system is called with t_high=500, t_low=300, t_step=50 and times_ns={"cool": 0.5}
    Then the cooling files are cool_500K.mdp, cool_450K.mdp, cool_400K.mdp, cool_350K.mdp, cool_300K.mdp
    And each has ref-t 500 (etc.), nsteps 250000 at dt 0.002, and pcoupl C-rescale with ref-p 1.0
    And protocol.yaml lists those temperatures

  @rq-9c761dab
  Scenario: Run script order
    Then run.sh calls grompp and mdrun for 01_em, 02_nvt_hot, 03_npt_compress, 04_npt_anneal and then each cooling step from high to low T, each starting from the previous step's output

  @rq-950969b7
  Scenario: Bilinear Tg fit
    Given synthetic data: density = 1.30 - 0.0002 (T - 400) below 400 K and 1.30 - 0.0006 (T - 400) above, at 600..200 K in 20 K steps
    When fit_tg is called
    Then Tg is 400 K within 1 K

  @rq-1e5b84fd
  Scenario: Tg analysis from energy files
    Given a bulk folder whose cooling steps have density .xvg data from a fake gmx energy
    When analyze_tg(folder) is called
    Then tg_density.csv has one row per temperature, tg_result.yaml has Tg and density_300, and the returned TgResult matches them

  @rq-a0aaaf7b
  Scenario: Missing energy file
    Given a bulk folder where cool_400K.edr is missing
    When analyze_tg(folder) is called
    Then SystemAnalysisError is raised naming cool_400K

  @rq-cca7d81d
  Scenario: Solvent force field mismatch
    Given an OPLS-AA polymer and solvent_ff pointing to amber99sb-ildn.ff
    When bulk_system is called
    Then SystemBuildError is raised about the [ defaults ] mismatch

  @rq-f75e9890
  Scenario: CLI bulk and analyze-tg
    Given the fake backend
    When an input "JOB = bulk / MONOMER = DMAPS / N = 3 / CHAINS = 2 / BOX = 3.0 / FF = fake / PACKER = builtin / OUTPUT = b" is run
    Then it exits 0 and prints the folder and the number of atoms
    And an input "JOB = analyze / OUTPUT = b" exits 1 with an error naming the first missing energy file

  @rq-da3d8837
  Scenario: Bulk density mode
    Given the fake backend
    When bulk_system("DMAPS", n=3, chains=4, density=1.0, initiator="KPS", ff="fake", packer="emc") is called
    Then protocol.yaml has a cubic box of edge (M / 1.0 g/cm3)^(1/3), system.gro is the 0.35 g/cm3 cube
    And 01_deform.mdp shrinks it to that edge, run.sh runs 01_em, 01_heat, 01_deform (from 01_heat.cpt), then 02_nvt_hot
    And giving both density and box_nm raises ValueError
```
