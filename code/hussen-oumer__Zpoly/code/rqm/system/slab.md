# Feature: LLPS Slab System (Direct Coexistence, Temperature Series) <!-- rq-0a2bbac9 -->

Builds the system for liquid–liquid phase separation by the **direct-coexistence (slab) method**:
the chains start as a dense layer in the middle of an elongated box, surrounded by water and salt,
so a dense (coacervate) phase and a dilute phase coexist along z. For UCST-type zwitterionic
polymers, the same equilibrated slab is run at a series of temperatures, and the analysis gives the
dense- and dilute-phase concentrations versus temperature (the phase diagram).

It builds on the shared pieces of `rqm/system/bulk.md` (chains with their own seeds, counterions,
water/ion force-field files, system topology, protocol files, packers).

```python
import zpoly as pp
pp.slab_system("DMAPS", n=12, initiator="KPS", dft="dft", chains=50, box_nm=(8, 8, 24), slab_nm=6,
               out="dmaps_slab")
# bash dmaps_slab/run.sh   (or sbatch dmaps_slab/run_slurm.sh)
pp.analyze_slab("dmaps_slab")     # dense / dilute concentration per temperature
```

## Box and Slab Geometry <!-- rq-5fbd6a40 -->

- **The user decides the system**: the number of chains (`chains`), the box (`box_nm = (Lx, Ly,
  Lz)`) and the thickness of the starting slab (`slab_nm`) are required; the code never chooses
  them. It only checks them (`ValueError`): **Lx = Ly**, **Lz ≥ 3 · Lx**, and 0 < `slab_nm` ≤ Lz.
- It reports the starting slab density (polymer + counterion mass / (Lx · Ly · `slab_nm`)) so the
  user can judge it; if the packer cannot fit the molecules, the error says to make the slab
  thicker, the box wider, or to use Packmol.
- The slab is centred at z = Lz / 2: every atom of the chains and counterions is placed with
  z ∈ [Lz/2 − t/2, Lz/2 + t/2] (both packers).

## Density Mode (`DENSITY`, `ASPECT`) <!-- rq-2c0d3c8f -->

Instead of the box and the slab, the user can give the **density** of the starting dense phase
(`density` in g/cm³, input key `DENSITY`) and `CHAINS`; the geometry follows, keeping the slab
criteria:

- The dense phase is a cube: Lx = Ly = slab thickness = (M / ρ)^(1/3), with M the mass of the chains
  and their counterions; Lz = `aspect` × Lx (`ASPECT`, default 3, at least 3).
- The chains and counterions are packed in a larger cube at `PACK_DENSITY` (0.35 g/cm³, or at ρ if
  that is lower) with the chosen packer (`emc` recommended), and written as `packed.gro`.
- `run.sh` first densifies and builds the slab with GROMACS:
  1. `00_em`: energy minimization of the dry cube;
  2. `00_heat`: 10 ps NVT at `t_equil` that generates the velocities; `00_compress`: NVT continuing
     from it, with `deform` shrinking the cube linearly to exactly Lx in `compress` ns (default 1)
     (no generated velocities with `deform`, so no `deform-init-flow`: any GROMACS version);
  3. `gmx trjconv -pbc mol` (whole molecules), then `gmx editconf -box Lx Lx Lz -center` (the dense
     cube in the middle of the elongated box);
  4. `gmx solvate` and `gmx genion` as for the built-in packer; then the usual protocol.
- `box_nm` and `slab_nm` must not be given with `density` (`ValueError`); `protocol.yaml` records the
  density, the aspect, the packing cube and the resulting box.

## Water and Salt <!-- rq-df0e695a -->

- Water model: TIP3P from the matching solvent force field (`rqm/system/bulk.md`), molecule `SOL`.
- Salt: `salt` = `"NaCl"` (default), `"KCl"` or `None`, at `salt_molar` (default 0.15 M). Ions come
  from the same solvent force field (`NA`, `K`, `CL`). The KPS counterions (K⁺) are added on top,
  so the system is neutral.
- With **`packer="packmol"`** (default), Packmol places everything: the chains and counterions in
  the slab region, then water and salt in the whole box. The number of waters is
  `round(33.4 nm⁻³ × (V_box − V_polymer))` with V_polymer = polymer mass / 1.3 g/cm³; salt pairs are
  `round(salt_molar × N_water / 55.5)`. `system.gro` is then complete.
- With **`packer="builtin"`**, the package packs the chains and counterions; water and salt are added
  by GROMACS in the first step of `run.sh` (`gmx solvate` with `spc216.gro`, then `gmx genion
  -conc salt_molar` with `-pname`/`-nname` for the salt, replacing `SOL`), which also updates
  `system.top`.
- `index.ndx` has a `Polymer` group (all chain atoms, which are always first in the system) for the
  analysis.

## TIP4P/2005 Water and Madrid-2019 Ions (`SOLVENT = tip4p2005`) <!-- rq-a57b1942 -->

`solvent="tip4p2005"` (input key `SOLVENT`, default `tip3p` = the water and ions of the matching
GROMACS force field as above) uses TIP4P/2005 water with the **Madrid-2019** ions (Zeron, Abascal,
Vega, J. Chem. Phys. 151, 134504 (2019)), built in (`system/data/madrid2019.yaml`):

- **Water** `SOL`: OW, HW1, HW2 and the virtual site MW; O–H 0.9572 Å, H–O–H 104.52°, O–M 0.1546 Å,
  q_H = 0.5564, q_M = −1.1128; σ_OO 0.315890 nm, ε_OO 0.774908 kJ/mol; rigid (`[ settles ]`, M from
  `[ virtual_sites3 ]`).
- **Ions** `NA`, `K`, `CL` with charges ±0.85 e and the paper's LJ parameters (Tables III and IV).
- **Ion–ion and ion–water pairs** are written in `[ nonbond_params ]` of the system's force field:
  the paper's σ_ij, ε_ij where it gives them, Lorentz–Berthelot where it says LB (Na–K). Pairs
  between the polymer and water or ions follow the polymer force field's combination rule.
- The ions are not rescaled again; the polymer's ionic groups must use the same scale:
  `CHARGE_SCALE = 0.85` is required (`ValueError` otherwise), so the box is exactly neutral.
- The built-in packer solvates with GROMACS' `tip4p.gro` box; Packmol places 4-site waters.

## Protocol (GROMACS) <!-- rq-0cd7642c -->

Same settings as `rqm/system/bulk.md` (Verlet, PME, 1.2 nm, h-bond constraints, 2 fs, V-rescale,
C-rescale):

| Step | File | What |
|---|---|---|
| (builtin packer only) | in `run.sh` | `gmx solvate`, `gmx genion` |
| 1 | `01_em.mdp` | steepest-descent minimization |
| 2 | `02_nvt.mdp` | NVT at `t_equil` (default 300 K), 1 ns, velocities generated |
| 3 | `03_npt_z.mdp` | NPT at `t_equil`, 1 bar, **semi-isotropic with zero compressibility in x/y**, so only Lz changes and Lx = Ly stays fixed; 5 ns |
| 4… | `prod_<T>K.mdp` | NVT production at each temperature (default 280, 300, 320, 340, 360 K), 200 ns each, every one starting from the end of step 3 with velocities generated at T |

- `run.sh` runs everything in order. `run_production.sh <T>` runs one production after step 3,
  so temperatures can run as separate jobs. `run_slurm.sh` runs `run.sh`.
- `times_ns` overrides `nvt`, `npt_z` and `production`; `temperatures` sets the series.
- `protocol.yaml` records the geometry, salt, temperatures and steps.

## Analysis <!-- rq-1dce521e -->

`analyze_slab(folder, *, skip_ns=50.0, bin_nm=0.1)`, for each `prod_<T>K`:

- Runs `gmx density` on `prod_<T>K.xtc` with the `Polymer` group, mass density, `-center` on that
  group, `-sl` = Lz / `bin_nm`, from `skip_ns` on, giving ρ(z) in kg/m³ (= g/L). A missing
  trajectory raises `SystemAnalysisError` naming the step.
- **Slab**: where ρ(z) is above half of its maximum (the longest such stretch around the centre);
  its length is the slab width.
- **Dense** concentration: the mean of ρ over the middle half of the slab. **Dilute**
  concentration: the mean of ρ over the bins more than 1.5 nm away from the slab, on both sides.
  If there are no such bins (the slab fills the box), `dilute` is `None` and a warning says the box
  is too short.
- Writes `density_profile_<T>K.csv` (z, ρ) and `phase_diagram.csv` (T, dense, dilute, width), and
  returns a `SlabResult` per temperature.

## Feature API <!-- rq-32d1505a -->

Importable from `zpoly` (and `zpoly.system`):

- `slab_system(monomer, n=None, *, chains, box_nm=None, slab_nm=None, density=None, aspect=3.0, <!-- rq-c2048297 -->
  initiator=None, ends=None, spacer=None, <!-- rq-e30b89e7 -->
  alkyl=None, dft=None, ff="gaff2", packer="packmol",
  salt="NaCl", salt_molar=0.15, t_equil=300.0, temperatures=(280, 300, 320, 340, 360),
  times_ns=None, seed=0, solvent_ff=None, out="slab", library=None) -> Path`
- `analyze_slab(folder, *, skip_ns=50.0, bin_nm=0.1, gmx="gmx") -> list[SlabResult]`; <!-- rq-59e62eef -->
  `SlabResult` has `temperature`, `dense`, `dilute` (g/L, `None` if undefined), `width_nm`,
  `profile` (tuple of (z nm, ρ)).
- `phase_concentrations(z, rho) -> (dense, dilute, width)` (the profile analysis alone). <!-- rq-840be29e -->

## Input File <!-- rq-72006b7b -->

`JOB = slab` with `CHAINS`, `BOX = Lx Ly Lz`, `SLAB` (required), `SALT = NaCl|KCl [molar] | none`,
`TEMPERATURES`, `TEQUIL`, `TIMES`, `PACKER`, `SOLVENT_FF`, `DFT`, `FF` and `OUTPUT`
(`rqm/input_file.md`). `JOB = analyze` with `OUTPUT = <folder>` runs `analyze_slab` on a slab folder.

## Out of Scope <!-- rq-4dc65b4e -->

- Free-energy or critical-point fitting of the phase diagram; non-rectangular boxes; coarse-grained
  slabs.

## Gherkin Scenarios <!-- rq-147b01b7 -->

```gherkin
Feature: LLPS slab system

  @rq-70bb800e
  Scenario: Box with Lz three times Lx
    Given the fake backend and the fixture solvent force field
    When slab_system("DMAPS", n=3, chains=6, box_nm=(3.4, 3.4, 10.2), slab_nm=3.4, initiator="KPS", ff="fake", packer="builtin", out=tmp/"s") is called
    Then the box in packed.gro is 3.4 x 3.4 x 10.2 nm

  @rq-185b551d
  Scenario: Box checks
    When slab_system is called with box_nm=(3.4, 3.4, 8.0) (Lz < 3 Lx), or box_nm=(3.4, 4.0, 12.0) (Lx ≠ Ly), or slab_nm=11
    Then ValueError is raised each time

  @rq-9b09c717
  Scenario: Chains start as a centred slab
    Given the slab system above
    Then every polymer and counterion atom has z within [Lz/2 - slab_nm/2, Lz/2 + slab_nm/2] (within the .gro precision)

  @rq-a6cd2763
  Scenario: Missing user inputs
    When slab_system is called without chains, box_nm or slab_nm
    Then TypeError is raised

  @rq-fee4e59e
  Scenario: Built-in packer leaves water and salt to GROMACS
    Given the slab system with packer="builtin"
    Then run.sh starts with gmx solvate and gmx genion -pname NA -nname CL -conc 0.15 before 01_em
    And system.top lists polymer 6 and K 12 before solvation

  @rq-9dfe556a
  Scenario: Packmol packs water and salt too
    Given packmol is on PATH
    When slab_system with packer="packmol" is called
    Then system.top lists polymer, K, SOL, NA and CL with the counts from the formulas
    And system.gro has all of those atoms and run.sh does not call solvate or genion

  @rq-1223ec9e
  Scenario: No salt
    When slab_system(..., salt=None, packer="builtin") is called
    Then run.sh solvates but does not add salt (genion is not called because the system is neutral)

  @rq-cb81b1ce
  Scenario: KCl salt
    When slab_system(..., salt="KCl", packer="builtin") is called
    Then genion uses -pname K -nname CL

  @rq-7bc134f2
  Scenario: Production temperatures
    When slab_system(..., temperatures=(290, 310)) is called
    Then prod_290K.mdp and prod_310K.mdp exist, are NVT with gen-vel yes at their temperature,
    And run.sh starts each from 03_npt_z.gro, and run_production.sh exists

  @rq-133ab185
  Scenario: z-only pressure coupling
    Then 03_npt_z.mdp has pcoupltype semiisotropic and compressibility "0 4.5e-5"

  @rq-03fc94cd
  Scenario: Index group
    Then index.ndx has a [ Polymer ] group with exactly the chain atoms (the first atoms of the system)

  @rq-d067b409
  Scenario: Phase concentrations from a profile
    Given a profile of 1000 g/L for |z - 15| < 3 nm and 20 g/L elsewhere on 0..30 nm, smoothed with tanh walls of 0.5 nm
    When phase_concentrations is called
    Then dense is 1000 within 2%, dilute is 20 within 2%, and the width is 6 nm within 0.3 nm

  @rq-4b2d638c
  Scenario: Slab analysis with gmx density
    Given a slab folder whose productions have trajectories and a fake gmx density
    When analyze_slab(folder) is called
    Then phase_diagram.csv has one row per temperature and a density_profile_<T>K.csv per temperature

  @rq-f0a89383
  Scenario: Missing trajectory
    Given a slab folder where prod_300K.xtc is missing
    When analyze_slab(folder) is called
    Then SystemAnalysisError is raised naming prod_300K

  @rq-42bf805c
  Scenario: CLI slab
    Given the fake backend
    When an input "JOB = slab / MONOMER = DMAPS / N = 3 / CHAINS = 4 / BOX = 3.4 3.4 10.2 / SLAB = 3.4 / FF = fake / PACKER = builtin / TEMPERATURES = 300 320 / OUTPUT = s" is run
    Then it exits 0 and s/prod_300K.mdp and s/prod_320K.mdp exist

  @rq-96c58164
  Scenario: TIP4P/2005 and Madrid-2019
    Given the fake backend
    When slab_system(..., solvent="tip4p2005", charge_scale=0.85, salt="NaCl", packer="builtin") is called
    Then the force field defines OW/HW/MW and NA/K/CL types, ions.itp has K +0.85 and CL -0.85, water.itp has a virtual site
    And [ nonbond_params ] lists K-OW 0.289040 1.400430, NA-CL 0.300512 1.438894 and NA-K by Lorentz-Berthelot
    And run.sh solvates with tip4p.gro, and the system's total charge is 0
    And with charge_scale=1.0 it raises ValueError

  @rq-7c268316
  Scenario: Density mode
    Given the fake backend
    When slab_system("DMAPS", n=3, chains=4, density=1.0, aspect=3, initiator="KPS", ff="fake", packer="emc") is called
    Then protocol.yaml has a box with Lx = Ly = (M / 1.0 g/cm3)^(1/3) and Lz = 3 Lx, and packed.gro is a cube at 0.35 g/cm3
    And 00_compress.mdp has a deform rate that shrinks the cube to Lx in 1 ns
    And run.sh runs 00_em, 00_compress, trjconv -pbc mol, editconf -box Lx Lx Lz -center, solvate and genion before 01_em
    And giving density together with box_nm, or aspect below 3, raises ValueError
```
