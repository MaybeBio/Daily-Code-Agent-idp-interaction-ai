# Sim Class & System Building

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [calvados/components.py](calvados/components.py)
- [calvados/sim.py](calvados/sim.py)

</details>



## Purpose and Scope

This page documents the `Sim` class and its `build_system` method, which serve as the central orchestration layer for CALVADOS simulations. The `Sim` class is responsible for:

- Reading configuration and component specifications from YAML files
- Instantiating molecular components based on their types
- Initializing force field parameters and OpenMM force objects
- Placing molecules in the simulation box according to specified topologies
- Assembling bonded and non-bonded interactions
- Building the complete OpenMM System ready for simulation

For information about configuring simulation parameters, see [Config Class - Simulation Parameters](#2.1). For details on molecular component definitions, see [Components Class - Molecular Definitions](#2.2). For the component class hierarchy itself, see [Component Class Hierarchy](#3.1). For molecule placement algorithms, see [Molecule Placement Strategies](#4.2).

---

## Sim Class Overview

The `Sim` class is defined in [calvados/sim.py:20-639]() and serves as the main simulation orchestrator. It coordinates the entire workflow from configuration parsing to system construction to MD execution.

### Key Attributes

| Attribute | Type | Description |
|-----------|------|-------------|
| `path` | str | Working directory for simulation outputs |
| `box` | np.ndarray | Box dimensions [x, y, z] in nm |
| `system` | openmm.System | OpenMM System object containing all forces and particles |
| `top` | md.Topology | MDTraj topology for trajectory I/O |
| `components` | np.ndarray | Array of Component objects (Protein, RNA, Lipid, etc.) |
| `ah`, `yu` | CustomNonbondedForce | Ashbaugh-Hatch and Yukawa force objects |
| `nparticles` | int | Total number of beads in the system |
| `pos` | list | Positions of all particles (built incrementally) |

### Key Methods

| Method | Purpose |
|--------|---------|
| `__init__` | Parse configuration, initialize attributes |
| `make_components` | Instantiate Component objects based on molecule_type |
| `build_system` | Main workflow: assemble complete OpenMM System |
| `add_forces_to_system` | Add all force objects to the OpenMM System |
| `add_interactions` | Add interactions for a single component |
| `place_molecule` | Place a molecule in the simulation box |
| `simulate` | Execute the MD simulation |

**Sources:** [calvados/sim.py:20-639]()

---

## Sim Class Architecture

```mermaid
graph TB
    subgraph "Sim Class Structure"
        INIT["__init__<br/>(config, components)"]
        MAKE["make_components()<br/>Instantiate Component objects"]
        BUILD["build_system()<br/>Main assembly workflow"]
        FORCES["add_forces_to_system()<br/>Add forces to OpenMM System"]
        SIM["simulate()<br/>Execute MD simulation"]
    end
    
    subgraph "Component Creation"
        MAKE_PROT["Protein instances"]
        MAKE_RNA["RNA instances"]
        MAKE_LIP["Lipid instances"]
        MAKE_CROWD["Crowder instances"]
    end
    
    subgraph "System Building Steps"
        INIT_TOP["Initialize md.Topology<br/>& openmm.System"]
        INIT_PARAMS["genParamsDH()<br/>Debye-Hückel parameters"]
        INIT_NB["init_nonbonded_interactions()<br/>AH & YU forces"]
        PLACE["place_molecule()<br/>Position molecules"]
        ADD_INT["add_interactions()<br/>Bonds, restraints, NB params"]
    end
    
    INIT --> MAKE
    INIT --> BUILD
    BUILD --> MAKE
    BUILD --> INIT_TOP
    BUILD --> INIT_PARAMS
    BUILD --> INIT_NB
    BUILD --> PLACE
    BUILD --> ADD_INT
    BUILD --> FORCES
    FORCES --> SIM
    
    MAKE --> MAKE_PROT
    MAKE --> MAKE_RNA
    MAKE --> MAKE_LIP
    MAKE --> MAKE_CROWD
```

**Sources:** [calvados/sim.py:20-649]()

---

## Initialization: `__init__` Method

The `Sim.__init__` method [calvados/sim.py:21-49]() parses the configuration and components dictionaries loaded from YAML files.

### Configuration Processing

```python
def __init__(self, path, config, components):
    self.path = path
    # Parse config dictionary into attributes
    for key, val in config.items():
        setattr(self, key, val)
    
    # Parse component defaults
    for key, val in components['defaults'].items():
        setattr(self, f'default_{key}', val)
    
    self.comp_dict = components['system']
    self.comp_defaults = components['defaults']
```

Key operations:
1. **Path Setup**: Sets working directory for all output files
2. **Config Parsing**: Converts all config.yaml parameters to `Sim` attributes
3. **Component Defaults**: Stores default component properties with `default_` prefix
4. **Unit Conversion**: Converts `eps_lj` from kcal/mol to kJ/mol (line 38)
5. **Equilibration Setup**: Initializes slab or external force restraints if enabled (lines 40-48)

### Restart Logic

| Condition | Behavior |
|-----------|----------|
| `restart == 'checkpoint'` and checkpoint file exists | Disable slab_eq and bilayer_eq (continue from checkpoint) |
| `slab_eq == True` | Initialize `rcent` force for slab equilibration |
| `ext_force == True` | Initialize custom external force with user expression |

**Sources:** [calvados/sim.py:21-49]()

---

## Component Creation: `make_components` Method

The `make_components` method [calvados/sim.py:50-98]() instantiates molecular components based on their `molecule_type` property.

### Component Type Decision Tree

```mermaid
graph TD
    START["For each component in comp_dict"]
    GET_TYPE["Get molecule_type<br/>(or use default_molecule_type)"]
    
    PROTEIN{"molecule_type<br/>== 'protein'?"}
    RNA{"molecule_type<br/>== 'rna'?"}
    LIPID{"molecule_type<br/>in ['lipid', 'cooke_lipid']?"}
    CROWDER{"molecule_type<br/>== 'crowder'?"}
    CYCLIC{"molecule_type<br/>== 'cyclic'?"}
    SEASTAR{"molecule_type<br/>== 'seastar'?"}
    PTM{"molecule_type<br/>== 'ptm_protein'?"}
    
    PROT_CLASS["Protein(name, properties, defaults)<br/>comp_setup = 'compact'"]
    RNA_CLASS["RNA(name, properties, defaults)<br/>comp_setup = 'spiral'"]
    LIP_CLASS["Lipid(name, properties, defaults)<br/>comp_setup = 'linear'"]
    CROWD_CLASS["Crowder(name, properties, defaults)<br/>comp_setup = 'compact'"]
    CYC_CLASS["Cyclic(name, properties, defaults)<br/>comp_setup = 'compact'"]
    SEAS_CLASS["Seastar(name, properties, defaults)<br/>comp_setup = 'compact'"]
    PTM_CLASS["PTMProtein(name, properties, defaults)<br/>comp_setup = 'compact'"]
    GEN_CLASS["Component(name, properties, defaults)<br/>comp_setup = 'linear'"]
    
    CALC_PROPS["comp.calc_properties(pH, verbose, comp_setup)"]
    INIT_RESTR["comp.init_restraint_force()<br/>(if restraint == True)"]
    APPEND["Append comp to self.components"]
    
    START --> GET_TYPE
    GET_TYPE --> PROTEIN
    PROTEIN -->|Yes| PROT_CLASS
    PROTEIN -->|No| RNA
    RNA -->|Yes| RNA_CLASS
    RNA -->|No| LIPID
    LIPID -->|Yes| LIP_CLASS
    LIPID -->|No| CROWDER
    CROWDER -->|Yes| CROWD_CLASS
    CROWDER -->|No| CYCLIC
    CYCLIC -->|Yes| CYC_CLASS
    CYCLIC -->|No| SEASTAR
    SEASTAR -->|Yes| SEAS_CLASS
    SEASTAR -->|No| PTM
    PTM -->|Yes| PTM_CLASS
    PTM -->|No| GEN_CLASS
    
    PROT_CLASS --> CALC_PROPS
    RNA_CLASS --> CALC_PROPS
    LIP_CLASS --> CALC_PROPS
    CROWD_CLASS --> CALC_PROPS
    CYC_CLASS --> CALC_PROPS
    SEAS_CLASS --> CALC_PROPS
    PTM_CLASS --> CALC_PROPS
    GEN_CLASS --> CALC_PROPS
    
    CALC_PROPS --> INIT_RESTR
    INIT_RESTR --> APPEND
```

### Component Property Calculation

After instantiation, each component's `calc_properties` method is called [calvados/sim.py:87]():

```python
comp.eps_lj = self.eps_lj  # Share LJ epsilon with component
comp.calc_properties(pH=self.pH, verbose=self.verbose, comp_setup=comp_setup)
```

This triggers:
1. Sequence reading from FASTA or PDB
2. Property calculation (sigmas, lambdas, charges, masses)
3. Initial coordinate generation (spiral, linear, or compact)
4. Bond force initialization
5. For proteins with restraints: PDB loading, distance map calculation, Go-model scaling

### Restraint Force Initialization

For components with `restraint == True` [calvados/sim.py:88-96]():

```python
if comp.restraint:
    if comp.restraint_type == 'go':
        comp.init_restraint_force(
            eps_lj=self.eps_lj, cutoff_lj=self.cutoff_lj,
            eps_yu=self.eps_yu, k_yu=self.k_yu
        )
    else:
        comp.init_restraint_force()
    self.use_restraints = True
```

Go-model restraints require additional force objects for scaled LJ and YU interactions (see [Restraints System](#3.5)).

**Sources:** [calvados/sim.py:50-98](), [calvados/components.py:11-108]()

---

## System Building Workflow: `build_system` Method

The `build_system` method [calvados/sim.py:130-223]() is the main orchestration function that assembles the complete OpenMM System.

### Build System Sequence Diagram

```mermaid
sequenceDiagram
    participant Sim as "Sim.build_system()"
    participant Top as "md.Topology"
    participant Sys as "openmm.System"
    participant Int as "interactions module"
    participant Comp as "Component objects"
    
    Sim->>Top: Initialize empty topology
    Sim->>Sys: Initialize empty System
    Sim->>Sys: setDefaultPeriodicBoxVectors(a, b, c)
    
    Sim->>Int: genParamsDH(temp, ionic)
    Int-->>Sim: eps_yu, k_yu
    
    Sim->>Sim: make_components()
    Sim->>Sim: count_components()
    
    Sim->>Int: init_nonbonded_interactions(eps_lj, cutoff_lj, eps_yu, k_yu, ...)
    Int-->>Sim: ah, yu (force objects)
    
    alt nlipids > 0 or ncookelipids > 0
        Sim->>Int: init_lipid_interactions(eps_lj, eps_yu, cutoff_yu, factor)
        Int-->>Sim: cos, cn (force objects)
    end
    
    loop For each component and molecule copy
        Sim->>Top: add_mdtraj_topol(comp)
        Sim->>Sys: add_particles_system(comp.mws)
        Sim->>Sim: place_molecule(comp) or place_bilayer(comp)
        Sim->>Sim: add_interactions(comp)
        
        alt slab_eq or ext_force
            Sim->>Sim: add_ext_restraints(comp)
        end
    end
    
    alt custom_restraints
        Sim->>Sim: map_custom_restraints()
        Sim->>Sim: add_custom_restraints()
    end
    
    Sim->>Top: Save trajectory topology to top.pdb
    Sim->>Sim: add_forces_to_system()
    Sim->>Sim: print_system_summary()
```

### Workflow Steps

| Step | Line Range | Description |
|------|------------|-------------|
| 1. Initialize Topology & System | 139-142 | Create empty `md.Topology` and `openmm.System`, set periodic box vectors |
| 2. Calculate DH Parameters | 146 | Call `genParamsDH(temp, ionic)` to get electrostatic parameters |
| 3. Make Components | 149 | Call `make_components()` to instantiate all Component objects |
| 4. Count Components | 150 | Organize components by type, validate topology compatibility |
| 5. Initialize Non-bonded Forces | 153-165 | Create AH and YU force objects; create lipid forces if needed |
| 6. Place Molecules | 172-211 | Loop over components, place molecules, add to topology and system |
| 7. Add Interactions | 207 | Add bonded/non-bonded parameters for each molecule |
| 8. Custom Restraints | 213-215 | Parse and add user-defined restraints if specified |
| 9. Save Topology | 217-220 | Save initial configuration as `top.pdb` |
| 10. Add Forces | 222 | Call `add_forces_to_system()` to attach all forces |
| 11. Summary | 223 | Print system statistics and save XML |

**Sources:** [calvados/sim.py:130-223]()

---

## Force Initialization

Force initialization occurs in two stages: global force object creation and per-component parameter addition.

### Global Force Object Creation

```mermaid
graph LR
    subgraph "Force Initialization Calls"
        GEN["genParamsDH(temp, ionic)"]
        INIT_NB["init_nonbonded_interactions(...)"]
        INIT_LIP["init_lipid_interactions(...)<br/>(if nlipids > 0)"]
    end
    
    subgraph "Force Objects Created"
        EPS_YU["eps_yu, k_yu<br/>(Debye-Hückel parameters)"]
        AH["self.ah<br/>(CustomNonbondedForce)<br/>Ashbaugh-Hatch"]
        YU["self.yu<br/>(CustomNonbondedForce)<br/>Yukawa"]
        COS["self.cos<br/>(CustomNonbondedForce)<br/>Cosine attraction"]
        CN["self.cn<br/>(CustomNonbondedForce)<br/>Charge-Nonpolar"]
    end
    
    GEN --> EPS_YU
    INIT_NB --> AH
    INIT_NB --> YU
    INIT_LIP --> COS
    INIT_LIP --> CN
```

#### Debye-Hückel Parameter Generation

[calvados/sim.py:146]() calls `interactions.genParamsDH(self.temp, self.ionic)`:

```python
self.eps_yu, self.k_yu = interactions.genParamsDH(self.temp, self.ionic)
```

This calculates:
- `eps_yu`: Prefactor for Yukawa potential based on temperature and dielectric constant
- `k_yu`: Inverse Debye screening length based on ionic strength

#### Non-bonded Force Initialization

[calvados/sim.py:153-154]() creates Ashbaugh-Hatch and Yukawa force objects:

```python
self.ah, self.yu = interactions.init_nonbonded_interactions(
    self.eps_lj, self.cutoff_lj, self.eps_yu, self.k_yu, self.cutoff_yu, self.fixed_lambda
)
```

These `CustomNonbondedForce` objects define the functional forms but contain no particles initially.

#### Lipid Force Initialization

If lipids are present [calvados/sim.py:156-165](), additional forces are created:

```python
if self.nlipids > 0:
    self.cos, self.cn = interactions.init_lipid_interactions(
        self.eps_lj, self.eps_yu, self.cutoff_yu, factor=1.9
    )
elif self.ncookelipids > 0:
    self.cos, self.cn = interactions.init_lipid_interactions(
        self.eps_lj, self.eps_yu, self.cutoff_yu, factor=3.0
    )
```

- `cos`: Cosine attraction between tail beads
- `cn`: Charge-nonpolar interaction between charged and nonpolar beads
- `factor`: Controls strength relative to eps_lj (1.9 for standard lipids, 3.0 for Cooke lipids)

**Sources:** [calvados/sim.py:146-165](), [calvados/interactions.py]()

---

## Particle and Topology Setup

For each molecule of each component, the system adds particles to both the MDTraj topology and the OpenMM System.

### Topology Grid Initialization

Before placing molecules, spatial grids are initialized based on topology type [calvados/sim.py:172-190]():

| Topology Type | Grid Initialization |
|---------------|---------------------|
| `'slab'` | `xyzgrid` for proteins/RNA in central slab region; additional grids for crowders in outer regions |
| `'grid'` | `xyzgrid` covering entire box for all molecules |
| Lipids present | `bilayergrid` for xy plane; `xyzgrid` for proteins/RNA in solution regions |

### Per-Molecule Loop

The main assembly loop [calvados/sim.py:192-211]():

```python
for cidx, comp in enumerate(self.components):
    for idx in range(comp.nmol):
        # 1. Add to MDTraj topology
        self.add_mdtraj_topol(comp)
        
        # 2. Add particles to OpenMM System
        self.add_particles_system(comp.mws)
        
        # 3. Place molecule in space
        if comp.molecule_type in ['protein','crowder','cyclic','seastar','ptm_protein']:
            xs = self.place_molecule(comp)
        elif comp.molecule_type in ['lipid','cooke_lipid']:
            xs = self.place_bilayer(comp)
        elif comp.molecule_type == 'rna':
            xs = self.place_molecule(comp)
        
        # 4. Add interactions
        self.add_interactions(comp)
        
        # 5. Add external restraints (if equilibrating)
        if (self.slab_eq or self.ext_force) and comp.ext_restraint:
            self.add_ext_restraints(comp)
```

### MDTraj Topology Construction

`add_mdtraj_topol(comp)` [calvados/sim.py:431-454]() creates a new chain in the topology with appropriate residues and bonds:

For **proteins/crowders**:
- Each residue → one CA atom
- Bonds added based on `comp.bond_check(i, i+1)`

For **RNA**:
- Each nucleotide → two atoms (phosphate "P" and base "N")
- Bonds added based on `comp.bond_check(i, j)`

### OpenMM Particle Addition

`add_particles_system(mws)` [calvados/sim.py:456-460]() adds particles with masses:

```python
for mw in mws:
    self.system.addParticle(mw*unit.amu)
```

Each particle is added sequentially, incrementing `self.nparticles`.

**Sources:** [calvados/sim.py:172-211, 292-336, 431-460]()

---

## Interaction Assembly: `add_interactions` Method

The `add_interactions` method [calvados/sim.py:378-422]() adds all interactions for a single molecule copy.

### Interaction Assembly Flow

```mermaid
graph TB
    START["add_interactions(comp)"]
    OFFSET["Calculate offset<br/>offset = nparticles - comp.nbeads"]
    
    subgraph "Non-bonded Parameters"
        ADD_AH["Add AH parameters<br/>self.ah.addParticle([sigma, lambda, type])"]
        ADD_YU["Add YU parameters<br/>self.yu.addParticle([q])"]
        ADD_LIP_COS["Add lipid cosine parameters<br/>self.cos.addParticle([sigma, lambda, type])"]
        ADD_LIP_CN["Add charge-nonpolar parameters<br/>self.cn.addParticle([sig^3, alpha, q, id])"]
    end
    
    subgraph "Bonded Interactions"
        ADD_BONDS["add_bonds(comp, offset)<br/>Add harmonic bonds<br/>Add exclusions"]
        ADD_ANGLES["add_angles(comp, offset)<br/>(RNA only)"]
        ADD_RESTR["add_restraints(comp, offset)<br/>(if comp.restraint)"]
    end
    
    WRITE["Write bond/restraint lists<br/>(if verbose)"]
    
    START --> OFFSET
    OFFSET --> ADD_AH
    ADD_AH --> ADD_YU
    
    ADD_YU --> ADD_LIP_COS
    ADD_YU --> ADD_LIP_CN
    ADD_LIP_COS --> ADD_BONDS
    ADD_LIP_CN --> ADD_BONDS
    
    ADD_BONDS --> ADD_ANGLES
    ADD_ANGLES --> ADD_RESTR
    ADD_RESTR --> WRITE
```

### Non-bonded Parameter Addition

#### Ashbaugh-Hatch Parameters

[calvados/sim.py:385-397]() loops over component beads:

```python
for sig, lam in zip(comp.sigmas, comp.lambdas):
    if comp.molecule_type in ['lipid', 'cooke_lipid']:
        self.ah.addParticle([sig*unit.nanometer, lam, 0])  # type=0
    elif comp.molecule_type == 'crowder':
        self.ah.addParticle([sig*unit.nanometer, lam, -1]) # type=-1
    else:  # protein, RNA
        self.ah.addParticle([sig*unit.nanometer, lam, 1])  # type=1
```

The `type` parameter controls interaction selectivity:
- `type=1`: Proteins/RNA (fully interacting)
- `type=0`: Lipids (selective interactions)
- `type=-1`: Crowders (repulsive only)

#### Yukawa Electrostatic Parameters

[calvados/sim.py:399-400]():

```python
for q in comp.qs:
    self.yu.addParticle([q])
```

Each bead's charge is added directly.

#### Lipid-Specific Parameters

If lipids are present [calvados/sim.py:403-406]():

```python
if self.nlipids > 0 or self.ncookelipids > 0:
    id_cn = 1 if comp.molecule_type == 'protein' else -1
    for sig, alpha, q in zip(comp.sigmas, comp.alphas, comp.qs):
        self.cn.addParticle([(sig/2)**3, alpha, q, id_cn])
```

### Bonded Interaction Addition

#### Bonds

`add_bonds(comp, offset)` [calvados/sim.py:338-342]() calls `comp.add_bonds(offset)`:

1. Component determines bonded pairs via `comp.bond_check(i, j)`
2. For each bonded pair, adds to `comp.hb` (HarmonicBondForce)
3. Returns exclusion map for non-bonded forces
4. `add_exclusions(exclusion_map)` excludes bonded pairs from AH and YU

#### Angles (RNA only)

[calvados/sim.py:411-412]():

```python
if comp.molecule_type == 'rna':
    self.add_angles(comp, offset)
```

RNA components have angle forces between phosphate beads.

#### Restraints

[calvados/sim.py:415-416]():

```python
if comp.restraint:
    self.add_restraints(comp, offset)
```

Calls `comp.add_restraints(offset)` which adds:
- **Harmonic restraints**: Between residue pairs in structured domains
- **Go-model restraints**: Scaled by AlphaFold confidence and PAE
- **Scaled LJ/YU**: For weakly restrained pairs in Go-model

**Sources:** [calvados/sim.py:378-422](), [calvados/components.py:85-96, 218-229, 240-272]()

---

## Adding Forces to System: `add_forces_to_system` Method

After all components are assembled, `add_forces_to_system()` [calvados/sim.py:225-271]() adds all force objects to the OpenMM System.

### Force Addition Order

```mermaid
graph TD
    START["add_forces_to_system()"]
    
    subgraph "Intermolecular Forces"
        ADD_YU["system.addForce(self.yu)"]
        ADD_AH["system.addForce(self.ah)"]
        ADD_LIP["system.addForce(self.cos)<br/>system.addForce(self.cn)<br/>(if lipids present)"]
    end
    
    subgraph "Intramolecular Forces"
        GET_FORCES["comp.get_forces()<br/>Populates comp.forces list"]
        ADD_COMP["system.addForce(force)<br/>for force in comp.forces"]
    end
    
    subgraph "External & Equilibration Forces"
        ADD_EXT["system.addForce(self.rcent)<br/>(if ext_force)"]
        ADD_SLAB["system.addForce(self.rcent)<br/>(if slab_eq)"]
        ADD_CUST["system.addForce(self.cres)<br/>(if custom_restraints)"]
        ADD_BARO["system.addForce(barostat)<br/>(if box_eq or bilayer_eq)"]
    end
    
    START --> ADD_YU
    ADD_YU --> ADD_AH
    ADD_AH --> ADD_LIP
    ADD_LIP --> GET_FORCES
    GET_FORCES --> ADD_COMP
    ADD_COMP --> ADD_EXT
    ADD_EXT --> ADD_SLAB
    ADD_SLAB --> ADD_CUST
    ADD_CUST --> ADD_BARO
```

### Component Forces

Each component's `get_forces()` method [calvados/components.py:98-99, 293-298]() populates its `forces` list:

**Protein**:
- `comp.hb` (HarmonicBondForce)
- `comp.cs` (restraint force, if `restraint=True`)
- `comp.scLJ`, `comp.scYU` (scaled forces for Go-model)

**RNA**:
- `comp.hb` (bonds between phosphate-phosphate and phosphate-base)
- `comp.scLJ_rna` (base-base interactions)
- `comp.ha` (HarmonicAngleForce for phosphate angles)
- `comp.cs` (restraints, if `restraint=True`)

**Lipid**:
- `comp.hb` (bonds or WCA-Fene bonds)
- `comp.ha` (angles for branching, if `molecule_type=='lipid'`)

### Barostat Forces

For pressure equilibration [calvados/sim.py:257-271]():

| Equilibration Type | Barostat |
|--------------------|----------|
| `box_eq=True` | `MonteCarloAnisotropicBarostat` with directional scaling |
| `bilayer_eq=True` | `MonteCarloMembraneBarostat` with XY isotropic, Z fixed |

**Sources:** [calvados/sim.py:225-271](), [calvados/components.py:98-99, 293-298, 363-366]()

---

## System Summary and Output

After force addition, `print_system_summary()` [calvados/sim.py:273-290]() outputs system statistics and optionally writes the System to XML format.

### System Statistics Output

```
{nparticles} particles in the system
---------- FORCES ----------
ah: {n_particles} particles, {n_exclusions} exclusions
yu: {n_particles} particles, {n_exclusions} exclusions
```

Additional output for special cases:
- Slab equilibration: Number of particles with external restraints
- Bilayer equilibration: Notification of zero lateral tension
- Box equilibration: Which axes (X/Y/Z) are scaling
- Custom restraints: Number of custom restraint bonds
- Component restraints: Number of restraints per component

### XML Serialization

[calvados/sim.py:277-278]():

```python
with open(f'{self.path}/{self.sysname}.xml', 'w') as output:
    output.write(openmm.XmlSerializer.serialize(self.system))
```

This saves the complete OpenMM System including all forces, particles, and parameters for inspection or reloading.

**Sources:** [calvados/sim.py:273-290]()

---

## Summary: Data Flow Through build_system

The following table summarizes the complete data flow:

| Stage | Input | Processing | Output |
|-------|-------|------------|--------|
| **Initialization** | config.yaml, components.yaml | Parse YAML → Sim attributes | Configured Sim object |
| **Component Creation** | comp_dict, comp_defaults | Instantiate by molecule_type, calc_properties() | Array of Component objects |
| **Force Initialization** | temp, ionic, eps_lj, cutoffs | genParamsDH(), init_*_interactions() | Force objects (ah, yu, cos, cn) |
| **Topology Setup** | Components | Initialize empty md.Topology and openmm.System | Empty topology and system |
| **Molecule Placement** | Component initial coords, topology type | place_molecule() or place_bilayer() | Positioned molecules in self.pos |
| **Interaction Assembly** | Component parameters | add_interactions() → addParticle() for forces | Populated force objects with parameters |
| **Bonded Terms** | Component bond_check() logic | add_bonds(), add_angles(), add_restraints() | Bonded forces and exclusions |
| **Force Addition** | All force objects | add_forces_to_system() | Complete OpenMM System ready for simulation |
| **Topology Output** | self.pos, self.top | MDTraj save to top.pdb | Initial configuration PDB file |

**Sources:** [calvados/sim.py:130-223]()

---