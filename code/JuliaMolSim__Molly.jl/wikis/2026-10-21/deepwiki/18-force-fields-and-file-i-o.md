# Force Fields and File I/O

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [src/setup.jl](src/setup.jl)

</details>



This page covers the loading and parsing of molecular force fields from XML files and the reading of coordinate data to set up molecular systems. It focuses on the file I/O capabilities that allow Molly.jl to interface with standard molecular simulation file formats, particularly OpenMM force field XML files and various coordinate formats supported by Chemfiles.

For information about the `System` data structure itself, see [System and Atoms](#3.1). For details about how interactions are calculated once loaded, see [Interaction System](#4). For unit handling in file I/O, see [Units and Physical Quantities](#7.2).

## Overview of Force Field and File I/O System

Molly.jl provides comprehensive support for reading molecular systems from files, with two main pathways for system setup:

```mermaid
graph TD
    XMLFiles["OpenMM XML Files<br/>(.xml)"] --> MFF["MolecularForceField"]
    CoordFiles["Coordinate Files<br/>(PDB, XYZ, etc.)"] --> System["System"]
    MFF --> System
    
    GromacsTOP["Gromacs Topology<br/>(.top)"] --> System
    GromacsGRO["Gromacs Coordinates<br/>(.gro)"] --> System
    
    RandomCoords["place_atoms()<br/>place_diatomics()"] --> System
    
    System --> Simulation["Ready for Simulation"]
    
    style MFF fill:#f9f9f9
    style System fill:#e1f5fe
```

The system supports two primary workflows: the OpenMM-compatible pathway using XML force fields with Chemfiles-readable coordinate formats, and the Gromacs pathway reading native Gromacs topology and coordinate files.

Sources: [src/setup.jl:1-1478]()

## Force Field Data Structures

### MolecularForceField

The `MolecularForceField` type serves as the central container for all force field parameters parsed from OpenMM XML files:

```mermaid
graph LR
    subgraph "MolecularForceField Components"
        AtomTypes["atom_types<br/>Dict{String, AtomType}"]
        ResidueTypes["residue_types<br/>Dict{String, ResidueType}"] 
        BondTypes["bond_types<br/>Dict{Tuple, HarmonicBond}"]
        AngleTypes["angle_types<br/>Dict{Tuple, HarmonicAngle}"]
        TorsionTypes["torsion_types<br/>Dict{Tuple, PeriodicTorsionType}"]
    end
    
    subgraph "Force Field Parameters"
        Weights["weight_14_coulomb<br/>weight_14_lj"]
        Order["torsion_order"]
        Attributes["attributes_from_residue"]
    end
    
    AtomTypes --> MFF["MolecularForceField"]
    ResidueTypes --> MFF
    BondTypes --> MFF
    AngleTypes --> MFF
    TorsionTypes --> MFF
    Weights --> MFF
    Order --> MFF
    Attributes --> MFF
```

The `AtomType` struct contains fundamental atomic parameters:

| Field | Type | Description |
|-------|------|-------------|
| `type` | `String` | Atom type identifier |
| `class` | `String` | Atom class (currently unused) |
| `element` | `String` | Chemical element symbol |
| `charge` | `Union{T, Missing}` | Partial charge |
| `mass` | `M` | Atomic mass |
| `σ` | `S` | Lennard-Jones sigma parameter |
| `ϵ` | `E` | Lennard-Jones epsilon parameter |

Sources: [src/setup.jl:136-144](), [src/setup.jl:181-191]()

### XML Parsing Process

The force field XML parsing follows OpenMM conventions and extracts parameters for different interaction types:

```mermaid
flowchart TD
    XMLFile["OpenMM XML File"] --> Parser["parsexml()"]
    Parser --> Root["Root Element"]
    
    Root --> AtomTypesSection["AtomTypes Section"]
    Root --> ResiduesSection["Residues Section"] 
    Root --> BondForceSection["HarmonicBondForce Section"]
    Root --> AngleForceSection["HarmonicAngleForce Section"]
    Root --> TorsionForceSection["PeriodicTorsionForce Section"]
    Root --> NonbondedSection["NonbondedForce Section"]
    
    AtomTypesSection --> ParseAtomTypes["Parse Atom Types<br/>(name, class, element, mass)"]
    ResiduesSection --> ParseResidues["Parse Residues<br/>(atom assignments, charges)"]
    BondForceSection --> ParseBonds["Parse Bond Parameters<br/>(k, length)"]
    AngleForceSection --> ParseAngles["Parse Angle Parameters<br/>(k, angle)"]
    TorsionForceSection --> ParseTorsions["Parse Torsion Parameters<br/>(periodicities, phases, ks)"]
    NonbondedSection --> ParseNonbonded["Parse LJ Parameters<br/>(sigma, epsilon, charges)"]
    
    ParseAtomTypes --> MFF["MolecularForceField"]
    ParseResidues --> MFF
    ParseBonds --> MFF
    ParseAngles --> MFF
    ParseTorsions --> MFF
    ParseNonbonded --> MFF
```

Sources: [src/setup.jl:204-356]()

## System Construction from Files

### OpenMM-Style System Construction

The primary `System` constructor reads coordinate files and applies force field parameters:

```mermaid
graph TD
    CoordFile["Coordinate File<br/>(PDB, XYZ, etc.)"] --> Chemfiles["Chemfiles.Trajectory"]
    FF["MolecularForceField"] --> SystemConstructor["System Constructor"]
    Chemfiles --> SystemConstructor
    
    SystemConstructor --> ParseTopology["Parse Topology<br/>(bonds, angles, dihedrals)"]
    SystemConstructor --> ParseCoords["Parse Coordinates<br/>(positions, boundary)"]
    SystemConstructor --> MatchTypes["Match Atom Types<br/>(residue templates)"]
    
    ParseTopology --> BuildInteractions["Build Interaction Lists<br/>(bonds, angles, torsions)"]
    MatchTypes --> BuildAtoms["Build Atom Arrays<br/>(charges, masses, LJ params)"]
    ParseCoords --> ProcessCoords["Process Coordinates<br/>(center, wrap)"]
    
    BuildInteractions --> FinalSystem["System"]
    BuildAtoms --> FinalSystem
    ProcessCoords --> FinalSystem
    
    FinalSystem --> SetupNeighbors["Setup Neighbor Finder"]
    FinalSystem --> SetupSolvent["Setup Implicit Solvent<br/>(optional)"]
```

The constructor performs several key operations:

1. **Topology Analysis**: Extracts bonds, angles, and dihedrals from the coordinate file using Chemfiles
2. **Residue Matching**: Matches atoms to force field residue templates by name
3. **Parameter Assignment**: Applies force field parameters to create interaction objects
4. **Coordinate Processing**: Centers and wraps coordinates according to boundary conditions

Sources: [src/setup.jl:471-1009]()

### Gromacs File Support

Molly.jl also supports reading Gromacs topology and coordinate files directly:

```mermaid
graph LR
    TopFile[".top File"] --> ParseSections["Parse Sections<br/>[atomtypes], [bonds], etc."]
    GroFile[".gro File"] --> ParseCoords["Parse Coordinates<br/>and Box Dimensions"]
    
    ParseSections --> BuildFF["Build Force Field<br/>Parameters"]
    ParseCoords --> BuildCoords["Build Coordinate<br/>Arrays"]
    
    BuildFF --> System["System"]
    BuildCoords --> System
```

The Gromacs reader processes topology sections sequentially:

| Section | Purpose | Data Extracted |
|---------|---------|----------------|
| `[atomtypes]` | Atom type definitions | Mass, charge, σ, ϵ |
| `[atoms]` | Atom assignments | Type mapping, charges |
| `[bonds]` | Bond definitions | Harmonic bond parameters |
| `[angles]` | Angle definitions | Harmonic angle parameters |
| `[dihedrals]` | Torsion definitions | RB torsion parameters |
| `[pairs]` | 1-4 interactions | Special pair listings |

Sources: [src/setup.jl:1011-1399]()

## Interaction List Construction

### Topology Processing

The system construction process builds interaction lists by matching atomic connectivity to force field parameters:

```mermaid
flowchart TD
    TopologyData["Topology from File<br/>(bonds, angles, dihedrals)"] --> MatchParameters["Match to Force Field Parameters"]
    
    subgraph "Bond Processing"
        BondPairs["Bond Pairs<br/>(atom1, atom2)"] --> LookupBondType["Lookup Bond Type<br/>in force_field.bond_types"]
        LookupBondType --> CreateBond["Create HarmonicBond<br/>(k, r0)"]
    end
    
    subgraph "Angle Processing"  
        AngleTriplets["Angle Triplets<br/>(atom1, atom2, atom3)"] --> LookupAngleType["Lookup Angle Type<br/>in force_field.angle_types"]
        LookupAngleType --> CreateAngle["Create HarmonicAngle<br/>(k, θ0)"]
    end
    
    subgraph "Torsion Processing"
        TorsionQuads["Torsion Quartets<br/>(atom1, atom2, atom3, atom4)"] --> WildcardSearch["Wildcard Search<br/>in force_field.torsion_types"]
        WildcardSearch --> CreateTorsion["Create PeriodicTorsion<br/>(periodicities, phases, ks)"]
    end
    
    MatchParameters --> BondPairs
    MatchParameters --> AngleTriplets  
    MatchParameters --> TorsionQuads
    
    CreateBond --> InteractionLists["InteractionList2Atoms<br/>InteractionList3Atoms<br/>InteractionList4Atoms"]
    CreateAngle --> InteractionLists
    CreateTorsion --> InteractionLists
```

The torsion parameter matching uses a sophisticated wildcard system where empty strings in the force field type keys act as wildcards, allowing flexible parameter assignment.

Sources: [src/setup.jl:604-818]()

### Special Topology Handling

The system handles several special cases during topology construction:

1. **N-terminal Residues**: Automatically adds missing hydrogen bonds and angles for N-terminal amino acids
2. **Water Molecules**: Adds missing O-H bonds and H-O-H angles for HOH residues  
3. **Terminal Residue Renaming**: Modifies residue names (e.g., "MET" → "NMET") to match appropriate force field templates

```mermaid
graph TD
    ResidueAnalysis["Analyze Residues"] --> CheckStandard{"Standard Residue?"}
    CheckStandard -->|Yes| CheckTerminal{"Terminal Position?"}
    CheckStandard -->|No| CheckWater{"Water (HOH)?"}
    
    CheckTerminal -->|N-terminal| AddNtermBonds["Add N-terminal H bonds<br/>and angles"]
    CheckTerminal -->|C-terminal| AddCtermBonds["Add C-terminal bonds"]
    CheckTerminal -->|Internal| StandardProcessing["Standard Processing"]
    
    CheckWater -->|Yes| AddWaterBonds["Add O-H bonds<br/>and H-O-H angle"]
    CheckWater -->|No| StandardProcessing
    
    AddNtermBonds --> ContinueProcessing["Continue Processing"]
    AddCtermBonds --> ContinueProcessing
    AddWaterBonds --> ContinueProcessing
    StandardProcessing --> ContinueProcessing
```

Sources: [src/setup.jl:521-565]()

## Coordinate and Random Placement

### File-based Coordinate Reading

Coordinate files are read using the Chemfiles library, which supports numerous formats:

```mermaid
graph LR
    ChemfilesFormats["Supported Formats<br/>PDB, XYZ, DCD, TRR, etc."] --> ChemfilesReader["Chemfiles.Trajectory"]
    ChemfilesReader --> ExtractData["Extract Positions<br/>and Unit Cell"]
    ExtractData --> ConvertUnits["Convert Å → nm"]
    ConvertUnits --> ProcessCoords["Center and Wrap<br/>Coordinates"]
    ProcessCoords --> SystemCoords["System.coords"]
```

The coordinate processing includes automatic unit conversion from Ångströms to nanometers and optional coordinate centering within the simulation box.

Sources: [src/setup.jl:917-926]()

### Random Coordinate Generation

For system initialization without coordinate files, Molly.jl provides functions to generate random atomic positions:

| Function | Purpose | Key Parameters |
|----------|---------|----------------|
| `place_atoms` | Random monatomic placement | `n_atoms`, `boundary`, `min_dist` |
| `place_diatomics` | Random diatomic molecules | `n_molecules`, `bond_length`, `aligned` |

Both functions ensure minimum distance constraints and account for periodic boundary conditions.

```mermaid
graph TD
    RandomPlacement["Random Placement Request"] --> CheckBoundary{"Finite Boundary?"}
    CheckBoundary -->|No| ThrowError["Throw ArgumentError"]
    CheckBoundary -->|Yes| CalculateMaxAtoms["Calculate Maximum<br/>Possible Atoms"]
    
    CalculateMaxAtoms --> CheckFeasibility{"n_atoms ≤ max_atoms?"}
    CheckFeasibility -->|No| ThrowError
    CheckFeasibility -->|Yes| PlacementLoop["Placement Loop"]
    
    PlacementLoop --> GenerateCoord["Generate Random<br/>Coordinate"]
    GenerateCoord --> CheckMinDist{"Min Distance<br/>Satisfied?"}
    CheckMinDist -->|Yes| AcceptCoord["Accept Coordinate"]
    CheckMinDist -->|No| IncrementAttempts["Increment Failed<br/>Attempts"]
    
    IncrementAttempts --> CheckMaxAttempts{"Max Attempts<br/>Reached?"}
    CheckMaxAttempts -->|Yes| ThrowError
    CheckMaxAttempts -->|No| PlacementLoop
    
    AcceptCoord --> CheckComplete{"All Atoms<br/>Placed?"}
    CheckComplete -->|No| PlacementLoop
    CheckComplete -->|Yes| ReturnCoords["Return Coordinates"]
```

Sources: [src/setup.jl:27-129]()

## Position Restraints and System Modification

### Adding Position Restraints

The `add_position_restraints` function demonstrates how to modify existing systems by adding new interaction types:

```mermaid
graph LR
    ExistingSystem["Existing System"] --> SelectAtoms["Select Atoms<br/>(atom_selector function)"]
    RestraintParams["Force Constants<br/>(k values)"] --> SelectAtoms
    ReferenceCoords["Reference Coordinates<br/>(restrain_coords)"] --> SelectAtoms
    
    SelectAtoms --> CreateRestraints["Create HarmonicPositionRestraint<br/>objects"]
    CreateRestraints --> BuildInteractionList["Build InteractionList1Atoms"]
    BuildInteractionList --> AddToSystem["Add to specific_inter_lists"]
    AddToSystem --> ModifiedSystem["Modified System"]
```

The function supports flexible atom selection through callback functions like `is_heavy_atom` and `is_any_atom`.

Sources: [src/setup.jl:1437-1477]()

---