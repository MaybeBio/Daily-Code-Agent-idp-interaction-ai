# Visualization and Output

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [docs/src/api.md](docs/src/api.md)
- [src/analysis.jl](src/analysis.jl)
- [src/loggers.jl](src/loggers.jl)
- [test/Project.toml](test/Project.toml)

</details>



This page covers the data export and visualization capabilities in Molly.jl, including writing trajectory files, creating structure outputs, and generating visual representations of molecular dynamics simulations. For information about data collection during simulations, see [Logging and Analysis](#8.1).

## File Output System

Molly.jl provides comprehensive file writing capabilities through two main approaches: direct structure writing and trajectory logging during simulations.

### Structure Writing

The `write_structure` function provides one-time structure output to various file formats:

```mermaid
graph TD
    WriteStructure["write_structure()"] --> FormatDetection["Format Detection"]
    FormatDetection --> PDBPath["PDB Format"]
    FormatDetection --> ChemfilesPath["Other Formats"]
    
    PDBPath --> BioStructures["BioStructures.jl"]
    ChemfilesPath --> Chemfiles["Chemfiles.jl"]
    
    BioStructures --> PDBOutput["PDB Files<br/>(.pdb)"]
    Chemfiles --> MultiFormat["Multiple Formats<br/>DCD, XTC, CIF, MOL2<br/>SDF, TRR, XYZ"]
    
    WriteStructure --> SystemData["System Data"]
    SystemData --> AtomsData["atoms_data"]
    SystemData --> Topology["topology"]
    SystemData --> Coordinates["coords"]
    SystemData --> Boundary["boundary"]
```

Sources: [src/loggers.jl:479-524]()

### Trajectory Writers

The `TrajectoryWriter` logger enables continuous output during simulations:

```mermaid
graph LR
    TrajectoryWriter --> Configuration["Configuration"]
    Configuration --> FilePath["filepath"]
    Configuration --> Format["format"]
    Configuration --> AtomSelection["atom_inds"]
    Configuration --> ExcludedRes["excluded_res"]
    Configuration --> WriteVelocities["write_velocities"]
    Configuration --> WriteBoundary["write_boundary"]
    
    TrajectoryWriter --> OutputMethods["Output Methods"]
    OutputMethods --> PDBMethod["PDB Method<br/>write_pdb_coords()"]
    OutputMethods --> ChemfilesMethod["Chemfiles Method<br/>write_chemfiles!()"]
    
    PDBMethod --> PDBTrajectory["Multi-model PDB<br/>MODEL/ENDMDL records"]
    ChemfilesMethod --> BinaryTrajectory["Binary Trajectories<br/>DCD, XTC, TRR"]
```

Sources: [src/loggers.jl:526-617]()

## Supported File Formats

Molly.jl supports a wide range of molecular file formats through integration with established Julia packages:

| Format | Extension | Library | Coordinates | Velocities | Topology | Notes |
|--------|-----------|---------|-------------|------------|----------|-------|
| PDB | `.pdb` | BioStructures.jl | ✓ | ✗ | ✓ | Special handling for residue info |
| DCD | `.dcd` | Chemfiles.jl | ✓ | ✓* | ✗ | Binary trajectory format |
| XTC | `.xtc` | Chemfiles.jl | ✓ | ✗ | ✗ | Compressed trajectory |
| TRR | `.trr` | Chemfiles.jl | ✓ | ✓* | ✗ | Full precision trajectory |
| CIF | `.cif` | Chemfiles.jl | ✓ | ✗ | ✓ | Crystallographic format |
| MOL2 | `.mol2` | Chemfiles.jl | ✓ | ✗ | ✓ | Tripos format |
| SDF | `.sdf` | Chemfiles.jl | ✓ | ✗ | ✓ | Structure data format |
| XYZ | `.xyz` | Chemfiles.jl | ✓ | ✗ | ✗ | Simple coordinate format |

*Velocity support depends on format capabilities

Sources: [src/loggers.jl:488-489](), [src/loggers.jl:533-537]()

## Coordinate and Unit Handling

The output system automatically handles unit conversions and coordinate transformations:

```mermaid
graph TD
    InputCoords["System Coordinates"] --> UnitCheck["Unit Detection"]
    UnitCheck --> NoUnits["NoUnits<br/>(assume nm)"]
    UnitCheck --> WithUnits["With Units"]
    
    NoUnits --> ConvertNm["Convert nm → Å<br/>(× 10)"]
    WithUnits --> ConvertUnit["ustrip(u\"Å\", coord)"]
    
    ConvertNm --> OutputCoords["Output Coordinates"]
    ConvertUnit --> OutputCoords
    
    InputVelocities["System Velocities"] --> VelUnitCheck["Velocity Unit Check"]
    VelUnitCheck --> NoVelUnits["NoUnits<br/>(assume nm/ps)"]
    VelUnitCheck --> WithVelUnits["With Units"]
    
    NoVelUnits --> ConvertVelNm["Convert nm/ps → Å/ps<br/>(× 10)"]
    WithVelUnits --> ConvertVelUnit["ustrip(u\"Å * ps^-1\", vel)"]
    
    ConvertVelNm --> OutputVelocities["Output Velocities"]
    ConvertVelUnit --> OutputVelocities
```

Sources: [src/loggers.jl:378-384](), [src/loggers.jl:448-470]()

## Visualization System

Molly.jl provides molecular visualization through the GLMakie extension. The visualization system is accessed via the `visualize` function:

```mermaid
graph TB
    VisualizeFunc["visualize()"] --> InputData["Input Data"]
    InputData --> CoordLogger["coord_logger"]
    InputData --> Boundary["boundary"]
    InputData --> OutputPath["out_filepath"]
    
    VisualizeFunc --> Parameters["Visualization Parameters"]
    Parameters --> Connections["connections<br/>(bond pairs)"]
    Parameters --> ConnectionFrames["connection_frames<br/>(frame-specific bonds)"]
    Parameters --> Trails["trails<br/>(motion trails)"]
    Parameters --> Colors["color, connection_color"]
    Parameters --> Sizes["markersize, linewidth"]
    Parameters --> BoundaryVis["show_boundary"]
    
    VisualizeFunc --> GLMakieExt["GLMakie Extension"]
    GLMakieExt --> Animation["3D Animation Output"]
    
    CoordLogger --> FrameData["Frame-by-frame<br/>Coordinate Data"]
    FrameData --> GLMakieExt
```

Sources: [src/analysis.jl:110-140]()

## Analysis and Visualization Extensions

Molly.jl uses package extensions to provide optional visualization and analysis capabilities:

```mermaid
graph LR
    MollyCore["Molly.jl Core"] --> Extensions["Package Extensions"]
    
    Extensions --> GLMakieExt["MollyGLMakieExt"]
    Extensions --> KernelDensityExt["KernelDensity Support"]
    
    GLMakieExt --> VisualizeImpl["visualize() Implementation"]
    VisualizeImpl --> Animation3D["3D Molecular Animation"]
    VisualizeImpl --> AtomRendering["Atom Sphere Rendering"]
    VisualizeImpl --> BondRendering["Bond Line Rendering"]
    VisualizeImpl --> TrailRendering["Motion Trail Effects"]
    
    KernelDensityExt --> RDFImpl["rdf() Implementation"]
    RDFImpl --> RadialDistribution["Radial Distribution<br/>Function Calculation"]
    
    UserCode["User Code"] --> ImportGLMakie["using GLMakie"]
    UserCode --> ImportKDE["using KernelDensity"]
    ImportGLMakie --> EnableVisualize["Enables visualize()"]
    ImportKDE --> EnableRDF["Enables rdf()"]
```

Sources: [docs/src/api.md:7-10](), [src/analysis.jl:142-152]()

## Boundary Handling in Output

The output system provides specific support for different boundary condition types:

```mermaid
graph TD
    BoundaryOutput["Boundary Output"] --> BoundaryType["Boundary Type Check"]
    
    BoundaryType --> CubicBoundary["CubicBoundary"]
    BoundaryType --> TriclinicBoundary["TriclinicBoundary"] 
    BoundaryType --> InfiniteBoundary["Infinite Boundary"]
    
    CubicBoundary --> CubicCRYST1["CRYST1 Record<br/>90° angles"]
    TriclinicBoundary --> TriclinicCRYST1["CRYST1 Record<br/>α, β, γ angles"]
    InfiniteBoundary --> NoCRYST1["No CRYST1 Record"]
    
    CubicCRYST1 --> PDBOutput["PDB File Output"]
    TriclinicCRYST1 --> PDBOutput
    NoCRYST1 --> PDBOutput
    
    BoundaryOutput --> ChemfilesBC["Chemfiles UnitCell"]
    ChemfilesBC --> BinaryFormats["Binary Format Output<br/>DCD, XTC, TRR"]
```

Sources: [src/loggers.jl:343-364](), [src/loggers.jl:441-443]()

## Integration with Logging System

The output system integrates seamlessly with Molly.jl's logging infrastructure:

```mermaid
graph TB
    SimulationLoop["Simulation Loop"] --> ApplyLoggers["apply_loggers!()"]
    
    ApplyLoggers --> LoggerTypes["Logger Types"]
    LoggerTypes --> PropertyLoggers["Property Loggers<br/>TemperatureLogger<br/>EnergyLogger"]
    LoggerTypes --> TrajectoryWriter["TrajectoryWriter"]
    LoggerTypes --> StructureWriter["StructureWriter"]
    
    PropertyLoggers --> InMemoryData["In-Memory Data<br/>Values accessible<br/>via values()"]
    
    TrajectoryWriter --> FileOutput["Continuous File Output"]
    StructureWriter --> FileOutput
    
    FileOutput --> MultiFormat["Multiple Format Support<br/>PDB, DCD, XTC, etc."]
    
    InMemoryData --> PostAnalysis["Post-Simulation Analysis"]
    PostAnalysis --> Visualize["visualize()"]
    PostAnalysis --> AnalysisFuncs["Analysis Functions<br/>rmsd(), rdf()"]
```

Sources: [src/loggers.jl:28-47](), [src/loggers.jl:558-617]()

---