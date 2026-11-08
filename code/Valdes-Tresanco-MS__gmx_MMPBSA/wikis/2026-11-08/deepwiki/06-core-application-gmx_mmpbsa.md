# Core Application (gmx_MMPBSA)

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [GMXMMPBSA/calculation.py](GMXMMPBSA/calculation.py)
- [GMXMMPBSA/createinput.py](GMXMMPBSA/createinput.py)
- [GMXMMPBSA/input_parser.py](GMXMMPBSA/input_parser.py)
- [GMXMMPBSA/main.py](GMXMMPBSA/main.py)
- [GMXMMPBSA/make_top.py](GMXMMPBSA/make_top.py)
- [GMXMMPBSA/utils.py](GMXMMPBSA/utils.py)
- [docs/changelog.md](docs/changelog.md)
- [docs/compatibility.md](docs/compatibility.md)
- [docs/gmx_MMPBSA_running.md](docs/gmx_MMPBSA_running.md)
- [docs/input_file.md](docs/input_file.md)
- [docs/news.json](docs/news.json)

</details>



## Purpose and Scope

This page documents the core calculation engine of `gmx_MMPBSA`, which orchestrates the end-to-end workflow from GROMACS input files to binding free energy results. The core application performs topology conversion, trajectory preparation, calculation execution using AmberTools engines, and statistical analysis of results.

For installation procedures, see [Installation and Setup](#2). For the graphical analyzer application, see [Analyzer Application (gmx_MMPBSA_ana)](#4). For input file configuration details, see [Input File Configuration](#3.2).

---

## Application Entry Point and Initialization

The `gmx_MMPBSA` application begins execution through the `gmxmmpbsa()` function, which serves as the main entry point after command-line invocation.

### Entry Flow

```mermaid
graph TB
    CLI["Command Line:<br/>gmx_MMPBSA"] --> gmxmmpbsa["gmxmmpbsa()<br/>app.py"]
    gmxmmpbsa --> setup["setup_run()<br/>main.py"]
    setup --> app_init["MMPBSA_App(MPI)<br/>main.py:65"]
    app_init --> parse_args["get_cl_args()<br/>main.py"]
    parse_args --> check_mode{Check Mode}
    
    check_mode -->|createinput| create["Generate mmpbsa.in<br/>Exit"]
    check_mode -->|clean| clean["Remove temp files<br/>Exit"]
    check_mode -->|infilehelp| help["Print input help<br/>Exit"]
    check_mode -->|normal| continue["Continue to<br/>Calculation"]
    
    continue --> read_input["read_input_file()<br/>main.py"]
    read_input --> check_parms["check_for_bad_input()<br/>main.py"]
    check_parms --> make_prmtops["make_prmtops()<br/>main.py"]
```

**Sources:** [GMXMMPBSA/main.py:65-110]()

### MMPBSA_App Class

The `MMPBSA_App` class in [GMXMMPBSA/main.py:65-110]() is the central orchestrator that manages:

- **MPI Communication**: Coordinates parallel execution across multiple processes via `mpi4py` or a `FakeMPI` fallback [GMXMMPBSA/main.py:58-60]().
- **Timing**: Tracks execution time for profiling via `Timer` objects [GMXMMPBSA/main.py:103-104]().
- **File Management**: Maintains `FILES` object with all input/output file paths.
- **Input Configuration**: Stores parsed input in `INPUT` dictionary [GMXMMPBSA/main.py:81]().

| Attribute | Type | Purpose |
|-----------|------|---------|
| `MPI` | MPI object | MPI communication handler [GMXMMPBSA/main.py:79]() |
| `mpi_rank` | int | Process rank (0 = master) [GMXMMPBSA/main.py:93]() |
| `mpi_size` | int | Total number of MPI processes [GMXMMPBSA/main.py:95]() |
| `INPUT` | dict | Parsed input file parameters [GMXMMPBSA/main.py:81]() |
| `master` | bool | True if rank == 0 [GMXMMPBSA/main.py:94]() |

**Sources:** [GMXMMPBSA/main.py:65-110]()

---

## Command-Line Interface

The command-line interface is defined by `GMXMMPBSA_ArgParser` in `commandlineparser.py`, which parses arguments into a `FILES` object.

### Required Arguments

```mermaid
graph LR
    subgraph "Complex Files"
        CS["-cs<br/>complex structure<br/>.tpr or .pdb"] 
        CT["-ct<br/>complex trajectory<br/>.xtc or .trr"]
        CI["-ci<br/>complex index<br/>.ndx"]
        CG["-cg<br/>complex groups<br/>rec & lig names"]
    end
    
    subgraph "Input File"
        IN["-i<br/>mmpbsa.in<br/>calculation params"]
    end
    
    subgraph "Optional Topology"
        CTOP["-cp<br/>complex topology<br/>.top"]
    end
    
    CS --> Parser["GMXMMPBSA_ArgParser"]
    CT --> Parser
    CI --> Parser
    CG --> Parser
    IN --> Parser
    CTOP --> Parser
    
    Parser --> FILES["FILES object"]
```

**Sources:** [GMXMMPBSA/main.py:68-69](), [GMXMMPBSA/make_top.py:114-118]()

---

## Calculation Workflow

The core application follows a linear pipeline from input processing to result generation, orchestrated by the `MMPBSA_App` class.

### Workflow Diagram

```mermaid
graph TB
    Start["Application Start"] --> ParseCLI["parser<br/>commandlineparser.py"]
    ParseCLI --> ReadInput["input_file<br/>input_parser.py:171"]
    ReadInput --> ValidateInput["check_for_bad_input()<br/>main.py"]
    
    ValidateInput --> BuildTop["buildTopology()<br/>CheckMakeTop.buildTopology()<br/>make_top.py:119"]
    
    BuildTop --> gmx2pdb["gmx2pdb()<br/>Extract PDB from TPR/XTC<br/>make_top.py:263"]
    gmx2pdb --> TopChoice{Topology<br/>Provided?}
    
    TopChoice -->|Yes| gmxtop["gmxtop2prmtop()<br/>Convert GROMACS→Amber<br/>make_top.py:562"]
    TopChoice -->|No| tleap["pdb2prmtop()+makeToptleap()<br/>Generate with tleap<br/>make_top.py"]
    
    gmxtop --> FileSetup["file_setup()<br/>main.py:111"]
    tleap --> FileSetup
    
    FileSetup --> CreateInputs["create_inputs()<br/>Generate .mdin files<br/>createinput.py:49"]
    CreateInputs --> MakeTrajs["make_trajectories()<br/>cpptraj processing<br/>make_trajs.py:46"]
    
    MakeTrajs --> LoadCalc["load_calc_list()<br/>main.py"]
    LoadCalc --> RunCalc["calc_list.run()<br/>Execute calculations<br/>calculation.py:114"]
    
    RunCalc --> ParseOut["Parse mdout files<br/>amber_outputs.py"]
    ParseOut --> CalcStats["Calculate statistics<br/>BindingStatistics<br/>amber_outputs.py:32"]
    CalcStats --> WriteOut["write_outputs()<br/>output_file.py:47"]
    
    WriteOut --> End["FINAL_RESULTS_MMPBSA.dat"]
    
    style BuildTop fill:#e3f2fd
    style MakeTrajs fill:#f3e5f5
    style RunCalc fill:#fff9c4
    style WriteOut fill:#e8f5e9
```

**Sources:** [GMXMMPBSA/main.py:111-150](), [GMXMMPBSA/make_top.py:119-261](), [GMXMMPBSA/calculation.py:88-148](), [GMXMMPBSA/createinput.py:49-61]()

---

## Input File Parsing

Input file parsing converts the AMBER-style namelist format into a Python dictionary structure.

### InputFile Class Architecture

```mermaid
graph TB
    InputFile["InputFile()<br/>input_parser.py:171"] --> addNamelist["addNamelist()<br/>Add namelist sections<br/>input_parser.py:143"]
    
    addNamelist --> general["Namelist('general')<br/>input_parser.py:106"]
    addNamelist --> gb["Namelist('gb')"]
    addNamelist --> pb["Namelist('pb')"]
    addNamelist --> decomp["Namelist('decomp')"]
    addNamelist --> ala["Namelist('ala')"]
    
    general --> vars1["Variables:<br/>temperature, startframe,<br/>endframe, interval, etc."]
    gb --> vars2["Variables:<br/>igb, saltcon,<br/>intdiel, etc."]
    
    InputFile --> Parse["Parse(filename)<br/>input_parser.py"]
    Parse --> OUTPUT["INPUT dictionary"]
```

The `Namelist` class [GMXMMPBSA/input_parser.py:106-160]() holds many different `Variable` instances. Each `Variable` object [GMXMMPBSA/input_parser.py:30-104]() has a name, datatype, default value, and description.

**Sources:** [GMXMMPBSA/input_parser.py:30-171](), [docs/input_file.md:116-142]()

---

## Topology Generation

The `CheckMakeTop` class [GMXMMPBSA/make_top.py:76-112]() manages the conversion from GROMACS topology files to AMBER prmtop format.

### Force Field Handling

The topology generation supports multiple force fields and assigns `PBRadii` [GMXMMPBSA/make_top.py:55]():

| ID | PBRadii Set | Source |
|----|-------------|--------|
| 1 | bondi | AmberTools |
| 2 | mbondi | AmberTools |
| 3 | mbondi2 | AmberTools |
| 7 | charmm_radii | CHARMM |

**Sources:** [GMXMMPBSA/make_top.py:55](), [GMXMMPBSA/make_top.py:562-606]()

---

## Calculation Execution

The `CalculationList` class [GMXMMPBSA/calculation.py:88-148]() manages sequential execution of all calculations.

### Calculation Type Hierarchy

```mermaid
graph TB
    Base["MultiCalculation<br/>calculation.py:151"] --> Energy["EnergyCalculation<br/>GB/PB/RISM<br/>calculation.py:36"]
    Energy --> PBEnergy["PBEnergyCalculation<br/>Poisson-Boltzmann<br/>calculation.py:36"]
    Base --> Nmode["NmodeCalc<br/>Normal mode entropy<br/>calculation.py:37"]
    Base --> QH["QuasiHarmCalc<br/>Quasi-harmonic entropy<br/>calculation.py:37"]
    Base --> IE["InteractionEntropyCalc<br/>IE calculation<br/>calculation.py:38"]
    Base --> C2["C2EntropyCalc<br/>C2 calculation<br/>calculation.py:38"]
```

**Sources:** [GMXMMPBSA/calculation.py:36-151]()

---

## Output Processing and Statistics

Output processing parses raw calculation results and computes binding statistics using the `EnergyVector` system.

### EnergyVector Class

The `EnergyVector` class [GMXMMPBSA/utils.py:52-151]() extends `np.ndarray` to store per-frame energies with statistical methods:

| Method | Description | Implementation |
|--------|-------------|---------|
| `avg()` | Average energy | `np.average(self)` [GMXMMPBSA/utils.py:98-99]() |
| `stdev()` | Propagated SD | `self.com_std or self.std()` [GMXMMPBSA/utils.py:86-87]() |
| `sem()` | Standard error | `self.std() / sqrt(len(self))` [GMXMMPBSA/utils.py:89-90]() |
| `corr_add()` | Correlated sum | `get_corrstd()` [GMXMMPBSA/utils.py:101-107]() |

**Sources:** [GMXMMPBSA/utils.py:52-159]()

### Final Output Generation

The `write_outputs()` function [GMXMMPBSA/output_file.py:47]() generates final result files, including `FINAL_RESULTS_MMPBSA.dat` and binary `.mmxsa` files for the GUI.

**Sources:** [GMXMMPBSA/output_file.py:47](), [GMXMMPBSA/main.py:47]()

---

## MPI Parallelization

The core application supports MPI-based parallelization for frame-parallel calculations.

### Frame Distribution Strategy

Each MPI rank processes a subset of trajectory frames. The `CalculationList.run()` method [GMXMMPBSA/calculation.py:114-148]() iterates through calculations, where each rank executes the `calc.run(rank)` command [GMXMMPBSA/calculation.py:142]().

The progress is tracked via `tqdm` in a separate thread for PB/GB calculations [GMXMMPBSA/calculation.py:53-87]().

**Sources:** [GMXMMPBSA/calculation.py:53-148](), [GMXMMPBSA/main.py:92-101]()

---

## Error Handling

Custom exceptions are defined in `exceptions.py` to handle various failure modes:

- `GMXMMPBSA_ERROR`: General failure exit [GMXMMPBSA/exceptions.py:45]().
- `InputError`: Raised during namelist/variable parsing [GMXMMPBSA/input_parser.py:25]().
- `CalcError`: Raised when a subprocess calculation fails [GMXMMPBSA/calculation.py:41]().

**Sources:** [GMXMMPBSA/exceptions.py:25-45](), [GMXMMPBSA/calculation.py:41]()

---