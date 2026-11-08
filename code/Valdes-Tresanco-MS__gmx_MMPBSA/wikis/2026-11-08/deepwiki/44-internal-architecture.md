# Internal Architecture

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [GMXMMPBSA/__init__.py](GMXMMPBSA/__init__.py)
- [GMXMMPBSA/app.py](GMXMMPBSA/app.py)
- [GMXMMPBSA/main.py](GMXMMPBSA/main.py)
- [GMXMMPBSA/make_top.py](GMXMMPBSA/make_top.py)
- [GMXMMPBSA/make_trajs.py](GMXMMPBSA/make_trajs.py)
- [GMXMMPBSA/parm_setup.py](GMXMMPBSA/parm_setup.py)
- [GMXMMPBSA/timer.py](GMXMMPBSA/timer.py)
- [GMXMMPBSA/utils.py](GMXMMPBSA/utils.py)
- [setup.py](setup.py)

</details>



This page provides a detailed technical reference of gmx_MMPBSA's internal architecture, covering the core application class, module organization, execution pipeline, and key data structures. This documentation is intended for developers who want to understand the codebase structure or contribute to the project.

**Scope**: This page focuses on the architectural design, class hierarchy, and module responsibilities. For specialized topics, see: [Error Handling and Exceptions](#7.2), [Parallel Processing and Performance](#7.3), [File Management and Temporary Files](#7.4), and [API Reference (MMPBSA_API)](#7.5).

---

## Main Application Class

The entire application is orchestrated by the `MMPBSA_App` class, which serves as the central coordinator for all calculation phases. This class manages state, coordinates subsystems, and enforces the execution pipeline.

### MMPBSA_App Structure

```mermaid
graph TB
    subgraph "MMPBSA_App Class (main.py)"
        INIT["__init__(MPI, stdout, stderr, size)<br/>Lines 72-110<br/>• Initialize MPI context<br/>• Setup I/O streams<br/>• Create timers<br/>• Set rank/size"]
        
        SETUP["file_setup()<br/>Lines 111-183<br/>• Create input files<br/>• Generate trajectories<br/>• Distribute frames to ranks"]
        
        RUN["run_mmpbsa(rank)<br/>Lines 184-210<br/>• Load calculation list<br/>• Execute calculations<br/>• Synchronize MPI"]
        
        LOADCALC["load_calc_list()<br/>Lines 212-232<br/>• Build calculation sequence<br/>• Append normal system calcs<br/>• Append mutant system calcs"]
        
        PARSEOUT["parse_output_files()<br/>(external)<br/>• Read mdout files<br/>• Extract energies<br/>• Compute statistics"]
        
        WRITEOUT["write_final_outputs()<br/>(external)<br/>• Generate FINAL_RESULTS<br/>• Generate decomp files<br/>• Write CSV data"]
    end
    
    subgraph "Key Attributes"
        INPUT["INPUT: dict<br/>Parsed input parameters"]
        FILES["FILES: Namespace<br/>File paths and flags"]
        SYSTEMS["normal_system: MMPBSA_System<br/>mutant_system: MMPBSA_System<br/>Topology systems"]
        CALCLIST["calc_list: CalculationList<br/>Ordered calculations"]
        TIMER["timer: Timer<br/>Performance profiling"]
    end
    
    INIT --> SETUP
    SETUP --> RUN
    RUN --> LOADCALC
    RUN --> PARSEOUT
    PARSEOUT --> WRITEOUT
    
    INIT -.creates.-> INPUT
    INIT -.creates.-> FILES
    INIT -.creates.-> TIMER
    LOADCALC -.builds.-> CALCLIST
    SETUP -.uses.-> SYSTEMS
```

**Sources**: [GMXMMPBSA/main.py:65-210](), [GMXMMPBSA/app.py:40-113]()

### Execution Entry Points

The application provides three entry points defined in `app.py`, each registered as console scripts in `setup.py`:

| Entry Point | Function | Purpose | Key Steps |
|------------|----------|---------|-----------|
| `gmx_MMPBSA` | `gmxmmpbsa()` | Main calculation engine | Setup → Parse args → Process → Calculate → Output |
| `gmx_MMPBSA_ana` | `gmxmmpbsa_ana()` | GUI analyzer | Initialize Qt → Load files → Display interface |
| `gmx_MMPBSA_test` | `gmxmmpbsa_test()` | Installation validation | Clone examples → Run tests → Report results |

**Sources**: [GMXMMPBSA/app.py:116-169](), [setup.py:56-61]()

---

## Module Organization

The codebase is organized into functional modules, each responsible for a specific phase of the calculation pipeline. This separation of concerns enables maintainability and parallel development.

### Core Module Responsibilities

```mermaid
graph LR
    subgraph "Input Layer"
        CLP["commandlineparser.py<br/>parser<br/>Parse CLI arguments"]
        INP["input_parser.py<br/>input_file<br/>Parse mmpbsa.in"]
    end
    
    subgraph "Topology Layer"
        MAKETOP["make_top.py<br/>CheckMakeTop<br/>GROMACS → Amber topology"]
        PARMSETUP["parm_setup.py<br/>MMPBSA_System<br/>Map complex/receptor/ligand"]
    end
    
    subgraph "Trajectory Layer"
        MAKETRAJ["make_trajs.py<br/>Trajectory<br/>Strip, image, distribute frames"]
    end
    
    subgraph "Calculation Layer"
        CREATEIN["createinput.py<br/>SanderInput subclasses<br/>Generate .mdin files"]
        CALC["calculation.py<br/>CalculationList<br/>EnergyCalculation<br/>Execute sander/cpptraj"]
    end
    
    subgraph "Output Layer"
        AMBEROUT["amber_outputs.py<br/>AmberOutput subclasses<br/>Parse mdout files"]
        OUTFILE["output_file.py<br/>write_outputs()<br/>Format statistics"]
    end
    
    subgraph "Support Layer"
        UTILS["utils.py<br/>EnergyVector<br/>Helper functions"]
        TIMER["timer.py<br/>Timer<br/>Performance tracking"]
        EXCEPT["exceptions.py<br/>Custom exceptions"]
    end
    
    CLP --> MAKETOP
    INP --> CREATEIN
    MAKETOP --> PARMSETUP
    PARMSETUP --> MAKETRAJ
    MAKETRAJ --> CREATEIN
    CREATEIN --> CALC
    CALC --> AMBEROUT
    AMBEROUT --> OUTFILE
    
    UTILS -.supports.-> AMBEROUT
    TIMER -.supports.-> CALC
    EXCEPT -.used by.-> MAKETOP
    EXCEPT -.used by.-> CALC
```

**Sources**: [GMXMMPBSA/main.py:31-50](), [GMXMMPBSA/make_top.py:76-128](), [GMXMMPBSA/parm_setup.py:62-113](), [GMXMMPBSA/make_trajs.py:43-108](), [GMXMMPBSA/utils.py:52-63]()

### Module-to-Phase Mapping

| Phase | Primary Module | Key Classes/Functions | Purpose |
|-------|---------------|----------------------|---------|
| **Configuration** | `commandlineparser.py` | `parser` | Parse command-line arguments |
| | `input_parser.py` | `input_file` | Parse input file namelists |
| **Topology Generation** | `make_top.py` | `CheckMakeTop` | Convert GROMACS files to Amber format |
| | `parm_setup.py` | `MMPBSA_System`, `Residue` | Map complex to receptor/ligand |
| **Trajectory Preparation** | `make_trajs.py` | `Trajectory`, `make_trajectories()` | Strip solvent, distribute frames |
| **Input Generation** | `createinput.py` | `create_inputs()` | Generate calculation-specific `.mdin` files |
| **Calculation Execution** | `calculation.py` | `CalculationList`, `EnergyCalculation`, `NmodeCalc` | Run external programs, manage parallel execution |
| **Output Processing** | `amber_outputs.py` | `GBout`, `PBout`, `NMODEout`, `BindingStatistics` | Parse output files, compute statistics |
| | `output_file.py` | `write_outputs()`, `data2pkl()` | Format final results, serialize data |
| **Utilities** | `utils.py` | `EnergyVector`, helper functions | Data structures, common operations |

**Sources**: [GMXMMPBSA/main.py:1-50](), [GMXMMPBSA/app.py:40-112]()

---

## Pipeline Architecture

gmx_MMPBSA follows a strict sequential pipeline where each phase must complete before the next begins. The `MMPBSA_App` class enforces this ordering and manages state transitions.

### Complete Execution Pipeline

```mermaid
graph TB
    START["Application Start<br/>gmxmmpbsa()"]
    
    subgraph "Phase 1: Initialization"
        P1A["setup_run()<br/>Signal handlers"]
        P1B["MMPBSA_App.__init__()<br/>MPI, I/O, timers"]
        P1C["get_cl_args()<br/>commandlineparser.py"]
    end
    
    subgraph "Phase 2: Input Processing"
        P2A["read_input_file()<br/>input_parser.py"]
        P2B["process_input()<br/>Validate parameters"]
        P2C["check_for_bad_input()<br/>Sanity checks"]
    end
    
    subgraph "Phase 3: Topology Generation"
        P3A["make_prmtops()<br/>CheckMakeTop.buildTopology()"]
        P3B["gmx2pdb()<br/>Extract PDB from TPR/trajs"]
        P3C["gmxtop2prmtop() or pdb2prmtop()<br/>Convert to Amber format"]
        P3D["loadcheck_prmtops()<br/>MMPBSA_System.Map()"]
    end
    
    subgraph "Phase 4: File Setup"
        P4A["file_setup()<br/>create_inputs()"]
        P4B["make_trajectories()<br/>Trajectory.Setup()"]
        P4C["cpptraj execution<br/>Strip, image, split frames"]
        P4D["make_mutant_trajectories()<br/>Alanine scanning"]
    end
    
    subgraph "Phase 5: Calculation"
        P5A["run_mmpbsa()<br/>load_calc_list()"]
        P5B["CalculationList.append()<br/>Build calculation sequence"]
        P5C["CalculationList.run()<br/>Execute all calculations"]
        P5D["Calculation.run()<br/>sander/cpptraj subprocess"]
    end
    
    subgraph "Phase 6: Output Processing"
        P6A["parse_output_files()<br/>AmberOutput classes"]
        P6B["Read mdout files<br/>Extract energies per frame"]
        P6C["Compute statistics<br/>Mean, SD, SEM"]
        P6D["write_final_outputs()<br/>FINAL_RESULTS_MMPBSA.dat"]
    end
    
    subgraph "Phase 7: Finalization"
        P7A["data2pkl()<br/>Serialize to .mmxsa"]
        P7B["finalize()<br/>Cleanup, timing report"]
    end
    
    START --> P1A
    P1A --> P1B
    P1B --> P1C
    P1C --> P2A
    P2A --> P2B
    P2B --> P2C
    P2C --> P3A
    P3A --> P3B
    P3B --> P3C
    P3C --> P3D
    P3D --> P4A
    P4A --> P4B
    P4B --> P4C
    P4C --> P4D
    P4D --> P5A
    P5A --> P5B
    P5B --> P5C
    P5C --> P5D
    P5D --> P6A
    P6A --> P6B
    P6B --> P6C
    P6C --> P6D
    P6D --> P7A
    P7A --> P7B
```

**Sources**: [GMXMMPBSA/app.py:40-112](), [GMXMMPBSA/main.py:65-210]()

### Phase-Specific Details

#### Phase 3: Topology Generation

This critical phase converts GROMACS topology files to Amber format:

```mermaid
graph TB
    TPR["GROMACS TPR/PDB files"]
    TOP["GROMACS TOP file (optional)"]
    
    GMX2PDB["gmx2pdb()<br/>Extract PDB structures<br/>Line 123"]
    
    PATH1["Has .top file?"]
    
    GMXTOP["gmxtop2prmtop()<br/>Direct topology conversion<br/>Line 125"]
    
    PDB2PRM["pdb2prmtop()<br/>Generate via tleap<br/>Line 127"]
    
    TLEAP["makeToptleap()<br/>Call tleap with force fields<br/>Line 128"]
    
    PRMTOP["COM.prmtop<br/>REC.prmtop<br/>LIG.prmtop"]
    
    MAPSYS["MMPBSA_System.Map()<br/>Map atoms to receptor/ligand<br/>Line 117"]
    
    TPR --> GMX2PDB
    TOP --> PATH1
    GMX2PDB --> PATH1
    
    PATH1 -->|Yes| GMXTOP
    PATH1 -->|No| PDB2PRM
    
    PDB2PRM --> TLEAP
    GMXTOP --> PRMTOP
    TLEAP --> PRMTOP
    
    PRMTOP --> MAPSYS
```

**Critical Classes**:
- `CheckMakeTop`: Main topology generation orchestrator [GMXMMPBSA/make_top.py:76-144]()
- `MMPBSA_System`: Maps complex topology to receptor/ligand subsystems [GMXMMPBSA/parm_setup.py:62-113]()

**Sources**: [GMXMMPBSA/make_top.py:76-144](), [GMXMMPBSA/parm_setup.py:62-113]()

#### Phase 4: Trajectory Preparation

The `Trajectory` class uses cpptraj to prepare trajectory files:

```mermaid
graph LR
    INPUT_TRAJ["Input Trajectories<br/>*.xtc/*.trr/*.pdb"]
    
    TRAJ_CLASS["Trajectory(prmtop, trajs, cpptraj)"]
    
    SETUP["Setup(start, end, interval)<br/>Configure frame range"]
    
    RMS["rms(mask)<br/>RMS fit to reference"]
    
    OUTTRAJ["Outtraj(filename, frames, filetype)<br/>Write output trajectory"]
    
    SPLIT["Split by MPI rank<br/>complex.nc.0, complex.nc.1, ..."]
    
    INPUT_TRAJ --> TRAJ_CLASS
    TRAJ_CLASS --> SETUP
    SETUP --> RMS
    RMS --> OUTTRAJ
    OUTTRAJ --> SPLIT
```

**Key Functions**:
- `make_trajectories()`: Main trajectory preparation [GMXMMPBSA/make_trajs.py:43-110]()
- `Trajectory.Setup()`: Configure frame selection [GMXMMPBSA/make_trajs.py:69]()

**Sources**: [GMXMMPBSA/make_trajs.py:43-110]()

---

## Key Data Structures

### EnergyVector Class

The `EnergyVector` is a NumPy array subclass that stores per-frame energies and provides automatic statistical methods:

```mermaid
graph TB
    subgraph "EnergyVector (utils.py:52-151)"
        VECTOR["EnergyVector(values, com_std)<br/>NumPy ndarray subclass"]
        
        METHODS["Statistical Methods:<br/>• mean() - Average<br/>• std() - Standard deviation<br/>• stdev() - Propagated SD<br/>• sem() - Standard error<br/>• semp() - Propagated SEM"]
        
        OPS["Arithmetic Operations:<br/>• __add__ - Propagate uncertainty<br/>• __sub__ - Propagate uncertainty<br/>• corr_add - Correlated addition<br/>• corr_sub - Correlated subtraction"]
    end
```

**Sources**: [GMXMMPBSA/utils.py:52-151]()

### MMPBSA_System Class

The `MMPBSA_System` class maintains the mapping between complex, receptor, and ligand topologies:

```mermaid
graph TB
    subgraph "MMPBSA_System (parm_setup.py:62-295)"
        INIT["__init__(complex_prmtop, receptor_prmtop, ligand_prmtop)<br/>Load all topology files"]
        
        MAP["Map(receptor_mask, ligand_mask)<br/>Build residue mappings"]
        
        RESLIST["res_list: List[Residue]<br/>All residues in complex with mappings"]
    end
    
    subgraph "Residue Class (parm_setup.py:33-60)"
        RES["Residue<br/>number: complex number<br/>receptor_number: receptor number<br/>ligand_number: ligand number<br/>selected: bool"]
    end
    
    INIT --> MAP
    MAP --> RESLIST
    RESLIST --> RES
```

**Sources**: [GMXMMPBSA/parm_setup.py:33-113]()

---

## Summary

The gmx_MMPBSA internal architecture follows a strict pipeline design:

1. **MMPBSA_App** orchestrates the workflow [GMXMMPBSA/main.py:65]().
2. **CheckMakeTop** handles topology conversion [GMXMMPBSA/make_top.py:76]().
3. **Trajectory** class prepares frames via cpptraj [GMXMMPBSA/make_trajs.py:68]().
4. **EnergyVector** manages data and statistical propagation [GMXMMPBSA/utils.py:52]().

**Sources**: [GMXMMPBSA/main.py:65](), [GMXMMPBSA/make_top.py:76](), [GMXMMPBSA/make_trajs.py:43](), [GMXMMPBSA/utils.py:52]()

---