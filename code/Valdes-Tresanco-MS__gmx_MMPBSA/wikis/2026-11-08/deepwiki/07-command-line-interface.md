# Command-Line Interface

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [GMXMMPBSA/GMXMMPBSA.sh](GMXMMPBSA/GMXMMPBSA.sh)
- [GMXMMPBSA/commandlineparser.py](GMXMMPBSA/commandlineparser.py)
- [GMXMMPBSA/tester.py](GMXMMPBSA/tester.py)
- [docs/env.yml](docs/env.yml)
- [docs/examples/gmx_MMPBSA_test.md](docs/examples/gmx_MMPBSA_test.md)
- [docs/gmx_MMPBSA_ana_command-line.md](docs/gmx_MMPBSA_ana_command-line.md)
- [docs/gmx_MMPBSA_command-line.md](docs/gmx_MMPBSA_command-line.md)
- [docs/installation.md](docs/installation.md)
- [docs/summary.md](docs/summary.md)
- [docs/versus.md](docs/versus.md)

</details>



This page documents the command-line interface for the core `gmx_MMPBSA` calculation application. It covers all available flags, required and optional input files, file format requirements, and execution modes (serial vs MPI parallel).

For information about the analyzer GUI application command-line, see [4.1](). For information about the testing utility command-line, see [2.2]().

## Overview

The `gmx_MMPBSA` command-line interface is defined in [GMXMMPBSA/commandlineparser.py:189-204]() using the `GMXMMPBSA_ArgParser` class, which extends Python's `ArgumentParser`. The application entry point in [GMXMMPBSA/app.py:40-112]() processes these arguments through the `MMPBSA_App` class.

Title: Command-Line Processing Architecture
```mermaid
graph TB
    CLI["gmx_MMPBSA command-line"]
    
    subgraph "Argument Groups"
        MISC["Miscellaneous Options<br/>-O, -prefix"]
        IO["Input/Output Files<br/>-i, -o, -do, -eo, -deo"]
        COMPLEX["Complex Group<br/>-cs, -ci, -cg, -ct, -cp, -cr"]
        RECEPTOR["Receptor Group<br/>-rs, -ri, -rg, -rt, -rp"]
        LIGAND["Ligand Group<br/>-lm, -ls, -li, -lg, -lt, -lp"]
        ACTIONS["Actions<br/>--rewrite-output, --clean"]
        SPECIAL["Special Modes<br/>--input-file-help, --create_input, -s"]
    end
    
    subgraph "Parser Classes"
        ARGPARSER["GMXMMPBSA_ArgParser"]
        OPTLIST["OptionList"]
    end
    
    subgraph "Application Flow (app.py)"
        APP["MMPBSA_App.get_cl_args"]
        INPUTPARSE["app.read_input_file"]
        TOPOGEN["app.make_prmtops"]
        TRAJPROC["make_trajs.py integration"]
        CALC["calculation.py integration"]
    end
    
    CLI --> MISC
    CLI --> IO
    CLI --> COMPLEX
    CLI --> RECEPTOR
    CLI --> LIGAND
    CLI --> ACTIONS
    CLI --> SPECIAL
    
    MISC --> ARGPARSER
    IO --> ARGPARSER
    COMPLEX --> ARGPARSER
    RECEPTOR --> ARGPARSER
    LIGAND --> ARGPARSER
    ACTIONS --> ARGPARSER
    SPECIAL --> ARGPARSER
    
    ARGPARSER --> OPTLIST
    OPTLIST --> APP
    APP --> INPUTPARSE
    INPUTPARSE --> TOPOGEN
    TOPOGEN --> TRAJPROC
    TRAJPROC --> CALC
```
**Sources:** [GMXMMPBSA/commandlineparser.py:189-204](), [GMXMMPBSA/app.py:40-112]()

## Argument Groups and Options

The command-line arguments are organized into logical groups. The following sections detail each group and its options.

### Miscellaneous Options

| Flag | Type | Default | Description |
|------|------|---------|-------------|
| `-O`, `--overwrite` | boolean | `False` | Allow output files to be overwritten |
| `-prefix` | string | `'_GMXMMPBSA_'` | Prefix for intermediate temporary files |
| `-sys_name` | string | `None` | System name. Overrides `sys_name` in the input file. |

**Sources:** [docs/gmx_MMPBSA_command-line.md:35-40]()

### Input and Output Files

| Flag | Type | Default | Required | Description |
|------|------|---------|----------|-------------|
| `-i` | file path | `None` | Yes | MM/PBSA input file containing calculation parameters |
| `-xvvfile` | file path | auto-detected | No | XVV solvent susceptibility file for 3D-RISM calculations |
| `-o` | file path | `'FINAL_RESULTS_MMPBSA.dat'` | No | Main output file with statistics |
| `-do` | file path | `'FINAL_DECOMP_MMPBSA.dat'` | No | Decomposition statistics output |
| `-eo` | file path | `None` | No | CSV-format per-frame energy output |
| `-deo` | file path | `None` | No | CSV-format per-residue decomposition output |
| `-nogui` | boolean | `True` | No | If True, does not open `gmx_MMPBSA_ana` after finish |
| `-s`, `--stability` | boolean | `False` | No | Perform stability calculation |

The `-xvvfile` default is determined at runtime in [GMXMMPBSA/commandlineparser.py:44-53]() by searching for the TIP3P XVV file in `$AMBERHOME` or falling back to the bundled file.

**Sources:** [GMXMMPBSA/commandlineparser.py:44-53](), [docs/gmx_MMPBSA_command-line.md:42-67]()

### Complex Group (Required)

The complex group contains the primary system files. At minimum, `-cs`, `-ci`, and `-cg` must be provided unless performing a stability calculation.

Title: Complex Input Validation Flow
```mermaid
graph LR
    subgraph "Complex Input Files"
        CS["-cs<br/>Structure File<br/>*.tpr, *.pdb, *.gro"]
        CI["-ci<br/>Index File<br/>*.ndx"]
        CG["-cg<br/>Receptor/Ligand Groups<br/>two integers/names"]
        CT["-ct<br/>Trajectory Files<br/>*.xtc, *.trr, *.pdb"]
        CP["-cp<br/>Topology File<br/>*.top"]
        CR["-cr<br/>Reference Structure<br/>*.pdb"]
    end
    
    subgraph "Validation Logic (commandlineparser.py)"
        STRUCT_CHECK["structure()<br/>@check_arg(['.tpr', '.pdb'])"]
        INDEX_CHECK["index()<br/>@check_arg(['.ndx'])"]
        TRAJ_CHECK["trajectory()<br/>@check_arg(['.xtc', '.trr', '.pdb'])"]
        TOP_CHECK["topology()<br/>@check_arg(['.top'])"]
    end
    
    CS --> STRUCT_CHECK
    CI --> INDEX_CHECK
    CT --> TRAJ_CHECK
    CP --> TOP_CHECK
    CR --> STRUCT_CHECK
```

| Flag | Format | Required | Description |
|------|--------|----------|-------------|
| `-cs` | `*.tpr`, `*.pdb`, `*.gro` | Yes | Complex structure file |
| `-ci` | `*.ndx` | Yes | GROMACS index file for the complex |
| `-cg` | `int int` | Yes | Receptor and ligand group numbers/names from index file |
| `-ct` | `*.xtc`, `*.trr`, `*.pdb` | Yes | Complex trajectory (multiple files allowed) |
| `-cp` | `*.top` | No | Complex Topology file. If defined, `-lm` is not needed. |
| `-cr` | `*.pdb` | No | Complex Reference Structure file for chain ID assignment. |

The argument type validators are implemented as decorators in [GMXMMPBSA/commandlineparser.py:56-66](). For example, the `structure()` validator at [GMXMMPBSA/commandlineparser.py:69-71]() ensures files exist and match required extensions.

**Sources:** [GMXMMPBSA/commandlineparser.py:56-71](), [docs/gmx_MMPBSA_command-line.md:69-95]()

### Receptor and Ligand Groups (Optional)

Required for multiple trajectory (MT) approach. If not provided, receptor and ligand info is generated from the complex.

| Flag | Format | Description |
|------|--------|-------------|
| `-rs` | `*.tpr`, `*.pdb` | Unbound receptor structure file |
| `-ri` | `*.ndx` | Receptor index file |
| `-rg` | `int` | Receptor group in receptor index file |
| `-rt` | `*.xtc`, `*.trr`, `*.pdb` | Receptor trajectory files |
| `-rp` | `*.top` | Receptor topology file |
| `-lm` | `*.mol2` | MOL2 file of unbound ligand (required for small molecules without topology) |
| `-ls` | `*.tpr`, `*.pdb` | Unbound ligand structure file |
| `-li` | `*.ndx` | Ligand index file |
| `-lg` | `int` | Ligand group in ligand index file |
| `-lt` | `*.xtc`, `*.trr`, `*.pdb` | Ligand trajectory files |
| `-lp` | `*.top` | Ligand topology file |

The `-lm` option is critical for small molecule ligands not parameterized in standard force fields [docs/gmx_MMPBSA_command-line.md:117-120]().

**Sources:** [docs/gmx_MMPBSA_command-line.md:97-133]()

### Special Actions and Modes

| Flag | Description |
|------|-------------|
| `--input-file-help` | Print all available options in the input file and exit |
| `--create_input` | Create a new input file with selected calculation type (`gb`, `pb`, `rism`, `ala`, `decomp`, `nmode`, `all`) |
| `--rewrite-output` | Skip calculations and reparse existing output files |
| `--clean` | Remove temporary files and exit |

**Sources:** [docs/gmx_MMPBSA_command-line.md:28-33](), [docs/gmx_MMPBSA_command-line.md:136-142]()

## Index Group Specification

The `-cg`, `-rg`, and `-lg` flags accept either **integer indices** or **group names** from the index file. This is handled by the `index_groups()` function:

Title: Group Index Parsing Logic
```mermaid
graph TB
    INPUT["Index Group Input<br/>-cg value"]
    
    TRY["Try int() conversion"]
    
    INT_SUCCESS["Integer: Group Number"]
    STR_RETURN["String: Group Name"]
    
    INPUT --> TRY
    TRY -->|success| INT_SUCCESS
    TRY -->|exception| STR_RETURN
```

The conversion logic is defined in [GMXMMPBSA/commandlineparser.py:112-117]():
```python
def index_groups(value):
    try:
        return int(value)
    except Exception:
        return value
```

**Sources:** [GMXMMPBSA/commandlineparser.py:112-117]()

## Execution Modes

### Serial Execution

The basic serial execution runs all preprocessing and calculations on a single process:

```bash
gmx_MMPBSA -O -i mmpbsa.in -cs com.tpr -ci index.ndx -cg 1 13 -ct com_traj.xtc
```

### MPI Parallel Execution

`gmx_MMPBSA` supports MPI for parallel frame processing. The application detects the MPI environment at startup in [GMXMMPBSA/app.py:53-55]().

Title: Parallel Execution Data Flow
```mermaid
graph TB
    subgraph "MPI Orchestration"
        START["mpirun -np N gmx_MMPBSA"]
        RANK0["Rank 0 (Master)"]
        RANKN["Rank 1..N (Workers)"]
    end
    
    subgraph "Workflow Stages"
        PRE["Preprocessing (Rank 0 only)<br/>Topology conversion, Traj setup"]
        DIST["Frame Distribution<br/>make_trajs.py"]
        CALC["Calculation<br/>calculation.py"]
        AGGR["Aggregation (Rank 0 only)<br/>Gathering EnergyVectors"]
    end
    
    START --> RANK0
    RANK0 --> PRE
    PRE --> DIST
    DIST --> RANK0
    DIST --> RANKN
    RANK0 --> CALC
    RANKN --> CALC
    CALC --> AGGR
```

**Key characteristics of MPI execution:**
- Preprocessing (topology generation) is performed by Rank 0.
- Trajectory frames are distributed across available MPI ranks.
- If the number of frames is less than the number of CPUs, the number of processors used will be equal to the number of frames [docs/examples/gmx_MMPBSA_test.md:61-64]().
- Normal mode calculations (`NMODE`) require significant memory per rank.

**Sources:** [GMXMMPBSA/app.py:53-55](), [docs/examples/gmx_MMPBSA_test.md:61-64]()

## Error Handling

Custom error handling is implemented through the `GMXMMPBSA_ArgParser.exit()` method at [GMXMMPBSA/commandlineparser.py:191-195]():

```python
def exit(self, status=0, message=None):
    if message:
        GMXMMPBSA_ERROR(message)
    sys.exit(status)
```

This ensures that command-line errors are logged via the `GMXMMPBSA_ERROR` exception before the application exits.

**Sources:** [GMXMMPBSA/commandlineparser.py:191-195](), [GMXMMPBSA/exceptions.py:28]()

## Testing Utility (gmx_MMPBSA_test)

The `gmx_MMPBSA_test` utility provides its own CLI for running validation examples.

| Flag | Description |
|------|-------------|
| `-t` | Test level (0=All, 1=Minimal, 2=Fast) or specific test numbers (3-25) |
| `-f` | Folder to store data |
| `-r` | Reuse existing test folder |
| `-n` | Number of processors for MPI per calculation |
| `-j` | Number of examples to run concurrently |

The mapping of test numbers to systems (e.g., 3: Protein-Ligand ST, 4: Protein-Protein) is defined in [GMXMMPBSA/tester.py:109-134]().

**Sources:** [GMXMMPBSA/tester.py:109-134](), [docs/examples/gmx_MMPBSA_test.md:10-68]()

---