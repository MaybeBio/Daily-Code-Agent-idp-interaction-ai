# File Management

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [docs/src/rmmol.md](docs/src/rmmol.md)
- [pacs/utils/rmfile.py](pacs/utils/rmfile.py)
- [pacs/utils/rmmol.py](pacs/utils/rmmol.py)

</details>



## Purpose and Scope

The PaCS-Toolkit provides two essential file management utilities that help manage disk space and improve performance during and after PaCS-MD simulations. This page documents the `rmfile` and `rmmol` utilities, which handle cleanup of temporary files and reduction of trajectory file sizes, respectively. For information about feature generation and trajectory analysis, see [Feature Generation](#6.2) and [Trajectory Fitting](#6.3).

## Overview of File Management Utilities

PaCS-MD simulations generate numerous files during execution, many of which are temporary or contain redundant information. The toolkit provides two key utilities to manage these files:

1. **rmfile**: Removes temporary and unnecessary files generated during simulation
2. **rmmol**: Reduces trajectory file sizes by removing unnecessary atoms (e.g., solvent molecules)

```mermaid
flowchart TD
    subgraph "PaCS-Toolkit File Management"
        rmfile["rmfile utility\n(Delete temporary files)"]
        rmmol["rmmol utility\n(Reduce trajectory file size)"]
    end
    
    MDRun["pacs mdrun"] --> |"Generates many files"| SimulationOutput["Simulation Files\n- Trajectory files\n- Checkpoint files\n- Topology files\n- Log files"]
    
    SimulationOutput --> rmfile
    SimulationOutput --> rmmol
    
    rmfile --> |"Removes temporary files"| CleanedOutput["Cleaned Output\n(No temporary files)"]
    rmmol --> |"Removes unnecessary atoms"| ReducedOutput["Reduced Output\n(Smaller trajectory files)"]
    
    CleanedOutput --> |"Takes less disk space"| Benefits["Benefits:\n- Reduced disk usage\n- Cleaner directories\n- Faster analysis"]
    ReducedOutput --> |"Faster to analyze"| Benefits
```

Sources: [pacs/utils/rmfile.py:1-117](), [pacs/utils/rmmol.py:1-316]()

## File Structure in PaCS-MD Simulations

Before discussing the file management utilities, it's helpful to understand how files are organized in a PaCS-MD simulation:

```mermaid
flowchart TD
    Trial["trial001/"] --> Cycle0["cycle000/"]
    Trial --> Cycle1["cycle001/"]
    Trial --> Cycle2["cycle002/"]
    Trial --> CycleN["..."]
    
    Cycle1 --> Summary["summary/"]
    Cycle1 --> Replica1["replica001/"]
    Cycle1 --> Replica2["replica002/"]
    Cycle1 --> ReplicaN["..."]
    
    Replica1 --> Traj["prd.xtc\n(Trajectory)"]
    Replica1 --> TPR["prd.tpr\n(Topology)"]
    Replica1 --> Log["*.log\n(Log files)"]
    Replica1 --> Other["Other files\n(.mdp, .cpt, etc.)"]
    
    Summary --> ProgressLog["progress.log\n(Status tracking)"]
    Summary --> AnalysisFiles["Analysis files"]
```

Sources: [pacs/utils/rmfile.py:10-31](), [pacs/utils/rmmol.py:11-33]()

## The rmfile Utility

The `rmfile` utility cleans up temporary files that are no longer needed after a simulation cycle has completed. It operates on a per-cycle basis and handles different files based on the molecular dynamics simulator used.

### How rmfile Works

```mermaid
sequenceDiagram
    participant User
    participant rmfile
    participant Settings as "MDsettings"
    participant FileSystem as "File System"
    
    User->>rmfile: "Execute rmfile utility"
    rmfile->>Settings: "Get cycle information"
    rmfile->>rmfile: "detect_n_replica(settings, cycle)"
    
    loop For each replica in cycle
        rmfile->>Settings: "Get replica directory"
        rmfile->>rmfile: "Determine simulator (gromacs/amber/namd)"
        
        alt gromacs
            rmfile->>FileSystem: "Remove mdout.mdp"
            rmfile->>FileSystem: "Remove *.cpt files"
            rmfile->>FileSystem: "Remove backup files"
            rmfile->>FileSystem: "Remove prd.gro"
            opt If not cycle 0 or max_cycle AND rmmol=true
                rmfile->>FileSystem: "Remove prd.tpr"
            end
        else amber
            rmfile->>FileSystem: "Remove prd.mdinfo"
            rmfile->>FileSystem: "Remove prd.log"
            rmfile->>FileSystem: "Remove prd.rst7"
        else namd
            Note over rmfile,FileSystem: Limited cleanup for NAMD files
        end
    end
    
    rmfile->>FileSystem: "Record completion in progress.log"
```

Sources: [pacs/utils/rmfile.py:42-109]()

### Files Removed by Simulator Type

The `rmfile` utility removes different files depending on the simulator being used:

| Simulator | Files Removed |
|-----------|---------------|
| GROMACS   | `.mdp` files (mdout.mdp), `.cpt` files, backup files (`#*`), output structure files (prd.gro), and optionally `.tpr` files |
| AMBER     | `.mdinfo` files, empty log files, output structure files (prd.rst7) |
| NAMD      | Limited cleanup (most NAMD files are preserved) |

Sources: [pacs/utils/rmfile.py:46-99]()

## The rmmol Utility

The `rmmol` utility reduces the size of trajectory files by removing unnecessary atoms, typically solvent molecules like water. This makes trajectory files smaller and faster to analyze.

### How rmmol Works

```mermaid
sequenceDiagram
    participant User
    participant rmmol
    participant Analyzer as "Analyzer (mdtraj/gromacs/cpptraj)"
    participant FileSystem as "File System"
    
    User->>rmmol: "Execute rmmol utility"
    
    alt First time
        rmmol->>rmmol: "make_top(settings)"
        rmmol->>Analyzer: "Create topology file without unnecessary atoms"
        Analyzer->>FileSystem: "Save rmmol_top.pdb"
    end
    
    rmmol->>rmmol: "detect_n_replica(settings, cycle)"
    
    par Process replicas in parallel
        loop For each replica
            rmmol->>FileSystem: "Check if trajectory exists"
            
            alt mdtraj analyzer
                rmmol->>Analyzer: "Load trajectory with mdtraj"
                rmmol->>Analyzer: "Select atoms to keep"
                rmmol->>Analyzer: "Create reduced trajectory"
                Analyzer->>FileSystem: "Save prd_rmmol.xtc"
            else gromacs analyzer
                rmmol->>Analyzer: "Run gmx trjconv with selection"
                Analyzer->>FileSystem: "Save prd_rmmol.xtc"
            else cpptraj analyzer
                rmmol->>Analyzer: "Create and run cpptraj script"
                Analyzer->>FileSystem: "Save prd_rmmol.nc"
            end
            
            opt If not last_cycle
                rmmol->>FileSystem: "Remove original trajectory file"
            end
        end
    end
    
    rmmol->>FileSystem: "Record completion in progress.log"
```

Sources: [pacs/utils/rmmol.py:134-280]()

### Analyzer-Specific Behavior

The `rmmol` utility works differently depending on the analyzer specified:

| Analyzer | Method |
|----------|--------|
| mdtraj   | Uses MDTraj to load the trajectory, select atoms to keep, and save a reduced trajectory |
| gromacs  | Uses GROMACS's `trjconv` tool to select atoms and convert the trajectory |
| cpptraj  | Uses a CPPTraj script to process the trajectory and remove unnecessary atoms |

Sources: [pacs/utils/rmmol.py:47-131](), [pacs/utils/rmmol.py:177-280]()

## Command-Line Usage

The PaCS-Toolkit provides command-line interfaces for both file management utilities.

### rmfile Command

The `rmfile` command is typically run automatically during PaCS-MD simulations but can also be run manually:

```
pacs rmfile -t TRIAL_NUMBER
```

This will clean up all unnecessary files for all cycles in the specified trial.

### rmmol Command

The `rmmol` command must be run manually after a simulation is complete. It has different options depending on the analyzer used:

#### For MDTraj:
```
pacs rmmol mdtraj -t TRIAL_NUMBER -k "not water" -e .xtc -m TOPOLOGY_FILE
```

#### For GROMACS:
```
pacs rmmol gmx -t TRIAL_NUMBER -k "not_water" -e .xtc -n index.ndx -g gmx --nojump
```

#### For CPPTraj:
```
pacs rmmol cpptraj -t TRIAL_NUMBER -k "not water" -e .nc -p TOPOLOGY_FILE
```

Sources: [docs/src/rmmol.md:1-76]()

## File Management Workflow in PaCS-MD

The following diagram shows how file management fits into the overall PaCS-MD workflow:

```mermaid
flowchart TD
    Start["Start PaCS-MD"] --> MDRun["pacs mdrun"]
    
    MDRun --> Cycle["Cycle Management"]
    Cycle --> Simulator["Run Simulator\n(gromacs/amber/namd)"]
    Simulator --> |"Generate trajectories"| Analyzer["Run Analyzer"]
    
    Analyzer --> |"Select replicas for next cycle"| NextCycle{"More\ncycles?"}
    NextCycle --> |"Yes"| RMFile["rmfile\n(Clean up current cycle)"]
    RMFile --> Cycle
    
    NextCycle --> |"No"| Complete["Simulation Complete"]
    Complete --> RMMol["rmmol\n(Reduce trajectory size)"]
    
    RMMol --> |"Smaller trajectories"| PostProcess["Post-Processing\n(Feature Generation, MSM, etc.)"]
```

Sources: [pacs/utils/rmfile.py:42-116](), [pacs/utils/rmmol.py:134-315]()

## Implementation Details

### File Detection Functions

Both utilities use similar functions to detect the number of replicas and cycles:

```python
def detect_n_replica(settings: MDsettings, cycle: int) -> int:
    # Detect number of replicas in a cycle
    for replica in range(1, 1000):
        replica_dir = Path(settings.each_replica(_cycle=cycle, _replica=replica))
        if not replica_dir.exists():
            n_replica = replica - 1
            break
    return n_replica

def detect_n_cycle(settings: MDsettings) -> int:
    # Detect maximum cycle number
    for cycle in range(1000):
        cycle_dir = Path(settings.each_cycle(_cycle=cycle))
        if not cycle_dir.exists():
            max_cycle = cycle - 1
            break
    return max_cycle
```

Sources: [pacs/utils/rmfile.py:10-31](), [pacs/utils/rmmol.py:11-33]()

### Parallel Processing in rmmol

The `rmmol` utility processes multiple replicas in parallel to improve performance:

```mermaid
flowchart TD
    rmmol["rmmol(settings, cycle, last_cycle)"] --> BatchProcessing["Batch replicas\ninto groups of n_parallel"]
    
    BatchProcessing --> |"Batch 1"| ParallelBatch1["Process replicas in parallel"]
    BatchProcessing --> |"Batch 2"| ParallelBatch2["Process replicas in parallel"]
    BatchProcessing --> |"..."| ParallelBatchN["..."]
    
    ParallelBatch1 --> |"For each replica"| ProcessReplica["Process single replica\nbased on analyzer type"]
    
    ProcessReplica --> |"mdtraj"| MDTrajProcess["rmmol_replica_mdtraj()"]
    ProcessReplica --> |"gromacs"| GmxProcess["rmmol_replica_gmx()"]
    ProcessReplica --> |"cpptraj"| CpptrajProcess["rmmol_replica_cpptraj()"]
    
    MDTrajProcess --> Wait["Wait for all processes to complete"]
    GmxProcess --> Wait
    CpptrajProcess --> Wait
    
    Wait --> NextBatch{"More\nbatches?"}
    NextBatch --> |"Yes"| BatchProcessing
    NextBatch --> |"No"| RecordFinished["Record completion in log"]
```

Sources: [pacs/utils/rmmol.py:134-174]()

## Best Practices and Considerations

1. **Run rmfile during simulations**: The `rmfile` utility is automatically run during PaCS-MD simulations to keep disk usage manageable.

2. **Run rmmol after simulation completion**: The `rmmol` utility should be run after a simulation is complete to reduce the size of trajectory files for analysis.

3. **Keep original files for the last cycle**: By default, `rmmol` keeps the original trajectory files for the last cycle as a precaution.

4. **Consider PBC handling**: When using `rmmol` with GROMACS, consider using the `--nojump` option if your trajectories were generated with PBC handling.

5. **Selection syntax**: Be careful with the atom selection syntax, which differs between analyzers:
   - MDTraj: `"not water"`
   - GROMACS: `"not_water"` (group name in index file)
   - CPPTraj: `"not water"` or `":!WAT"`

Sources: [pacs/utils/rmmol.py:218-222](), [docs/src/rmmol.md:6-9]()

## Summary

The file management utilities in PaCS-Toolkit provide essential functionality for managing the files generated during PaCS-MD simulations. The `rmfile` utility removes temporary and unnecessary files, while the `rmmol` utility reduces the size of trajectory files by removing unnecessary atoms. These utilities help keep disk usage manageable and improve performance during analysis.

---