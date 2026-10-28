# GPU Acceleration

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [ext/MollyCUDAExt.jl](ext/MollyCUDAExt.jl)
- [ext/MollyGLMakieExt.jl](ext/MollyGLMakieExt.jl)
- [src/energy.jl](src/energy.jl)
- [src/force.jl](src/force.jl)
- [src/kernels.jl](src/kernels.jl)
- [test/gradients.jl](test/gradients.jl)

</details>



This section covers GPU acceleration capabilities in Molly.jl, including array type abstraction, GPU kernels for force and energy calculations, and hardware-specific optimizations. GPU acceleration enables high-performance molecular dynamics simulations by leveraging parallel computing architectures.

For information about the core force calculation algorithms, see [Force and Energy Calculations](#4.3). For details about neighbor finding algorithms, see [Neighbor Finding](#3.3).

## Array Type System and Backend Abstraction

Molly.jl uses an array type abstraction system that allows the same simulation code to run on different hardware backends without modification. The system automatically detects GPU arrays and dispatches to appropriate GPU kernels.

```mermaid
graph TB
    subgraph "Array Types"
        CPU_Array["Array<br/>(CPU)"]
        CUDA_Array["CuArray<br/>(NVIDIA GPU)"]
        ROC_Array["ROCArray<br/>(AMD GPU)"]
        Metal_Array["MtlArray<br/>(Apple GPU)"]
        OneAPI_Array["oneArray<br/>(Intel GPU)"]
    end
    
    subgraph "Backend Detection"
        get_backend["get_backend()"]
        Backend_CPU["CPU Backend"]
        Backend_CUDA["CUDA Backend"]
        Backend_ROCm["ROCm Backend"]
        Backend_Metal["Metal Backend"]
    end
    
    subgraph "Kernel Dispatch"
        KA_Kernels["KernelAbstractions.jl<br/>Generic Kernels"]
        CUDA_Kernels["MollyCUDAExt<br/>Optimized CUDA Kernels"]
        Fallback_CPU["CPU Fallback"]
    end
    
    CPU_Array --> Backend_CPU
    CUDA_Array --> Backend_CUDA
    ROC_Array --> Backend_ROCm
    Metal_Array --> Backend_Metal
    
    get_backend --> Backend_CPU
    get_backend --> Backend_CUDA
    get_backend --> Backend_ROCm
    get_backend --> Backend_Metal
    
    Backend_CPU --> Fallback_CPU
    Backend_CUDA --> CUDA_Kernels
    Backend_ROCm --> KA_Kernels
    Backend_Metal --> KA_Kernels
    
    CUDA_Kernels --> KA_Kernels
```

**Sources:** [src/force.jl:135-151](), [src/kernels.jl:35-43](), [ext/MollyCUDAExt.jl:13]()

The system determines GPU capability through the `uses_gpu_neighbor_finder` function, which returns `true` for CUDA arrays and enables GPU-specific neighbor finding algorithms.

## GPU Force Calculation Architecture

GPU force calculations use a two-level dispatch system: generic KernelAbstractions.jl kernels for broad GPU compatibility, and specialized CUDA kernels for optimal NVIDIA GPU performance.

```mermaid
graph TD
    subgraph "Force Calculation Entry Points"
        forces_nounits["forces_nounits!()"]
        GPU_Detection{"`GPU Array?`"}
        CPU_Path["CPU Force Calculation"]
        GPU_Path["GPU Force Calculation"]
    end
    
    subgraph "GPU Buffer Management"
        ForcesBuffer["ForcesBuffer"]
        fs_mat["fs_mat: Force Matrix"]
        box_mins["box_mins: Spatial Bounds"]
        box_maxs["box_maxs: Spatial Bounds"] 
        morton_seq["morton_seq: Ordering"]
    end
    
    subgraph "GPU Kernel Types"
        Pairwise_GPU["pairwise_force_gpu!()"]
        Specific_GPU["specific_force_gpu!()"]
        KA_Pairwise["pairwise_force_kernel_nl!"]
        KA_Specific["specific_force_*_atoms_kernel!"]
    end
    
    subgraph "CUDA Optimizations"
        CUDA_Nonl["pairwise_force_kernel_nonl!"]
        CUDA_Nl["force_kernel!"]
        Morton_Order["Morton Space-Filling Curve"]
        Warp_Shuffles["Warp-Level Shuffles"]
    end
    
    forces_nounits --> GPU_Detection
    GPU_Detection -->|CPU| CPU_Path
    GPU_Detection -->|GPU| GPU_Path
    
    GPU_Path --> ForcesBuffer
    ForcesBuffer --> fs_mat
    ForcesBuffer --> box_mins
    ForcesBuffer --> box_maxs
    ForcesBuffer --> morton_seq
    
    GPU_Path --> Pairwise_GPU
    GPU_Path --> Specific_GPU
    
    Pairwise_GPU --> KA_Pairwise
    Specific_GPU --> KA_Specific
    
    Pairwise_GPU -->|CUDA| CUDA_Nonl
    Pairwise_GPU -->|CUDA| CUDA_Nl
    
    CUDA_Nl --> Morton_Order
    CUDA_Nl --> Warp_Shuffles
```

**Sources:** [src/force.jl:372-405](), [src/kernels.jl:45-64](), [ext/MollyCUDAExt.jl:45-86]()

### GPU Buffer Initialization

GPU force calculations require specialized buffers managed by the `ForcesBuffer` struct:

- `fs_mat`: Force accumulation matrix with dimensions `(D, N)` where `D` is spatial dimensions and `N` is number of atoms
- `box_mins` and `box_maxs`: Spatial bounding boxes for neighbor finding optimization
- `morton_seq`: Morton-ordered atom indices for spatial locality
- `compressed_eligible` and `compressed_special`: Compressed boolean matrices for interaction eligibility

**Sources:** [src/force.jl:126-151]()

### Pairwise Force Kernels

The `pairwise_force_kernel_nl!` kernel processes neighbor lists in parallel, with each GPU thread handling one atom pair interaction:

```mermaid
graph LR
    subgraph "Kernel Execution"
        Thread_ID["Thread ID"]
        Neighbor_Pair["neighbors[thread_id]<br/>(i, j, special)"]
        Force_Calc["sum_pairwise_forces()"]
        Atomic_Add["Atomic Force Accumulation"]
    end
    
    Thread_ID --> Neighbor_Pair
    Neighbor_Pair --> Force_Calc
    Force_Calc --> Atomic_Add
```

**Sources:** [src/kernels.jl:66-83]()

### Specific Interaction Kernels

Specific interactions (bonds, angles, torsions) use separate kernels for each interaction type:

- `specific_force_1_atoms_kernel!`: Position restraints and external forces
- `specific_force_2_atoms_kernel!`: Bond interactions  
- `specific_force_3_atoms_kernel!`: Angle interactions
- `specific_force_4_atoms_kernel!`: Torsion interactions

Each kernel processes interaction lists in parallel with atomic force accumulation.

**Sources:** [src/kernels.jl:85-220]()

## GPU Energy Calculation Architecture

GPU potential energy calculations follow a similar pattern to force calculations but accumulate scalar energy values instead of vector forces.

```mermaid
graph TD
    subgraph "Energy Calculation Flow"
        PE_Entry["potential_energy()"]
        GPU_PE_Check{"`GPU Array?`"}
        CPU_PE["CPU Energy Calculation"]
        GPU_PE["GPU Energy Calculation"]
    end
    
    subgraph "GPU Energy Kernels"
        Pairwise_PE["pairwise_pe_gpu!()"]
        Specific_PE["specific_pe_gpu!()"]
        PE_Vector["pe_vec_nounits"]
        KA_PE_Kernels["KernelAbstractions.jl<br/>Energy Kernels"]
    end
    
    subgraph "Energy Accumulation"
        Atomic_Energy["Atomic Energy Addition"]
        Scalar_Sum["Final Energy Sum"]
        Energy_Units["Apply Energy Units"]
    end
    
    PE_Entry --> GPU_PE_Check
    GPU_PE_Check -->|CPU| CPU_PE
    GPU_PE_Check -->|GPU| GPU_PE
    
    GPU_PE --> Pairwise_PE
    GPU_PE --> Specific_PE
    
    Pairwise_PE --> PE_Vector
    Specific_PE --> PE_Vector
    
    PE_Vector --> KA_PE_Kernels
    KA_PE_Kernels --> Atomic_Energy
    Atomic_Energy --> Scalar_Sum
    Scalar_Sum --> Energy_Units
```

**Sources:** [src/energy.jl:256-289](), [src/kernels.jl:222-365]()

Energy kernels use `Atomix.@atomic` operations to safely accumulate energy contributions from multiple threads into a single scalar value.

## CUDA-Specific Optimizations

The `MollyCUDAExt` extension provides highly optimized CUDA kernels that leverage NVIDIA GPU-specific features for maximum performance.

### Morton Space-Filling Curve Ordering

CUDA kernels use Morton space-filling curves to improve spatial locality and memory access patterns:

```mermaid
graph TB
    subgraph "Morton Ordering Process"
        Coords["3D Coordinates"]
        Scale["Scale by Cutoff Distance"]
        Morton_Code["generalized_morton_code()"]
        Sort["Sort by Morton Code"]
        Reordered["Reordered Atom Indices"]
    end
    
    subgraph "Spatial Benefits"
        Cache_Locality["Improved Cache Locality"]
        Memory_Coalescing["Memory Coalescing"]
        Branch_Coherence["Branch Coherence"]
    end
    
    Coords --> Scale
    Scale --> Morton_Code  
    Morton_Code --> Sort
    Sort --> Reordered
    
    Reordered --> Cache_Locality
    Reordered --> Memory_Coalescing
    Reordered --> Branch_Coherence
```

**Sources:** [ext/MollyCUDAExt.jl:111-131]()

### Warp-Level Optimizations

CUDA kernels use warp-level primitives for efficient data sharing and reduced memory traffic:

- `CUDA.shfl_sync()`: Share data between threads in a warp without shared memory
- Diagonal force calculation patterns to maximize data reuse
- Atomic operations for safe force accumulation across warps

**Sources:** [ext/MollyCUDAExt.jl:325-606]()

### No-Neighbor-List Algorithm

For systems with dense interactions, CUDA kernels can bypass neighbor lists entirely using an optimized O(N²) algorithm with spatial partitioning:

```mermaid
graph LR
    subgraph "Spatial Partitioning"
        Atoms["N Atoms"]
        Blocks["32×32 Atom Blocks"]
        Bounds["Bounding Box Culling"]
        Interactions["Force Calculations"]
    end
    
    subgraph "Warp Execution"
        Warp_Thread["Warp Thread"]
        Row_Forces["Calculate Row Forces"]
        Shuffle["Shuffle Atom Data"]
        Accumulate["Accumulate Forces"]
    end
    
    Atoms --> Blocks
    Blocks --> Bounds
    Bounds --> Interactions
    
    Interactions --> Warp_Thread
    Warp_Thread --> Row_Forces
    Row_Forces --> Shuffle
    Shuffle --> Accumulate
```

**Sources:** [ext/MollyCUDAExt.jl:840-906]()

## GPU Neighbor Finding

GPU systems use specialized neighbor finding strategies depending on the hardware backend:

| Backend | Neighbor Finder | Algorithm |
|---------|----------------|-----------|
| CUDA | `GPUNeighborFinder` | Optimized spatial hashing |
| Other GPUs | `DistanceNeighborFinder` | Distance-based culling |
| CPU Fallback | Any neighbor finder | Standard CPU algorithms |

The `GPUNeighborFinder` is only available for CUDA arrays and provides the best performance for NVIDIA GPUs.

**Sources:** [src/kernels.jl:47-50](), [ext/MollyCUDAExt.jl:13]()

## Configuration and Performance Tuning

GPU acceleration can be configured through environment variables:

| Variable | Default | Purpose |
|----------|---------|---------|
| `MOLLY_GPUNTHREADS_PAIRWISE` | 512 | GPU threads for pairwise interactions |
| `MOLLY_GPUNTHREADS_SPECIFIC` | 128 | GPU threads for specific interactions |

The system automatically selects appropriate thread counts and block sizes based on GPU capabilities and interaction counts.

**Sources:** [src/kernels.jl:35-43](), [ext/MollyCUDAExt.jl:33-43]()

### Usage Example

To use GPU acceleration, simply pass GPU arrays to the system constructor:

```julia
# CPU simulation
system = System(atoms=atoms, coords=coords, ...)

# GPU simulation  
system = System(atoms=CuArray(atoms), coords=CuArray(coords), ...)
```

The simulation algorithms automatically detect GPU arrays and dispatch to appropriate GPU kernels without code changes.

**Sources:** [test/gradients.jl:89-90]()

---