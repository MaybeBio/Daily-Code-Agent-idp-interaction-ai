# Force and Energy Calculations

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [src/energy.jl](src/energy.jl)
- [src/force.jl](src/force.jl)

</details>



This page covers how forces and potential energies are computed from molecular interactions in Molly.jl. It details the computational engine that takes interaction definitions and calculates the resulting forces on atoms and system energies.

For information about the specific interaction types themselves, see [Pairwise Interactions](#4.1) and [Bonded Interactions](#4.2). For how these calculated forces are used in simulation algorithms, see [Integrators and Simulators](#5.1).

## Overview

Force and energy calculations in Molly.jl follow a modular architecture that separates the computational engine from the interaction definitions. The system supports three categories of interactions:

- **Pairwise interactions**: Non-bonded interactions like Lennard-Jones and Coulomb
- **Specific interactions**: Bonded interactions with fixed atom lists (bonds, angles, torsions)  
- **General interactions**: Complex interactions like implicit solvent models

All calculations support both CPU execution (with multi-threading) and GPU acceleration, with automatic unit handling throughout.

## Force Calculation Architecture

### Main Calculation Flow

The primary entry point for force calculations is the `forces()` function, which orchestrates the computation across all interaction types:

```mermaid
flowchart TD
    forces["forces()"] --> find_neighbors["find_neighbors()"]
    forces --> forces_nounits["forces_nounits!()"]
    
    forces_nounits --> filter_pairwise["Filter pairwise interactions"]
    filter_pairwise --> pairwise_nonl["Pairwise (no neighbor list)"]
    filter_pairwise --> pairwise_nl["Pairwise (with neighbor list)"]
    
    forces_nounits --> filter_specific["Filter specific interactions"]
    filter_specific --> sils_1["InteractionList1Atoms"]
    filter_specific --> sils_2["InteractionList2Atoms"] 
    filter_specific --> sils_3["InteractionList3Atoms"]
    filter_specific --> sils_4["InteractionList4Atoms"]
    
    forces_nounits --> general_inters["General interactions"]
    
    pairwise_nonl --> pairwise_forces["pairwise_forces!()"]
    pairwise_nl --> pairwise_forces
    
    sils_1 --> specific_forces["specific_forces!()"]
    sils_2 --> specific_forces
    sils_3 --> specific_forces  
    sils_4 --> specific_forces
    
    general_inters --> atoms_calc["AtomsCalculators.forces()"]
    
    pairwise_forces --> accumulate["Force accumulation"]
    specific_forces --> accumulate
    atoms_calc --> accumulate
    
    accumulate --> result["Final forces with units"]
```

Sources: [src/force.jl:159-168](), [src/force.jl:170-206]()

### Force Calculation Methods

Individual interactions implement the `force()` function with multiple dispatch based on the number of atoms involved:

| Function Signature | Interaction Type | Returns |
|-------------------|------------------|---------|
| `force(inter, vec_ij, atom_i, atom_j, ...)` | Pairwise | Single force vector |
| `force(inter, coord_i, boundary, atom_i, ...)` | 1-atom specific | `SpecificForce1Atoms` |
| `force(inter, coord_i, coord_j, boundary, ...)` | 2-atom specific | `SpecificForce2Atoms` |
| `force(inter, coord_i, coord_j, coord_k, ...)` | 3-atom specific | `SpecificForce3Atoms` |
| `force(inter, coord_i, coord_j, coord_k, coord_l, ...)` | 4-atom specific | `SpecificForce4Atoms` |

The `SpecificForceNAtoms` structs contain force vectors for each participating atom:

```mermaid
classDiagram
    class SpecificForce2Atoms {
        +SVector f1
        +SVector f2
    }
    
    class SpecificForce3Atoms {
        +SVector f1  
        +SVector f2
        +SVector f3
    }
    
    class SpecificForce4Atoms {
        +SVector f1
        +SVector f2 
        +SVector f3
        +SVector f4
    }
```

Sources: [src/force.jl:44-51](), [src/force.jl:54-93]()

### CPU vs GPU Execution Paths

The system automatically selects appropriate calculation methods based on the array type:

```mermaid
flowchart LR
    system_coords["System.coords"] --> check_type{Array type?}
    
    check_type -->|"CPU Array"| cpu_path["CPU calculation path"]
    check_type -->|"AbstractGPUArray"| gpu_path["GPU calculation path"]
    
    cpu_path --> threading{"n_threads > 1?"}
    threading -->|"Yes"| pairwise_forces_threads["pairwise_forces_threads!()"]
    threading -->|"No"| pairwise_forces_single["pairwise_forces!()"]
    
    gpu_path --> forces_buffer["ForcesBuffer allocation"]
    forces_buffer --> pairwise_force_gpu["pairwise_force_gpu!()"]
    forces_buffer --> specific_force_gpu["specific_force_gpu!()"]
    
    pairwise_forces_threads --> cpu_result["CPU forces"]
    pairwise_forces_single --> cpu_result
    pairwise_force_gpu --> gpu_result["GPU forces"]
    specific_force_gpu --> gpu_result
```

Sources: [src/force.jl:372-405](), [src/force.jl:253-316]()

## Energy Calculation Architecture

### Energy Components

Total energy calculation combines kinetic and potential energy components:

```mermaid
flowchart TD
    total_energy["total_energy()"] --> kinetic_energy["kinetic_energy()"]
    total_energy --> potential_energy["potential_energy()"]
    
    kinetic_energy --> ke_formula["E_k = 1/2 * Σ(m_i * v_i²)"]
    ke_formula --> temperature["temperature()"]
    
    potential_energy --> pairwise_pe["Pairwise PE"]
    potential_energy --> specific_pe["Specific PE"] 
    potential_energy --> general_pe["General PE"]
    
    pairwise_pe --> filter_pe_interactions["Filter interactions by neighbor usage"]
    specific_pe --> filter_pe_lists["Filter by interaction list type"]
    general_pe --> atoms_calc_pe["AtomsCalculators.potential_energy()"]
    
    filter_pe_interactions --> pe_calculation["pairwise_pe() / pairwise_pe_threads()"]
    filter_pe_lists --> pe_calculation_specific["specific_pe()"]
    
    pe_calculation --> sum_pe["Sum all PE components"]
    pe_calculation_specific --> sum_pe
    atoms_calc_pe --> sum_pe
    
    sum_pe --> total_pe["Total potential energy"]
```

Sources: [src/energy.jl:15-21](), [src/energy.jl:77-118]()

### Potential Energy Calculation Methods

Similar to force calculations, interactions implement `potential_energy()` methods:

```mermaid
flowchart LR
    pe_dispatch["potential_energy() dispatch"] --> pairwise["Pairwise interactions"]
    pe_dispatch --> specific_1["1-atom specific"]
    pe_dispatch --> specific_2["2-atom specific"] 
    pe_dispatch --> specific_3["3-atom specific"]
    pe_dispatch --> specific_4["4-atom specific"]
    
    pairwise --> pairwise_pe_calc["pairwise_pe() / pairwise_pe_threads()"]
    specific_1 --> specific_pe_calc["specific_pe()"]
    specific_2 --> specific_pe_calc
    specific_3 --> specific_pe_calc
    specific_4 --> specific_pe_calc
    
    pairwise_pe_calc --> energy_sum["Energy accumulation"]
    specific_pe_calc --> energy_sum
```

Sources: [src/energy.jl:120-160](), [src/energy.jl:209-254]()

## GPU Acceleration

### GPU-Specific Data Structures

GPU calculations use specialized buffers and data structures for efficient memory management:

```mermaid
classDiagram
    class ForcesBuffer {
        +fs_mat : GPU matrix
        +box_mins : GPU array  
        +box_maxs : GPU array
        +morton_seq : GPU array
        +compressed_eligible : GPU array
        +compressed_special : GPU array
    }
    
    class GPUNeighborFinder {
        +initialized : Bool
    }
    
    ForcesBuffer --> "manages" : GPUNeighborFinder
```

Sources: [src/force.jl:126-151]()

### GPU Function Variants

The system provides GPU-specific variants of force and energy functions that can be overridden for optimization:

```mermaid
flowchart TD
    interaction["Interaction type"] --> cpu_default["Default CPU implementation"]
    interaction --> gpu_override{"GPU override available?"}
    
    gpu_override -->|"Yes"| force_gpu["force_gpu()"]
    gpu_override -->|"No"| force_fallback["force() fallback"]
    
    gpu_override -->|"Yes"| pe_gpu["potential_energy_gpu()"] 
    gpu_override -->|"No"| pe_fallback["potential_energy() fallback"]
    
    force_gpu --> gpu_kernels["GPU kernel execution"]
    pe_gpu --> gpu_kernels
```

Sources: [src/force.jl:46-51](), [src/energy.jl:291-296]()

## Threading and Performance

### Multi-Threading Strategy

CPU calculations support multi-threading through work distribution across threads:

```mermaid
flowchart TD
    n_threads{"n_threads > 1?"} -->|"Yes"| thread_setup["Initialize force chunks"]
    n_threads -->|"No"| single_thread["Single-threaded execution"]
    
    thread_setup --> thread_dispatch["Threads.@threads dispatch"]
    thread_dispatch --> chunk_work["Each thread works on chunk_i:n_threads:n_atoms"]
    
    chunk_work --> accumulate_chunks["Accumulate results from all chunks"]
    single_thread --> direct_calc["Direct calculation"]
    
    accumulate_chunks --> final_forces["Final forces"]
    direct_calc --> final_forces
```

Sources: [src/force.jl:253-316](), [src/energy.jl:162-207]()

### Performance Considerations

| Aspect | Implementation | Performance Impact |
|--------|---------------|-------------------|
| Neighbor lists | Used by most pairwise interactions | O(N) vs O(N²) scaling |
| Threading | Work distribution by atom index | Near-linear speedup with cores |
| GPU acceleration | Kernel-based computation | High throughput for large systems |
| Unit handling | Stripped during calculation, restored at end | Minimal overhead |
| Memory allocation | Buffer reuse on GPU, chunks on CPU | Reduced GC pressure |

Sources: [src/force.jl:118-124](), [src/force.jl:135-151]()

## Unit Handling

Forces and energies maintain proper units throughout the calculation process:

```mermaid
flowchart LR
    input["Input coordinates with units"] --> strip["ustrip_vec.()"]
    strip --> calc_nounits["Calculation without units"]
    calc_nounits --> check["check_force_units() / check_energy_units()"]
    check --> restore["Restore units"]
    restore --> output["Output with proper units"]
    
    system_units["System.force_units / System.energy_units"] --> restore
```

The system ensures unit consistency by:
- Stripping units before calculation for performance
- Validating unit compatibility during accumulation  
- Restoring proper units to final results
- Supporting both dimensional and dimensionless calculations

Sources: [src/force.jl:164-167](), [src/force.jl:222-224](), [src/energy.jl:36-39]()

---