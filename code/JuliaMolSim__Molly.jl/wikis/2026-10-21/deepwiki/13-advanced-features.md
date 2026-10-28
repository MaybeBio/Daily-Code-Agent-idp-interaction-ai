# Advanced Features

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [README.md](README.md)
- [docs/src/differentiable.md](docs/src/differentiable.md)
- [docs/src/examples.md](docs/src/examples.md)
- [src/interactions/implicit_solvent.jl](src/interactions/implicit_solvent.jl)

</details>



This page covers sophisticated capabilities that extend basic molecular dynamics functionality in Molly.jl. These features enable high-performance computing, machine learning integration, and advanced simulation techniques that go beyond standard MD simulations.

For basic simulation setup and core molecular dynamics concepts, see [Core Architecture](#3) and [Simulation Engine](#5). For system configuration and force field setup, see [System Setup and Configuration](#7).

## Overview

Molly.jl provides several advanced capabilities that distinguish it from traditional molecular dynamics packages:

```mermaid
graph TB
    subgraph "Core MD System"
        CoreSim["simulate!()"]
        System["System"]
        Forces["Force Calculation"]
    end
    
    subgraph "Advanced Features"
        GPU["GPU Acceleration<br/>KernelAbstractions.jl"]
        Diff["Differentiable Simulation<br/>Enzyme.jl"]
        ImplicitSolv["Implicit Solvent<br/>ImplicitSolventOBC<br/>ImplicitSolventGBN2"]
        MLPot["ML Potentials<br/>Neural Networks"]
        ASECalc["ASE Integration<br/>ASECalculator"]
        HighPerf["High Performance<br/>Multithreading<br/>Float32 Support"]
    end
    
    subgraph "External Integrations"
        Python["Python Ecosystem<br/>PythonCall.jl"]
        Visualization["Advanced Visualization<br/>GLMakie.jl"]
        AtomsBases["Standard Interfaces<br/>AtomsBase.jl<br/>AtomsCalculators.jl"]
    end
    
    CoreSim --> GPU
    CoreSim --> Diff
    Forces --> ImplicitSolv
    System --> MLPot
    System --> ASECalc
    
    GPU --> HighPerf
    Diff --> MLPot
    ASECalc --> Python
    System --> Visualization
    System --> AtomsBases
    
    style CoreSim fill:#e1f5fe
    style GPU fill:#fff3e0
    style Diff fill:#e8f5e8
    style ImplicitSolv fill:#f3e5f5
```

**Advanced Features Architecture Overview**

Sources: [README.md:17-46](), [docs/src/differentiable.md:1-20]()

## GPU Acceleration

Molly.jl supports GPU acceleration across all major GPU backends through `KernelAbstractions.jl`. This includes NVIDIA CUDA, AMD ROCm, Apple Metal, and Intel oneAPI devices. GPU acceleration is automatically enabled when appropriate array types are used.

```mermaid
graph LR
    subgraph "GPU Backend Support"
        CUDA["MollyCUDAExt<br/>NVIDIA GPUs"]
        ROCm["GPU Support<br/>AMD ROCm"]
        Metal["GPU Support<br/>Apple Metal"]
        OneAPI["GPU Support<br/>Intel oneAPI"]
    end
    
    subgraph "Core Abstractions"
        KA["KernelAbstractions.jl<br/>Backend Abstraction"]
        ArrayTypes["Array Types<br/>CuArray, ROCArray"]
        Kernels["GPU Kernels<br/>Force Calculations"]
    end
    
    subgraph "Accelerated Components"
        PairwiseForces["Pairwise Interactions<br/>LennardJones, Coulomb"]
        NeighborLists["Neighbor Finding<br/>GPU Optimized"]
        ImplicitSolvent["Implicit Solvent<br/>GBSA Calculations"]
    end
    
    KA --> CUDA
    KA --> ROCm
    KA --> Metal
    KA --> OneAPI
    
    ArrayTypes --> Kernels
    Kernels --> PairwiseForces
    Kernels --> NeighborLists
    Kernels --> ImplicitSolvent
    
    KA --> ArrayTypes
```

**GPU Acceleration Architecture**

For detailed information on GPU setup, performance optimization, and backend-specific considerations, see [GPU Acceleration](#6.1).

Sources: [README.md:38-39]()

## Differentiable Simulation

Differentiable molecular dynamics enables automatic differentiation through entire simulations using `Enzyme.jl`. This unique capability allows for gradient-based optimization of force field parameters, neural network potentials, and other simulation parameters.

```mermaid
graph TB
    subgraph "Differentiable Components"
        EnzymeAD["Enzyme.jl<br/>Automatic Differentiation"]
        ForwardMode["Forward Mode AD<br/>Small Parameter Sets"]
        ReverseMode["Reverse Mode AD<br/>Large Parameter Sets"]
    end
    
    subgraph "Gradient-Enabled Simulations"
        ForceFieldOpt["Force Field<br/>Optimization"]
        NeuralPotentials["Neural Network<br/>Potentials"]
        ParameterTuning["Parameter<br/>Tuning"]
    end
    
    subgraph "Applications"
        GB99dms["GB99dms<br/>Force Field Training"]
        LossFunc["Custom Loss<br/>Functions"]
        MLTraining["Machine Learning<br/>Integration"]
    end
    
    EnzymeAD --> ForwardMode
    EnzymeAD --> ReverseMode
    
    ForwardMode --> ForceFieldOpt
    ReverseMode --> NeuralPotentials
    ReverseMode --> ParameterTuning
    
    ForceFieldOpt --> GB99dms
    NeuralPotentials --> MLTraining
    ParameterTuning --> LossFunc
```

**Differentiable Simulation Capabilities**

Key features include gradient calculation through:
- Pairwise and specific interactions
- Neighbor list updates  
- Stochastic thermostats
- Complex loss functions

For implementation details, training examples, and performance considerations, see [Differentiable Simulation](#6.2).

Sources: [README.md:45](), [docs/src/differentiable.md:9-21]()

## Implicit Solvent Models

Advanced implicit solvent models approximate solvent effects without explicit water molecules, significantly reducing computational cost while maintaining physical accuracy.

```mermaid
graph TB
    subgraph "GBSA Implementations"
        AbstractGBSA["AbstractGBSA<br/>Base Interface"]
        OBC["ImplicitSolventOBC<br/>Onufriev-Bashford-Case"]
        GBN2["ImplicitSolventGBN2<br/>GBn2 Model"]
    end
    
    subgraph "Core Calculations"
        BornRadii["born_radii_and_grad()<br/>Born Radii Calculation"]
        SurfaceArea["Surface Area<br/>ACE Approximation"]
        Electrostatics["Electrostatic<br/>Interactions"]
    end
    
    subgraph "GPU Optimization"
        BornGPU["gbsa_born_gpu()<br/>GPU Born Radii"]
        ForceGPU["gbsa_force_*_gpu()<br/>GPU Force Kernels"]
        KernelOpt["Optimized Kernels<br/>Memory Efficient"]
    end
    
    AbstractGBSA --> OBC
    AbstractGBSA --> GBN2
    
    OBC --> BornRadii
    GBN2 --> BornRadii
    BornRadii --> SurfaceArea
    BornRadii --> Electrostatics
    
    BornRadii --> BornGPU
    Electrostatics --> ForceGPU
    BornGPU --> KernelOpt
    ForceGPU --> KernelOpt
```

**Implicit Solvent Model Architecture**

Available models:
- **OBC (Onufriev-Bashford-Case)**: Parameter sets I and II with ACE surface area approximation
- **GBN2**: Advanced model with neck correction and optimized parameters

For model selection, parameterization, and usage examples, see [Implicit Solvent Models](#6.3).

Sources: [src/interactions/implicit_solvent.jl:1-600]()

## Machine Learning Integration

Molly.jl supports integration with machine learning frameworks and external calculators for advanced potential energy surfaces and force calculations.

| Integration Type | Implementation | Use Case |
|------------------|----------------|----------|
| Neural Network Potentials | `Flux.jl` models as general interactions | Custom ML potentials |
| Python ASE Calculators | `ASECalculator` with `PythonCall.jl` | MACE, psi4, other ASE-compatible tools |
| AtomsCalculators.jl | `MollyCalculator` interface | Cross-package compatibility |
| Custom Interactions | User-defined force functions | Domain-specific physics |

Example neural network integration pattern:
```julia
# Neural network as general interaction
struct NeuralNetworkPotential{T}
    model::T
end

function AtomsCalculators.forces(sys, inter::NeuralNetworkPotential; kwargs...)
    # Custom force calculation using neural network
    # Returns force vector for each atom
end
```

Sources: [docs/src/examples.md:442-535](), [README.md:42-44]()

## Performance and Hardware Optimization

Advanced performance features enable efficient simulations across different hardware configurations:

```mermaid
graph LR
    subgraph "Precision Control"
        Float64["Float64<br/>High Precision"]
        Float32["Float32<br/>Memory Efficient"]
        MixedPrec["Mixed Precision<br/>Optimized Performance"]
    end
    
    subgraph "Parallelization"
        MultiThread["Automatic<br/>Multithreading"]
        GPU["GPU Backends<br/>Massive Parallelism"]
        SIMD["SIMD Instructions<br/>Vectorization"]
    end
    
    subgraph "Memory Optimization"
        NeighborOpt["Optimized Neighbor Lists<br/>CellListMap.jl"]
        ArrayReuse["Array Reuse<br/>Minimal Allocation"]
        Streaming["GPU Memory<br/>Streaming"]
    end
    
    Float64 --> MultiThread
    Float32 --> GPU
    MixedPrec --> SIMD
    
    MultiThread --> NeighborOpt
    GPU --> ArrayReuse
    SIMD --> Streaming
```

**Performance Optimization Strategies**

Key optimization features:
- **Precision flexibility**: Support for `Float64`, `Float32`, and custom numeric types
- **Automatic multithreading**: Thread-safe parallel force calculations
- **Memory management**: Efficient GPU memory usage and CPU cache optimization
- **Neighbor list optimization**: Integration with `CellListMap.jl` for efficient spatial searches

Sources: [README.md:37-40]()

## External Ecosystem Integration

Molly.jl integrates with the broader Julia and scientific computing ecosystems:

| Package | Integration | Purpose |
|---------|-------------|---------|
| `AtomsBase.jl` | Standard atomic system interface | Interoperability |
| `AtomsCalculators.jl` | Calculator interface | Force/energy computation |
| `Chemfiles.jl` | File I/O support | Trajectory reading/writing |
| `Unitful.jl` | Physical units | Type-safe calculations |
| `GLMakie.jl` | Visualization | Interactive simulation viewing |
| `SimpleCrystals.jl` | Crystal structure generation | Solid-state systems |

This ecosystem integration enables Molly.jl to serve as both a standalone simulation package and a component in larger computational workflows.

Sources: [README.md:21-22](), [README.md:35-36](), [README.md:41-42]()

---