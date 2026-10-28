# Testing Framework

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [test/basic.jl](test/basic.jl)
- [test/energy_conservation.jl](test/energy_conservation.jl)
- [test/protein.jl](test/protein.jl)
- [test/runtests.jl](test/runtests.jl)

</details>



## Purpose and Scope

The testing framework in Molly.jl provides comprehensive validation of molecular dynamics functionality across different hardware platforms, numerical precision levels, and use cases. This framework ensures correctness, performance, and compatibility of simulations on CPU and GPU architectures while maintaining numerical accuracy and energy conservation properties.

For information about performance benchmarking and optimization strategies, see [Benchmarking and Performance](#9.2).

## Test Organization and Structure

The test suite is organized into specialized test files that cover different aspects of the molecular dynamics package:

```mermaid
graph TD
    runtests["runtests.jl<br/>Main test runner"]
    basic["basic.jl<br/>Core functionality"]
    protein["protein.jl<br/>Protein systems"]
    interactions["interactions.jl<br/>Force calculations"]
    simulation["simulation.jl<br/>Integration methods"]
    gradients["gradients.jl<br/>Differentiable MD"]
    agent["agent.jl<br/>Agent-based modeling"]
    minimization["minimization.jl<br/>Energy minimization"]
    energy_conservation["energy_conservation.jl<br/>Energy conservation"]
    
    runtests --> basic
    runtests --> protein
    runtests --> interactions
    runtests --> simulation
    runtests --> gradients
    runtests --> agent
    runtests --> minimization
    runtests --> energy_conservation
    
    basic --> spatial["Spatial calculations"]
    basic --> neighbors["Neighbor finding"]
    basic --> analysis["Analysis functions"]
    basic --> replica["Replica systems"]
    
    protein --> peptide["Peptide simulations"]
    protein --> openmm["OpenMM comparisons"]
    protein --> implicit["Implicit solvent"]
```

The main test runner coordinates execution based on environment variables and available hardware capabilities.

Sources: [test/runtests.jl:1-136](), [test/basic.jl:1-524](), [test/protein.jl:1-298]()

## Test Groups and Conditional Execution

Tests are organized into groups that can be selectively executed based on the `GROUP` environment variable:

| Group | Description | Files Included |
|-------|-------------|----------------|
| `All` | Complete test suite | All test files |
| `Protein` | Protein-specific tests | `protein.jl` |
| `Gradients` | Differentiable simulation tests | `gradients.jl` |
| `NotGradients` | All tests except gradients | `basic.jl`, `interactions.jl`, `simulation.jl`, `agent.jl`, `minimization.jl`, `protein.jl` |

```mermaid
flowchart TD
    start["Test Execution Start"]
    group_check{"GROUP environment<br/>variable check"}
    all_tests["All Tests<br/>GROUP = 'All'"]
    protein_tests["Protein Tests<br/>GROUP = 'Protein'"]
    gradient_tests["Gradient Tests<br/>GROUP = 'Gradients'"]
    no_gradient_tests["Non-Gradient Tests<br/>GROUP = 'NotGradients'"]
    error["Error:<br/>Unrecognized GROUP"]
    
    start --> group_check
    group_check --> |"All"| all_tests
    group_check --> |"Protein"| protein_tests
    group_check --> |"Gradients"| gradient_tests
    group_check --> |"NotGradients"| no_gradient_tests
    group_check --> |"Other"| error
    
    all_tests --> aqua["Aqua.jl quality checks"]
    all_tests --> basic_tests["basic.jl"]
    all_tests --> other_tests["Other test files"]
```

Sources: [test/runtests.jl:29-35](), [test/runtests.jl:113-136]()

## Hardware Platform Testing

The testing framework supports multiple hardware platforms and array types through conditional execution based on hardware availability:

### CPU Testing
- **Single-threaded**: Always available with `Array` type
- **Multi-threaded**: Conditional on `Threads.nthreads() > 1`
- **Parallel tests**: Executed when Julia runs with multiple threads

### GPU Testing
GPU tests are controlled by the `GPUTESTS` environment variable and hardware availability:

```mermaid
graph TD
    gpu_check["GPU Tests Enabled?<br/>GPUTESTS != '0'"]
    
    subgraph "GPU Platform Detection"
        cuda_check["CUDA.functional()"]
        rocm_check["AMDGPU.functional()"]
        oneapi_check["oneAPI.functional()"]
        metal_check["Metal.functional()"]
    end
    
    subgraph "Array Types"
        cpu_array["Array<br/>(CPU)"]
        cu_array["CuArray<br/>(NVIDIA)"]
        roc_array["ROCArray<br/>(AMD)"]
        one_array["oneArray<br/>(Intel)"]
        metal_array["MtlArray<br/>(Apple)"]
    end
    
    gpu_check --> cuda_check
    gpu_check --> rocm_check
    gpu_check --> oneapi_check
    gpu_check --> metal_check
    
    cuda_check --> |"Available"| cu_array
    rocm_check --> |"Available"| roc_array
    oneapi_check --> |"Available"| one_array
    metal_check --> |"Available"| metal_array
    
    cpu_array --> test_execution["Test Execution<br/>on Platform"]
    cu_array --> test_execution
    roc_array --> test_execution
    one_array --> test_execution
    metal_array --> test_execution
```

The framework maintains an `array_list` tuple that contains all available array types for testing across platforms.

Sources: [test/runtests.jl:58-98](), [test/runtests.jl:67-98]()

## Test Environment Configuration

### Environment Variables
Key environment variables control test behavior:

| Variable | Default | Purpose |
|----------|---------|---------|
| `GROUP` | `"All"` | Selects test group to run |
| `GPUTESTS` | `"1"` | Enables/disables GPU tests |
| `DEVICE` | `"0"` | Specifies GPU device number |
| `VISTESTS` | `"1"` | Enables/disables visualization tests |
| `CI` | - | Indicates CI environment (affects memory-intensive tests) |

### CI-Specific Behavior
Tests adapt to CI environments to handle memory constraints:

```julia
const running_CI = haskey(ENV, "CI")
if running_CI
    @warn "Some CPU gradient tests will not be run as this is CI"
end
```

Sources: [test/runtests.jl:37-41](), [test/runtests.jl:29-35](), [test/runtests.jl:58-61]()

## Core Test Categories

### Spatial and Boundary Tests
Tests for coordinate systems, periodic boundaries, and spatial utility functions:

```mermaid
graph LR
    subgraph "Spatial Tests"
        vector_1d["vector_1D<br/>Periodic distance"]
        vector_3d["vector<br/>3D periodic vectors"]
        wrap_coord["wrap_coord_1D<br/>Coordinate wrapping"]
        boundaries["Boundary Types<br/>Cubic/Rectangular/Triclinic"]
    end
    
    subgraph "Boundary Tests"
        cubic["CubicBoundary<br/>Tests"]
        rectangular["RectangularBoundary<br/>Tests"]
        triclinic["TriclinicBoundary<br/>Tests"]
    end
    
    vector_1d --> boundaries
    vector_3d --> boundaries
    wrap_coord --> boundaries
    boundaries --> cubic
    boundaries --> rectangular
    boundaries --> triclinic
```

Sources: [test/basic.jl:1-216]()

### Neighbor Finding Tests
Comprehensive testing of neighbor finding algorithms:

- **DistanceNeighborFinder**: Brute force O(N²) algorithm
- **TreeNeighborFinder**: Tree-based spatial partitioning
- **CellListMapNeighborFinder**: Cell list algorithm
- **GPUNeighborFinder**: GPU-optimized neighbor finding

Tests verify consistency across different algorithms and hardware platforms.

Sources: [test/basic.jl:218-338]()

### Protein System Tests
Complex molecular system validation including:

- **Peptide simulations**: Small protein systems with full force fields
- **OpenMM comparisons**: Force and energy validation against OpenMM reference
- **Implicit solvent**: OBC2 and GBN2 implicit solvent models
- **Multi-precision**: Float32 and Float64 testing

Sources: [test/protein.jl:1-298]()

## Energy Conservation Testing

Dedicated energy conservation tests ensure numerical stability:

```mermaid
graph TD
    energy_test["Energy Conservation Test"]
    
    subgraph "Test Parameters"
        system_size["2000 atoms<br/>Lennard-Jones"]
        cutoffs["Multiple Cutoff Types<br/>Distance/Shifted/Spline"]
        platforms["CPU/GPU Platforms"]
    end
    
    subgraph "Validation Criteria"
        initial_energy["E₀ = total_energy(sys)"]
        final_energies["E_final after simulation"]
        tolerance["max|ΔE| < 5e-4 kJ/mol"]
    end
    
    energy_test --> system_size
    energy_test --> cutoffs
    energy_test --> platforms
    
    system_size --> initial_energy
    cutoffs --> initial_energy
    platforms --> initial_energy
    
    initial_energy --> final_energies
    final_energies --> tolerance
```

Sources: [test/energy_conservation.jl:1-95]()

## External Integration Testing

### AtomsBase.jl Integration
Tests conversion between Molly systems and AtomsBase standard:

```julia
# AtomsBase conversion testing
ab_sys = make_test_system().system
molly_sys = System(ab_sys; energy_units=u"kJ", force_units=u"kJ/Å")
test_approx_eq(ab_sys, molly_sys; common_only=true)
```

### AtomsCalculators.jl Interface
Validates calculator interface compliance:

```julia
calc = MollyCalculator(pairwise_inters=(coul,), force_units=u"kJ/Å", energy_units=u"kJ")
pe = AtomsCalculators.potential_energy(ab_sys, calc)
fs = AtomsCalculators.forces(ab_sys, calc)
```

Sources: [test/basic.jl:483-523]()

## Test Data and Fixtures

### Test Data Organization
```mermaid
graph TD
    data_dir["data/<br/>Test data directory"]
    
    subgraph "Force Fields"
        ff_dir["force_fields/<br/>XML force field files"]
        ff99["ff99SBildn.xml"]
        tip3p["tip3p_standard.xml"]
        his["his.xml"]
    end
    
    subgraph "Reference Data"
        openmm_dir["openmm_6mrr/<br/>OpenMM reference"]
        coords_ref["Coordinate references"]
        energy_ref["Energy references"]
        force_ref["Force references"]
    end
    
    subgraph "Temporary Files"
        temp_files["Temporary test files<br/>DCD, PDB, XYZ, etc."]
    end
    
    data_dir --> ff_dir
    data_dir --> openmm_dir
    data_dir --> temp_files
    
    ff_dir --> ff99
    ff_dir --> tip3p
    ff_dir --> his
    
    openmm_dir --> coords_ref
    openmm_dir --> energy_ref
    openmm_dir --> force_ref
```

Sources: [test/runtests.jl:100-109]()

---