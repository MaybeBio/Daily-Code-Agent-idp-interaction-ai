using Droplet
using Test
using Random

@testset "Droplet.jl" begin
    include("test_init.jl")
    include("test_hamiltonians.jl")
    include("test_connected_components.jl")
    include("test_density_droplet.jl")
    include("test_interaction_potentials.jl")
    include("test_polymonomer.jl")
    include("test_heterogeneous_phase0.jl")
    include("test_heterogeneous_phase2.jl")
    include("test_heterogeneous_phase3.jl")
    include("test_heterogeneous_phase5.jl")
    include("test_heterogeneous_phase6.jl")
    include("test_updateMonomers.jl")
    include("test_pattern_stats.jl")
    include("test_centerMassPattern.jl")
    include("test_contact_matrix.jl")
    include("test_find_patterns.jl")
    include("test_grow_pattern.jl")
    include("test_checkpoint.jl")
    include("test_restore_species.jl")
    include("test_rebuild_lattice.jl")
    include("test_restore.jl")
    include("test_integration_checkpoint.jl")
    # Three-way interaction tests
    include("test_utils.jl")
    include("test_couplingPot.jl")
    include("test_ham_coupling.jl")
    include("test_deltaHam_coupling.jl")
    include("test_integration_coupling.jl")
    # Time-varying Gillespie
    include("test_generate_tv.jl")
    # Lattice KMC (Metropolis-rate Gillespie) v1: BulkMonomer
    include("test_kmc_phase1.jl")
    include("test_kmc_phase2.jl")
    include("test_kmc_phase3.jl")
    include("test_kmc_phase4.jl")
    include("test_kmc_phase5.jl")
    # 2D axisymmetric two-sphere BVP
    include("test_bvp_two_spheres.jl")
    # Pseudo-arclength continuation of that BVP in phi_inf (through the folds)
    include("test_bvp_arclength.jl")
    # 1D planar slab BVP (parallel plates)
    include("test_bvp_parallel_plates.jl")
    # 1D isolated-sphere BVP vs. the prewetting graphical construction
    # (~7 min; set BVP_SCAN_TRIALS to shorten the robustness scan)
    include("test_bvp_solver.jl")
    # Rouse polymer overdamped Langevin solver
    include("test_rouse.jl")
    include("test_rouse_sojourn.jl")
    # 1D backward-FPE first-passage-time solver
    include("test_fpt_solver.jl")
end
