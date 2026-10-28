module Droplet

    using Random
    using Statistics
    using StatsBase
    using Distributions
    using Serialization
    using LinearAlgebra
    using SparseArrays
    using JSON
    using DataStructures
    using OrderedCollections
    using TimerOutputs
    using Printf
    using NonlinearSolve
    using NonlinearSolve.ADTypes: AutoFiniteDiff
    using NonlinearSolve.LinearSolve: KrylovJL_GMRES
    using NonlinearSolve.SciMLBase

    include("types.jl")
    include("interactions.jl")
    include("hamiltonian.jl")
    include("metropolis.jl")
    include("move.jl")
    include("port.jl")
    include("init.jl")
    include("utils.jl")
    include("gillespie.jl")
    include("stats.jl")
    include("rouse.jl")
    include("ctmc.jl")
    include("generate.jl")
    include("get_residence_time.jl")

    # 2D axisymmetric Flory-Huggins BVP (two spheres) + 1D isolated-sphere
    # reference solver.  Independent of the Metropolis core.
    include("bvp_solver.jl")
    include("bvp_two_spheres.jl")
    # Pseudo-arclength continuation of the two-sphere BVP in phi_inf: tracks the
    # solution branch THROUGH the capillary-bridging limit points, which the
    # prescribed-parameter solves in bvp_two_spheres.jl cannot.
    include("bvp_arclength.jl")

    # General 1D backward Fokker–Planck first-passage-time solver.  Independent
    # of the Metropolis core and of how the potential is built.
    include("fpt_solver.jl")

    export Parameters, with, fb, dfb_dphi, d2fb_dphi2, psi_linear,
           compute_mu_inf, compute_spinodal, is_spherical,
           IsolatedSphere1D, solve_isolated_sphere,
           PlanarSlab1D, solve_parallel_plates, psi_linear_plates,
           excess_free_energy, plates_pair_potential,
           Geometry, build_geometry, decay_length,
           scatter, scatter!, laplacian, laplacian!,
           make_residual, make_residual_inplace, build_linear_operator,
           initial_guess, solve_two_spheres,
           interp_phi, cell_fluid_fraction,
           phi_on_sphere, excess_free_energy_two_spheres,
           quadratic_ghost_coefficients, find_clean_stencil,
           ArcTrace, arclength_phi_inf, polish_at_phi_inf, gap_center_phi,
           FPTProblem, solve_fpt, boltzmann_ic, fpt_mfpt_quadrature

end # module Droplet
