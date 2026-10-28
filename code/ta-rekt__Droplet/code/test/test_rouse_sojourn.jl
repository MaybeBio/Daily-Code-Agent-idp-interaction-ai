using Droplet
using Test
using Random

# Tests for A/B sojourn-time tracking in the Rouse Langevin solver
# (src/rouse.jl), ported from viz_util/overdamped_langevin.py and the
# scalar_reference of test/test_sojourn_tracking.py.
#
# Geometry: the tracked-bead distance d is partitioned into
#     A = {d ≤ a}        B = {d ≥ b}        buffer (a, b),  a ≤ b
# with Schmitt-trigger (hysteresis) assignment. A sojourn that opens at step s
# and is closed by a crossing detected at step i has duration i - s; the first
# sojourn of each region enters at step 0. Durations are in integration steps.

# Faithful Julia transcription of the Python `scalar_reference` double loop,
# operating on a 1-D distance trace (single walker). `positions[1]` is step 0
# (matching distanceTraj[1] from init!); `positions[idx]` is step i = idx - 1.
function scalar_reference(positions::AbstractVector, a, b)
    n = length(positions)
    A = Int[]
    B = Int[]
    in_A = positions[1] <= a
    in_B = positions[1] > a
    A_entry = 0
    B_entry = 0
    for idx in 2:n
        i = idx - 1
        r = positions[idx]
        if r >= b && in_A
            B_entry = i
            push!(A, B_entry - A_entry)
            in_A = false
            in_B = true
        end
        if r <= a && in_B
            A_entry = i
            push!(B, A_entry - B_entry)
            in_B = false
            in_A = true
        end
    end
    return A, B
end

@testset verbose=true "Rouse sojourn tracking" begin

    # -- Handcrafted single-walker cases pin down every hysteresis rule.
    # (These exercise the reference logic; the real solver is validated against
    # the same reference on a live trajectory in the equivalence testset below.)
    @testset "handcrafted reference cases" begin
        a, b = -1.0, 1.0

        # A at 0, B at 3, A at 5, B at 7
        A, B = scalar_reference([-2.0, -0.5, 0.5, 1.5, 0.3, -1.5, 0.5, 1.2], a, b)
        @test A == [3, 2]
        @test B == [2]

        # never reaches b -> nothing recorded
        A, B = scalar_reference([-2.0, -0.9, -0.5, -3.0, -1.5], a, b)
        @test A == Int[]
        @test B == Int[]

        # buffer excursion does not reset the A entry
        A, B = scalar_reference([-2.0, 0.0, -1.5, 0.0, 1.5], a, b)
        @test A == [4]
        @test B == Int[]

        # touching the boundary (>= b, <= a) counts as a crossing
        A, B = scalar_reference([-1.0, 1.0, -1.0], a, b)
        @test A == [1]
        @test B == [1]

        # initial position in the buffer (r0 > a) is classified as B
        A, B = scalar_reference([0.0, -2.0], a, b)
        @test A == Int[]
        @test B == [1]
    end

    # -- The real solver's inlined tracker must equal the reference applied to
    # the trajectory it actually produces (bit-exact, on a live stochastic run).
    @testset "solver tracker ≡ reference on live trajectory" begin
        a, b = 2.5, 3.2
        r = RousePolymer(; nbeads = Int32(10), nsteps = Int32(50_000),
                           dt = 0.01 / 3, saveMode = :distance_only, seed = 23,
                           trackedBead1 = Int32(1), trackedBead2 = Int32(10),
                           a = a, b = b)
        init!(r); simulate!(r)

        A_ref, B_ref = scalar_reference(r.distanceTraj, a, b)
        @test r.A_intervals == A_ref
        @test r.B_intervals == B_ref

        # the run is meaningful (many crossings) and structurally sane
        @test length(r.A_intervals) + length(r.B_intervals) > 10
        @test all(r.A_intervals .>= 0)
        @test all(r.B_intervals .>= 0)
        # A/B sojourns strictly alternate, so the counts differ by at most one
        @test abs(length(r.A_intervals) - length(r.B_intervals)) <= 1
    end

    # -- :sojourn_only must reproduce the exact same sojourns as :distance_only
    # (tracking is mode-independent) while dropping the distance trace.
    @testset "sojourn_only ≡ distance_only sojourns, no distance trace" begin
        kwargs = (; nbeads = Int32(10), nsteps = Int32(20_000), dt = 0.01 / 3,
                   seed = 71, trackedBead1 = Int32(1), trackedBead2 = Int32(10),
                   a = 2.5, b = 3.2)
        rd = RousePolymer(; kwargs..., saveMode = :distance_only)
        rs = RousePolymer(; kwargs..., saveMode = :sojourn_only)
        init!(rd); simulate!(rd)
        init!(rs); simulate!(rs)

        @test rs.distanceTraj === nothing
        @test rd.distanceTraj !== nothing
        @test rs.A_intervals == rd.A_intervals
        @test rs.B_intervals == rd.B_intervals
    end

    @testset "no bounds ⇒ no tracking" begin
        r = RousePolymer(; nbeads = Int32(10), nsteps = Int32(500), dt = 1e-3,
                           seed = 5, trackedBead1 = Int32(1), trackedBead2 = Int32(10),
                           saveMode = :distance_only)
        init!(r); simulate!(r)
        @test r.A_intervals === nothing
        @test r.B_intervals === nothing
    end

    @testset "validation errors" begin
        # :sojourn_only requires both bounds
        r = RousePolymer(; nbeads = Int32(5), nsteps = Int32(10), dt = 1e-3,
                           saveMode = :sojourn_only)
        @test_throws ErrorException init!(r)

        # a <= b is required
        r2 = RousePolymer(; nbeads = Int32(5), nsteps = Int32(10), dt = 1e-3,
                            saveMode = :distance_only, a = 3.0, b = 1.0)
        @test_throws AssertionError init!(r2)
    end
end
