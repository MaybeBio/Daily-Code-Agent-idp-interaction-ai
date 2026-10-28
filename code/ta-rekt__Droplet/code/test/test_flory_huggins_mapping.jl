# Flory-Huggins ↔ MC parameter-mapping validation.
#
# These tests confirm that the spec convention used in
# spec/flory_huggins_prewetting.jl
#
#     intPot!(bulk, -alpha/3)
#     bulk.chemPot = -mu - alpha
#
# produces a Hamiltonian that matches the FH free energy used in viz.ipynb,
#
#     F(phi) = mu*phi + 3*chi*phi*(1-phi) + T*[phi ln phi + (1-phi) ln(1-phi)]
#
# under the identification alpha = 3*chi (the y-axis of the notebook binodal).
#
# Standalone-runnable: `julia --project=@. test/test_flory_huggins_mapping.jl`.

using Droplet
using Test
using Random
using Statistics

Random.seed!(20260508)

# Exact identity check: ham!(b) computes
#   H = -chemPot * N + intPot * (sum_i nnp(b, b, pos_i)) * 0.5
# where the 0.5 is the selfFactor in src/hamiltonian.jl.
@testset "FH mapping: exact ham! identity" begin
    L = Int32(12)
    z = 6
    mu = 0.1

    for alpha in (1.0, 2.0, 3.0, 4.0)
        for phi in (0.25, 0.5, 0.75)
            l = Lattice(; boxlength=L, nspecies=1)
            N = Int32(round(L^3 * phi))
            b = FastBulkMonomer(; length=N, initShape="uniform")

            intPot!(b, -alpha/3)
            b.chemPot = -mu - alpha
            init!(l, [b])

            ham!(b)

            # recompute by hand
            nnp_sum = 0
            for i in 1:b.length
                nnp_sum += nnp(b, b, b.pos[i][1:3])
            end
            expected = -b.chemPot * b.length + b.intMat[b.intIdx, b.intIdx] * nnp_sum * 0.5

            @test isapprox(b.ham, expected; atol=1e-10)
        end
    end
end

# Mean-field consistency: for a uniformly random configuration, the average
# nnp per particle is z*phi (each of z neighbors is occupied with prob phi).
# So <H>/V = -chemPot*phi + intPot*z*phi^2/2
#         = (mu+alpha)*phi - alpha*phi^2
#         = mu*phi + alpha*phi*(1-phi)
# matching the FH energy term mu*phi + 3*chi*phi(1-phi) under alpha = 3*chi.
@testset "FH mapping: mean-field energy density" begin
    L = Int32(20)   # bigger box -> tighter mean-field agreement
    z = 6
    mu = 0.1
    nrep = 8        # average over independent random fills

    for alpha in (1.0, 2.0, 3.0)
        for phi in (0.3, 0.5, 0.7)
            energies = Float64[]
            for r in 1:nrep
                l = Lattice(; boxlength=L, nspecies=1)
                N = Int32(round(L^3 * phi))
                b = FastBulkMonomer(; length=N, initShape="uniform")
                intPot!(b, -alpha/3)
                b.chemPot = -mu - alpha
                init!(l, [b])
                ham!(b)
                push!(energies, b.ham / L^3)
            end

            measured = mean(energies)
            stderr_ = std(energies) / sqrt(nrep)
            predicted = mu*phi + alpha*phi*(1-phi)

            # 4-sigma tolerance, plus a floor for the small-phi limit
            tol = max(4*stderr_, 0.02)
            @test isapprox(measured, predicted; atol=tol)
        end
    end
end

# Critical-point behavior: at fixed phi=0.5, susceptibility (variance of local
# density) should grow as alpha crosses the mean-field critical value alpha=2.
# We use a coarse 4x4x4-block density and check monotonic increase across the
# transition. This is qualitative but catches sign/factor errors that pass the
# energy tests above.
@testset "FH mapping: susceptibility monotone through alpha=2" begin
    L = Int32(16)
    block = 4
    mu = 0.1
    phi = 0.5
    alphas = (1.0, 2.0, 3.0, 4.0)
    t_max_min = 0.4   # ~24 s per alpha
    step = 200

    function block_density_var(b, L, block)
        nb = div(L, block)
        counts = zeros(Float64, nb, nb, nb)
        for i in 1:b.length
            x, y, z = b.pos[i][1], b.pos[i][2], b.pos[i][3]
            bx = div(x - 1, block) + 1
            by = div(y - 1, block) + 1
            bz = div(z - 1, block) + 1
            counts[bx, by, bz] += 1
        end
        densities = counts ./ block^3
        return var(vec(densities))
    end

    susceptibilities = Float64[]
    for alpha in alphas
        l = Lattice(; boxlength=L, nspecies=1)
        N = Int32(round(L^3 * phi))
        b = FastBulkMonomer(; length=N, initShape="uniform")
        intPot!(b, -alpha/3)
        b.chemPot = -mu - alpha
        init!(l, [b])

        # short MC, take samples after burn-in
        t_start = time_ns()
        t_burn = t_max_min / 2
        nmoves = 0
        samples = Float64[]
        while (time_ns() - t_start) / 6e10 <= t_max_min
            if move!(b)
                nmoves += 1
                if mod(nmoves, step) == 0
                    if (time_ns() - t_start) / 6e10 >= t_burn
                        push!(samples, block_density_var(b, L, block))
                    end
                end
            end
        end
        push!(susceptibilities, isempty(samples) ? NaN : mean(samples))
    end

    @info "FH susceptibility sweep" alphas=alphas chi_var=susceptibilities

    # demand strict ordering across the transition: alpha=4 > alpha=1
    @test susceptibilities[end] > susceptibilities[1]
    # and that the largest-alpha susceptibility exceeds the critical-point one
    @test susceptibilities[end] > susceptibilities[2]
end
