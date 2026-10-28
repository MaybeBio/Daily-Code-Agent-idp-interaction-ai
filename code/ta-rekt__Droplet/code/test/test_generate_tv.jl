using Test
using Droplet
using StatsBase
using JSON

@testset "generate_tv" begin
    # setup: 1bs circuit
    circuit = JSON.parsefile(joinpath(@__DIR__, "..", "spec", "circuits.json"))["1bs"]
    circuitId = "1bs"
    nstates = circuit["nstates"]
    edges = [Tuple(e) for e in circuit["edges"]]
    input_edges = [Tuple(e) for e in circuit["input_edges"]]
    off_indicator = circuit["off_indicator"]
    on_ind = [x == 0 for x in off_indicator]

    K = 0.5; eta_UA = -0.5; eta_BI = 0.0; epsilon = -7.75; phi_R = 0.02
    kboff = K; kuoff = K*exp(eta_BI); konb = K*phi_R; koffb = K*exp(eta_UA+epsilon)
    kuon = K*exp(eta_BI+epsilon); kbon = K; koffu = K*exp(eta_UA); konu = K*phi_R
    rates = [kboff, kuoff, konb, koffb, kuon, kbon, koffu, konu]

    N = nstates
    idist = ones(N) / N
    cr = 0.01

    function _make_Q(c)
        weighted_edges = [(e[1], e[2], rates[i]) for (i, e) in enumerate(edges)]
        build_mat(weighted_edges, input_edges, c)
    end

    @testset "constant c recovers generate" begin
        nruns = 500

        frac_const = zeros(nruns)
        Q0 = _make_Q(cr)
        for i in 1:nruns
            s, t = generate(100000, Q0, idist)
            frac_const[i] = sum(t[j] for j in eachindex(s) if on_ind[s[j]]) / sum(t)
        end

        frac_tv = zeros(nruns)
        for i in 1:nruns
            Q2 = _make_Q(cr)
            s, t = generate_tv(Q2, input_edges, N, fill(cr, 2000), 50.0, idist)
            frac_tv[i] = sum(t[j] for j in eachindex(s) if on_ind[s[j]]) / sum(t)
        end

        @test abs(mean(frac_tv) - mean(frac_const)) / mean(frac_const) < 0.05
    end

    @testset "total time matches c_series * T" begin
        T = 50.0
        c_series = fill(cr, 200)
        Q2 = _make_Q(cr)
        s, t = generate_tv(Q2, input_edges, N, c_series, T, idist)
        @test sum(t) <= length(c_series) * T + 1e-10
        @test sum(t) > (length(c_series) - 1) * T
    end

    @testset "length(s) == length(t)" begin
        Q2 = _make_Q(cr)
        s, t = generate_tv(Q2, input_edges, N, fill(cr, 100), 10.0, idist)
        @test length(s) == length(t)
    end

    @testset "rate change is reflected" begin
        nruns = 200
        T = 50.0
        n_half = 500
        c_series = vcat(fill(0.003, n_half), fill(0.018, n_half))

        frac_first = zeros(nruns)
        frac_second = zeros(nruns)
        for i in 1:nruns
            Q2 = _make_Q(c_series[1])
            s, t = generate_tv(Q2, input_edges, N, c_series, T, idist)

            # split trajectory at midpoint
            mid = n_half * T
            wall = 0.0
            on1 = 0.0; tot1 = 0.0; on2 = 0.0; tot2 = 0.0
            for j in eachindex(s)
                dt = t[j]
                is_on = on_ind[s[j]]
                if wall + dt <= mid
                    tot1 += dt; is_on && (on1 += dt)
                elseif wall >= mid
                    tot2 += dt; is_on && (on2 += dt)
                else
                    dt1 = mid - wall; dt2 = dt - dt1
                    tot1 += dt1; is_on && (on1 += dt1)
                    tot2 += dt2; is_on && (on2 += dt2)
                end
                wall += dt
            end
            frac_first[i] = on1 / max(tot1, 1e-30)
            frac_second[i] = on2 / max(tot2, 1e-30)
        end

        # higher concentration → higher ON fraction
        @test mean(frac_second) > mean(frac_first)
    end
end
