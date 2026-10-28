using Test
using LinearAlgebra

"""
Compute closed-form rates from energy parameters.
(Standalone copy for testing — must stay in sync with dev.jl)
"""
function rates_closed_form(phi_R, phi_T, eta_BI, eta_UA, epsilon, K=1.0)
    a = exp(-eta_BI)
    b = exp(-eta_UA)
    g = exp(-epsilon)

    k_u_off = K
    k_off_b = K
    k_u_on  = K
    k_off_u = K

    k_b_off = a * K
    k_on_b  = b * g * K
    k_b_on  = a * g * K
    k_on_u  = b * K

    return [k_b_off, k_u_off, k_on_b, k_off_b, k_u_on, k_b_on, k_off_u, k_on_u]
end

"""
Build the 4×4 rate matrix Q for states {UI, BI, BA, UA}.
Full transition rates include concentration prefactors ϕ_T, ϕ_R.
"""
function build_Q(phi_R, phi_T, rates)
    k_b_off, k_u_off, k_on_b, k_off_b, k_u_on, k_b_on, k_off_u, k_on_u = rates

    # States: 1=UI, 2=BI, 3=BA, 4=UA
    Q = zeros(4, 4)
    Q[1, 2] = phi_T * k_b_off   # UI → BI
    Q[2, 1] = k_u_off            # BI → UI
    Q[2, 3] = phi_R * k_on_b     # BI → BA
    Q[3, 2] = k_off_b            # BA → BI
    Q[3, 4] = k_u_on             # BA → UA
    Q[4, 3] = phi_T * k_b_on     # UA → BA
    Q[4, 1] = k_off_u            # UA → UI
    Q[1, 4] = phi_R * k_on_u     # UI → UA

    # Diagonal: negative sum of outgoing rates
    for i in 1:4
        Q[i, i] = -sum(Q[i, :])
    end
    return Q
end

"""
Solve for steady-state distribution from rate matrix Q.
Convention: Q[i,j] = rate from state i to state j.
Steady state satisfies πQ = 0, i.e. Q^T π = 0 (column vector).
"""
function steady_state(Q)
    A = copy(transpose(Q))
    A[end, :] .= 1.0
    rhs = zeros(size(A, 1))
    rhs[end] = 1.0
    return A \ rhs
end

# ──────────────────────────────────────────────────────────
# Test parameters
# ──────────────────────────────────────────────────────────
test_params = [
    # (phi_R, phi_T, eta_BI, eta_UA, epsilon, description)
    (0.02,  0.5,  0.0,  -0.5, -7.75,  "dev.jl defaults"),
    (0.1,   0.1,  -1.0, -1.0, -2.0,   "moderate energies"),
    (0.5,   0.5,  0.0,   0.0,  0.0,   "zero energies"),
    (0.01,  2.0,  -3.0, -0.1, -5.0,   "strong η_BI"),
    (1.0,   1.0,  -1.0, -1.0, -1.0,   "unit concentrations"),
    (0.001, 0.001, 0.5,  0.5,  0.5,   "positive energies"),
]

@testset "rates_closed_form" begin

    @testset "Positivity — $desc" for (phi_R, phi_T, eta_BI, eta_UA, eps_, desc) in test_params
        for K in [0.1, 1.0, 10.0]
            rates = rates_closed_form(phi_R, phi_T, eta_BI, eta_UA, eps_, K)
            for (i, r) in enumerate(rates)
                @test r > 0
            end
        end
    end

    @testset "Detailed balance — $desc" for (phi_R, phi_T, eta_BI, eta_UA, eps_, desc) in test_params
        a = exp(-eta_BI)
        b = exp(-eta_UA)
        g = exp(-eps_)

        rates = rates_closed_form(phi_R, phi_T, eta_BI, eta_UA, eps_)
        k_b_off, k_u_off, k_on_b, k_off_b, k_u_on, k_b_on, k_off_u, k_on_u = rates

        @test k_b_off / k_u_off ≈ a       atol=1e-14   # DB1
        @test k_on_b  / k_off_b ≈ b * g   atol=1e-14   # DB2
        @test k_b_on  / k_u_on  ≈ a * g   atol=1e-14   # DB3
        @test k_on_u  / k_off_u ≈ b       atol=1e-14   # DB4
    end

    @testset "Kolmogorov cycle condition — $desc" for (phi_R, phi_T, eta_BI, eta_UA, eps_, desc) in test_params
        rates = rates_closed_form(phi_R, phi_T, eta_BI, eta_UA, eps_)
        k_b_off, k_u_off, k_on_b, k_off_b, k_u_on, k_b_on, k_off_u, k_on_u = rates

        # Product of full rates around the cycle UI→BI→BA→UA→UI
        forward  = (phi_T * k_b_off) * (phi_R * k_on_b) * k_u_on * k_off_u
        backward = (phi_R * k_on_u) * (phi_T * k_b_on) * k_off_b * k_u_off

        @test forward ≈ backward  rtol=1e-12
    end

    @testset "Steady state = Boltzmann — $desc" for (phi_R, phi_T, eta_BI, eta_UA, eps_, desc) in test_params
        a = exp(-eta_BI)
        b = exp(-eta_UA)
        g = exp(-eps_)

        rates = rates_closed_form(phi_R, phi_T, eta_BI, eta_UA, eps_)
        Q = build_Q(phi_R, phi_T, rates)
        pi_ss = steady_state(Q)

        # Expected Boltzmann weights
        w = [1.0, phi_T * a, phi_T * phi_R * a * b * g, phi_R * b]
        # States: 1=UI, 2=BI, 3=BA, 4=UA
        Z = sum(w)
        pi_expected = w ./ Z

        for i in 1:4
            @test pi_ss[i] ≈ pi_expected[i]  rtol=1e-10
        end
    end

    @testset "Zero entropy production — $desc" for (phi_R, phi_T, eta_BI, eta_UA, eps_, desc) in test_params
        rates = rates_closed_form(phi_R, phi_T, eta_BI, eta_UA, eps_)
        Q = build_Q(phi_R, phi_T, rates)
        pi_ss = steady_state(Q)

        # EP = Σ_{i≠j} J_ij ln(π_i Q_ij / π_j Q_ji) where J_ij = π_i Q_ij - π_j Q_ji
        EP = 0.0
        for i in 1:4, j in 1:4
            if i != j && Q[i, j] > 0 && Q[j, i] > 0
                Jij = pi_ss[i] * Q[i, j] - pi_ss[j] * Q[j, i]
                EP += Jij * log(pi_ss[i] * Q[i, j] / (pi_ss[j] * Q[j, i]))
            end
        end
        @test abs(EP) < 1e-12
    end

    @testset "Linear scaling with K — $desc" for (phi_R, phi_T, eta_BI, eta_UA, eps_, desc) in test_params
        rates1 = rates_closed_form(phi_R, phi_T, eta_BI, eta_UA, eps_, 1.0)
        rates2 = rates_closed_form(phi_R, phi_T, eta_BI, eta_UA, eps_, 3.7)

        for i in eachindex(rates1)
            @test rates2[i] / rates1[i] ≈ 3.7  rtol=1e-14
        end
    end

    @testset "Flux balance (master equation) — $desc" for (phi_R, phi_T, eta_BI, eta_UA, eps_, desc) in test_params
        rates = rates_closed_form(phi_R, phi_T, eta_BI, eta_UA, eps_)
        Q = build_Q(phi_R, phi_T, rates)
        pi_ss = steady_state(Q)

        # πQ = 0 (row vector convention)
        residual = transpose(pi_ss) * Q
        for j in 1:4
            @test abs(residual[j]) < 1e-12
        end
    end

    @testset "Default K=1" begin
        rates_default = rates_closed_form(0.1, 0.5, -1.0, -0.5, -2.0)
        rates_explicit = rates_closed_form(0.1, 0.5, -1.0, -0.5, -2.0, 1.0)
        @test rates_default == rates_explicit
    end

end
