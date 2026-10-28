# Validation figures for the 1D backward-FPE first-passage-time solver
# (src/fpt_solver.jl).  For a few representative potentials it overlays the
# solver's FPT density against a histogram of first-passage times from a direct
# Euler–Maruyama Brownian-dynamics simulation of the *same* overdamped Langevin
# equation, confirming the two agree.
#
#   julia --project=@. viz_util/fpt_solver_validation.jl
#
# Writes PNG + PDF to media/figures/<today>/.

ENV["GKSwstype"] = "100"   # headless GR backend

using Droplet
using Plots
using Random
using Statistics
using Printf
using Dates

# --- Euler–Maruyama first-passage sampler in the same 1D potential -----------
# Overdamped: dx = -D β U'(x) dt + sqrt(2 D dt) ξ.  Reflecting at xmin,
# absorbing at xmax; initial position drawn Boltzmann on [xmin, xmax].
function bd_first_passage(Ufun, Uprime, xmin, xmax; D, kBT, dt, nsamples, seed)
    rng = MersenneTwister(seed)
    β = 1.0 / kBT
    c1 = D * β * dt
    c2 = sqrt(2 * D * dt)
    umin = Ufun(xmin)   # potential increases away from the reflecting wall here
    times = Float64[]
    maxstep = 5_000_000
    for _ in 1:nsamples
        x0 = xmin
        while true
            xp = xmin + (xmax - xmin) * rand(rng)
            rand(rng) < exp(-β * (Ufun(xp) - umin)) && (x0 = xp; break)
        end
        xc = x0
        for step in 1:maxstep
            xc += -c1 * Uprime(xc) + c2 * randn(rng)
            xc <= xmin && (xc = 2 * xmin - xc)        # reflect
            if xc >= xmax                              # absorb
                push!(times, step * dt); break
            end
        end
    end
    return times
end

# --- Cases (mirroring test/test_fpt_solver.jl) -------------------------------
cases = [
    (name = "Free diffusion",  latex = "U(x) = 0",
     xmin = 0.0, xmax = 2.0, N = 401,
     Ufun = x -> 0.0,            Uprime = x -> 0.0),
    (name = "Constant force",  latex = "U(x) = 1.5\\,x",
     xmin = 0.0, xmax = 2.0, N = 801,
     Ufun = x -> 1.5 * x,        Uprime = x -> 1.5),
    (name = "Double well",     latex = "U(x) = 2(x^2-1)^2",
     xmin = -1.5, xmax = 1.5, N = 1001,
     Ufun = x -> 2.0 * (x^2 - 1)^2, Uprime = x -> 8.0 * x * (x^2 - 1)),
]

D = 1.0; kBT = 1.0; dt = 1e-3; nsamples = 20_000

pot_panels = Plots.Plot[]
fpt_panels = Plots.Plot[]

for c in cases
    x = collect(range(c.xmin, c.xmax; length = c.N))
    U = c.Ufun.(x)
    prob = FPTProblem(x, U; D = D, kBT = kBT,
                      left_bc = :reflecting, right_bc = :absorbing)
    res = solve_fpt(prob)

    times = bd_first_passage(c.Ufun, c.Uprime, c.xmin, c.xmax;
                             D = D, kBT = kBT, dt = dt,
                             nsamples = nsamples, seed = 20260701)
    mfpt_bd = mean(times)
    @printf("%-16s  MFPT solver = %.4f   MFPT Langevin = %.4f   (%.1f%%)\n",
            c.name, res.mfpt, mfpt_bd, 100 * (mfpt_bd - res.mfpt) / res.mfpt)

    # potential panel
    pp = plot(x, U; lw = 2, color = :black, legend = false,
              xlabel = "x", ylabel = "U(x) / k_BT",
              title = "$(c.name)\n\$$(c.latex)\$")
    vline!(pp, [c.xmin]; color = :dodgerblue, ls = :dash, lw = 2)  # reflecting
    vline!(pp, [c.xmax]; color = :crimson,    ls = :solid, lw = 2) # absorbing
    push!(pot_panels, pp)

    # FPT panel: Langevin histogram + solver density
    tmax = quantile(times, 0.99)
    edges = range(0, tmax; length = 60)
    pf = histogram(times; bins = edges, normalize = :pdf,
                   color = :dodgerblue, alpha = 0.45, linecolor = :white,
                   label = "Langevin (n=$nsamples)",
                   xlabel = "first-passage time t", ylabel = "PDF  f(t)",
                   title = @sprintf("MFPT: solver %.3f / BD %.3f", res.mfpt, mfpt_bd))
    keep = res.t .<= tmax
    plot!(pf, res.t[keep], res.pdf[keep]; lw = 3, color = :crimson,
          label = "backward-FPE solver")
    vline!(pf, [res.mfpt]; color = :crimson, ls = :dot, lw = 2, label = "MFPT (solver)")
    xlims!(pf, 0, tmax)
    push!(fpt_panels, pf)
end

fig = plot(pot_panels..., fpt_panels...;
           layout = (2, 3), size = (1500, 850), left_margin = 6Plots.mm,
           bottom_margin = 6Plots.mm, dpi = 150,
           plot_title = "Backward-FPE FPT solver vs Langevin Brownian dynamics")

outdir = joinpath("media", "figures", Dates.format(today(), "yyyymmdd"))
mkpath(outdir)
png = joinpath(outdir, "fpt_solver_validation.png")
pdf = joinpath(outdir, "fpt_solver_validation.pdf")
savefig(fig, png); savefig(fig, pdf)
println("saved -> ", png)
println("saved -> ", pdf)
