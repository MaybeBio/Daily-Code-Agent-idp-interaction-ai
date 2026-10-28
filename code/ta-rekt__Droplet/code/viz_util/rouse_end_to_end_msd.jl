# Is the end-to-end coordinate of the Rouse chain actually diffusing with
# D_eff = 2 kBT/gamma, the constant the 1D backward-FPE reduction assumes?
#
# Runs a FREE chain at the production parameters of spec/rouse_mfpt_attractive.jl
# (no wall, no pair potential -- the question is about the chain's own dynamics) and
# writes the mean-square displacement of the tracked-bead distance on log-spaced lags
# to a JSON next to the figure, so viz_util/rouse_end_to_end_msd.py can re-plot it
# against exact discrete-Rouse theory without re-simulating.
#
#   julia --project=@. viz_util/rouse_end_to_end_msd.jl [dt/tau_bond] [nsteps] [seed]

using Droplet, JSON, Statistics

const S        = 300                    # bonds between the tracked beads
const NBEADS   = S + 1
const KUHN     = 20.0
const KBT      = 1.0
const FRICTION = 1.3477088948787062e-08
const STIFF    = 3 * KBT / KUHN^2
const TAU_BOND = FRICTION / STIFF
const TAU_ROUSE = TAU_BOND * S^2 / pi^2
const D_EFF    = 2 * KBT / FRICTION     # the reduction's constant
const OUT      = "media/figures/20260914/rouse_end_to_end_msd.json"

dtf    = parse(Float64, get(ARGS, 1, "0.01"))
nsteps = parse(Int,     get(ARGS, 2, "2000000"))
seed   = parse(Int,     get(ARGS, 3, "20260914"))
DT     = dtf * TAU_BOND

r = RousePolymer(; boltzmannConst = 1.0, temperature = KBT, friction = FRICTION,
                   kuhnLength = KUHN, nbeads = Int32(NBEADS), nsteps = Int32(nsteps),
                   dt = DT, trackedBead1 = 1, trackedBead2 = NBEADS,
                   saveMode = :distance_only, seed = seed)
init!(r); t0 = time(); simulate!(r)
@info "simulate!" seconds=round(time() - t0, digits=1) t_total=nsteps*DT tau_rouse=TAU_ROUSE

x    = r.distanceTraj
lags = unique(round.(Int, exp10.(range(0, log10(nsteps ÷ 20); length = 40))))
msd  = [mean((@view(x[1+L:end]) .- @view(x[1:end-L])) .^ 2) for L in lags]

mkpath(dirname(OUT))
open(OUT, "w") do io
    JSON.print(io, Dict("lags_steps" => lags, "t" => lags .* DT, "msd" => msd,
                        "dt" => DT, "D_eff" => D_EFF, "tau_bond" => TAU_BOND,
                        "tau_rouse" => TAU_ROUSE, "s" => S, "a" => KUHN,
                        "nsteps" => nsteps, "seed" => seed,
                        "mean_x" => mean(x), "rms_x" => sqrt(mean(x .^ 2))))
end
println("wrote $OUT")
