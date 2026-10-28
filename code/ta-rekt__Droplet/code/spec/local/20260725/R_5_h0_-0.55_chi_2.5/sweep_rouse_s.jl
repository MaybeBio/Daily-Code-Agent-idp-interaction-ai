# Fisher information rate vs the Rouse bead count PARAMS.rouse.s, at FIXED
# concentration.  Thin driver -- no new physics: it just reruns the existing two
# stages (bvp_potential_interp.jl -> fisher_rate_from_fpt_solver.jl) for each s in
# `svals`, harvests the fisher rate (total / esc / ret) at every phi_inf in the
# (small) concentration grid, and writes them all to one JSON for plotting.
#
# Everything except s is taken from param.jl as-is: the concentration grid is
# PARAMS.sweep.phi_inf_values (edit that ONCE below in param.jl -- see note), the
# potential is PARAMS.fpt.potential, etc.  Only rouse.s is varied, by rewriting the
# `const PARAMS = FisherParams(...)` line each iteration; param.jl is restored on
# exit.  Both stages are include-d into throwaway modules so Droplet compiles ONCE
# and the 50 iterations reuse it (far faster than 100 cold `julia` starts).
#
# The fisher stage already differentiates wrt phi_inf via its 3-point phi stencil,
# so to read the rate AT a concentration c you only need c bracketed by one
# uniformly-spaced neighbour on each side in phi_inf_values (the middle fisher_rate
# entry is "at c").  Spacing must be uniform.
#
#   julia --project=@. spec/jobs/20260726/sweep_rouse_s.jl

using Droplet
import JSON

const HERE       = @__DIR__
const PARAM_FILE = joinpath(HERE, "param.jl")
const INTERP     = joinpath(HERE, "bvp_potential_interp.jl")
const FISHER     = joinpath(HERE, "fisher_rate_from_fpt_solver.jl")

svals = collect(range(10, 300; length = 64))   # 50 bead counts, 10 .. 300

orig = read(PARAM_FILE, String)                 # restore this verbatim at the end
include(PARAM_FILE)                             # Main's PARAMS -> the (s-independent) output paths
fisher_out = joinpath(PARAMS.sim.output_dir, PARAMS.sim.fisher_file)

# Run one stage in a throwaway module.  Anonymous Module()s don't get a single-arg
# `include`, so define one (-> Base.include into this module) for the stage's own
# `include("param.jl")` to resolve; each stage thus rebuilds a fresh const PARAMS.
function run_stage(path)
    m = Module(:Stage)
    Core.eval(m, :(include(p) = Base.include($m, p)))
    Base.include(m, path)
end

results = Vector{Dict{String,Any}}()
try
    for (i, s) in enumerate(svals)
        # inject this s, leaving every other default in param.jl untouched
        patched = replace(orig,
            r"const PARAMS = FisherParams\([^\n]*\)" =>
            "const PARAMS = FisherParams(rouse = RouseParams(s = $s))")
        write(PARAM_FILE, patched)

        # rerun both stages, each in a fresh module (fresh const PARAMS per stage)
        run_stage(INTERP)
        run_stage(FISHER)

        per = JSON.parsefile(fisher_out)["per_phi"]     # esc/ret live here, not top level
        push!(results, Dict(
            "s"               => s,
            "phi_infs"        => [r["phi_inf"]         for r in per],
            "fisher_rate"     => [r["fisher_rate"]     for r in per],
            "fisher_rate_esc" => [r["fisher_rate_esc"] for r in per],
            "fisher_rate_ret" => [r["fisher_rate_ret"] for r in per]))
        println(">>> [$i/$(length(svals))] s=$(round(s; digits = 2))  nphi_out=$(length(per))")
    end

    outpath = joinpath(PARAMS.sim.output_dir, "fisher_vs_rouse_s.json")
    open(outpath, "w") do io
        JSON.print(io, Dict("svals" => svals, "runs" => results))
    end
    println("saved -> ", outpath)
finally
    write(PARAM_FILE, orig)                      # always restore param.jl
end
