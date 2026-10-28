# Fisher information rate vs the Rouse bead count PARAMS.rouse.s, at FIXED
# concentration.  Thin driver -- no new physics: it just reruns the existing two
# stages (bvp_potential_interp.jl -> fisher_rate_from_fpt_solver.jl) for each s in
# `svals`, harvests the fisher rate (total / esc / ret) at every phi_inf in the
# (small) concentration grid, and writes them all to one JSON for plotting.
#
# Everything except s and the concentration grid is taken from param.jl as-is (the
# potential is PARAMS.fpt.potential, etc).  Only rouse.s and sweep.phi_inf_values are
# varied, by rewriting the `const PARAMS = FisherParams(...)` line each iteration;
# param.jl is restored on exit.  Both stages are include-d into throwaway modules so
# Droplet compiles ONCE and the 50 iterations reuse it (far faster than 100 cold
# `julia` starts).
#
# The fisher stage already differentiates wrt phi_inf via its 3-point phi stencil, so
# reading the rate AT a concentration c needs only c plus its next two uniformly-spaced
# neighbours: the stage emits nphi-2 = 1 record, labelled with the FIRST of the three
# and weighted by P(.|c) and c^2 (the score is a central difference over the outer pair,
# so it carries a half-step offset -- that is the convention the full sweeps use too).
# Spacing must be uniform.
#
# Under :bvp those three values cannot be typed in freely: bvp_potential_interp
# reconciles phi_inf_grid(PARAMS) against the phi_inf the BVP sweep actually computed,
# and it cannot interpolate BVP rows in phi_inf.  So the triple is PICKED FROM the
# sweep's own grid, around `phi_target` below.  That is what lets a big pre-existing
# sweep sit in output_dir and be reused, instead of having to move the output_cpu*.json
# out of the way to stop the reconciliation from firing.
#
# This copy drives the ARCLENGTH (d-major) job in this directory: the BVP sweep here
# parallelises over the GAP, so each record carries a whole phi_inf axis rather than one
# phi_inf.  The interp stage in this directory reads both layouts; `bracket_phis` below
# likewise takes the phi_inf axis from either shape.
#
#   julia --project=@. spec/jobs/20260807/R_5_h0_-0.4_chi_3.0_bcr_0.04/sweep_rouse_s.jl

using Droplet
import JSON

const HERE       = @__DIR__
const PARAM_FILE = joinpath(HERE, "param.jl")
const INTERP     = joinpath(HERE, "bvp_potential_interp.jl")
const FISHER     = joinpath(HERE, "fisher_rate_from_fpt_solver.jl")

svals = collect(range(5, 1000; length = 64))   # 50 bead counts, 10 .. 300

orig = read(PARAM_FILE, String)                 # restore this verbatim at the end
include(PARAM_FILE)                             # Main's PARAMS -> the (s-independent) output paths
fisher_out = joinpath(PARAMS.sim.output_dir, PARAMS.sim.fisher_file)
phi_target = (PARAMS.sweep.phi_frac_max - PARAMS.sweep.phi_frac_min)/2   # reported AT

# The three phi_inf to run: the sweep node nearest `target` plus its two neighbours,
# taken from the BVP sweep's OWN grid so bvp_potential_interp can always match them.
# Only phi_inf is read out of the shards here (the stage itself does the real
# aggregation); :linear reads no BVP files at all, so there param.jl still decides.
function bracket_phis(target)
    PARAMS.fpt.potential == :bvp || return PARAMS.sweep.phi_inf_values
    files = filter(f -> occursin(r"output_cpu\d+\.json$", f),
                   readdir(PARAMS.sim.output_dir; join = true))
    isempty(files) && error("no output_cpu*.json under $(PARAMS.sim.output_dir)")
    phis = Float64[]
    for f in files, (_, recs) in JSON.parsefile(f), (_, rec) in recs
        v = rec["phi_inf"]
        v === nothing && continue
        # d-major (arclength, this job): one record per GAP carrying the whole phi_inf
        # axis -- every record shares it, which the interp stage verifies.
        # phi-major: one record per phi_inf.  Either way this collects the axis.
        if v isa AbstractVector
            append!(phis, Float64[x for x in v if x !== nothing])
        else
            push!(phis, Float64(v))
        end
    end
    phis = sort(unique(phis))
    length(phis) >= 3 || error("BVP sweep has $(length(phis)) phi_inf value(s); need >= 3")
    # Anchor the triple's FIRST node on phi_target.  The fisher stage labels record k
    # with phi_infs[k] and weights it by P(.|phi_infs[k]) and phi_infs[k]^2 (only the
    # score itself is the central difference over the outer pair), so the rate belongs
    # to phi_1, NOT to the middle value -- same convention as the full sweeps, which
    # keeps this driver's point directly comparable to those curves.
    i = clamp(argmin(abs.(phis .- target)), 1, length(phis) - 2)
    return phis[i:i+2]
end

phi_run = bracket_phis(phi_target)
println("phi_target=$phi_target -> running phi_inf = $phi_run  (rate reported at $(phi_run[1]))")

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
        # inject this s and the bracketing phi triple, leaving every other default in
        # param.jl untouched (the rewrite DROPS any kwargs that line already had, so
        # anything else you want to change belongs in the kwdef defaults, not here)
        patched = replace(orig,
            r"const PARAMS = FisherParams\([^\n]*\)" =>
            "const PARAMS = FisherParams(rouse = RouseParams(s = $s), " *
            "sweep = SweepParams(phi_inf_values = $phi_run))")
        write(PARAM_FILE, patched)

        # rerun both stages, each in a fresh module (fresh const PARAMS per stage)
        run_stage(INTERP)
        run_stage(FISHER)

        per = JSON.parsefile(fisher_out)["per_phi"]     # esc/ret live here, not top level
        push!(results, Dict(
            "s"               => s,
            "phi_infs"        => [r["phi_inf"]         for r in per],
            "fisher"          => [r["fisher"]          for r in per],
            "fisher_esc"      => [r["fisher_esc"]      for r in per],
            "fisher_ret"      => [r["fisher_ret"]      for r in per],
            "escape_mfpt"     => [r["escape_mfpt"]     for r in per],
            "return_mfpt"     => [r["return_mfpt"]     for r in per],
            "phi_local"       => [r["phi_local"]       for r in per],
            "dphi_local"      => [r["dphi_local"]      for r in per]))
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
