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
# The fisher stage already differentiates wrt phi_inf via its 3-point phi stencil,
# so to read the rate AT a concentration c you only need c bracketed by one
# uniformly-spaced neighbour on each side in phi_inf_values (the middle -- and here
# only -- fisher_rate entry is "at c").  Spacing must be uniform.
#
# Under :bvp those three values cannot be typed in freely: bvp_potential_interp
# reconciles phi_inf_grid(PARAMS) against the phi_inf the BVP sweep actually computed,
# and it cannot interpolate BVP rows in phi_inf.  So the triple is PICKED FROM the
# sweep's own grid, around `phi_target` below.  That is what lets a big pre-existing
# sweep sit in output_dir and be reused, instead of having to move the output_cpu*.json
# out of the way to stop the reconciliation from firing.
#
#   julia --project=@. spec/jobs/20260726/sweep_rouse_s.jl

using Droplet
import JSON

const HERE       = @__DIR__
const PARAM_FILE = joinpath(HERE, "param.jl")
const INTERP     = joinpath(HERE, "bvp_potential_interp.jl")
const FISHER     = joinpath(HERE, "fisher_rate_from_fpt_solver.jl")

svals = collect(range(10, 300; length = 64))   # 50 bead counts, 10 .. 300
phi_target = 0.012                              # concentration the rate is reported AT

orig = read(PARAM_FILE, String)                 # restore this verbatim at the end
include(PARAM_FILE)                             # Main's PARAMS -> the (s-independent) output paths
fisher_out = joinpath(PARAMS.sim.output_dir, PARAMS.sim.fisher_file)

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
        # One scalar phi_inf per record is the layout this driver's interp stage reads.
        # Arclength sweeps (spec/jobs/20260807/...) parallelise over the GAP and store a
        # whole phi_inf branch per record; that data needs its own interp stage, not this
        # one, so say so here instead of dying in a Float64 conversion.
        v isa AbstractVector &&
            error("$f stores phi_inf as a vector (arclength layout: parallel axis = gap).\n" *
                  "bvp_potential_interp.jl here expects one scalar phi_inf per record -- " *
                  "point PARAMS.sim.output_dir at a sweep whose parallel axis IS phi_inf.")
        push!(phis, Float64(v))
    end
    phis = sort(unique(phis))
    length(phis) >= 3 || error("BVP sweep has $(length(phis)) phi_inf value(s); need >= 3")
    # Anchor the triple's FIRST node on phi_target: the fisher stage labels record k with
    # phi_infs[k] and weights it by P(.|phi_infs[k]) and phi_infs[k]^2, so the rate
    # belongs to phi_1, not to the middle value (same convention as the full sweeps).
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
