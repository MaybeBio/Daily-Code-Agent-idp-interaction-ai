import Pkg
using JSON
using Statistics
using Droplet
using Base.Filesystem
using LinearAlgebra

BLAS.set_num_threads(1)

let
    ### command-line arguments
    specfile = ARGS[1]
    cpukey = parse(Int32, ARGS[2])
    ncpus = parse(Int32, ARGS[3])
    outdir = ARGS[4]

    ### output directory
    Base.Filesystem.mkpath(outdir)

    ### specification file
    specpath = string("../", specfile)
    include(specpath)

    ### stamp run metadata into the spec's PARAMS (if it defines one), so every
    ### record's `params` dict records where this data came from and how it ran.
    if isdefined(Main, :PARAMS) && hasproperty(Main.PARAMS, :sim)
        Main.PARAMS.sim.spec_path  = specfile
        Main.PARAMS.sim.output_dir = outdir
        Main.PARAMS.sim.ncpus      = Int(ncpus)
    end

    nparams_per_cpu = div(nparams, ncpus)

    ### output data
    results = Dict{Int64, Dict{Int64, Dict{String, Any}}}()

    for k in 1:nparams_per_cpu

        keya = (cpukey - 1) * nparams_per_cpu + k
        alpha = alphas[keya]
        results[keya] = Dict{Int64, Dict{String, Any}}()

        for (keyb, beta) in enumerate(betas)

            res = Dict{String, Any}()
            args = [slist, keya, alpha, beta, blurb]

            results[keya][keyb] = target!(res, args...)
        end

        ### save checkpoint
        checkpoint_at_param = string(checkpoint_path, "_keya", keya)
        try
            checkpoint!(checkpoint_at_param, slist)
        catch e
            # @warn "checkpoint! skipped for this slist type: $(typeof(slist))"
        end

    end
    
    ### write data
    file = string(outdir, "/output_cpu$(cpukey).json")

    open(file, "w") do f
        JSON.print(f, results, 4)
    end
end
