let 
    using Droplet
    using StatsBase
    using Distributions
    using Plots
    using Profile
    using OrderedCollections
    using Base.Filesystem
end

# ## 1010011111110000100000101010111100010100001100101110000
# # Testing Hamiltonian - PolyMonomer implementation
# @time let
#     date = "20251208"
#     outdir = "data/$(date)/"
#     Base.Filesystemkpath(outdir)
#     outfile = "test"
#     outpath = string(outdir, outfile)

#     l = Lattice(;
#         boxlength = 4,
#         nspecies = 1
#     )

#     p = Polymer(;
#         length = 4,
#         maxlength = 8,
#         nMonomerTypes = 1,
#         initPos = [1,2,2],
#         pattern=:uniform
#     )
    
#     initMonomers!(p)
#     a = p.monomers[1]

#     # self and interaction potentials
#     intPot!(a, -1.0)

#     slist = [p]

#     # initialize lattice and species
#     init!(l, slist)

#     #### iterating multiple moves
#     nmoves = [0]
#     frame = [1]
#     frame_joint = 1
#     nframes = 10000
#     step = 1
#     deltaH = [0.0]
#     niter = 0
#     moveMade = [false]

#     hamsA = Vector{Float64}(undef, 0)
#     hams_ = [hamsA]
#     hams = []
#     accDeltaHA = Vector{Float64}(undef, 0)
#     accDeltaH_ = [accDeltaHA]
#     accDeltaH = []
    
#     # first frame, initial energy
#     push!(hams_[1], getHam(a))
#     push!(accDeltaH_[1], getHam(a))

#     bulkSpecies = Vector{Any}(undef, 0)

#     for (sp, idx) in a.intSpecies
#         push!(hams, hams_[idx])
#         push!(accDeltaH, accDeltaH_[idx])

#         if typeof(sp) != PolyMonomer
#             push!(bulkSpecies, (sp, idx))
#         end
#     end

#     ham!(a)
#     p.deltaH = 0.0

#     color = ["blue (A)"]
#     exportSpecies!(slist, string(outpath, "_frame$(1)"))

#     movetypes = ["add", "kink", "remove"]
#     bad_frame = 1

#     while frame_joint <= nframes - 2
#         for (s, i) in a.intSpecies

#             s_prev = deepcopy(p)
#             slist_prev = [s_prev]

#             if move!(s)
#                 nmoves[i] += 1
#                 deltaH[i] += getDeltaHam(s)

#                 push!(accDeltaH[i], accDeltaH[i][end] + deltaH[i])
#                 deltaH[i] = 0
                
#                 ham!(s)
#                 push!(hams[i], getHam(s))

#                 frame[i] += 1
#                 frame_joint += 1

#                 if hams[1] != accDeltaH[1]
#                     bad_frame = frame_joint
#                     println("Inconsistency at frame: ", frame_joint)
#                     println("Move type: ", movetypes[p.moveType])
#                     println("Bead: ", p.moveIdx)

#                     exportSpecies!(slist_prev, string(outpath, "_frame$(2)"))
#                     exportSpecies!(slist, string(outpath, "_frame$(3)"))
#                     frame_joint = nframes
#                 end
#             end
#         end
#     end

#     # plots = []
#     # colors = ["blue"]
#     # hamy = hams[1][bad_frame-1:bad_frame+9]
#     # accy = accDeltaH[1][bad_frame-1:bad_frame+9]

#     # pl = plot(eachindex(hamy), hamy;
#     # label = "ham",
#     # xlabel = "frame",
#     # ylabel = "Energy",
#     # color = colors[1],
#     # legend = :topright)

#     # plot!(pl, eachindex(accy), accy;
#     # label = "accDeltaH",
#     # color = colors[1],
#     # linestyle = :dash)

#     # push!(plots, pl)

#     # plot(plots[1]; layout = (1, 1), size = (600, 300))
# end

# ## 1010011111110000100000101010111100010100001100101110000
# # Testing Hamiltonian and its increment
# @time let
#     date = "20251208"
#     outdir = "data/$(date)/"
#     Base.Filesystemkpath(outdir)
#     outfile = "test"
#     outpath = string(outdir, outfile)

#     l = Lattice(;
#         boxlength = 16,
#         nspecies = 3
#     )

#     p = Polymer(;
#         length = 16,
#         maxlength = 64,
#         nMonomerTypes = 1,
#         initPos = [1,8,8],
#         pattern=:uniform
#     )
    
#     initMonomers!(p)
#     a = p.monomers[1]

#     b = FastBulkMonomer(;
#         length = 64,
#         initShape = "uniform"
#     )

#     c = FastBulkMonomer(;
#         length = 64,
#         initShape = "uniform"
#     )

#     # self and interaction potentials
#     intPot!(a, -1.0)
#     intPot!(b, -1.0)
#     intPot!(c, -1.0)
#     intPot!(a, b, -1.0)
#     intPot!(a, c, -1.0)

#     slist = [p, b, c]
#     # slist = [a, b]

#     # initialize lattice and species
#     init!(l, slist)

#     #### iterating multiple moves
#     nmoves = [0,0,0]
#     frame = [1,1,1]
#     frame_joint = 1
#     nframes = 1000
#     step = 1
#     deltaH = [0.0,0.0,0.0]
#     niter = 0
#     moveMade = [false, false, false]

#     hamsA = Vector{Float64}(undef, 0)
#     hamsB = Vector{Float64}(undef, 0)
#     hamsC = Vector{Float64}(undef, 0)
#     hams_ = [hamsA, hamsB, hamsC]
#     hams = []
#     accDeltaHA = Vector{Float64}(undef, 0)
#     accDeltaHB = Vector{Float64}(undef, 0)
#     accDeltaHC = Vector{Float64}(undef, 0)
#     accDeltaH_ = [accDeltaHA, accDeltaHB, accDeltaHC]
#     accDeltaH = []
    
#     # first frame, initial energy
#     push!(hams_[1], getHam(a))
#     push!(hams_[2], getHam(b))
#     push!(hams_[3], getHam(c))
#     push!(accDeltaH_[1], getHam(a))
#     push!(accDeltaH_[2], getHam(b))
#     push!(accDeltaH_[3], getHam(c))

#     bulkSpecies = Vector{Any}(undef, 0)

#     for (sp, idx) in a.intSpecies
#         push!(hams, hams_[idx])
#         push!(accDeltaH, accDeltaH_[idx])

#         if typeof(sp) != PolyMonomer
#             push!(bulkSpecies, (sp, idx))
#         end
#     end

#     ham!(a)
#     ham!(b)
#     ham!(c)
#     p.deltaH = 0.0
#     b.deltaH = 0.0
#     c.deltaH = 0.0
#     p.deltaHP = 0.0
#     b.deltaHP = 0.0
#     c.deltaHP = 0.0

#     color = ["blue (A)", "orange (B)", "green (C)"]
#     exportSpecies!(slist, string(outpath, "_frame$(1)"))

#     while frame_joint <= nframes - 2
#         for (s, i) in a.intSpecies
#             if move!(s)
#                 nmoves[i] += 1
#                 deltaH[i] += getDeltaHam(s)

#                 # update deltaH for the partners
#                 for (partner, j) in a.intSpecies
#                     if j != i
#                         deltaH[j] += getDeltaHP(partner)
#                         # @show p.name, deltaH[j]
#                     end
#                 end

#                 # @show s.name, deltaH[i] 
#                 push!(accDeltaH[i], accDeltaH[i][end] + deltaH[i])
#                 deltaH[i] = 0
                
#                 ham!(s)
#                 push!(hams[i], getHam(s))

#                 frame[i] += 1
#                 frame_joint += 1

#                 # println("frame ", frame_joint, ", ", color[i], " moved (", frame[i], ")")
#                 # exportSpecies!(slist, string(outpath, "_frame$(frame_joint)"))

#                 if hams[1] != accDeltaH[1] || hams[2] != accDeltaH[2] || hams[3] != accDeltaH[3]
#                     println("inconsistency at frame ", frame_joint)
#                     break
#                 end
#             end
#         end
#     end

#     # gyration radii of a straight line and a sphere of volume 'length'
#     # exportSpecies!(slist, string(outpath, "_frame$(2)"))

#     println("A acceptance rate: ", nmoves[1] / nframes)
#     println("B acceptance rate: ", nmoves[2] / nframes)

#     @assert hams[1] == accDeltaH[1]
#     @assert hams[2] == accDeltaH[2]
#     @assert hams[3] == accDeltaH[3]

#     # plots = []
#     # colors = ["blue", "orange", "green"]

#     # for i in eachindex(hams)
#     #     pl = plot(eachindex(hams[i]), hams[i][:];
#     #     label = "ham",
#     #     xlabel = "frame",
#     #     ylabel = "Energy",
#     #     color = colors[i],
#     #     legend = :topright)

#     #     plot!(pl, eachindex(accDeltaH[i]), accDeltaH[i][:];
#     #     label = "accDeltaH",
#     #     color = colors[i],
#     #     linestyle = :dash)

#     #     push!(plots, pl)
#     # end

#     # plot(plots[1], plots[2], plots[3]; layout = (3, 1), size = (600, 600))
# end

# @time let
#     date = "20251208"
#     outdir = "data/$(date)/"
#     Base.Filesystemkpath(outdir)
#     outfile = "test"
#     outpath = string(outdir, outfile)

#     ############ droplet boundary ############

#     l = Lattice(;
#         boxlength = 20,
#         nspecies = 3
#     )

#     tf = FastBulkMonomer(;
#         length = 32,
#         initShape = "uniform"
#     )

#     # rnap = FastBulkMonomer(;
#     #     length = 32,
#     #     initShape = "uniform"
#     # )

#     ### rod 1
#     dna1 = FastBulkMonomer(;
#         length = 7,
#         initShape = "line",
#         initParams = 1,
#         initPos = [1,5,8]
#     )

#     enhancer = FastBulkMonomer(;
#         length = 6,
#         initShape = "line",
#         initParams = 1,
#         initPos = [8,5,8]
#     )

#     dna2 = FastBulkMonomer(;
#         length = 7,
#         initShape = "line",
#         initParams = 1,
#         initPos = [14,5,8]
#     )


#     ### rod 1
#     dna1_l = FastBulkMonomer(;
#         length = 7,
#         initShape = "line",
#         initParams = 1,
#         initPos = [1,10,8]
#     )

#     enhancer_l = FastBulkMonomer(;
#         length = 6,
#         initShape = "line",
#         initParams = 1,
#         initPos = [8,10,8]
#     )

#     dna2_l = FastBulkMonomer(;
#         length = 7,
#         initShape = "line",
#         initParams = 1,
#         initPos = [14,10,8]
#     )

#     # # self-interactions
#     # intPot!(tf, -1.0)
#     # intPot!(rnap, -1.0)

#     # # multi-species interactions
#     # intPot!(tf, dna1, -0.5)  
#     # intPot!(tf, dna2, -0.5)  
#     # intPot!(tf, enhancer, -2.0)

#     # intPot!(rnap, dna1, -0.5)  
#     # intPot!(rnap, dna2, -0.5)  
#     # intPot!(rnap, enhancer, -0.5)

#     # intPot!(tf, rnap, -2.0)

#     slist = [dna1, dna2, enhancer, dna1_l, dna2_l, enhancer_l, tf]
    
#     # initialize lattice and species
#     init!(l, slist)

#     exportSpecies!(slist, string(outpath, "_frame$(1)"))
# end

# ## 1010011111110000100000101010111100010100001100101110000
# # Multi-component polymer

# @time let

#     date = "20251208"
#     outdir = "data/$(date)/"
#     Base.Filesystemkpath(outdir)
#     outfile = "test"
#     outpath = string(outdir, outfile)

#     ############ droplet boundary ############

#     l = Lattice(;
#         boxlength = 42,
#         nspecies = 3
#     )

#     # heterogeneous polymer with patterned regions
#     poly = Polymer(;
#         length = 32,
#         maxlength = 256,
#         initShape = "line",
#         nMonomerTypes = 2,  # base, pattern1
#         initPos = [5, 12, 12],
#         pattern = :custom
#     )

#     # define pattern: 256 beads with two patterned regions
#     pattern = zeros(Int, 32)
#     p1_start, p1_end = 12, 16
#     pattern[1:p1_start-1] .= 1         # base type (40 beads)
#     pattern[p1_start:p1_end] .= 2      # pattern 1 (8 beads)
#     pattern[p1_end+1:32] .= 1

#     poly.customPattern = pattern
    
#     initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)

#     init!(l, [poly])

#     find_patterns(poly)

# end

## 1010011111110000100000101010111100010100001100101110000
# Generate sphere and test histogram

@time let
    date = "20251208"
    outdir = "data/$(date)/"
    Base.Filesystem.mkpath(outdir)
    outfile = "sphere_test"
    outpath = string(outdir, outfile)

    # create lattice
    l = Lattice(;
        boxlength = 40,
        nspecies = 1
    )

    # create sphere with specified radius
    sphere_radius = 5
    sphere_length = Int32(floor(4/3 * π * sphere_radius^3))
    normalHist = zeros(100)

    sphere = FastBulkMonomer(;
        length = 64,
        initShape = "sphere",
        initParams = 3,
        initPos = [20, 20, 20]
    )

    # initialize
    init!(l, [sphere])

    @show sphere.length

    # compute center of mass
    cm = centerMass(sphere)
    maxDistance = 10.0
    total = 0.0
    
    for _ in range(1, stop=10000)
        cm .+= rand(Normal(0.0, 0.005), 3)
        total += updateNormalHistogram!(sphere, normalHist, cm, maxDistance)
    end

    ## create bin edges and centers for plotting
    # nbins = length(normalHist)
    # bin_edges = range(0, maxDistance, length=nbins+1)
    # bin_centers = [(bin_edges[i] + bin_edges[i+1])/2 for i in 1:nbins]
    # normalHist ./= total
    # normalHist .*= sphere.length
    
    # # create histogram plot
    # pl = plot(
    #     bin_centers, normalHist,
    #     label="Radial density",
    #     xlabel="Distance from center (lattice units)",
    #     ylabel="Number of particles",
    #     title="Sphere radial density profile",
    #     legend=:topright,
    #     color=:blue,
    #     linewidth=2
    # )

    # display(pl)
    
    # also export to the standard format
    exportSpecies!([sphere], string(outpath, "_frame1"))
end

# ## 1010011111110000100000101010111100010100001100101110000
# # Testing interpolations

# let
#     using JSON
#     using Plots
#     using Interpolations

#     # Load the JSON data
#     json_path = "data/20260218/spherical_phi_vs_phi_inf_run_interface1.json"
#     data = JSON.parsefile(json_path)
    
#     # Extract the arrays
#     phi_inf_values = data["phi_inf_values"]
#     phi_at_probe = data["phi_at_probe"]
#     phi_func = LinearInterpolation(phi_inf_values, phi_at_probe, extrapolation_bc=Line())

#     # @show phi_inf_values
#     @show phi_func(0.237)

# end

# ## 1010011111110000100000101010111100010100001100101110000
# # Compute Fisher information rate from a concentration profile

# let
#     using JSON
#     using Plots
#     using Interpolations

#     # Load the JSON data
#     # json_path1 = "data/20260218/spherical_phi_vs_phi_inf_run_enhancer2.json"
#     # json_path2 = "data/20260218/spherical_phi_vs_phi_inf_run_enhancer2_collapsed.json"
#     # data1 = JSON.parsefile(json_path1)
#     # data2 = JSON.parsefile(json_path2)
    
#     # # Extract the arrays and convert to Float64
#     # phi_inf_values1 = Float64.(data1["phi_inf_values"])
#     # phi_inf_values2 = Float64.(data2["phi_inf_values"])
#     # phi_at_probe1 = Float64.(data1["phi_at_probe"])
#     # phi_at_probe2 = Float64.(data2["phi_at_probe"])
#     # phi_at_probe = similar(phi_at_probe1)

#     json_path1 = "data/20260304/empirical_phi_local_proximal_promoter_d_100.json"
#     data1 = JSON.parsefile(json_path1)
    
#     phi_inf_values1 = Float64.(data1["concentration"])
#     phi_at_probe1 = Float64.(data1["local_concentration_smoothed"])
#     phi_at_probe = phi_at_probe1

#     # for (idx, (phi1, phi2)) in enumerate(zip(phi_at_probe1, phi_at_probe2))
#     #     phi_infty = phi_inf_values1[idx]

#     #     if phi_infty < 0.2
#     #         phi_at_probe[idx] = phi1
#     #     else
#     #         phi_at_probe[idx] = phi2
#     #     end
#     # end

#     # Create monotonic spline interpolation with extrapolation
#     # Using Steffen monotonic interpolation to avoid oscillations at the transition
#     itp = interpolate(phi_inf_values1, phi_at_probe, SteffenMonotonicInterpolation())
#     phi_func = extrapolate(itp, Line())

#     fisher_at_probe = zeros(length(phi_at_probe))
#     fisher_at_probe_control = zeros(length(phi_at_probe))
#     activity_at_probe = zeros(length(phi_at_probe))
#     activity_at_probe_control = zeros(length(phi_at_probe))

#     ## load circuit specifications
#     circuits_file = "/home/tohme/.julia/dev/Droplet/spec/circuits.json"
#     circuit_specs = Dict()

#     open(circuits_file, "r") do f
#         circuit_specs = JSON.parse(f)
#     end

#     circuitId = "1bs"
#     circuit = circuit_specs[circuitId]
#     edges = circuit["edges"]
#     input_edges = circuit["input_edges"]
#     nstates = circuit["nstates"]
#     off_ind = circuit["off_indicator"]
#     on_ind = BitArray(1 .- off_ind)
#     start_state = circuit["start_state"]
#     hitting_set = circuit["hitting_set"]

#     reltol = 1e-10
#     eps = 1e-15

#     d0 = OneBurstDecoyModel(;
#         circuitId = circuitId,
#         nstates = nstates,
#         edges = edges,
#         input_edges = input_edges,
#         logrates = nothing,
#         cr = 1,
#         cw = nothing,
#         gamma = nothing,
#         ind_set = on_ind,
#         start_state = start_state,
#         hitting_set = hitting_set,
#         eps = eps,
#         reltol = reltol,
#         maxevals = nothing,
#     )

#     # kboff, kuoff, konb, koffb, kuon, kbon, koffu = 1, 1e-5, 1, 1e-5, 1, 1e-5, 1

#     eta_UA = -0.5
#     eta_BI = 0.0
#     epsilon = -7.75

#     K = 1
#     koffu = K
#     kuoff = K
#     koffb = K
#     kuon = K
#     konu = exp(-eta_UA)*K
#     kboff = exp(-eta_BI)*K
#     konb = exp(-eta_UA-epsilon)*K
#     kbon = exp(-eta_BI-epsilon)*K

#     rates = [kboff, kuoff, konb, koffb, kuon, kbon, koffu]
#     kmax = maximum(rates)

#     d0.logrates = log10.(rates)

#     d1 = deepcopy(d0)
#     d2 = deepcopy(d0)

#     for (idx, phi) in enumerate(phi_at_probe)    

#         d0.cr = phi
#         d0.cw = phi_inf_values1[idx]
        
#         # take a step in phi_infty
#         h = 0.02 * d0.cw
#         d1.cw = d0.cw - h
#         d2.cw = d0.cw - 2*h

#         d1.cr = phi_func(d1.cw)
#         d2.cr = phi_func(d2.cw)

#         @show d0.cr, d1.cr, h

#         initAcc!(d0)
#         initAcc!(d1)
#         initAcc!(d2)

#         activity_at_probe[idx] = pi_on(d0)
#         fisher_at_probe[idx] = FisherRateAcc(d0, d1, d2)
#     end

#     for (idx, phi) in enumerate(phi_inf_values1)    

#         d0.cr = phi
#         h = 0.02 * d0.cr
#         d1.cr = d0.cr - h
#         d2.cr = d0.cr - 2*h

#         d0.cw = d0.cr
#         d1.cw = d1.cr
#         d2.cw = d2.cr

#         initAcc!(d0)
#         initAcc!(d1)
#         initAcc!(d2)

#         activity_at_probe_control[idx] = pi_on(d0)
#         fisher_at_probe_control[idx] = FisherRateAcc(d0, d1, d2)
#     end
    
#     # Create the plot
#     pl1 = plot(
#         phi_inf_values1, phi_at_probe,
#         xlabel="\$\\phi_\\infty\$",
#         ylabel="\$\\phi(z)\$",
#         title="Probe concentration vs bulk concentration",
#         legend=false,
#         linewidth=1,
#         color=:blue
#     )

#     pl2 = plot(
#         phi_inf_values1, activity_at_probe,
#         xlabel="\$\\phi_\\infty\$",
#         ylabel="\$p_{ss}(\\phi(r_0))\$",
#         title="Steady-state activity vs bulk concentration",
#         legend=false,
#         linewidth=1,
#         color=:red
#     )

#     plot!(
#         pl2,
#         phi_inf_values1, activity_at_probe_control,
#         linewidth=1,
#         color=:red,
#         linestyle=:dash
#     )

#     pl3 = plot(
#         phi_inf_values1, fisher_at_probe,
#         xlabel="\$\\phi_\\infty\$",
#         ylabel="\$\\dot{F}(\\phi_\\infty)\$",
#         title="Fisher rate vs bulk concentration",
#         legend=false,
#         label="using \$\\phi(r_0)\$",
#         linewidth=1,
#         color=:green
#     )

#     plot!(
#         pl3,
#         phi_inf_values1, fisher_at_probe_control,
#         linewidth=1,
#         color=:green,
#         linestyle=:dash,
#         label="control",
#         legend=:bottomright,
#         yscale=:log10
#     )

#     hline!([0.16*kmax], label="\$0.16k_{max}\$")
    
#     display(plot(pl1, pl2, pl3; layout = (3, 1), size = (600, 800)))

#     savefig("media/figures/20260304/fisher_rate_local_proximal_promoter_d_100.png")
# end


# ## 1010011111110000100000101010111100010100001100101110000
# # Generate ON/OFF runs

# let
#     using JSON
#     using Plots

#     # json_path1 = "data/20260304/empirical_phi_local_proximal_promoter_d_100.json"
#     # data1 = JSON.parsefile(json_path1)

#     ## load circuit specifications
#     circuits_file = "/home/tohme/.julia/dev/Droplet/spec/circuits.json"
#     circuit_specs = Dict()

#     open(circuits_file, "r") do f
#         circuit_specs = JSON.parse(f)
#     end

#     circuitId = "1bs"
#     circuit = circuit_specs[circuitId]
#     edges = circuit["edges"]
#     input_edges = circuit["input_edges"]
#     nstates = circuit["nstates"]
#     off_ind = circuit["off_indicator"]
#     on_ind = BitArray(1 .- off_ind)
#     start_state = circuit["start_state"]
#     hitting_set = circuit["hitting_set"]

#     reltol = 1e-10
#     eps = 1e-15

#     d0 = OneBurstDecoyModel(;
#         circuitId = circuitId,
#         nstates = nstates,
#         edges = edges,
#         input_edges = input_edges,
#         logrates = nothing,
#         cr = 1,
#         cw = nothing,
#         gamma = nothing,
#         ind_set = on_ind,
#         start_state = start_state,
#         hitting_set = hitting_set,
#         eps = eps,
#         reltol = reltol,
#         maxevals = nothing,
#     )

#     # energies
#     eta_UA = -0.5
#     eta_BI = 0.0
#     epsilon = -7.75

#     # rnap concentration
#     phi_R = 0.02

#     # transition rates
#     K = 0.5
#     kuoff = K*exp(eta_BI)
#     koffb = K*exp(eta_UA+epsilon)
#     kuon = K*exp(eta_BI+epsilon)
#     koffu = K*exp(eta_UA)
#     konu = K*phi_R
#     kboff = K
#     konb = K*phi_R
#     kbon = K

#     rates = [kboff, kuoff, konb, koffb, kuon, kbon, koffu, konu]
#     kmax = maximum(rates)

#     d0.logrates = log10.(rates)

#     nruns = 1
#     nsteps = 100000
#     idist = ones(4) ./ 4

#     function _idist0(idist, set_ind)
#         idist0 = idist[findall(set_ind)]
#         idist0 ./= sum(idist0)
#         return idist0
#     end

#     # initialize data containers
#     state_trajs = zeros((nruns, nsteps))
#     time_trajs = zeros((nruns, nsteps))
#     state_traj = zeros(nsteps)
#     time_traj = zeros(nsteps)
#     on_times = Dict()
#     off_times = Dict()
#     llr = Dict()
#     res_times = Dict()

#     phi_global = range(0.003, stop=0.018, length=16)

#     _idist0_on = _idist0(idist, on_ind)
#     _idist0_off = _idist0(idist, off_ind)

#     # loop across concentrations
#     for (idx, phi) in enumerate(phi_global)    
#         d0.cr = phi

#         # equivalent to initFull without rate constraints
#         initFreeWeights!(d0)
#         initOneBurst!(d0)

#         # Gillespie simulation
#         state_traj, time_traj = generate(nsteps, d0.Q, idist)

#         # residence times
#         s_in_on, rt, onrt, offrt = get_residence_time(state_traj, time_traj, on_ind)

#         res_times[idx] = rt
#         on_times[idx] = onrt
#         off_times[idx] = offrt
#         # llr[idx] = llr_traj
#     end

#     # Generate complementary CDFs from on_times and off_times
#     phi_key = 1
    
#     on_times_sorted = sort(on_times[phi_key])
#     off_times_sorted = sort(off_times[phi_key])
    
#     ccdf_on = 1.0 .- collect(1:length(on_times_sorted)) ./ length(on_times_sorted)
#     ccdf_off = 1.0 .- collect(1:length(off_times_sorted)) ./ length(off_times_sorted)
    
#     # pl1 = plot(
#     #     on_times_sorted, log10.(ccdf_on),
#     #     xlabel="Duration (time units)",
#     #     ylabel="Complementary CDF",
#     #     title="ON state residence times",
#     #     legend=false,
#     #     color=:green,
#     #     linewidth=1,
#     #     # yscale=:log10,
#     #     xlim=(1, 100000),
#     #     xscale=:log10
#     # )

#     # plot!(pl1,
#     #     on_times_sorted[1:50000], -on_times_sorted[1:50000].*0.12,
#     #     xlabel="Duration (time units)",
#     #     ylabel="Complementary CDF",
#     #     title="ON state residence times",
#     #     legend=false,
#     #     color=:black,
#     #     linewidth=1,
#     #     # yscale=:log10,
#     #     xlim=(1, 100000),
#     #     xscale=:log10
#     # )
    
#     # pl2 = plot(
#     #     off_times_sorted, log10.(ccdf_off),
#     #     xlabel="Duration (time units)",
#     #     ylabel="Complementary CDF",
#     #     title="OFF state residence times",
#     #     legend=false,
#     #     color=:red,
#     #     linewidth=1,
#     #     # yscale=:log10,
#     #     xlim=(1, 10000),
#     #     xscale=:log10
#     # )
    
#     # display(plot(pl1, pl2; layout = (1, 2), size = (600, 300)))

#     # savefig("media/figures/20260415/on_off_times_test.png")
    
#     # Save llr trajectories to JSON
#     Base.Filesystem.mkpath("data/20260415/")
#     llr_data = Dict(
#         "phi_global" => collect(phi_global),
#         "llr" => Dict(string(k) => v for (k, v) in llr)
#     )
#     open("data/20260415/llr_trajs.json", "w") do f
#         JSON.print(f, llr_data, 4)
#     end
# end
