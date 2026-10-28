### Homogeneous-bulk Flory-Huggins binodal validation.
###
### A periodic, single-species lattice gas (no surface) at fixed volume fraction
### phi=0.5, sweeping the FH coupling alpha across the mean-field critical point
### alpha_c = 2. Designed to be small and fast so the binodal can be measured
### without waiting for hours.
###
### Mapping (matches viz.ipynb FH formula F = mu*phi + 3*chi*phi(1-phi) + T*S):
###     alpha = 3 * chi_notebook
###     intPot!(bulk, -alpha/3)
###     bulk.chemPot = -mu - alpha
### Validated by test/test_flory_huggins_mapping.jl.

checkpoint_path = string("data/checkpoints/", splitext(basename(@__FILE__))[1])

l = Lattice(;
    boxlength = 24,
    nspecies = 1
)

bulk = FastBulkMonomer(;
    length = 1
)

slist = [bulk]

# minimum cluster size to count as a "droplet"
minsize = 10

### simulation parameters

blurb = """
Homogeneous bulk binodal: no surface, periodic BCs, phi=0.5, alpha sweep above
the mean-field critical point. Measures dilute and dense phase densities for
direct comparison against the Flory-Huggins binodal.
"""

t_max = 1                  # runtime in minutes — start small per CLAUDE.md
t_burn_in = t_max / 2      # discard first half as burn-in
step = 50                  # sample every <step> accepted moves

alphas = [2.5, 3.0, 3.5, 4.0]   # above critical alpha=2 for clear two-phase signal
betas = [0.5]                    # symmetric volume fraction
mu = 0.1
nparams = length(alphas)

### functions

function target!(res, args...)
    t_start = time_ns()
    t_current = 0.0

    slist, keya, alpha, beta, blurb = args

    phi_bar = beta
    chi = alpha     # spec convention: alpha = 3*chi_notebook

    bulk = slist[1]
    bulk.length = Int32(round(l.boxlength^3 * phi_bar))
    bulk.chemPot = -mu - chi
    intPot!(bulk, -chi/3)

    init!(l, slist)

    nmoves = 0
    iter = 0
    frame = 1

    ndroplets = Int[]
    largest_droplet_size = Int[]
    largest_droplet_density = Float64[]
    dilute_density = Float64[]
    dense_density_samples = Float64[]
    dilute_density_samples = Float64[]
    hams = Float64[]
    accs = Float64[]

    droplets = filter(x -> length(x) >= minsize, cc(bulk))
    push!(ndroplets, length(droplets))
    push!(largest_droplet_size, 1)
    push!(largest_droplet_density, 0.0)
    push!(dilute_density, diluteDensity(bulk, droplets))
    push!(hams, bulk.ham)
    push!(accs, 0.0)

    while t_current <= t_max
        t_current = (time_ns() - t_start) / 6e10
        iter += 1

        if move!(bulk)
            nmoves += 1

            if mod(nmoves, step) == 0
                frame += 1

                ccs = cc(bulk)
                droplets = filter(x -> length(x) >= minsize, ccs)

                if isempty(droplets)
                    maxsize = 1
                    dens = 0.0
                else
                    sizes = map(length, droplets)
                    maxsize, idx = findmax(sizes)
                    dmax = droplets[idx]
                    dens = maxsize / dropletVolume(bulk, dmax)
                end

                dilute = diluteDensity(bulk, droplets)

                push!(ndroplets, length(droplets))
                push!(largest_droplet_size, maxsize)
                push!(largest_droplet_density, dens)
                push!(dilute_density, dilute)
                push!(hams, hams[end] + bulk.deltaH)
                push!(accs, iter == 0 ? 0.0 : nmoves / iter)

                bulk.deltaH = 0

                if t_current >= t_burn_in
                    push!(dilute_density_samples, dilute)
                    push!(dense_density_samples, dens)
                end
            end
        end
    end

    coords = Matrix{Int32}(undef, bulk.length, 4)

    res["blurb"] = blurb
    res["species_name"] = bulk.name
    res["nframes"] = frame
    res["step"] = step
    res["concentration"] = alpha
    res["temperature"] = beta
    res["boxlength"] = bulk.boxlength
    res["length"] = bulk.length
    res["mu"] = mu
    res["t_burn_in"] = t_burn_in
    res["ndroplets"] = ndroplets
    res["largest_droplet_size"] = largest_droplet_size
    res["largest_droplet_density"] = largest_droplet_density
    res["dilute_density"] = dilute_density
    res["dilute_density_mean"] = isempty(dilute_density_samples) ? NaN : sum(dilute_density_samples) / length(dilute_density_samples)
    res["dense_density_mean"] = isempty(dense_density_samples) ? NaN : sum(dense_density_samples) / length(dense_density_samples)
    res["nsamples_post_burn"] = length(dilute_density_samples)
    res["hamiltonian"] = hams
    res["acceptance_rate"] = accs
    res["final_coordinates"] = extractCoords!(bulk, coords)
    res["runtime_min"] = (time_ns() - t_start) / 6e10

    return res
end
