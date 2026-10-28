## load circuit specifications

circuits_file = "/home/tohme/.julia/dev/Droplet/spec/circuits.json"
circuit_specs = JSON.parsefile(circuits_file)

circuitId = "1bs"
circuit = circuit_specs[circuitId]
edges = circuit["edges"]
input_edges = circuit["input_edges"]
nstates = circuit["nstates"]
off_ind = circuit["off_indicator"]
on_ind = BitArray(1 .- off_ind)
start_state = circuit["start_state"]
hitting_set = circuit["hitting_set"]

reltol = 1e-10
eps = 1e-15

d0 = OneBurstDecoyModel(;
    circuitId = circuitId,
    nstates = nstates,
    edges = edges,
    input_edges = input_edges,
    logrates = nothing,
    cr = 1,
    cw = nothing,
    gamma = nothing,
    ind_set = on_ind,
    start_state = start_state,
    hitting_set = hitting_set,
    eps = eps,
    reltol = reltol,
    maxevals = nothing,
)

# fixed parameters
eta_UA = -0.5
eta_BI = 0.0
epsilon = -7.75
phi_R = 0.02

K = 0.5
kuoff = K*exp(eta_BI)
koffb = K*exp(eta_UA+epsilon)
kuon = K*exp(eta_BI+epsilon)
koffu = K*exp(eta_UA)
konu = K*phi_R
kboff = K
konb = K*phi_R
kbon = K

rates = [kboff, kuoff, konb, koffb, kuon, kbon, koffu, konu]
kmax = maximum(rates)

d0.logrates = log10.(rates)

slist = [d0]

checkpoint_path = string("data/checkpoints/", splitext(basename(@__FILE__))[1])

nruns = 1000
nsteps = 100000
t_max = 1000000.0             # total continuous time window
snapshot_interval = 1000.0   # snapshot every dt continuous time units
nbins = 50

# varied parameters
alphas = collect(range(0.003, stop=0.018, length=16))  # tf concentration
# alphas = [0.010]
betas = [0]
nparams = length(alphas)

blurb = """
Empirical PDF of accumulated transcription (ON-time) for a 4-state CTMC promoter model.
Sweep over TF concentration. For each concentration, run nruns Gillespie trajectories,
snapshot accumulated ON time at fixed continuous-time points separated by snapshot_interval,
build leave-one-out histograms, and store log p(m_t | c) traces.
"""

function target!(res, args...)
    t_start = time_ns()

    # species list and sweep parameters
    slist, keya, alpha, beta, blurb = args

    d0 = slist[1]

    idist = ones(4) ./ 4
    on_ind_vec = Vector{Bool}(on_ind)

    # fixed continuous-time snapshot grid
    t_grid = collect(snapshot_interval:snapshot_interval:t_max)
    nsnapshots = length(t_grid)

    # determine alpha and alpha_next
    alpha_next = (keya < length(alphas)) ? alphas[keya + 1] : alpha
    dc = alpha_next - alpha

    # --- helper: generate nruns trajectories and return acc_on matrix ---
    function _generate_acc_on!(acc_on, Q)
        for k in 1:nruns
            state_traj, time_traj = generate(nsteps, Q, idist)

            cumulative = 0.0
            wall_time = 0.0
            snap_idx = 1

            for i in 1:nsteps
                if snap_idx > nsnapshots
                    break
                end

                dt = time_traj[i]
                in_on = on_ind_vec[state_traj[i]]
                next_wall = wall_time + dt

                while snap_idx <= nsnapshots && t_grid[snap_idx] <= next_wall
                    if in_on
                        acc_on[k, snap_idx] = cumulative + (t_grid[snap_idx] - wall_time)
                    else
                        acc_on[k, snap_idx] = cumulative
                    end
                    snap_idx += 1
                end

                if in_on
                    cumulative += dt
                end
                wall_time = next_wall
            end

            while snap_idx <= nsnapshots
                acc_on[k, snap_idx] = cumulative
                snap_idx += 1
            end
        end
    end

    # --- helper: build histogram from a column and return (lo, bin_width, counts) ---
    function _build_hist(col)
        lo, hi = minimum(col), maximum(col)
        if hi - lo < 1e-30
            return (lo, 0.0, zeros(Int, nbins))
        end
        bw = (hi - lo) / nbins
        inv_bw = 1.0 / bw
        counts = zeros(Int, nbins)
        for k in eachindex(col)
            b = clamp(floor(Int, (col[k] - lo) * inv_bw) + 1, 1, nbins)
            counts[b] += 1
        end
        return (lo, bw, counts)
    end

    # --- helper: evaluate log-density of value x against a histogram (no leave-one-out) ---
    function _log_density(x, lo, bw, counts, ntotal)
        if bw == 0.0
            return 0.0
        end
        b = clamp(floor(Int, (x - lo) / bw) + 1, 1, nbins)
        c = counts[b]
        if c > 0
            return log(c / (ntotal * bw))
        else
            return log(1.0 / (2 * ntotal * bw))
        end
    end

    # --- 1. generate trajectories under alpha ---
    acc_on = zeros(nruns, nsnapshots)
    d0.cr = alpha
    initFreeWeights!(d0)
    initOneBurst!(d0)
    _generate_acc_on!(acc_on, d0.Q)

    _pi_on = pi_on(d0)
    _variance = variance(d0)
    _fisher_rate_acc = FisherRateAcc(d0)

    # --- 2. generate trajectories under alpha_next (for the alternative histogram) ---
    acc_on_next = zeros(nruns, nsnapshots)
    d0.cr = alpha_next
    initFreeWeights!(d0)
    initOneBurst!(d0)
    _generate_acc_on!(acc_on_next, d0.Q)

    # --- 3. compute LLR traces using leave-one-out for alpha, full hist for alpha_next ---
    # llr[k, j] = log p(m_{t_j}^{(k)} | c) - log p(m_{t_j}^{(k)} | c + dc)
    llr = zeros(nruns, nsnapshots)

    for j in 1:nsnapshots
        col      = @view acc_on[:, j]
        col_next = @view acc_on_next[:, j]

        # histogram under alpha (for leave-one-out evaluation of the traces)
        lo_a, bw_a, counts_a = _build_hist(col)

        # histogram under alpha_next
        lo_b, bw_b, counts_b = _build_hist(col_next)

        if bw_a == 0.0 && bw_b == 0.0
            continue
        end

        # precompute bin indices for leave-one-out under alpha
        bin_indices_a = zeros(Int, nruns)
        if bw_a > 0.0
            inv_bw_a = 1.0 / bw_a
            for k in 1:nruns
                bin_indices_a[k] = clamp(floor(Int, (col[k] - lo_a) * inv_bw_a) + 1, 1, nbins)
            end
        end

        for k in 1:nruns
            # log p(m | c) with leave-one-out
            if bw_a > 0.0
                c_loo = counts_a[bin_indices_a[k]] - 1
                if c_loo > 0
                    log_p_alpha = log(c_loo / ((nruns - 1) * bw_a))
                else
                    log_p_alpha = log(1.0 / (2 * nruns * bw_a))
                end
            else
                log_p_alpha = 0.0
            end

            # log p(m | c + dc) — full histogram (trace k is independent of this histogram)
            log_p_next = _log_density(col[k], lo_b, bw_b, counts_b, nruns)

            llr[k, j] = log_p_alpha - log_p_next
        end
    end

    # --- 4. store results ---
    res["blurb"] = blurb
    res["circuit_id"] = circuitId
    res["logrates"] = d0.logrates
    res["concentration"] = alpha
    res["concentration_next"] = alpha_next
    res["dc"] = dc
    res["nruns"] = nruns
    res["nsteps"] = nsteps
    res["t_max"] = t_max
    res["snapshot_interval"] = snapshot_interval
    res["t_grid"] = t_grid
    res["nbins"] = nbins
    res["p_on"] = _pi_on
    res["sigma_sq"] = _variance
    res["fisher_rate_acc"] = _fisher_rate_acc

    # store LLR traces: run index -> vector of snapshots
    res["llr"] = Dict(k => llr[k, :] for k in 1:nruns)

    # store accumulated transcription traces
    res["acc_on"] = Dict(k => acc_on[k, :] for k in 1:nruns)

    t_current = (time_ns() - t_start) / 6e10
    res["runtime_min"] = t_current

    return res
end