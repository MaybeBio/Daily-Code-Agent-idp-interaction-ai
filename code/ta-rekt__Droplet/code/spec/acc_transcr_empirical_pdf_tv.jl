## load circuit specifications

using Printf

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
epsilon = -6.4
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
d0.logrates = log10.(rates)

slist = [d0]

checkpoint_path = string("data/checkpoints/", splitext(basename(@__FILE__))[1])

# --- time-varying parameters ---
nruns_per_trace = 100        # Gillespie runs per MC trace
Tc = 1.0                     # concentration update interval (CTMC time units per MC frame)
snapshot_interval = 100.0   # snapshot every dt continuous time units
nbins = 100
nsplit = 12

# mc_data_dir = "/home/tohme/.julia/dev/Droplet/data/20260503"
# mc_data_dir = "/home/tohme/.julia/dev/Droplet/data/20260423"
mc_data_dir = "/home/tohme/.julia/dev/Droplet/data/20260624"

# varied parameters
alphas = collect(range(0.003, stop=0.014, length=12))
# alphas = [0.003]
betas = [0]
nparams = length(alphas)

blurb = """
Time-varying concentration Gillespie. MC polymer traces drive local TF concentration.
For each global concentration alpha, load 16 MC traces, run $(nruns_per_trace) Gillespie
trajectories per trace via generate_tv, snapshot accumulated ON-time, compute LLR.
"""

# --- moving-average smoother ---
function _moving_average(v::Vector{Float64}, w::Int)
    n = length(v)
    out = similar(v)
    half = w ÷ 2
    cumsum_v = cumsum(v)
    for i in 1:n
        lo = max(1, i - half)
        hi = min(n, i + half)
        out[i] = (cumsum_v[hi] - (lo > 1 ? cumsum_v[lo-1] : 0.0)) / (hi - lo + 1)
    end
    return out
end


const _MA_WINDOW = 1000   # moving-average window for MC concentration traces

# --- load MC concentration traces for a given alpha ---
function _load_mc_traces(alpha)
    folder = joinpath(mc_data_dir, @sprintf("proximal_promoter_transc_noen_long_d_30_c_%.3f", alpha))
    traces_raw = Vector{Float64}[]
    for f in sort(readdir(folder))
        endswith(f, ".json") || continue
        data = JSON.parsefile(joinpath(folder, f))
        for ka in keys(data)
            lc = Float64.(data[ka]["1"]["local_concentration"])
            push!(traces_raw, _moving_average(lc, _MA_WINDOW))
            # push!(traces, lc)
        end
    end

    traces = Vector{Float64}[]
    for trace in traces_raw
        trace_len = length(trace)
        for i in range(1, stop=nsplit)
            st=1+(i-1)*floor(Int, trace_len/nsplit)
            fn=i*floor(Int, trace_len/nsplit)
            push!(traces, collect(trace[st:fn]))
        end
    end

    return traces
end

function target!(res, args...)
    t_start = time_ns()

    slist, keya, alpha, beta, blurb = args
    d0 = slist[1]

    idist = ones(4) ./ 4
    on_ind_vec = Vector{Bool}(on_ind)
    N = d0.nstates

    # determine alpha_next for LLR
    alpha_next = (keya < length(alphas)) ? alphas[keya + 1] : alpha
    dc = alpha_next - alpha

    # load MC traces
    traces = _load_mc_traces(alpha)
    traces_next = _load_mc_traces(alpha_next)
    ntraces = length(traces)
    min_len = minimum(length.(vcat(traces, traces_next)))

    # traces, traces_next = Vector(), Vector()
    # for k in range(1, ntraces)
    #     phi_local_avg = sum(traces_raw[k])/length(traces_raw[k])
    #     phi_local_next_avg = sum(traces_raw_next[k])/length(traces_raw_next[k])

    #     push!(traces, ones(min_len) .* phi_local_avg)
    #     push!(traces_next, ones(min_len) .* phi_local_next_avg)
    # end

    # use shortest trace to define t_max
    t_max_local = min_len * Tc
    t_grid = collect(snapshot_interval:snapshot_interval:t_max_local)
    nsnapshots = length(t_grid)
    nruns_total = ntraces * nruns_per_trace

    # --- helper: accumulate ON-time from (s, t) onto snapshot grid ---
    function _acc_on_from_trace!(row, s, t)
        cumulative = 0.0
        wall_time = 0.0
        snap_idx = 1
        for i in eachindex(s)
            snap_idx > nsnapshots && break
            dt = t[i]
            in_on = on_ind_vec[s[i]]
            next_wall = wall_time + dt
            while snap_idx <= nsnapshots && t_grid[snap_idx] <= next_wall
                row[snap_idx] = in_on ? cumulative + (t_grid[snap_idx] - wall_time) : cumulative
                snap_idx += 1
            end
            in_on && (cumulative += dt)
            wall_time = next_wall
        end
        while snap_idx <= nsnapshots
            row[snap_idx] = cumulative
            snap_idx += 1
        end
    end

    # --- helper: generate trajectories from MC traces ---
    function _generate_acc_on_tv!(acc_on, mc_traces)
        for (ti, trace) in enumerate(mc_traces)
            c_series = trace[1:min_len]
            for ri in 1:nruns_per_trace
                k = (ti - 1) * nruns_per_trace + ri

                # build fresh Q at c_series[1]
                weighted_edges = [(e[1], e[2], rates[i]) for (i, e) in enumerate(d0.edges)]
                Q = build_mat(weighted_edges, d0.input_edges, max(c_series[1], 1e-30))

                s, t = generate_tv(Q, d0.input_edges, N, c_series, Tc, idist)
                _acc_on_from_trace!(@view(acc_on[k, :]), s, t)
            end
        end
    end

    # --- helper: build histogram ---
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

    function _log_density(x, lo, bw, counts, ntotal)
        bw == 0.0 && return 0.0
        b = clamp(floor(Int, (x - lo) / bw) + 1, 1, nbins)
        c = counts[b]
        c > 0 ? log(c / (ntotal * bw)) : log(1.0 / (2 * ntotal * bw))
    end

    # --- 1. generate acc_on under alpha ---
    acc_on = zeros(nruns_total, nsnapshots)
    _generate_acc_on_tv!(acc_on, traces)

    # --- 2. generate acc_on under alpha_next ---
    acc_on_next = zeros(nruns_total, nsnapshots)
    _generate_acc_on_tv!(acc_on_next, traces_next)

    # --- 3. compute LLR traces ---
    llr = zeros(nruns_total, nsnapshots)

    for j in 1:nsnapshots
        col      = @view acc_on[:, j]
        col_next = @view acc_on_next[:, j]

        lo_a, bw_a, counts_a = _build_hist(col)
        lo_b, bw_b, counts_b = _build_hist(col_next)

        (bw_a == 0.0 && bw_b == 0.0) && continue

        bin_indices_a = zeros(Int, nruns_total)
        if bw_a > 0.0
            inv_bw_a = 1.0 / bw_a
            for k in 1:nruns_total
                bin_indices_a[k] = clamp(floor(Int, (col[k] - lo_a) * inv_bw_a) + 1, 1, nbins)
            end
        end

        for k in 1:nruns_total
            if bw_a > 0.0
                c_loo = counts_a[bin_indices_a[k]] - 1
                log_p_alpha = (c_loo > 0 ? log(c_loo / ((nruns_total - 1) * bw_a)) 
                                         : log(1.0 / (2 * nruns_total * bw_a)))
            else
                log_p_alpha = 0.0
            end
            log_p_next = _log_density(col[k], lo_b, bw_b, counts_b, nruns_total)
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
    res["nruns_per_trace"] = nruns_per_trace
    res["ntraces"] = ntraces
    res["nruns_total"] = nruns_total
    res["T"] = Tc
    res["t_max"] = t_max_local
    res["snapshot_interval"] = snapshot_interval
    res["t_grid"] = t_grid
    res["nbins"] = nbins
    res["llr_drift"] = sum(llr[:,end]) / (nruns_total * t_grid[end])

    # --- variability of the final LLR drift, split by source (law of total variance) ---
    # rows are grouped by MC trace: trace ti occupies rows
    #   (ti-1)*nruns_per_trace+1 : ti*nruns_per_trace
    # the runs within a trace share one concentration path, so they differ only by
    # transcription noise; different traces are independent concentration realizations.
    endpoint = llr[:, end]
    trace_means = [mean(@view endpoint[(ti-1)*nruns_per_trace+1 : ti*nruns_per_trace]) for ti in 1:ntraces]
    trace_vars  = [var(@view endpoint[(ti-1)*nruns_per_trace+1 : ti*nruns_per_trace]) for ti in 1:ntraces]

    var_within  = mean(trace_vars)          # transcription noise, conditional on the path
    var_between = ntraces > 1 ? var(trace_means) : 0.0  # concentration fluctuations across traces
    var_total   = var_within + var_between  # full spread of one LLR trace (≈ var(endpoint))

    inv_T = 1.0 / t_grid[end]
    res["llr_drift_stdev"]      = sqrt(var_total)   * inv_T   # total variability (error bar)
    res["llr_drift_std_transc"] = sqrt(var_within)  * inv_T   # from transcription noise
    res["llr_drift_std_conc"]   = sqrt(var_between) * inv_T   # from concentration fluctuations
    # standard error of the *mean* drift, clustered on the ntraces independent paths:
    res["llr_drift_sem"] = sqrt(var_between / ntraces + var_within / nruns_total) * inv_T

    # res["llr"] = Dict(k => llr[k, :] for k in 1:nruns_total)
    # res["acc_on"] = Dict(k => acc_on[k, :] for k in 1:nruns_total)

    t_current = (time_ns() - t_start) / 6e10
    res["runtime_min"] = t_current

    return res
end
