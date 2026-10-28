# Time-Varying Concentration Gillespie Generator

## Overview

Implement a Gillespie simulation driven by a time-varying concentration signal from Monte Carlo polymer simulations. The MC simulation produces a time series of local TF concentration at a promoter binding site, sampled every `step=10000` accepted moves. This concentration trace is fed into the CTMC promoter model, replacing the constant concentration used in `generate`.

**Physical setup:**
- 4-state CTMC promoter model (`1bs` circuit), with states: Unbound-OFF, Bound-OFF, Bound-ON, Unbound-ON
- `input_edges = [[1,2], [4,3]]` — TF binding transitions, scaled by concentration
- Rate constant `K` is chosen so 1 CTMC time unit = 1 MC frame
- Concentration is piecewise-constant: updated every `T` CTMC time units from the MC trace
- `T = 50` by default (smoothing window), so MC frame `k` covers CTMC time `[(k-1)*T, k*T)`

**Data source:** `data/20260415/proximal_promoter_transc_d_50_c_{alpha}/output_cpu{i}.json`
- 16 concentration folders (`alpha` ∈ 0.003:0.001:0.018)
- 16 CPU files per folder, each containing 1 MC trajectory
- Each trajectory has a `local_concentration` array of length ~70k frames
- Total: 16 traces per global concentration → 1600 Gillespie runs per concentration (100 per MC trace)

---

## Phase 1: `generate_tv` in `src/generate.jl`

### Description

Add a new exported function that runs a Gillespie simulation with piecewise-constant concentration, updating Q every `T` time units.

### Interface

```julia
export generate_tv
function generate_tv(Q, input_edges, N, c_series, T, idist)
    # Q         : initial rate matrix (will be mutated in-place via update_mat!)
    # input_edges : concentration-dependent edge list
    # N         : number of states (size(Q, 1))
    # c_series  : Vector{Float64}, concentration at each interval
    # T         : Float64, duration of each concentration interval
    # idist     : initial state distribution
    #
    # Returns: (s::Vector{Int}, t::Vector{Float64})
    #   variable-length state and sojourn-time vectors
    #   total simulation time ≈ length(c_series) * T
end
```

### Algorithm

```
t_max = length(c_series) * T
wall_time = 0.0
interval_idx = 1
c_prev = c_series[1]

# build initial Q (already done by caller), compute _ddist
DD = _ddist(Q)

s = Int[]; t = Float64[]
push!(s, sample(eachindex(idist), Weights(idist)))

while wall_time < t_max
    si = s[end]
    q = -Q[si, si]
    dt = rand(Exponential(1/q))

    # check if this step crosses an interval boundary
    next_boundary = interval_idx * T
    if wall_time + dt >= next_boundary && interval_idx < length(c_series)
        # clamp sojourn to boundary, record partial step
        dt_partial = next_boundary - wall_time
        push!(t, dt_partial)
        wall_time = next_boundary

        # advance to next interval
        interval_idx += 1
        c_new = c_series[interval_idx]
        update_mat!(Q, N, input_edges, c_prev, c_new)
        DD = _ddist(Q)
        c_prev = c_new

        # re-draw remaining sojourn from NEW rate (memoryless property)
        # no state transition occurred — stay in same state, draw fresh
        continue
    end

    # normal step within interval
    push!(t, dt)
    wall_time += dt

    if wall_time >= t_max
        break
    end

    # transition
    dd = @view DD[si, :]
    push!(s, sample(eachindex(idist), Weights(dd)))
end
```

**Key design decision — boundary handling:** When a sojourn crosses a concentration boundary, we split it: record the partial time up to the boundary, update Q, and redraw the remaining sojourn from the new exponential rate. This is exact because the exponential distribution is memoryless — conditioning on "no jump before boundary" and then restarting with the new rate is equivalent to a non-homogeneous Poisson process with piecewise-constant intensity.

**Key design decision — pre-allocation:** Use `sizehint!` on `s` and `t` with an estimate based on the mean total rate. This avoids excessive re-allocation while keeping the interface clean.

### Changes

**File:** `src/generate.jl` — append the new function after the existing `generate!`.

### Tests (`test/test_generate_tv.jl`)

1. **Constant concentration recovers `generate`:**
   ```julia
   # Build Q at concentration c=0.01, run generate(10000, Q, idist)
   # Build same Q, run generate_tv(Q, input_edges, N, fill(0.01, 1000), 10.0, idist)
   # Compare: mean sojourn time per state, steady-state occupancy fractions
   # Tolerance: 5% relative error (statistical)
   @test abs(mean_on_frac_tv - mean_on_frac_const) / mean_on_frac_const < 0.05
   ```

2. **Output length scales with trace length:**
   ```julia
   c_short = fill(0.01, 100)
   c_long  = fill(0.01, 1000)
   s1, t1 = generate_tv(Q1, input_edges, N, c_short, 50.0, idist)
   s2, t2 = generate_tv(Q2, input_edges, N, c_long,  50.0, idist)
   @test sum(t2) / sum(t1) ≈ 10.0 atol=1.0
   ```

3. **Total simulation time matches `length(c_series) * T`:**
   ```julia
   s, t = generate_tv(Q, input_edges, N, c_series, T, idist)
   @test sum(t) ≤ length(c_series) * T
   @test sum(t) > (length(c_series) - 1) * T  # at least got close
   ```

4. **Q is correctly updated — rate changes are reflected:**
   ```julia
   # c_series = [0.001, 0.001, ..., 0.1, 0.1, ...]  (low then high)
   # After the switch, transitions on input_edges should be ~100x faster
   # Measure mean sojourn in first vs second half
   # First half should have longer sojourns (lower binding rate)
   ```

5. **`s` and `t` have same length:**
   ```julia
   @test length(s) == length(t)
   ```

---

## Phase 2: Accumulated ON-time from variable-length traces

### Description

Factor the snapshot accumulation logic from `_generate_acc_on!` into a reusable function that works with variable-length `(s, t)` output from `generate_tv`. This stays in the spec file (not in `src/`), since it's experiment-specific.

### Interface (inline in spec)

```julia
function _acc_on_from_trace!(acc_on_row, s, t, on_ind_vec, t_grid)
    # acc_on_row : pre-allocated Vector{Float64} of length nsnapshots
    # s, t       : state and sojourn vectors from generate or generate_tv
    # on_ind_vec : Vector{Bool} indicating ON states
    # t_grid     : snapshot time points (sorted, increasing)
    #
    # Fills acc_on_row[j] = accumulated ON time at t_grid[j]
end
```

The logic is identical to what's currently inlined in `_generate_acc_on!`:

```julia
cumulative = 0.0
wall_time = 0.0
snap_idx = 1
nsnapshots = length(t_grid)

for i in eachindex(s)
    snap_idx > nsnapshots && break
    dt = t[i]
    in_on = on_ind_vec[s[i]]
    next_wall = wall_time + dt

    while snap_idx <= nsnapshots && t_grid[snap_idx] <= next_wall
        acc_on_row[snap_idx] = in_on ? cumulative + (t_grid[snap_idx] - wall_time) : cumulative
        snap_idx += 1
    end

    in_on && (cumulative += dt)
    wall_time = next_wall
end

# hold last value for remaining snapshots
while snap_idx <= nsnapshots
    acc_on_row[snap_idx] = cumulative
    snap_idx += 1
end
```

### Tests (inline in spec, print-based)

1. **Constant concentration via `generate_tv` matches `generate` pipeline:**
   ```julia
   # Run 200 trajectories each way, compare mean accumulated ON at last snapshot
   # Should agree within ~10% (statistical)
   ```

2. **Monotonicity:** `acc_on_row` is non-decreasing for every run.

3. **Final value bounded:** `acc_on_row[end] ≤ t_max` (can't spend more time ON than total time).

---

## Phase 3: Load MC concentration traces

### Description

Write a loading function (in the spec file) that reads all 16 CPU files for a given global concentration `alpha` and returns a flat vector of `local_concentration` arrays.

### Interface (inline in spec)

```julia
function _load_mc_traces(alpha)
    # alpha : Float64, global concentration (e.g., 0.010)
    # Returns: Vector{Vector{Float64}} — one concentration trace per MC run
    
    folder = @sprintf("data/20260415/proximal_promoter_transc_d_50_c_%.3f", alpha)
    traces = Vector{Float64}[]
    
    for cpufile in sort(readdir(folder, join=true))
        endswith(cpufile, ".json") || continue
        data = JSON.parsefile(cpufile)
        for keya in keys(data)
            lc = Float64.(data[keya]["1"]["local_concentration"])
            push!(traces, lc)
        end
    end
    return traces
end
```

Each CPU file contains 1 keya with 1 beta, so this yields exactly 16 traces per concentration.

### Tests (inline, print-based)

1. `length(_load_mc_traces(0.010)) == 16`
2. All traces have length > 50000
3. `mean(mean.(traces))` is positive and correlated with `alpha`

---

## Phase 4: New spec file `acc_transcr_empirical_pdf_tv.jl`

### Description

Create the sweep spec that combines Phases 1–3. For each `alpha`:
1. Load 16 MC traces
2. For each trace, run `nruns_per_trace = 100` Gillespie trajectories via `generate_tv`
3. Snapshot accumulated ON-time at regular intervals
4. Build LOO histograms and compute LLR vs next alpha

### Key parameters

```julia
nruns_per_trace = 100    # Gillespie runs per MC trace
T = 50.0                 # concentration update interval (CTMC time units per MC frame)
snapshot_interval = 1000.0
nbins = 50
```

### Structure

The spec file reuses the circuit setup from `acc_transcr_empirical_pdf.jl` (copy the preamble). The `target!` function:

```julia
function target!(res, args...)
    slist, keya, alpha, beta, blurb = args
    d0 = slist[1]
    
    # load MC traces for this alpha
    traces = _load_mc_traces(alpha)
    
    # determine alpha_next for LLR
    alpha_next = (keya < length(alphas)) ? alphas[keya + 1] : alpha
    
    # load MC traces for alpha_next
    traces_next = _load_mc_traces(alpha_next)
    
    # compute t_max from shortest trace (conservative)
    min_trace_len = minimum(length.(traces))
    t_max_local = min_trace_len * T
    t_grid = collect(snapshot_interval:snapshot_interval:t_max_local)
    nsnapshots = length(t_grid)
    
    ntraces = length(traces)
    nruns_total = ntraces * nruns_per_trace  # 16 * 100 = 1600
    
    # --- generate acc_on under alpha ---
    acc_on = zeros(nruns_total, nsnapshots)
    for (ti, trace) in enumerate(traces)
        # truncate trace to min_trace_len
        c_series = trace[1:min_trace_len]
        for ri in 1:nruns_per_trace
            k = (ti - 1) * nruns_per_trace + ri
            
            # build fresh Q for this trace's initial concentration
            d0.cr = c_series[1]
            initFreeWeights!(d0)
            initQ!(d0)
            Q = copy(d0.Q)
            
            s, t = generate_tv(Q, d0.input_edges, d0.nstates, c_series, T, idist)
            _acc_on_from_trace!(@view(acc_on[k, :]), s, t, on_ind_vec, t_grid)
        end
    end
    
    # --- generate acc_on_next under alpha_next ---
    # (same structure, using traces_next)
    
    # --- compute LLR (reuse _build_hist, _log_density from constant spec) ---
    # (identical to acc_transcr_empirical_pdf.jl)
end
```

**Important detail — fresh Q per run:** `generate_tv` mutates Q in-place via `update_mat!`. Each Gillespie run needs its own copy of the initial Q. We `copy(d0.Q)` before each call.

**Important detail — initial concentration for Q:** The initial Q should be built with `c_series[1]` (the first MC frame's local concentration), not with `alpha` (the global/bulk concentration). This ensures `update_mat!` ratios are correct from the start.

### Parallelization

The sweep runs 1 `alpha` per CPU (16 CPUs for 16 concentrations). Within each CPU, the 1600 trajectories run sequentially. This matches the existing `runJL` + `sweep.jl` workflow:

```bash
./runJL sim/sweep.jl spec/acc_transcr_empirical_pdf_tv.jl 16 data/20260417/acc_transcr_tv
```

### Output (`res` dict)

```julia
res["blurb"] = blurb
res["circuit_id"] = circuitId
res["logrates"] = d0.logrates
res["concentration"] = alpha
res["concentration_next"] = alpha_next
res["dc"] = alpha_next - alpha
res["nruns_per_trace"] = nruns_per_trace
res["ntraces"] = ntraces
res["nruns_total"] = nruns_total
res["T"] = T
res["t_max"] = t_max_local
res["snapshot_interval"] = snapshot_interval
res["t_grid"] = t_grid
res["nbins"] = nbins
res["llr"] = Dict(k => llr[k, :] for k in 1:nruns_total)
res["runtime_min"] = t_current
```

### Tests

1. **Smoke test with 1 alpha, 1 trace, 10 runs:**
   ```bash
   # Set alphas = [0.010], nruns_per_trace = 10
   ./runJL sim/sweep.jl spec/acc_transcr_empirical_pdf_tv.jl 1 data/20260417/acc_transcr_tv_test
   ```
   Verify output JSON has correct structure and `nruns_total = 160` (16 traces × 10).

2. **Sanity check — constant vs time-varying:** Run the constant-concentration spec and the TV spec side by side for the same `alpha`. The TV spec should produce *different* LLR drift values, confirming that concentration fluctuations matter.

---

## Phase 5: End-to-End Validation

### Steps

1. Run the TV spec on a single concentration with small parameters, verify JSON output.

2. Visualize LLR drift curve (notebook cell): load TV results, fit slope through origin, plot drift/dc vs concentration. Overlay with the constant-concentration curve from the earlier experiment.

3. Add `test_generate_tv.jl` to `runtests.jl`:
   ```julia
   include("test_generate_tv.jl")
   ```

4. Run the full test suite:
   ```bash
   cd /home/tohme/.julia/dev/Droplet && julia --project=@. test/runtests.jl
   ```

---

## Risk Assessment

| Risk | Likelihood | Impact | Mitigation |
|------|-----------|--------|------------|
| `update_mat!` division by zero when `c_prev = 0` | High | Critical | Guard: if `c_prev < eps`, rebuild Q from scratch via `build_mat` instead of ratio update |
| MC traces with `local_concentration = 0` at start | High | Critical | Same guard as above; initial Q uses `max(c_series[1], eps)` |
| Gillespie steps too fine relative to T | Low | Medium | `T = 50` ensures ~50 Gillespie steps per interval at typical rates — adequate |
| Memory from 1600 variable-length trajectories | Low | Low | Only `acc_on` matrix is stored (1600 × nsnapshots), trajectories are discarded |
| `generate_tv` slower than `generate` due to boundary checks | Low | Low | Inner loop adds one comparison per step; negligible vs exponential RNG |

---

## Dependency Graph

```
Phase 1 (generate_tv in src/generate.jl)
    │
    ├── Phase 2 (acc_on accumulation helper, inline in spec)
    │
    └── Phase 3 (MC trace loading, inline in spec)
            │
            └── Phase 4 (new spec file, combines 1+2+3)
                    │
                    └── Phase 5 (validation, test suite)
```

Phases 2 and 3 are independent of each other but both depend on Phase 1 being tested.

---

## Files Modified / Created

| File | Action |
|------|--------|
| `src/generate.jl` | Add `generate_tv` function |
| `test/test_generate_tv.jl` | New test file for Phase 1 |
| `spec/acc_transcr_empirical_pdf_tv.jl` | New spec file (Phase 4) |
| `test/runtests.jl` | Add `include("test_generate_tv.jl")` (Phase 5) |
