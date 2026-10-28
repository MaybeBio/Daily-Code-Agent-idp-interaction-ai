
export init!, simulate!

# Rouse connectivity matrix A: tridiagonal with diag (1,2,2,...,2,1), off-diag -1.
# Stored as SymTridiagonal so mul!(drift, A, pos) is O(N): the dense N×N product
# dominated the step (24x slower at N = 301).
function _build_adjacency_matrix(N::Integer)
    dv = fill(2.0, N)
    dv[1] = 1.0
    dv[N] = 1.0
    return SymTridiagonal(dv, fill(-1.0, N - 1))
end

# Characteristic scales: timeScale = γ/κ, posScale = √(2 kBT/κ); κ = 3 kBT/b².
function _set_scaling_params!(r::RousePolymer)
    r.stiffness = 3 * r.boltzmannConst * r.temperature / r.kuhnLength^2
    r.timeScale = r.friction / r.stiffness
    r.posScale  = sqrt(2 * r.boltzmannConst * r.temperature / r.stiffness)
    r.dtDimless = r.dt / r.timeScale
    return r
end

# Allocate working buffers and trajectory storage.
function _set_simulation_vars!(r::RousePolymer)
    N = Int(r.nbeads)
    r.adjacency = _build_adjacency_matrix(N)
    r.pos       = zeros(Float64, N, 3)
    r.drift     = zeros(Float64, N, 3)
    r.noise     = zeros(Float64, N, 3)
    # :sojourn_only skips the (potentially huge) distance trace entirely.
    r.distanceTraj = r.saveMode === :sojourn_only ? nothing : zeros(Float64, r.nsteps)
    if r.saveMode === :full
        r.posTraj = zeros(Float64, r.nsteps, N, 3)
    else
        r.posTraj = nothing
    end
    # Sojourn buffers, allocated only when both bounds are set.
    if r.a !== nothing && r.b !== nothing
        @assert r.a <= r.b "sojourn bounds require a <= b (got a=$(r.a), b=$(r.b))"
        r.A_intervals = Int[]
        r.B_intervals = Int[]
    end
    return r
end

# Overwrite the N×3 noise buffer in-place from rng.
function _set_noise!(r::RousePolymer, rng::AbstractRNG)
    randn!(rng, r.noise)
    return r.noise
end

# Sample an equilibrium Gaussian chain: bond vectors ~ N(0, b²/3) per axis.
function _sample_equilibrium_chain(rng::AbstractRNG, N::Integer, b::Float64)
    pos = zeros(Float64, N, 3)
    sigma = b / sqrt(3)
    @inbounds for i in 2:N
        pos[i, 1] = pos[i - 1, 1] + sigma * randn(rng)
        pos[i, 2] = pos[i - 1, 2] + sigma * randn(rng)
        pos[i, 3] = pos[i - 1, 3] + sigma * randn(rng)
    end
    return pos
end

# Initialize scales, buffers, initial positions, and the first sample of distanceTraj.
function init!(r::RousePolymer)
    if r.saveMode === :sojourn_only && (r.a === nothing || r.b === nothing)
        error("saveMode :sojourn_only requires both bounds a and b to be set")
    end
    _set_scaling_params!(r)
    _set_simulation_vars!(r)

    rng = r.seed === nothing ? Random.default_rng() : Random.MersenneTwister(r.seed)

    if r.initPos === nothing
        r.initPos = _sample_equilibrium_chain(rng, Int(r.nbeads), r.kuhnLength)
    end
    r.pos .= r.initPos ./ r.posScale

    b1, b2 = Int(r.trackedBead1), Int(r.trackedBead2)
    dx = (r.pos[b1, 1] - r.pos[b2, 1]) * r.posScale
    dy = (r.pos[b1, 2] - r.pos[b2, 2]) * r.posScale
    dz = (r.pos[b1, 3] - r.pos[b2, 3]) * r.posScale
    if r.distanceTraj !== nothing
        r.distanceTraj[1] = sqrt(dx * dx + dy * dy + dz * dz)
    end

    if r.saveMode === :full
        @views r.posTraj[1, :, :] .= r.pos .* r.posScale
    end
    return r
end

# Dimensionless Euler-Maruyama:
#   x ← x − (A·x)·dτ + ξ·√dτ,   ξ ~ N(0,1) per bead per dim.
# The tracked-bead distance is recorded each step (unless saveMode == :sojourn_only);
# positions only when saveMode == :full. When bounds a, b are set, A/B sojourn
# durations are accumulated regardless of saveMode.
function simulate!(r::RousePolymer)
    rng = r.seed === nothing ? Random.default_rng() : Random.MersenneTwister(r.seed)
    # consume the same draws that init! used so the seed produces reproducible streams
    # (init! constructs its own rng from the same seed -- safe to re-seed here)

    sdt = sqrt(r.dtDimless)
    ps  = r.posScale
    b1  = Int(r.trackedBead1)
    b2  = Int(r.trackedBead2)
    fullMode = r.saveMode === :full
    f   = r.attractiveForce
    invKappaPos = 1 / (r.stiffness * ps)
    rMin = r.r_min === nothing ? nothing : r.r_min / ps   # hard repulsion in dimensionless units

    dtraj = r.distanceTraj                  # nothing in :sojourn_only mode

    # Sojourn tracking state (Schmitt-trigger). Initial classification matches
    # the scalar reference: in_A iff the initial distance ≤ a, else in_B (so a
    # starting position inside the buffer is assigned to B). The first sojourn of
    # each region is taken to enter at step 0.
    track = r.a !== nothing && r.b !== nothing
    if track
        aa = r.a::Float64
        bb = r.b::Float64
        Aint = r.A_intervals::Vector{Int}
        Bint = r.B_intervals::Vector{Int}
        empty!(Aint); empty!(Bint)
        d0x = (r.pos[b1, 1] - r.pos[b2, 1]) * ps
        d0y = (r.pos[b1, 2] - r.pos[b2, 2]) * ps
        d0z = (r.pos[b1, 3] - r.pos[b2, 3]) * ps
        d0  = sqrt(d0x*d0x + d0y*d0y + d0z*d0z)
        in_A = d0 <= aa
        in_B = !in_A
        A_entry = 0
        B_entry = 0
    end

    @inbounds for t in 2:r.nsteps
        mul!(r.drift, r.adjacency, r.pos)

        # Add the radial pair force between the two tracked beads to drift.
        # Convention: f(d_phys) = −dU/dr; the dimensionless force on i1 is
        #   F/(κ·posScale) = f(d_phys) · (x_i1 − x_i2) / (κ·posScale·|Δx|).
        # We *subtract* it from drift so the existing update
        #   pos ← pos − drift·dτ + noise·√dτ
        # adds the force contribution with the correct sign.
        if f !== nothing
            dxv  = r.pos[b1, 1] - r.pos[b2, 1]
            dyv  = r.pos[b1, 2] - r.pos[b2, 2]
            dzv  = r.pos[b1, 3] - r.pos[b2, 3]
            dmag = sqrt(dxv*dxv + dyv*dyv + dzv*dzv)
            if dmag > 0
                coeff = -invKappaPos * f(dmag * ps) / dmag
                r.drift[b1, 1] += coeff * dxv
                r.drift[b1, 2] += coeff * dyv
                r.drift[b1, 3] += coeff * dzv
                r.drift[b2, 1] -= coeff * dxv
                r.drift[b2, 2] -= coeff * dyv
                r.drift[b2, 3] -= coeff * dzv
            end
        end

        randn!(rng, r.noise)
        @. r.pos = r.pos - r.drift * r.dtDimless + r.noise * sdt

        # Hard repulsion: reflect the tracked-bead separation off r_min.
        # Mirrors r_curr ← 2 r_min − r_curr by rescaling the separation vector
        # about the pair midpoint so the new distance is 2 r_min − d.
        if rMin !== nothing
            sx = r.pos[b1, 1] - r.pos[b2, 1]
            sy = r.pos[b1, 2] - r.pos[b2, 2]
            sz = r.pos[b1, 3] - r.pos[b2, 3]
            smag = sqrt(sx*sx + sy*sy + sz*sz)
            if smag < rMin && smag > 0
                scale = (2 * rMin - smag) / smag
                mx = 0.5 * (r.pos[b1, 1] + r.pos[b2, 1])
                my = 0.5 * (r.pos[b1, 2] + r.pos[b2, 2])
                mz = 0.5 * (r.pos[b1, 3] + r.pos[b2, 3])
                r.pos[b1, 1] = mx + 0.5 * scale * sx
                r.pos[b1, 2] = my + 0.5 * scale * sy
                r.pos[b1, 3] = mz + 0.5 * scale * sz
                r.pos[b2, 1] = mx - 0.5 * scale * sx
                r.pos[b2, 2] = my - 0.5 * scale * sy
                r.pos[b2, 3] = mz - 0.5 * scale * sz
            end
        end

        dx = (r.pos[b1, 1] - r.pos[b2, 1]) * ps
        dy = (r.pos[b1, 2] - r.pos[b2, 2]) * ps
        dz = (r.pos[b1, 3] - r.pos[b2, 3]) * ps
        dist = sqrt(dx * dx + dy * dy + dz * dz)
        if dtraj !== nothing
            dtraj[t] = dist
        end

        # Sojourn update. i is the scalar-reference step index: distanceTraj[1]
        # is step 0, so the distance computed at loop index t is step i = t-1.
        # A→B is evaluated before B→A, matching the two-`if` order of the reference.
        if track
            i = t - 1
            if dist >= bb && in_A
                push!(Aint, i - A_entry)
                B_entry = i
                in_A = false
                in_B = true
            end
            if dist <= aa && in_B
                push!(Bint, i - B_entry)
                A_entry = i
                in_B = false
                in_A = true
            end
        end

        if fullMode
            @views r.posTraj[t, :, :] .= r.pos .* ps
        end
    end
    return r
end
