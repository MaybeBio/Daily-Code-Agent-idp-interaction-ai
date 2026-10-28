# . - -- .- .-.- . .- .-- -.- .- .- --.-.... - -.-
# TODO: Finish this shit
export generate
function generate(n, Q, idist)

    s, t = zeros(Int, n), zeros(n)
    
    s[1] = sample(eachindex(idist), Weights(idist))
    q = -Q[s[1], s[1]]

    _DD = _ddist(Q)

    for i in 1:n-1
        si = s[i]
        q = -Q[si, si]

        # t
        e = Exponential(1/q)
        t[i] = rand(e)

        # s        
        _dd = _DD[si, :]
        s[i+1] = sample(eachindex(_dd), Weights(_dd))
    end

    # t
    si = s[end]
    q = -Q[si, si]
    e = Exponential(1/q)
    t[end] = rand(e)

    return s, t
end

export generate_tv
function generate_tv(Q, input_edges, N, c_series, T, idist)
    n_intervals = length(c_series)
    t_max = n_intervals * T

    # base rates (at c=1) for input edges — allows direct assignment, no ratio issues
    c0 = c_series[1]
    base_rates = [Q[n1, n2] / max(c0, 1e-10) for (n1, n2) in input_edges]

    function _update_Q!(c)
        for (j, (n1, n2)) in enumerate(input_edges)
            Q[n1, n2] = base_rates[j] * c
        end
        for ri in 1:N
            Q[ri, ri] = 0.0
            Q[ri, ri] = -sum(@view Q[ri, :])
        end
    end

    DD = _ddist(Q)
    s = Int[]; t = Float64[]
    sizehint!(s, 10000); sizehint!(t, 10000)

    si = sample(eachindex(idist), Weights(idist))
    push!(s, si)
    wall = 0.0       # absolute time
    idx = 1          # current interval index
    sojourn = 0.0    # accumulated sojourn in current state across boundary splits

    while true
        q = -Q[si, si]
        dt = rand(Exponential(1/q))

        boundary = idx * T
        if wall + dt >= boundary && idx < n_intervals
            sojourn += boundary - wall
            wall = boundary
            idx += 1
            _update_Q!(c_series[idx])
            DD = _ddist(Q)
            # memoryless: redraw from new rate, same state
            continue
        end

        # transition happens
        wall += dt
        sojourn += dt

        if wall >= t_max
            push!(t, sojourn - (wall - t_max))  # clamp to t_max
            break
        end

        push!(t, sojourn)
        sojourn = 0.0

        dd = @view DD[si, :]
        si = sample(eachindex(idist), Weights(dd))
        push!(s, si)
    end

    return s, t
end

export generate!
function generate!(m::LLROptModel)
    # TODO: inline
    s, t = generate(m.n, m.Q_c, m.idist)
    
    m.s = s
    m.t = t
    return m
end