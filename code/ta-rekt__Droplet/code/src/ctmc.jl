## . - -- .- .-.- . .- .-- -.- .- .- --.-.... - -.-
# weighted_edges [(n1, n2, w), ...]
# input_edges [(n1, n2), ...]
# c
export build_mat
function build_mat(weighted_edges, input_edges, c::T) where T<:Any
    # max node index
    N = -1
    for (n1, n2, w) in weighted_edges
        N = max(N, n1, n2)
    end
    @assert N > 0
    
    # weights matrix
    Q = zeros(Float64, N, N)
    for (n1, n2, w) in weighted_edges
        Q[n1, n2] = w
    end

    # up input edges
    for (n1, n2) in input_edges
        Q[n1, n2] = Q[n1, n2] * c
    end

    # normalize
    for ri in 1:N
        Z = sum(view(Q, ri, :))
        Q[ri, ri] = -Z
    end

    return Q
end

export update_mat!
function update_mat!(Q, N, input_edges, c_prev, c_new)

    # up input edges
    for (n1, n2) in input_edges
        Q[n1, n2] = Q[n1, n2] * (c_new / c_prev)
    end

    # normalize
    for ri in 1:N
        Z = sum(view(Q, ri, :))
        Q[ri, ri] = -Z
    end

    return Q
end

export init_Qs!
function init_Qs!(m::LLROptModel)
    m.Q_c = build_mat(m.weighted_edges, m.input_edges, m.c)
    m.Q_cdc = build_mat(m.weighted_edges, m.input_edges, m.c + m.dc)
    return m
end

## . - -- .- .-.- . .- .-- -.- .- .- --.-.... - -.-
export P_ss
function P_ss(d::OneBurstDecoyModel)
    A = d.Q .* 1.0
    A[:, end] .= 1.0
    b = zeros(size(A, 1))
    b[end] = 1.0
    return clamp.(A' \ b, 1e-10, 1)
end

## . - -- .- .-.- . .- .-- -.- .- .- --.-.... - -.-
export EP_ss
function EP_ss(d::OneBurstDecoyModel)
    nstates = size(d.Q, 1)
    pss = P_ss(d)
    EP = 0.0

    for i in 1:nstates, j in 1:nstates
        if d.Q[j, i] != 0 && pss[j] > 0
            Jij = pss[i] * d.Q[i,j] - pss[j] * d.Q[j, i]
            EP += Jij * log(pss[i] * d.Q[i,j] / (pss[j] * d.Q[j, i]))
        end
    end
    return 0.5 * EP
end

## . - -- .- .-.- . .- .-- -.- .- .- --.-.... - -.-
export Z
function Z(d::OneBurstDecoyModel)
    pss = P_ss(d)
    Pi = d.Q .* 0.0
    for ri in 1:size(d.Q, 1)
        Pi[ri, :] .= pss
    end
    return inv(Pi - d.Q) - Pi
end 

## . - -- .- .-.- . .- .-- -.- .- .- --.-.... - -.-
export pi_on
# on_ind :: vector
function pi_on(d::OneBurstDecoyModel)
    return P_ss(d)' * d.ind_set
end

## . - -- .- .-.- . .- .-- -.- .- .- --.-.... - -.-
export variance
function variance(d::OneBurstDecoyModel)
    Z_ = Z(d)
    pss = P_ss(d)

    return 2 * (d.ind_set .* pss)' * (Z_ * d.ind_set)
end

## . - -- .- .-.- . .- .-- -.- .- .- --.-.... - -.-
export sharpness
function sharpness(d::OneBurstDecoyModel, h::Float64=0.01)

    d1 = deepcopy(d)
    d2 = deepcopy(d)
    d1.cr = d.cr - h
    d2.cr = d.cr - 2*h
    
    initAcc!(d1)
    initAcc!(d2)

    pi0 = pi_on(d)
    pi1 = pi_on(d1)
    pi2 = pi_on(d2)

    # using the three-point backward finite difference

    dpion = pi2 - 4*pi1 + 3*pi0
    dpion /= (2 * h)

    return abs(d.cr * dpion)
end

## . - -- .- .-.- . .- .-- -.- .- .- - -.-.... - -.-
export FisherRateAcc
function FisherRateAcc(d::OneBurstDecoyModel, h::Float64=0.01)

    sharp = sharpness(d, h)
    var = variance(d)

    return sharp^2 / var
end

## . - -- .- .-.- . .- .-- -.- .- .- --.-.... - -.-
function _select(M, _set, _comp)
    # TODO: check views
    _M = zeros(eltype(M), length(_set), length(_comp))
    _select!(M, _M, _set, _comp)
    _M
end

export _select!
function _select!(Msrc, Mdest, _set, _comp)
    # TODO: check views
    for (i, I) in enumerate(_set), (j, J) in enumerate(_comp)
        Mdest[i, j] = Msrc[I, J]
    end
    return Mdest
end

## . - -- .- .-.- . .- .-- -.- .- .- --.-.... - -.-
export _ddist
function _ddist(Q)
    k = size(Q, 1)
    _DD = zeros(k, k)
    for ri in 1:k, ci in 1:k
        ri == ci && continue
        _DD[ri, ci] = - Q[ri, ci]/ Q[ri, ri]
    end
    return _DD
end

## . - -- .- .-.- . .- .-- -.- .- .- --.-.... - -.-
export  _jump_count
function _jump_count(s, on_ind)
    _jump_count = 1
    _last_is_on = on_ind[s[1]]
    for si in s
        (on_ind[si] ⊻ _last_is_on) || continue
        _jump_count += 1
        _last_is_on = on_ind[si]
    end
    _jump_count
end

