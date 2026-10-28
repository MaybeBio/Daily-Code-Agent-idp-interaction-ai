# . - -- .- .-.- . .- .-- -.- .- .- --.-.... - -.-
export get_residence_time
function get_residence_time(s, t, on_ind)
    
    jc = _jump_count(s, on_ind)
    rt = zeros(Float64, jc)

    _last_is_on = on_ind[s[1]]

    ji = 1 # jump index
    for (si, ti) in zip(s, t)
        if (on_ind[si] ⊻ _last_is_on)
            _last_is_on = on_ind[si]
            ji += 1
        end
        rt[ji] += ti
    end
    # postprocess
    s_in_on = on_ind[s[1]]
    onrt = s_in_on ? rt[1:2:end] : rt[2:2:end]
    offrt = s_in_on ? rt[2:2:end] : rt[1:2:end]
    return s_in_on, rt, onrt, offrt
end

# . - -- .- .-.- . .- .-- -.- .- .- --.-.... - -.-
export get_residence_time!
function get_residence_time!(m::LLROptModel)
    s_in_on, rt, onrt, offrt = get_residence_time(m.s, m.t, m.on_ind)
    m.s_in_on = s_in_on
    m.rt = rt
    m.onrt = onrt
    m.offrt = offrt
    return m
end