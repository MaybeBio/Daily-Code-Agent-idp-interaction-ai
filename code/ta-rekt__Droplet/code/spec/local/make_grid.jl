# make_grid.jl -- generate one job folder per point of a parameter grid.
#
# For each combination of the swept axes it creates a folder under OUTROOT,
# copies the verbatim spec file(s), and writes a copy of param.jl with the
# per-cell + fixed overrides substituted IN PLACE (comments/formatting kept).
#
# Reuse for any sweep by editing only the CONFIG block: name each param.jl
# knob by its `section.field` path, where `section` is a field of FisherParams
# (surface / bvp / rouse / fpt / sweep) and `field` is a field of that struct.
# The engine finds the field's default -- whether it is written inline in the
# `Parameters( ... )` call or as a `struct XxxParams` default -- and rewrites it.
#
#   julia spec/jobs/make_grid.jl

# ============================ CONFIG (edit me) ==============================
const SRCDIR  = joinpath(@__DIR__, "20260721")   # template param.jl + spec live here
const OUTROOT = joinpath(@__DIR__, "20260722")   # one folder per grid cell goes here

# Swept axes: name => scalar values.  Their cartesian product is the grid.
const AXES = [
    "R"  => [2, 3, 4, 5],
    "h0" => [-0.4, -0.5, -0.6, -0.7],
]

# Per-cell overrides: axis combo `c` -> ["section.field" => value, ...].
overrides(c) = [
    "surface.R"            => float(c["R"]),   # sphere radius
    "sweep.omega_b_values" => [c["h0"]],       # the run sweeps this axis; 1 value/folder
]

# Fixed overrides applied to EVERY folder.
const FIXED = [
    "bvp.Nrho"             => 200,
    "bvp.dz_target"        => 0.05,
    "bvp.n_d"              => 50,
    "bvp.d_max_factor"     => 1.0,
    "sweep.n_phi"          => 16,
    "sweep.phi_inf_values" => Float64[],       # empty => the n_phi grid is used
]

# Folder name from an axis combo.  `_g` drops a trailing `.0` on integral reals.
_g(x) = (x isa Real && isinteger(x)) ? string(Int(x)) : string(x)
label(c) = "R_$(_g(c["R"]))_h0_$(_g(c["h0"]))"

# Files copied verbatim into every folder (besides the generated param.jl).
const COPY_VERBATIM = ["fisher_two_spheres.jl"]
# ============================================================================


# ================================ ENGINE ====================================
# Render a value as the Julia literal to splice into the source.
_lit(x::Symbol)         = repr(x)          # :bvp
_lit(x::AbstractString) = repr(x)          # "foo"
_lit(x)                 = string(x)        # 200, 0.05, 1.0, [-0.4], Float64[]

# Index of the ')' matching the '(' at src[open].
function _match_paren(src::AbstractString, open::Int)
    depth = 0
    i = open
    while i <= lastindex(src)
        c = src[i]
        c == '(' && (depth += 1)
        c == ')' && (depth -= 1; depth == 0 && return i)
        i = nextind(src, i)
    end
    error("unbalanced parens from index $open")
end

# Apply one "section.field" => value override to the param.jl source text.
function apply_override(src::String, path::AbstractString, value)::String
    section, field = split(path, '.'; limit = 2)
    lit = _lit(value)

    # Locate `section::Type` inside FisherParams to learn the section's type.
    sm = match(Regex("\\b$(section)\\s*::\\s*(\\w+)\\s*="), src)
    sm === nothing && error("section `$section` not found (path `$path`)")
    typ = sm.captures[1]

    # The section's default constructor call `Type( ... )`.
    call = findnext("$typ(", src, sm.offset)
    call === nothing && error("constructor `$typ(` not found for `$section`")
    popen  = last(call)                       # the '('
    pclose = _match_paren(src, popen)         # its matching ')'
    expr   = src[popen:pclose]                # "( ... )"

    fm = match(Regex("(\\b$(field)\\s*=\\s*)([^,\\)\\n]*)"), expr)
    if fm !== nothing
        # Field is set inline in the constructor call -> rewrite it there.
        v0, v1 = fm.offsets[2], fm.offsets[2] + lastindex(fm.captures[2]) - 1
        new_expr = expr[1:prevind(expr, v0)] * lit * expr[nextind(expr, v1):end]
        return src[1:prevind(src, popen)] * new_expr * src[nextind(src, pclose):end]
    end

    # Otherwise the field's default lives in `struct Type ... end`.
    sblk = findfirst("struct $typ", src)
    sblk === nothing && error("`struct $typ` not found for field `$field`")
    eend = match(r"(?m)^end", src[last(sblk):end])
    eend === nothing && error("no closing `end` for `struct $typ`")
    lo, hi = first(sblk), last(sblk) + eend.offset + 1   # struct header .. end
    block  = src[lo:hi]

    lm = match(Regex("(?m)^(\\s*$(field)\\s*::[^=\\n]*=\\s*)([^\\n]*?)(\\s*(?:#.*)?)\$"),
               block)
    lm === nothing && error("field `$field` not found in `struct $typ`")
    v0, v1 = lm.offsets[2], lm.offsets[2] + lastindex(lm.captures[2]) - 1
    new_block = block[1:prevind(block, v0)] * lit * block[nextind(block, v1):end]
    return src[1:prevind(src, lo)] * new_block * src[nextind(src, hi):end]
end

apply_all(src, ovs) = foldl((s, kv) -> apply_override(s, kv[1], kv[2]), ovs; init = src)

# Cartesian product of AXES -> Vector{Dict(name => value)}.
function _combos(axes)
    combos = [Dict{String,Any}()]
    for (name, vals) in axes
        combos = [merge(c, Dict(name => v)) for c in combos for v in vals]
    end
    return combos
end

function main()
    template = read(joinpath(SRCDIR, "param.jl"), String)
    mkpath(OUTROOT)
    for c in _combos(AXES)
        dir = joinpath(OUTROOT, label(c))
        mkpath(dir)
        for f in COPY_VERBATIM
            cp(joinpath(SRCDIR, f), joinpath(dir, f); force = true)
        end
        write(joinpath(dir, "param.jl"), apply_all(template, vcat(FIXED, overrides(c))))
        println("wrote ", joinpath(basename(OUTROOT), label(c)))
    end
end

main()
# ============================================================================
