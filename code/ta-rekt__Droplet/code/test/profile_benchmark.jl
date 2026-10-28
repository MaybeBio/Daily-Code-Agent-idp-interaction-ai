#!/usr/bin/env julia
"""
Profiling script for the multi_enhancer_search_benchmark.
Measures wall clock time in each major function using TimerOutputs.
Also checks for type instabilities and allocations.
"""

using Droplet
using TimerOutputs
using Profile
using Random
using InteractiveUtils

# Create a global timer
const to = TimerOutput()

# Re-export timer macros for easy use
macro timed_call(name, expr)
    quote
        @timeit to $(name) $(esc(expr))
    end
end

println("="^60)
println("PROFILING MULTI-ENHANCER SEARCH BENCHMARK")
println("="^60)

#=============================================================================
# PHASE 1: Setup (same as benchmark)
=============================================================================#

Random.seed!(42)  # For reproducibility

l = Lattice(;
    boxlength = 42,
    nspecies = 2
)

poly = Polymer(;
    length = 32,
    maxlength = 1250,
    initShape = "line",
    nMonomerTypes = 12,
    initPos = [5, 21, 21],
    pattern = :uniform
)

initMonomers!(poly, defaultChemPot = 1.0)

tf = BulkMonomer(;
    length = 1,
    initShape = "uniform"
)

slist = [poly, tf]

# Set up interactions
base = poly.monomers[1]
mainPattern1 = poly.monomers[3]
mainPattern2 = poly.monomers[5]
mainPattern3 = poly.monomers[7]
mainPattern4 = poly.monomers[9]
mainPattern5 = poly.monomers[11]

concentration = 0.01
temperature = 0.66
beta = -0.4

for (key, pat) in poly.monomers
    intPot!(pat, 0.5)
    intPot!(base, pat, 0.5)
end

intPot!(mainPattern1, tf, beta)
intPot!(mainPattern2, tf, beta)
intPot!(mainPattern3, tf, beta)
intPot!(mainPattern4, tf, beta)
intPot!(mainPattern5, tf, beta)
intPot!(tf, -1.0 / temperature)

tf.length = Int32(round(l.boxlength^3 * concentration))
updateMonomers!(poly, base.intIdx)
init!(l, slist)
poly.deltaHP = 0.0

# Build polymer patterns
patternVec = zeros(Int32, poly.maxlength)
p1_start, p1_end = 201, 210
p2_start, p2_end = 411, 420
p3_start, p3_end = 621, 630
p4_start, p4_end = 831, 840
p5_start, p5_end = 1041, 1050

headPattern1 = poly.monomers[2]
tailPattern1 = poly.monomers[4]
tailPattern2 = poly.monomers[6]
tailPattern3 = poly.monomers[8]
tailPattern4 = poly.monomers[10]
tailPattern5 = poly.monomers[12]

patternVec[1:p1_start-1] .= headPattern1.intIdx
patternVec[p1_start:p1_end] .= mainPattern1.intIdx
patternVec[p1_end+1:p2_start-1] .= tailPattern1.intIdx
patternVec[p2_start:p2_end] .= mainPattern2.intIdx
patternVec[p2_end+1:p3_start-1] .= tailPattern2.intIdx
patternVec[p3_start:p3_end] .= mainPattern3.intIdx
patternVec[p3_end+1:p4_start-1] .= tailPattern3.intIdx
patternVec[p4_start:p4_end] .= mainPattern4.intIdx
patternVec[p4_end+1:p5_start-1] .= tailPattern4.intIdx
patternVec[p5_start:p5_end] .= mainPattern5.intIdx
patternVec[p5_end+1:poly.maxlength-7] .= tailPattern5.intIdx
patternVec[poly.maxlength-6:poly.maxlength] .= base.intIdx

growPolymer!(poly, minAttempts = 900000)
updatePattern!(l, poly, patternVec)

println("\nSetup complete:")
println("  Polymer length: $(poly.length)")
println("  TF count: $(tf.length)")
println("  Box length: $(l.boxlength)")

#=============================================================================
# PHASE 2: Type stability check
=============================================================================#

println("\n" * "="^60)
println("TYPE STABILITY CHECK")
println("="^60)

println("\n--- move!(::Polymer) ---")
@code_warntype move!(poly)

println("\n--- move!(::BulkMonomer) ---")
@code_warntype move!(tf)

#=============================================================================
# PHASE 3: Micro-benchmarks of hot functions
=============================================================================#

println("\n" * "="^60)
println("MICRO-BENCHMARKS (10,000 iterations each)")
println("="^60)

# Warm up
for _ in 1:100
    move!(poly)
    move!(tf)
end

# Benchmark move!(poly)
n_iters = 10000
t_poly_moves = @elapsed begin
    for _ in 1:n_iters
        move!(poly)
    end
end
println("\nmove!(poly): $(round(t_poly_moves * 1e6 / n_iters, digits=2)) μs/call")

# Benchmark move!(tf)
t_tf_moves = @elapsed begin
    for _ in 1:n_iters
        move!(tf)
    end
end
println("move!(tf): $(round(t_tf_moves * 1e6 / n_iters, digits=2)) μs/call")

#=============================================================================
# PHASE 4: Allocation tracking
=============================================================================#

println("\n" * "="^60)
println("ALLOCATION TRACKING (1000 iterations)")
println("="^60)

# Force GC before measuring
GC.gc()

n_alloc_test = 1000

# Measure allocations for move!(poly)
allocs_poly = @allocated begin
    for _ in 1:n_alloc_test
        move!(poly)
    end
end
println("\nmove!(poly) allocations: $(allocs_poly ÷ n_alloc_test) bytes/call")

# Measure allocations for move!(tf)
allocs_tf = @allocated begin
    for _ in 1:n_alloc_test
        move!(tf)
    end
end
println("move!(tf) allocations: $(allocs_tf ÷ n_alloc_test) bytes/call")

#=============================================================================
# PHASE 5: Timed simulation run
=============================================================================#

println("\n" * "="^60)
println("TIMED SIMULATION (equivalent to benchmark)")
println("="^60)

# Reset timer
reset_timer!(to)

movingsps = [(poly, tf.intIdx + 1), (tf, tf.intIdx)]

function run_timed_simulation!(poly, tf, movingsps, to, max_frames, step)
    nmoves = 0
    niter = 0
    nframes = 1
    
    @timeit to "total_simulation" begin
        while nframes <= max_frames
            niter += 1

            for (s, i) in movingsps
                @timeit to "move_attempt" begin
                    if s isa Polymer
                        @timeit to "move_poly" begin
                            accepted = move!(s)
                        end
                    else
                        @timeit to "move_tf" begin
                            accepted = move!(s)
                        end
                    end
                end
                
                if accepted
                    nmoves += 1
                    if mod(nmoves, step) == 0
                        nframes += 1
                    end
                end
            end
        end
    end
    
    return niter, nmoves
end

niter_result, nmoves_result = run_timed_simulation!(poly, tf, movingsps, to, 50, 10000)

println("\nSimulation stats:")
println("  Total iterations: $niter_result")
println("  Accepted moves: $nmoves_result")
println("  Acceptance rate: $(round(nmoves_result / (niter_result * length(movingsps)) * 100, digits=2))%")

println("\nTimer breakdown:")
show(to)
println()

#=============================================================================
# PHASE 6: CPU Profile (sampling)
=============================================================================#

println("\n" * "="^60)
println("CPU PROFILE (sampling, 100k iterations)")
println("="^60)

# Clear previous profile data
Profile.clear()

# Profile with sampling
@profile begin
    for _ in 1:100000
        move!(poly)
        move!(tf)
    end
end

println("\nTop functions by sample count:")
Profile.print(maxdepth=10, mincount=100, noisefloor=2)

println("\n" * "="^60)
println("PROFILING COMPLETE")
println("="^60)
