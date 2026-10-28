#!/usr/bin/env julia
"""
Minimal allocation tracking script.
Run with: julia --track-allocation=user test/alloc_track.jl
Then look at *.mem files in src/
"""

using Droplet
using Random
using Profile

Random.seed!(42)

# Setup (same as benchmark)
l = Lattice(; boxlength = 42, nspecies = 2)
poly = Polymer(; length = 32, maxlength = 1250, initShape = "line",
               nMonomerTypes = 12, initPos = [5, 21, 21], pattern = :uniform)
initMonomers!(poly, defaultChemPot = 1.0)
tf = BulkMonomer(; length = 1, initShape = "uniform")
slist = [poly, tf]

base = poly.monomers[1]
for (key, pat) in poly.monomers
    intPot!(pat, 0.5)
    intPot!(base, pat, 0.5)
end
intPot!(poly.monomers[3], tf, -0.4)
intPot!(poly.monomers[5], tf, -0.4)
intPot!(poly.monomers[7], tf, -0.4)
intPot!(poly.monomers[9], tf, -0.4)
intPot!(poly.monomers[11], tf, -0.4)
intPot!(tf, -1.0 / 0.66)

tf.length = Int32(round(l.boxlength^3 * 0.01))
updateMonomers!(poly, base.intIdx)
init!(l, slist)
poly.deltaHP = 0.0

# Pattern setup
patternVec = zeros(Int32, poly.maxlength)
headPattern1 = poly.monomers[2]
tailPattern1 = poly.monomers[4]
tailPattern2 = poly.monomers[6]
tailPattern3 = poly.monomers[8]
tailPattern4 = poly.monomers[10]
tailPattern5 = poly.monomers[12]

patternVec[1:200] .= headPattern1.intIdx
patternVec[201:210] .= poly.monomers[3].intIdx
patternVec[211:410] .= tailPattern1.intIdx
patternVec[411:420] .= poly.monomers[5].intIdx
patternVec[421:620] .= tailPattern2.intIdx
patternVec[621:630] .= poly.monomers[7].intIdx
patternVec[631:830] .= tailPattern3.intIdx
patternVec[831:840] .= poly.monomers[9].intIdx
patternVec[841:1040] .= tailPattern4.intIdx
patternVec[1041:1050] .= poly.monomers[11].intIdx
patternVec[1051:1243] .= tailPattern5.intIdx
patternVec[1244:1250] .= base.intIdx

growPolymer!(poly, minAttempts = 900000)
updatePattern!(l, poly, patternVec)

# Warm up
for _ in 1:1000
    move!(poly)
    move!(tf)
end

# Clear allocation counters
Profile.clear_malloc_data()

# Run the hot loop
println("Running allocation tracking...")
for _ in 1:10000
    move!(poly)
    move!(tf)
end

println("Done. Check *.mem files in src/ directory.")
