using Droplet
using Random

function countBeadTypes(p::Polymer)
    counts = Dict{Int32, Int}()
    key = p.head
    while key != 0
        intIdx = p.pos[key][6]
        counts[intIdx] = get(counts, intIdx, 0) + 1
        key = p.pos[key][5]
    end
    return counts
end

function getPatternVector(p::Polymer)
    result = Int32[]
    key = p.head
    while key != 0
        push!(result, p.pos[key][6])
        key = p.pos[key][5]
    end
    return result
end

function isPatternContiguous(p::Polymer)
    if p.patternTails === nothing || isempty(p.patternTails)
        return true
    end
    patternVec = getPatternVector(p)
    for intIdx in keys(p.patternTails)
        if intIdx == p.baseIntIdx
            continue
        end
        indices = findall(x -> x == intIdx, patternVec)
        if isempty(indices)
            continue
        end
        if indices != collect(indices[1]:indices[end])
            return false
        end
    end
    return true
end

Random.seed!(42)

l = Lattice(boxlength = 64, nspecies = 1)
customPattern = vcat(fill(1, 16), fill(2, 16), fill(1, 32))
poly = Polymer(
    length = 64,
    maxlength = 128,
    nMonomerTypes = 2,
    pattern = :custom,
    customPattern = customPattern
)

initMonomers!(poly, defaultChemPot=0.0, defaultIntPot=0.0)
init!(l, [poly])

m2 = poly.monomers[2]
initial_count = countBeadTypes(poly)[m2.intIdx]
println("Initial pattern count: ", initial_count)
println("Initial patternTails: ", poly.patternTails)

for iter in 1:5000
    lengthBefore = poly.length
    countBefore = countBeadTypes(poly)[m2.intIdx]
    
    # Use the real move! function
    result = move!(poly)
    
    lengthAfter = poly.length
    countAfter = countBeadTypes(poly)[m2.intIdx]
    
    if !isPatternContiguous(poly)
        println("FRAGMENTED at move ", iter, "!")
        println("  Before: length=", lengthBefore, ", count=", countBefore)
        println("  After:  length=", lengthAfter, ", count=", countAfter)
        println("  Result: ", result)
        println("  Pattern vector: ", getPatternVector(poly))
        break
    end
    
    if countAfter != initial_count
        println("COUNT CHANGED at move ", iter, ": ", initial_count, " -> ", countAfter)
        break
    end
end

println("Completed - if no errors above, tests passed for 5000 moves")
