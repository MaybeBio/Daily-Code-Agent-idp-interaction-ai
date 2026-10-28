using Droplet

# Minimal test to debug BulkMonomer deltaHam!
l = Lattice(boxlength=16, nspecies=1)
b = BulkMonomer(length=10)
intPot!(b, -1.0)
init!(l, [b])

# Initial ham
ham!(b)
println("Initial ham: ", b.ham)

# Make 5 moves and track both incremental and full recomputation
for i in 1:5
    # Store current ham before move
    ham_before = b.ham
    
    # Try a move
    moveIdx = rand(1:b.length)
    b.moveIdx = moveIdx
    identifyBond!(b)
    proposeMove!(b)
    
    if !overlap(b)
        # Compute deltaH
        old_pos = copy(b.bi)
        new_pos = copy(b.newPos)
        
        # Check neighbors before deltaHam!
        old_neighbors = nnp(b, b, old_pos)
        new_neighbors = nnp(b, b, new_pos)
        
        # Check if old and new positions are neighbors
        are_neighbors = (abs(old_pos[1] - new_pos[1]) + abs(old_pos[2] - new_pos[2]) + abs(old_pos[3] - new_pos[3])) == 1
        
        deltaHam!(b)
        delta = b.deltaH
        println("\nMove $i:")
        println("  Old pos: $old_pos, neighbors: $old_neighbors")
        println("  New pos: $new_pos, neighbors: $new_neighbors")
        println("  Are old/new neighbors? $are_neighbors")
        if are_neighbors
            println("  Adjusted new_neighbors (excluding self): $(new_neighbors - 1)")
        end
        println("  Expected delta (old formula): $((new_neighbors - old_neighbors - 1) * -1.0)")
        println("  Expected delta (new formula): $((new_neighbors - old_neighbors) * -1.0 * 0.5)")
        println("  Actual deltaH: $delta")
        println("  ham before + deltaH: $(ham_before + delta)")
        
        # Make the move
        makeMove!(b)
        
        # Recompute from scratch
        ham!(b)
        println("  ham after (recomputed): $(b.ham)")
        
        # Manually count interactions
        total_neighbors = 0.0
        for i in 1:b.length
            pos = b.pos[i][1:3]
            neighbors = nnp(b, b, pos)
            total_neighbors += neighbors
            if neighbors > 0
                println("    Particle $i at $pos has $neighbors neighbors")
            end
        end
        println("  Total neighbor count: $total_neighbors")
        println("  Expected ham: $(- b.chemPot * b.length + total_neighbors * -1.0 * 0.5)")
        println("  Difference from old formula: $(b.ham - (ham_before + (new_neighbors - old_neighbors - 1) * -1.0))")
        println("  Difference from new formula: $(b.ham - (ham_before + delta))")
    else
        println("\nMove $i: overlap, skipped")
    end
end
