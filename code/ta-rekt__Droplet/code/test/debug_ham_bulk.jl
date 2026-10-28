using Droplet
using Random

# Debug script to find the rare bug where b_hams diverges from b_accDeltaH
# in heterogeneous polymer with bulk interaction test

function run_debug_test(seed::Int)
    Random.seed!(seed)
    
    l = Lattice(;
        boxlength = 20,
        nspecies = 4
    )

    p = Polymer(;
        length = 10,
        maxlength = 30,
        nMonomerTypes = 2,
        initPos = [1,10,10],
        pattern = :custom,
        customPattern = [1,1,1,1,2,2,2,2,1,1]
    )
    
    initMonomers!(p, defaultChemPot=1.0, defaultIntPot=-0.3)
    
    m1 = p.monomers[1]
    m2 = p.monomers[2]

    b = FastBulkMonomer(;
        length = 64,
        initShape = "uniform"
    )

    c = FastBulkMonomer(;
        length = 64,
        initShape = "uniform"
    )

    intPot!(m1, -0.3)
    intPot!(m2, -0.5)
    intPot!(b, -0.2)
    intPot!(c, -0.1)
    
    intPot!(m1, b, -0.1)
    intPot!(m1, c, -0.05)
    intPot!(m2, b, -0.4)
    intPot!(m2, c, -0.2)

    slist = [p, b, c]
    init!(l, slist)

    # Initial hamiltonians
    ham!(p)
    ham!(b)
    ham!(c)
    
    p_accDeltaH = p.ham
    b_accDeltaH = b.ham
    c_accDeltaH = c.ham
    
    p.deltaH = 0.0
    b.deltaH = 0.0
    c.deltaH = 0.0
    p.deltaHP = 0.0
    b.deltaHP = 0.0
    c.deltaHP = 0.0

    movable = [p, b, c]
    names = ["polymer", "bulk_b", "bulk_c"]
    
    frame = 0
    max_frames = 10000
    
    # Track last successful polymer move info
    last_poly_move = nothing
    
    while frame < max_frames
        chooser = rand(1:3)
        s = movable[chooser]
        
        if move!(s)
            frame += 1
            
            # Capture polymer move info AFTER successful move
            if s isa Polymer
                last_poly_move = (
                    frame = frame,
                    moveType = s.moveType,
                    ki = s.ki,
                    kj = s.kj,
                    moveIdx = s.moveIdx,
                    newPosi = s.newPosi,
                    newPosj = s.newPosj,
                    polymer_length = s.length,
                    deltaH = s.deltaH,
                    deltaHP = s.deltaHP
                )
            end
            
            # Before accumulating, record the deltas
            mover_deltaH = getDeltaHam(s)
            partner_deltaHPs = [getDeltaHP(movable[j]) for j in 1:3 if j != chooser]
            
            # Accumulate deltaH for mover
            if chooser == 1
                p_accDeltaH += mover_deltaH
            elseif chooser == 2
                b_accDeltaH += mover_deltaH
            else
                c_accDeltaH += mover_deltaH
            end
            
            # Accumulate deltaHP for partners
            partner_idx = 1
            for j in 1:3
                if j != chooser
                    if j == 1
                        p_accDeltaH += getDeltaHP(movable[j])
                    elseif j == 2
                        b_accDeltaH += getDeltaHP(movable[j])
                    else
                        c_accDeltaH += getDeltaHP(movable[j])
                    end
                    partner_idx += 1
                end
            end
            
            # Compute fresh ham for all species
            ham!(p)
            ham!(b)
            ham!(c)
            
            # Check for discrepancy
            p_diff = abs(p.ham - p_accDeltaH)
            b_diff = abs(b.ham - b_accDeltaH)
            c_diff = abs(c.ham - c_accDeltaH)
            
            if p_diff > 1e-10 || b_diff > 1e-10 || c_diff > 1e-10
                println("\n=== DISCREPANCY FOUND at frame $frame (seed=$seed) ===")
                println("Mover: $(names[chooser])")
                
                if chooser == 1 && last_poly_move !== nothing
                    info = last_poly_move
                    println("Move type: $(info.moveType)")
                    println("  ki=$(info.ki), kj=$(info.kj)")
                    println("  moveIdx=$(info.moveIdx)")
                    println("  newPosi=$(info.newPosi), newPosj=$(info.newPosj)")
                    println("  polymer_length: $(info.polymer_length)")
                    println("  polymer deltaH: $(info.deltaH)")
                    println("  polymer deltaHP: $(info.deltaHP)")
                else
                    println("Move type: bulk move")
                end
                
                println("\nMover deltaH: $mover_deltaH")
                println("Partner deltaHPs: $partner_deltaHPs")
                
                println("\nPolymer:")
                println("  fresh ham: $(p.ham)")
                println("  accumulated: $p_accDeltaH")
                println("  diff: $p_diff")
                
                println("\nBulk b:")
                println("  fresh ham: $(b.ham)")
                println("  accumulated: $b_accDeltaH")  
                println("  diff: $b_diff")
                
                println("\nBulk c:")
                println("  fresh ham: $(c.ham)")
                println("  accumulated: $c_accDeltaH")
                println("  diff: $c_diff")
                
                return (false, frame, chooser, chooser == 1 && last_poly_move !== nothing ? last_poly_move.moveType : 0)
            end
        end
    end
    
    return (true, frame, 0, 0)
end

# Run multiple times with different seeds
println("Searching for the rare bug...")
for seed in 1:1000
    result, frame, chooser, moveType = run_debug_test(seed)
    if !result
        println("\nBug found with seed $seed at frame $frame")
        println("Mover index: $chooser, Move type: $moveType")
        break
    end
    if seed % 100 == 0
        println("Tested $seed seeds, no bug found yet...")
    end
end
println("Done searching")
