using TimerOutputs
using Droplet  # your package/module

const global to = TimerOutput()

# Redefine run() to instrument function calls
function run_profiled(niter)
    for i in 1:niter
        # Wrap the move! call with a timer
        @timeit to "move!" begin
            move!(s1)
        end
    end
end

let 
    l = Lattice(;
        boxlength = 32,
        nspecies = 2
    )

    # polymer
    global s1 = Species(;
        length = 32,
        maxlength = 256,
        boxlength = l.boxlength,
        intPot = -1,
        chemPot = 1
    )

    # dilute monomers
    global s2 = Species(;
        length = 256,
        maxlength = 256,
        boxlength = l.boxlength,
        intPot = 0.0,
        chemPot = 0.0
    )

    # initializing lattice and species
    initPoly!(s1)
    initDilute!(s2)
    initLattice!(l, [s1, s2])
    
    # initialize hamiltonians
    initHamPoly!(s1)
    hamInitial = s1.ham

    # Initialize niter and s1 appropriately before profiling
    niter = 1_000_000

    # Run the instrumented version and print the timing report
    run_profiled(niter)
    println("Profiling report:")
    show(to; allocations=false)

end