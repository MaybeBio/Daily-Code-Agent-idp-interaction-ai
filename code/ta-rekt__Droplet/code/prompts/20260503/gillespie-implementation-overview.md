# Gillespie implementation overview

### Building the list of rates and moves

- each possible move is associated with a particle, a move type, a direction and a rate. there could be a list of 5-tuples which records all this (the fifth element being the cumulative sum of rates up to that entry), and must be sorted by rate before each new move
- list is built up once at initialization, and updated after each local move
- list is built by looping through each particle in the lattice, calling buildMoves. this function lists all moves from the given particle into non-occupied neighboring lattice points, for BulkMonomer and analogously for Polymer (to be defined in more detail)
- for each acceptable move, an entry is pushed to the list of the form (particleId, moveType, direction, rate, cumulativeRate). particleId is some global index to be specified, there is only one moveType for BulkMonomer, and three for Polymer: add, kink, remove.
- rate is calculated as k_ij = k0 * exp(-T*deltaH(i,j)), where k0 is a global timescale that does not matter for the purposes of this algorithm, it gets normalized away during the move selection. deltaH is computed by the deltaHam function.
- cumulativeRate of entry i is the sum of rates of all entries 1 < j < i 


### Updating the list of rates at each move

- each time a move happens, the only rates that can be affected are those in the union of the neighborhood of the moved particle's old position and that of its new position (for polymer, the union of both particles that were moved). these can be obtained by the nnIdxs function in utils.jl
- after a move, the affected neighborhood is looped over, and the entries in the rateList are updated for each element
- the rateList might need to be locally re-sorted when that happens. this needs to be done efficiently. we need to think of the right data structure for this
- the value of the sum of rates is updated (without re-summing over all rates, just add the increment from the local updates)


### Sampling a new move

- before each move, a uniform random number is sampled between 1 and Q, Q being the sum of all rates in the rateList (the cumulativeRate of the last element)
- a binary search over the cumulativeRates determines the position of the particle to be moved

