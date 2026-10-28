Right now there is a type called PolyMonomer, which is meant to be a constituent of the Polymer structure as far as interaction potential and chemical potential is concerned. It is meant to fill in these attributes on behalf of the Polymer class in order to make multicomponent interaction possible with bead-dependent interaction energies.

They do not store their own position (latS and pos), and do not posess any of the runtime and initialization containers relating to position, as for the BulkMonomers and FastBulkMonomers types. They only behave like these other types as far as the interaction matrix and intSpecies machinery is concerned, as specified in the interaction.jl file.

I would like to be able to initialize PolyMonomers in the following way:

l = Lattice(;
    boxlength = 4,
    nspecies = 3
)

tf = FastBulkMonomer(;
    length = 12,
    initShape = "uniform"
)

poly = Polymer(;
    length = l.boxlength,
    maxlength = 8,
    nMonomerTypes = 3,
    chemPot = 0.0,
    intPot = 0.0
)

initMonomers!(poly)

dna = poly.monomers[1]
promoter = poly.monomers[2]

dna.name = "dna"
promoter.name = "promoter"

intPot!(tf, dna, 0.3)
intPot!(tf, promoter, 2.0)
species = [tf, dna, promoter]

init!(l, species)

For this to be possible, hamiltonian needs to be set up to use a monomer's position and lattice in two different ways:
1. if it is a bulk monomer, use its own position and lattice
2. if it is a polymonomer, use it's parent polymer's position and lattice.

This will come up when initializing PolyMonomer species. I have begun setting this up, I need you to go over all the relevant parts of the code and complete the process, up to the point that the init! function works as above. After that is done, please write and run a comprehensive set of tests or the init and related functions for the multicomponent polymer functionalities. 

Don't modify the hamiltonian functions for now beyond what is required for init!. Don't start implementing the modifications required for move and other functionalities. This will be done later.