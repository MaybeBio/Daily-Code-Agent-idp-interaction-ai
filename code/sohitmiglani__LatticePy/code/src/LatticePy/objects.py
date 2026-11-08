class amino_acid():
    """One bead (residue) of a lattice polymer.

    Attributes
    ----------
    polarity : int
        +1 or -1; selects the contact energies in ``lattice.bond_energies``.
    coordinates : list of int
        Current lattice site [x, y, z].
    next, previous : list of int or None
        Coordinates of the bonded neighbours along the chain; None at the chain ends.
    polymer : int
        Id of the chain this bead belongs to.
    neighbors : dict
        Currently unused.
    """
    def __init__(self, polarity, coordinates, polymer_id=None):
        """Create a bead with the given polarity at ``coordinates``, belonging to chain ``polymer_id``.

        The bead starts unbonded (``next`` and ``previous`` are None); the caller links it into a chain.
        """
        self.polarity = polarity
        self.coordinates = coordinates
        self.next = None
        self.previous = None
        self.polymer = polymer_id
        self.neighbors = dict()