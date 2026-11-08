import sys
import copy
import math
import plotly
import random
from datetime import datetime
import numpy as np
import pandas as pd
import networkx as nx
import plotly.io as pio
from .objects import amino_acid
import matplotlib.pyplot as plt
import plotly.graph_objects as go

system_random = random.SystemRandom()
pio.renderers.default = 'iframe_connected'

class lattice():
    """Lattice polymers on a periodic simple cubic lattice, sampled by Metropolis Monte Carlo.

    Beads are stored in ``space``, a dict keyed by ``str([x, y, z])``. ``start`` and ``last``
    map the site of each chain's first and last bead to the chain's polymer id. Every
    non-bonded nearest-neighbour pair of beads contributes ``bond_energies[p1][p2]``. By default
    that is -2.3 - E_c for (-1, -1), -1 - E_c for (-1, +1) and -E_c for (+1, +1); a custom table
    can be passed to the constructor.

    Running totals kept up to date by the moves:
    ``energy``, ``native_contacts`` (number of non-bonded contacts),
    ``non_covalent_hydrophobic_contacts`` (number of (-1, -1) contacts), and the per-chain
    lists ``individual_energy``, ``individual_ncs`` and ``individual_nchcs`` (intramolecular
    contacts only), indexed by polymer id.
    """
    def __init__(self, bound, E_c, bond_energies=None, beta=0, lattice_type='simple_cubic', record_moves=False, allow_cluster_move=True):
        """Create an empty lattice.

        Parameters
        ----------
        bound : int
            Largest coordinate. Sites run from 0 to ``bound``, so the box is ``bound + 1`` sites wide
            and coordinates wrap modulo ``bound + 1``.
        E_c : float
            Constant subtracted from every contact energy in the default table. Ignored when
            ``bond_energies`` is given.
        bond_energies : dict of dict, optional
            Contact energy for each pair of polarities, ``bond_energies[p1][p2]``, e.g.
            ``{-1: {-1: -2.3, 1: -1}, 1: {-1: -1, 1: 0}}``. Used as given. If None, the default
            table is built from ``E_c``.
        beta : float, optional
            Inverse temperature. Note that ``simulate`` sets it to ``beta_lower_bound`` at the start of a run.
        lattice_type : str, optional
            Only 'simple_cubic' is implemented.
        record_moves : bool, optional
            If True, record the name of every proposed move (``move_records``) and whether it was
            accepted (``acceptance_records``).
        allow_cluster_move : bool, optional
            If True, ``simulate`` includes ``cluster_move`` when there is more than one chain.
        """
        self.lattice_type = lattice_type
        self.record_moves = record_moves
        self.bound = bound
        self.space = dict()
        self.start = dict()
        self.last = dict()
        self.energy = 0
        self.native_contacts = 0
        self.non_covalent_hydrophobic_contacts = 0
        self.beta = beta
        self.n_mcmc = 0
        self.length_of_polymer = 0
        self.energy_records = []
        self.beta_records = []
        self.native_contacts_records = []
        self.nchc_records = []
        self.move_records = []
        self.acceptance_records = []
        self.records = []
        self.n_polymers = 0
        self.plateau_time = 0
        self.plateaued = False
        self.individual_energy = []
        self.individual_energy_records = []
        self.individual_ncs = []
        self.individual_ncs_records = []
        self.individual_nchcs = []
        self.individual_nchcs_records = []
        self.allow_cluster_move = allow_cluster_move
        self.cluster_records = []
        self.cluster_size_records = []

        if bond_energies is None:
            bond_energies = {1: {}, -1: {}}
            bonds = [-1,-1,1,-1,1,1]
            energies = [-2.3-E_c, -1-E_c, -E_c]
            for i in range(3):
                bond_energies[bonds[2*i]][bonds[2*i+1]] = energies[i]
                bond_energies[bonds[2*i+1]][bonds[2*i]] = energies[i]
        elif not isinstance(bond_energies, dict):
            raise TypeError('bond_energies must be a dict of dicts such as {-1: {-1: -2.3, 1: -1}, 1: {-1: -1, 1: 0}}, '
                            'got ' + repr(bond_energies) + '. Pass beta by keyword (beta=...).')
        self.bond_energies = bond_energies

    def periodic_coordinate(self, coordinates):
        """Wrap a position into the periodic box.

        Parameters
        ----------
        coordinates : sequence of 3 numbers
            Position, possibly outside the box.

        Returns
        -------
        list of int
            Each component taken modulo ``bound + 1`` and rounded.
        """
        return [round(coordinates[0]%(self.bound +1)), 
                round(coordinates[1]%(self.bound +1)), 
                round(coordinates[2]%(self.bound +1))]
                
    def find_valid_neighbors(self, aa):
        """Return the lattice sites a bead could have a non-bonded contact with.

        These are the 6 nearest-neighbour sites of ``aa`` (wrapped into the box), minus the sites
        of its bonded neighbours ``aa.next`` and ``aa.previous``. Occupancy is not checked.
        """
        if self.lattice_type == 'simple_cubic':
            neighbors = []
            for i in range(3):
                neighbor = aa.coordinates.copy()
                neighbor[i] += 1
                neighbor_rev = aa.coordinates.copy()
                neighbor_rev[i] += -1
                neighbors.append(neighbor)
                neighbors.append(neighbor_rev)
            neighbors = [self.periodic_coordinate(neighbor) for neighbor in neighbors]
            neighbors = [neighbor for neighbor in neighbors if str(neighbor) not in [str(aa.next),str(aa.previous)]]
            return neighbors
    
    def update_system_energy(self, nearest_neighbors=True):
        """Recount the energy and all contact counts from scratch and store them on the lattice.

        Every non-bonded nearest-neighbour pair is counted once. Sets ``energy``, ``native_contacts``,
        ``non_covalent_hydrophobic_contacts`` and the per-chain lists ``individual_energy``,
        ``individual_ncs`` and ``individual_nchcs`` (indexed by polymer id).

        Parameters
        ----------
        nearest_neighbors : bool, optional
            If False, no contacts are counted and everything is set to zero.
        """
        original_connections=nx.Graph()
        energy = 0
        hydrophobic_contacts = 0
        individual_ncs = [0]*len(self.last)
        individual_nchcs = [0]*len(self.last)
        individual_energies = [0]*len(self.last)

        for last in list(self.last.keys()):
            polymer = self.last[last]
            indiv_ncs = 0
            indiv_nchcs = 0
            indiv_energy = 0
            aa = self.space[str(last)]
            while True:
                if nearest_neighbors:
                    neighbors = self.find_valid_neighbors(aa)
                    neighbors = [neighbor for neighbor in neighbors if (str(aa.coordinates), str(neighbor)) not in original_connections.edges and str(neighbor) in self.space]
                    for neighbor in neighbors:
                        energy += round(self.bond_energies[aa.polarity][self.space[str(neighbor)].polarity],2)
                        original_connections.add_edge(str(aa.coordinates), str(neighbor))
                        if aa.polarity == -1 and self.space[str(neighbor)].polarity == -1:
                            hydrophobic_contacts += 1
                            if aa.polymer == self.space[str(neighbor)].polymer:
                                indiv_nchcs += 1
                        if aa.polymer == self.space[str(neighbor)].polymer:
                            indiv_ncs += 1
                            indiv_energy += round(self.bond_energies[aa.polarity][self.space[str(neighbor)].polarity],2)
                if aa.previous is not None:
                    aa = self.space[str(aa.previous)]
                else:
                    break
            individual_ncs[polymer] = indiv_ncs
            individual_nchcs[polymer] = indiv_nchcs
            individual_energies[polymer] = round(indiv_energy,2)
        
        self.individual_ncs = individual_ncs
        self.individual_nchcs = individual_nchcs
        self.individual_energy = individual_energies
        self.native_contacts = len(original_connections.edges)
        self.non_covalent_hydrophobic_contacts = hydrophobic_contacts
        self.energy = round(energy,2)

    def measure_system_energy(self, new_space, new_starts, new_lasts, nearest_neighbors=True):
        """Count the energy and contacts of a proposed configuration without changing the lattice.

        Same counting as ``update_system_energy``, applied to ``new_space`` with chain ends
        ``new_starts`` and ``new_lasts``. Used by ``cluster_move``.

        Returns
        -------
        tuple
            (individual_ncs, individual_nchcs, individual_energies, native_contacts,
            hydrophobic_contacts, energy); the per-chain lists are indexed by polymer id.
        """
        original_connections=nx.Graph()
        energy = 0
        hydrophobic_contacts = 0
        individual_ncs = [0]*len(new_lasts)
        individual_nchcs = [0]*len(new_lasts)
        individual_energies = [0]*len(new_lasts)

        for last in list(new_lasts.keys()):
            polymer = new_lasts[last]
            indiv_ncs = 0
            indiv_nchcs = 0
            indiv_energy = 0
            aa = new_space[str(last)]
            while True:
                if nearest_neighbors:
                    neighbors = self.find_valid_neighbors(aa)
                    neighbors = [neighbor for neighbor in neighbors if (str(aa.coordinates), str(neighbor)) not in original_connections.edges and str(neighbor) in new_space]
                    for neighbor in neighbors:
                        energy += round(self.bond_energies[aa.polarity][new_space[str(neighbor)].polarity],2)
                        original_connections.add_edge(str(aa.coordinates), str(neighbor))
                        if aa.polarity == -1 and new_space[str(neighbor)].polarity == -1:
                            hydrophobic_contacts += 1
                            if aa.polymer == new_space[str(neighbor)].polymer:
                                indiv_nchcs += 1
                        if aa.polymer == new_space[str(neighbor)].polymer:
                            indiv_ncs += 1
                            indiv_energy += round(self.bond_energies[aa.polarity][new_space[str(neighbor)].polarity],2)
                if aa.previous != None:
                    aa = new_space[str(aa.previous)]
                else:
                    break
            individual_ncs[polymer] = indiv_ncs
            individual_nchcs[polymer] = indiv_nchcs
            individual_energies[polymer] = round(indiv_energy,2)

        return individual_ncs, individual_nchcs, individual_energies, len(original_connections.edges), hydrophobic_contacts, round(energy,2)

    def find_subsystem_energy(self, originals, replacements, start_or_end):
        """Compute the change in energy and contacts if the beads at ``originals`` move to ``replacements``.

        Only contacts involving the moved beads are counted, before and after the move. In the new
        configuration, beads adjacent in ``replacements`` are bonded and so excluded, and sites
        vacated by the move are treated as empty.

        Parameters
        ----------
        originals : list of str
            Current sites of the moving beads, consecutive along one chain.
        replacements : list of str
            New sites, aligned with ``originals``.
        start_or_end : int
            Not used here; kept so the call matches ``move_chain``.

        Returns
        -------
        tuple
            (delta_energy, delta_contacts, delta_NCHC, delta_individual_ncs,
            delta_individual_nchcs, delta_individual_energy). The individual deltas refer to the
            chain of the moved beads.
        """
        original_energy = 0
        original_indiv_energy = 0
        non_existing = [i for i in originals if i not in replacements]
        original_connections=nx.Graph()
        original_NCHC = 0
        original_indiv_ncs = 0
        original_indiv_nchcs = 0

        for coord in originals:
            aa = self.space[coord]
            neighbors = self.find_valid_neighbors(aa)
            neighbors = [neighbor for neighbor in neighbors if str(neighbor) in self.space and (coord, str(neighbor)) not in original_connections.edges ]
            for neighbor in neighbors:
                original_energy += self.bond_energies[self.space[str(aa.coordinates)].polarity][self.space[str(neighbor)].polarity]
                original_connections.add_edge(str(aa.coordinates),str(neighbor))
                if aa.polymer == self.space[str(neighbor)].polymer:
                    original_indiv_ncs += 1
                    original_indiv_energy += self.bond_energies[self.space[str(aa.coordinates)].polarity][self.space[str(neighbor)].polarity]
                if self.space[str(aa.coordinates)].polarity == -1 and self.space[str(neighbor)].polarity == -1:
                    original_NCHC += 1
                    if aa.polymer == self.space[str(neighbor)].polymer:
                        original_indiv_nchcs += 1
        
        new_connections=nx.Graph()
        new_NCHC = 0
        new_energy = 0
        new_indiv_energy = 0
        new_indiv_ncs = 0
        new_indiv_nchcs = 0
        
        i = 0
        for (original, aa) in zip(originals, replacements):
            original_aa = copy.copy(self.space[original])
            original_aa.coordinates = [int(i) for i in aa.strip('][').split(', ')].copy()
            neighbors = self.find_valid_neighbors(original_aa)
            exclude = copy.copy(non_existing)
            
            if len(originals) > 1:
                if i == 0:
                    exclude.append(replacements[i+1])
                elif i == len(originals) - 1:
                    exclude.append(replacements[i-1])
                else:
                    exclude.append(replacements[i+1])
                    exclude.append(replacements[i-1])
            
            neighbors = [neighbor for neighbor in neighbors if str(neighbor) not in exclude and (aa, str(neighbor)) not in new_connections.edges]
            for neighbor in neighbors:
                if str(neighbor) in replacements:
                    index = replacements.index(str(neighbor))
                    new_aa = self.space[originals[index]]
                    new_polarity = new_aa.polarity
                    new_energy += self.bond_energies[self.space[original].polarity][new_polarity]
                    new_connections.add_edge(aa, str(neighbor))
                    
                    if original_aa.polymer == new_aa.polymer:
                        new_indiv_ncs += 1
                        new_indiv_energy += self.bond_energies[self.space[original].polarity][new_polarity]
                    
                    if self.space[original].polarity == -1 and new_polarity == -1:
                        new_NCHC += 1
                        if original_aa.polymer == new_aa.polymer:
                            new_indiv_nchcs += 1
                    new_connections.add_edge(aa,str(neighbor))
                    
                elif str(neighbor) in self.space:
                    new_energy += self.bond_energies[self.space[original].polarity][self.space[str(neighbor)].polarity]
                    new_connections.add_edge(aa,str(neighbor))
                    
                    if original_aa.polymer == self.space[str(neighbor)].polymer:
                        new_indiv_ncs += 1
                        new_indiv_energy += self.bond_energies[self.space[original].polarity][self.space[str(neighbor)].polarity]
                    
                    if self.space[original].polarity == -1 and self.space[str(neighbor)].polarity == -1:
                        new_NCHC += 1
                        if original_aa.polymer == self.space[str(neighbor)].polymer:
                            new_indiv_nchcs += 1
            i += 1
        
        return round(new_energy - original_energy,2), len(new_connections.edges) - len(original_connections.edges), new_NCHC - original_NCHC, new_indiv_ncs - original_indiv_ncs, new_indiv_nchcs - original_indiv_nchcs, round(new_indiv_energy - original_indiv_energy,2)

    def validate_chain(self):

        """Check that every chain is intact and the lattice is consistent; raise RuntimeError if not.

        Checks that polymer ids 0 to ``n_polymers - 1`` all appear in ``start`` and ``last``, that every
        coordinate is inside the box, that each chain has ``length_of_polymer`` beads linked through
        ``next``, that no chain visits a site twice, that every bond has length 1 (or ``bound`` across
        the periodic boundary), and that the number of occupied sites equals
        ``n_polymers * length_of_polymer``.

        Returns
        -------
        bool
            True if all checks pass.
        """
        for i in range(self.n_polymers):
            if i not in list(self.start.values()):
                raise RuntimeError('Polymer numbers not distinct in start')
            if i not in list(self.last.values()):
                raise RuntimeError('Polymer numbers not distinct in last')
        
        for start in list(self.start.keys()):
            i = 0
            records = []
            start = self.space[start]
            while i < self.length_of_polymer:
                if max(start.coordinates) > self.bound or min(start.coordinates) < 0:
                    raise RuntimeError('Periodic Boundary Conditions have been violated.')
                records.append(str(start.coordinates))
                i+= 1
                next_coordinate = copy.copy(start.next)
                if next_coordinate is None:
                    if i==self.length_of_polymer:
                        break
                    else:
                        raise RuntimeError('The chain has broken. The polymer does not have the specified number of residues.')
                if str(next_coordinate) in records:
                    raise RuntimeError('The chain has two residues on the same site')

                if math.dist(start.coordinates, start.next) not in [1, self.bound]:
                    raise RuntimeError('The distance between neighbors is not valid. The chain has broken.')

                start = self.space[str(next_coordinate)]
                
        if len(self.space) != self.n_polymers*self.length_of_polymer:
            raise RuntimeError('There is an overlap between amino acids of different polymers.')
        return True

    def move_success(self, deltaE):
        """Metropolis acceptance test.

        Returns True with probability min(1, exp(-beta * deltaE)).
        """
        if system_random.random() < np.exp(-deltaE*self.beta):
            return True
        else:
            return False
    
    def move_chain(self, originals, replacements, inflection_point, start_or_end):
        """Propose moving a run of beads; apply the move if the Metropolis test accepts it.

        Computes the change in energy with ``find_subsystem_energy``. If accepted, updates the running
        energy and contact totals, the per-chain stats of the moved chain, the bead coordinates and
        ``next``/``previous`` links, the ``start``/``last`` keys and ``space``. A replacement site may
        coincide with another moved bead's original site.

        Parameters
        ----------
        originals : list of str
            Current sites of the moving beads, consecutive along one chain, starting from the bead
            bonded to ``inflection_point``.
        replacements : list of str
            New sites, aligned with ``originals``.
        inflection_point : list of int or None
            Site of the unmoved bead bonded to ``originals[0]``; its link to the moved segment is
            rewired. None when the whole chain moves.
        start_or_end : int
            0 if the moved beads lie towards the chain's first bead from ``inflection_point``,
            1 if they lie towards its last bead.

        Returns
        -------
        bool
            True if the move was accepted.
        """
        deltaE, deltaNC, deltaNCHC, delta_indiv_ncs, delta_indiv_nchcs, delta_indiv_energy = self.find_subsystem_energy(originals, replacements, start_or_end)

        if self.move_success(deltaE):
            self.energy += round(deltaE,2)
            self.native_contacts += deltaNC
            self.non_covalent_hydrophobic_contacts += deltaNCHC
            polymer_id = self.space[originals[0]].polymer
            self.individual_ncs[polymer_id] += delta_indiv_ncs
            self.individual_nchcs[polymer_id] += delta_indiv_nchcs
            self.individual_energy[polymer_id] += round(delta_indiv_energy,2)
            all_objects = []
            for aa_step in range(len(originals)):
                original = copy.copy(originals[aa_step])
                replacement = copy.copy(replacements[aa_step])
                replacement_int = [int(i) for i in replacement.strip('][').split(', ')].copy()
                original_aa = copy.copy(self.space[original])
                original_aa.coordinates = replacement_int.copy()
                if start_or_end == 0:
                    if aa_step != 0: 
                        original_aa.next = [int(i) for i in replacements[aa_step-1].strip('][').split(', ')].copy()
                    if original_aa.previous is not None and len(replacements) > 1 and aa_step < len(originals)-1:
                        original_aa.previous = [int(i) for i in replacements[aa_step+1].strip('][').split(', ')].copy()
                elif start_or_end == 1:
                    if aa_step != 0:
                        original_aa.previous = [int(i) for i in replacements[aa_step-1].strip('][').split(', ')].copy()
                    if original_aa.next is not None and len(replacements) > 1 and aa_step < len(originals)-1:
                        original_aa.next = [int(i) for i in replacements[aa_step+1].strip('][').split(', ')].copy()
                all_objects.append(copy.copy(original_aa))
                
            for aa_step in range(len(originals)):
                original = originals[aa_step]
                original_int = [int(i) for i in original.strip('][').split(', ')].copy()
                replacement = replacements[aa_step]
                replacement_int = [int(i) for i in replacement.strip('][').split(', ')].copy()
                replacement_object = copy.copy(all_objects[aa_step])
                if aa_step ==0:
                    if original in self.last.keys():
                        del self.last[copy.copy(original)]
                        self.last[replacement] = polymer_id
                    elif original in self.start.keys():
                        del self.start[copy.copy(original)]
                        self.start[replacement] = polymer_id
                    
                    if len(originals) != self.length_of_polymer:
                        if start_or_end == 0:
                            self.space[str(inflection_point)].previous = replacement_int.copy()
                        else:
                            self.space[str(inflection_point)].next = replacement_int.copy()
                if aa_step == len(originals)-1:
                    if original in self.last.keys():
                        del self.last[copy.copy(original)]
                        self.last[replacement] = polymer_id
                    elif original in self.start.keys():
                        del self.start[copy.copy(original)]
                        self.start[replacement] = polymer_id
                    if start_or_end == 0:
                        if replacement_object.previous is not None:
                            self.space[str(replacement_object.previous)].next = replacement_int
                    else:
                        if replacement_object.next is not None:
                            self.space[str(replacement_object.next)].previous = replacement_int
                self.space[replacement] = replacement_object
                if original not in replacements:
                    del self.space[original]
            if self.record_moves:
                self.acceptance_records.append(1)
            return True
        else:
            if self.record_moves:
                self.acceptance_records.append(0)
            return False

    def add_protein(self, sequence=None, placement='straight', n_polymers=1):
        """Place chains given as an amino-acid sequence.

        Letters in GAVCPLIMWFKRH become +1 and letters in STYNQDE become -1; the resulting
        polarities are passed to ``add_polymer``.

        Parameters
        ----------
        sequence : str
            One-letter amino-acid sequence.
        placement : {'straight', 'randomly'}, optional
            Passed to ``add_polymer``.
        n_polymers : int, optional
            Number of copies to place.
        """
        polymer = []

        for aa in list(sequence):
            if aa in 'GAVCPLIMWFKRH':
                polymer.append(+1)
            elif aa in 'STYNQDE':
                polymer.append(-1)
            else:
                raise InputError('Unrecognized amino acid in sequence: {}'.format(aa)) # noqa: F821

        if placement in ['straight','randomly']:
            self.add_polymer(polymer, n_polymers=n_polymers, placement=placement)
        else:
            raise InputError('Unrecognized type of polymer placement: {}'.format(type)) # noqa: F821
            
    def add_polymer(self, polymer, n_polymers=1, placement='straight'):
        """Place ``n_polymers`` copies of a chain with the given polarities.

        Polymer ids continue from the chains already on the lattice, so ``add_polymer`` can be
        called more than once. With 'straight' placement, the first chain starts at a random site and
        each later chain starts 1 or 2 sites, along a random axis, from the start of a randomly chosen
        chain already placed. With 'randomly' placement every chain starts at a random site. A
        placement that fails is retried from a new start. Each placed chain is checked with
        ``validate_chain``; the energy is recounted once at the end. Sets ``length_of_polymer`` and
        ``n_polymers``.

        Parameters
        ----------
        polymer : list of int
            Polarity (+1 or -1) of each bead, from first to last.
        n_polymers : int, optional
            Number of copies to place.
        placement : {'straight', 'randomly'}, optional
            'straight' lays each chain along one random axis and direction (needs
            ``len(polymer) <= bound + 1``, otherwise ValueError). 'randomly' grows a self-avoiding
            random walk, choosing each step uniformly among the free neighbouring sites.
        """
        length = len(polymer)
        if placement == 'straight' and length > self.bound + 1:
            raise ValueError('Straight placement needs length <= bound + 1 ({}), got {}.'.format(self.bound + 1, length))
        polymers_placed = 0

        while polymers_placed < n_polymers:
            polymer_id = len(self.start)
            new_coords = None
            if len(self.start) > 0 and placement=='straight':
                new_coords = system_random.choice(list(self.start.keys()))
                new_coords = [int(i) for i in new_coords.strip('][').split(', ')].copy()
                axis_start = system_random.choice([0,1,2])
                change = system_random.choice([1,2])
                new_coords[axis_start] += system_random.choice([+change, -change])
            else:
                x = system_random.randint(0, self.bound)
                y = system_random.randint(0, self.bound)
                z = system_random.randint(0, self.bound)
                new_coords = [x, y, z]

            new_coords = self.periodic_coordinate(new_coords)
            if str(new_coords) in self.space:
                continue

            all_placements = []
            all_placements.append(new_coords)
            if placement == 'straight':
                axis = system_random.choice([0,1,2])
                direction = system_random.choice([1, -1])
                for i in range(1, length):
                    tries = 0
                    while tries < 5:
                        tries += 1
                        new = new_coords.copy()
                        new[axis] += direction*i
                        new = self.periodic_coordinate(new)
                        if str(new) not in self.space:
                            all_placements.append(new.copy())
                            break
            elif placement == 'randomly':
                coordinates = new_coords.copy()
                for i in range(length-1):
                    neighbors = []
                    for axis in range(3):
                        neighbor = coordinates.copy()
                        neighbor[axis] += 1
                        neighbor_rev = coordinates.copy()
                        neighbor_rev[axis] += -1
                        neighbors.append(self.periodic_coordinate(neighbor))
                        neighbors.append(self.periodic_coordinate(neighbor_rev))
                    neighbors = [neighbor for neighbor in neighbors if neighbor not in all_placements and str(neighbor) not in self.space]
                    if len(neighbors) == 0:
                        break
                    coordinates = system_random.choice(neighbors)
                    all_placements.append(coordinates)

            if len(all_placements) == length:
                self.start[str(new_coords)] = polymer_id
                self.space[str(new_coords)] = amino_acid(polymer[0], new_coords, polymer_id)
                current_aa = self.space[str(new_coords)]
                for i in range(1, length):
                    next_coordinates = all_placements[i]
                    self.space[str(current_aa.coordinates)].next = next_coordinates
                    self.space[str(next_coordinates)] = amino_acid(polymer[i], next_coordinates, polymer_id)
                    self.space[str(next_coordinates)].previous = current_aa.coordinates.copy()
                    current_aa = self.space[str(next_coordinates)]
                polymers_placed += 1
                self.length_of_polymer = length
                self.n_polymers += 1
                self.last[str(next_coordinates)] = polymer_id
                self.validate_chain()

        self.update_system_energy()


    def end_move(self):
        """End move: rotate a chain's end bead by 90 degrees about its bonded neighbour.

        Picks a random chain end and one of the 4 sites around the neighbour perpendicular to the
        current bond. Rejected if that site is occupied; otherwise passed to ``move_chain``.

        Returns
        -------
        bool or None
            True if the move was accepted; a falsy value otherwise.
        """
        if self.record_moves:
            self.move_records.append('end_move')
        start_or_end=system_random.choice([0,1])
        coordinates = system_random.choice([ list(self.start.keys()), list(self.last.keys())][start_or_end])
        coordinates = self.periodic_coordinate([int(i) for i in coordinates.strip('][').split(', ')])
        current_aa = copy.copy(self.space[str(coordinates)])

        next_coordinates = coordinates.copy()
        back = None
        if start_or_end == 0:
            back = current_aa.next.copy()
        else:
            back = current_aa.previous.copy()
        axis = None
        for i in range(3):
            if back[i] !=  next_coordinates[i]:
                next_coordinates[i] = copy.copy(back[i])
                axis = i
        next_axis = system_random.choice([i for i in [0,1,2] if i != axis])
        direction = system_random.choice([1, -1])
        next_coordinates[next_axis] += direction
        next_coordinates = self.periodic_coordinate(next_coordinates)
        if str(next_coordinates) not in self.space:
            return self.move_chain([str(coordinates)], [str(next_coordinates)], back, start_or_end)

        if self.record_moves:
            self.acceptance_records.append(0)

    def crankshaft_move(self):
        """Crankshaft move: rotate a U-shaped pair of beads by 90 degrees.

        Picks a random chain and a random bead i (0 to L-4). If beads i to i+3 form a U (beads i and
        i+3 adjacent), beads i+1 and i+2 rotate by +90 or -90 degrees, chosen at random, about the
        axis through beads i and i+3. U shapes across the periodic boundary are recognised. Rejected
        if beads i to i+3 are not a U or a target site is occupied.

        Returns
        -------
        bool or None
            True if the move was accepted; a falsy value otherwise.
        """
        if self.record_moves:
            self.move_records.append('crankshaft')
        aa = self.space[system_random.choice(list(self.start.keys()))]
        for _ in range(system_random.randint(0, self.length_of_polymer-4)):
            aa = self.space[str(aa.next.copy())]

        axis_1 = None
        axis_2 = None
        axis_3 = None
        coordinates = aa.coordinates.copy()
        first_coordinates = aa.next.copy()
        second_coordinates = self.space[str(first_coordinates)].next.copy()
        third_coordinates = self.space[str(second_coordinates)].next.copy()
        first_direction = None
        back_direction = None

        for i in range(3):
            if coordinates[i] != first_coordinates[i]:
                first_direction = first_coordinates[i] - coordinates[i]
                axis_1 = i
        for j in range(3):
            if second_coordinates[j] != third_coordinates[j]:
                back_direction = third_coordinates[j] - second_coordinates[j]
                axis_3 = j
        for k in range(3):
            if first_coordinates[k] != second_coordinates[k]:
                axis_2 = k

        if axis_1 == axis_3 and axis_1 != axis_2 and (first_direction + back_direction) % (self.bound + 1) == 0:
            axis_of_replacement = [i for i in [0,1,2] if i not in [axis_1, axis_2]][0]
            direction = system_random.choice([1, -1])
            first_try = first_coordinates.copy()
            second_try = second_coordinates.copy()
            first_try[axis_1] += back_direction
            second_try[axis_1] += back_direction
            first_try[axis_of_replacement] += direction
            second_try[axis_of_replacement] += direction
            first_try = self.periodic_coordinate(first_try)
            second_try = self.periodic_coordinate(second_try)

            if str(first_try) not in self.space and str(second_try) not in self.space:
                return self.move_chain([str(second_coordinates), str(first_coordinates)], [str(second_try), str(first_try)], third_coordinates, 0)

        if self.record_moves:
            self.acceptance_records.append(0)

    def corner_move(self):
        """Corner move: move a chain's last two beads together.

        Picks a random chain end. The second bead jumps to one of the 4 sites next to the third bead
        perpendicular to their bond, chosen at random, and the end bead keeps its bond direction
        relative to the second bead. Rejected if either target site is occupied.

        Returns
        -------
        bool or None
            True if the move was accepted; a falsy value otherwise.
        """
        if self.record_moves:
            self.move_records.append('corner_move')
        polymer_id=system_random.choice(range(len(self.start.keys())))
        start_or_end=system_random.choice([0,1])
        coordinates = [ list(self.start.keys()), list(self.last.keys())][start_or_end][polymer_id]
        coordinates = [int(i) for i in coordinates.strip('][').split(', ')]
        current_aa = self.space[str(coordinates)]
        if start_or_end == 0:
            next_coordinates = current_aa.next.copy()
        else:
            next_coordinates = current_aa.previous.copy()

        second_rep = next_coordinates.copy()
        if start_or_end == 0:
            third = self.space[str(second_rep)].next.copy()
        else:
            third = self.space[str(second_rep)].previous.copy()
        first_direction = 0
        axis = None
        for i in range(3):
            if second_rep[i] !=  third[i]:
                first_direction = third[i] - second_rep[i]
                second_rep[i] = third[i]
                axis = i
        next_axis = system_random.choice([i for i in [0,1,2] if i != axis])
        next_direction = system_random.choice([1, -1])
        second_rep[next_axis] += next_direction
        first_rep = coordinates.copy()
        first_rep[axis] += first_direction
        first_rep[next_axis] += next_direction

        first_rep = self.periodic_coordinate(first_rep)
        second_rep = self.periodic_coordinate(second_rep)

        if str(first_rep) not in self.space and str(second_rep) not in self.space:
            originals = [str(next_coordinates), str(coordinates)]
            replacements = [str(second_rep), str(first_rep)]
            return self.move_chain(originals, replacements, third, start_or_end)

        if self.record_moves:
            self.acceptance_records.append(0)
            return False

    def corner_move_anywhere(self):
        """Rigidly translate the segment between a random pivot bead and a chain end.

        Picks a random chain end and a pivot bead 1 to L-2 beads in from it (shorter distances are
        more likely). The whole segment between the pivot and the end is translated so that the bead
        next to the pivot moves to one of the 4 sites at 90 degrees to their old bond, chosen at
        random. Rejected if any target site is occupied, including the segment's own old sites.

        Returns
        -------
        bool
            True if the move was accepted.
        """
        if self.record_moves:
            self.move_records.append('corner_anywhere')
        polymer_id=system_random.choice(range(len(self.start.keys())))
        start_or_end=system_random.choice([0,1])
        coordinates = [ list(self.start.keys()), list(self.last.keys())][start_or_end][polymer_id]
        coordinates = [int(i) for i in coordinates.strip('][').split(', ')]
        steps = system_random.randint(2, self.length_of_polymer-2)
        all_selections = []

        if start_or_end == 0:
            for step in range(steps):
                coordinates = self.space[str(coordinates)].next.copy()
                all_selections.append(coordinates)
        else:
            for step in range(steps):
                coordinates = self.space[str(coordinates)].previous.copy()
                all_selections.append(coordinates)

        inflection_point = system_random.choice(all_selections)
        if start_or_end == 0:
            original = self.space[str(inflection_point)].previous.copy()
        else:
            original = self.space[str(inflection_point)].next.copy()

        first_axis = None
        first_direction = 0
        for i in range(3):
            if original[i] != inflection_point[i]:
                first_direction = inflection_point[i] - original[i]
                first_axis = i
        next_axis = system_random.choice([i for i in [0,1,2] if i != first_axis])
        next_direction = system_random.choice([1, -1])

        to_be_replaced = []
        replacements = []
        while original is not None:
            new_position = original.copy()
            new_position[first_axis] += first_direction
            new_position[next_axis] += next_direction
            new_position = self.periodic_coordinate(new_position)
            if str(new_position) in self.space:
                if self.record_moves:
                    self.acceptance_records.append(0)
                return False
            to_be_replaced.append(str(original))
            replacements.append(str(new_position))
            if start_or_end == 0:
                original = copy.copy(self.space[str(original)].previous)
            else:
                original = copy.copy(self.space[str(original)].next)

        return self.move_chain(to_be_replaced, replacements, inflection_point, start_or_end)

    def corner_flip(self):
        """Corner flip (kink jump): move one interior bead across the corner it sits on.

        Picks a random chain and a random interior bead (1 to L-2). If the bead's two bonds are
        perpendicular, it moves to the opposite corner of the square they span. Rejected if the bead
        is not at a corner or the target site is occupied.

        Returns
        -------
        bool or None
            True if the move was accepted; a falsy value otherwise.
        """
        if self.record_moves:
            self.move_records.append('corner_flip')
        start = system_random.choice(list(self.start.keys()))
        center_coordinates = self.space[start].next.copy()
        for _ in range(system_random.randint(0, self.length_of_polymer-3)):
            center_coordinates = self.space[str(center_coordinates)].next.copy()

        before = self.space[str(center_coordinates)].previous
        after = self.space[str(center_coordinates)].next

        before_axis = None
        first_direction = 0
        for i in range(3):
            if before[i] != center_coordinates[i]:
                first_direction = before[i] - center_coordinates[i]
                before_axis = i

        next_axis = None
        next_direction = 0
        for j in range(3):
            if after[j] != center_coordinates[j]:
                next_direction = after[j] - center_coordinates[j]
                next_axis = j

        if before_axis != next_axis:
            replacement = center_coordinates.copy()
            replacement[before_axis] += first_direction
            replacement[next_axis] += next_direction
            replacement = self.periodic_coordinate(replacement)
            if str(replacement) not in self.space:
                return self.move_chain([str(center_coordinates)], [str(replacement)], after, 0)

        if self.record_moves:
            self.acceptance_records.append(0)
        
    def reptation_move(self):
        """Reptation (slithering snake): move a chain forward by one bead along its own contour.

        Picks a random chain end and one of the 6 directions. If that site is free, the end bead moves
        into it and every other bead moves to the site of the bead ahead of it, so the far end gives
        up its site. Rejected if the target site is occupied.

        Returns
        -------
        bool
            True if the move was accepted.
        """
        if self.record_moves:
            self.move_records.append('reptation')
        polymer_id=system_random.choice(range(len(self.start.keys())))
        start_or_end=system_random.choice([0,1])
        coordinates = [ list(self.start.keys()), list(self.last.keys())][start_or_end][polymer_id]
        originals = []
        replacements = []
        originals.append(copy.copy(coordinates))
        coordinates = [int(i) for i in coordinates.strip('][').split(', ')]

        main_axis = system_random.choice([0, 1, 2])
        main_direction = system_random.choice([1, -1])
        rep = coordinates.copy()
        rep[main_axis] += main_direction
        rep = self.periodic_coordinate(rep)
        if str(rep) not in self.space:
            replacements.append(str(rep))
            main = coordinates.copy()
            while True:
                if start_or_end == 0:
                    main = self.space[str(main)].next
                else:
                    main = self.space[str(main)].previous
                if main is None:
                    break
                replacements.append(originals[-1])
                originals.append(str(main))

            originals.reverse()
            replacements.reverse()
            return self.move_chain(originals, replacements, None, start_or_end)

        if self.record_moves:
            self.acceptance_records.append(0)
        return False
    
    def rotation_move(self):
        """Rotate a whole chain by 90 degrees about its middle bead.

        Picks a random chain, one of the 3 coordinate planes and +90 or -90 degrees. The pivot is the
        bead at index L // 2, which stays in place, so the opposite rotation exactly undoes the move.
        Beads may land on sites the chain itself occupied; the move is rejected if any bead lands on
        another chain.

        Returns
        -------
        bool
            True if the move was accepted.
        """
        if self.record_moves:
            self.move_records.append('rotation_move')
        
        def rotate(center, point, angle):
            """Rotate the 2-D point ``point`` about ``center`` by ``angle`` radians; return rounded integers."""
            c1, c2 = center
            p1, p2 = point
            nx = c1 + math.cos(angle) * (p1 - c1) - math.sin(angle) * (p2 - c2)
            ny = c2 + math.sin(angle) * (p1 - c1) + math.cos(angle) * (p2 - c2)
            return round(nx), round(ny)
    
        start_point = system_random.choice(list(self.start.keys()))
        all_records = []
        point = [int(i) for i in start_point.strip('][').split(', ')]
        while point is not None:
            aa = self.space[str(point)]
            all_records.append(aa.coordinates)
            point = copy.copy(aa.next)
        center = all_records[len(all_records)//2]
        first_axis = system_random.choice([0,1,2])
        second_axis = system_random.choice([i for i in [0,1,2] if i != first_axis])
        angle = math.radians(system_random.choice([-90, 90]))
        all_replacements = []
        for coordinates in all_records:
            rep = coordinates.copy()
            n1, n2 = rotate([center[first_axis], center[second_axis]], [rep[first_axis], rep[second_axis]], angle)
            rep[first_axis] = n1
            rep[second_axis] = n2
            rep = self.periodic_coordinate(rep)
            if str(rep) not in self.space or rep in all_records:
                all_replacements.append(str(rep))
        all_records = [str(i) for i in all_records]
        if len(all_records) == len(all_replacements):
            return self.move_chain(all_records, all_replacements, None, 1)
        else:
            if self.record_moves:
                self.acceptance_records.append(0)
            return False

    def transform_move(self):
        """Translate a whole chain by one site.

        Picks a random chain and one of the 6 directions. Beads may land on sites the chain itself
        occupied; the move is rejected if any bead would land on another chain.

        Returns
        -------
        bool
            True if the move was accepted.
        """
        if self.record_moves:
            self.move_records.append('transform_move')
        point = system_random.choice(list(self.start.keys()))
        point = [int(i) for i in point.strip('][').split(', ')]
        axis = system_random.choice([0, 1, 2])
        direction = system_random.choice([1, -1])

        originals = []
        while point is not None:
            originals.append(str(point))
            point = self.space[str(point)].next

        replacements = []
        for original in originals:
            rep = [int(i) for i in original.strip('][').split(', ')]
            rep[axis] += direction
            rep = self.periodic_coordinate(rep)
            if str(rep) in self.space and str(rep) not in originals:
                if self.record_moves:
                    self.acceptance_records.append(0)
                return False
            replacements.append(str(rep))

        return self.move_chain(originals, replacements, None, 1)

    def complexes(self):
        """Build the graph of which chains touch which.

        Nodes are the polymer ids of all chains, including chains with no contacts. Two chains are
        joined by an edge if any of their beads form a non-bonded contact; contacts within a chain
        give a self-loop.

        Returns
        -------
        networkx.Graph
        """
        complexes = nx.Graph()
        complexes.add_nodes_from(self.start.values())
        for first in list(self.start.keys()):
            next = copy.copy(first)
            while next is not None:
                aa = self.space[str(next)]
                neighbors = self.find_valid_neighbors(aa)
                neighbors = [neighbor for neighbor in neighbors if str(neighbor) in self.space]
                for neighbor in neighbors:
                    if (aa.polymer, self.space[str(neighbor)].polymer) not in complexes.edges:
                        complexes.add_edge(aa.polymer, self.space[str(neighbor)].polymer)
                next = aa.next
        return complexes

    def calculate_number_of_clusters(self):
        """Return the number of clusters (connected components of ``complexes()``), counting isolated chains."""
        clusters = []
        for cluster in nx.connected_components(self.complexes()):
            clusters.append(list(cluster))
        return len(clusters)

    def collect_networks(self):
        """Return the clusters as a list of lists of polymer ids (connected components of ``complexes()``)."""
        clusters = []
        for cluster in nx.connected_components(self.complexes()):
            clusters.append(list(cluster))
        return clusters

    def cluster_move(self):
        """Translate a whole cluster of touching chains by one site.

        Clusters are the connected components of ``complexes()``. Picks a cluster uniformly and one of
        the 6 directions, builds the shifted configuration, recounts its energy with
        ``measure_system_energy`` and applies the Metropolis test. Rejected if there is at most one
        cluster. Shifts that bring the cluster into contact with another cluster are allowed (the
        clusters merge); this move never proposes the reverse split.

        Returns
        -------
        bool
            True if the move was accepted.
        """
        if self.record_moves:
            self.move_records.append('cluster move')
        clusters = []
        for cluster in nx.connected_components(self.complexes()):
            clusters.append(list(cluster))

        if len(clusters) <= 1:
            if self.record_moves:
                self.acceptance_records.append(0)
            return False

        cluster_choice = system_random.choice(clusters)
        direction = system_random.choice([0,1,2])
        step = system_random.choice([1, -1])
        new_space = dict()
        new_starts = copy.deepcopy(self.start)
        new_lasts = copy.deepcopy(self.last)
    
        for polymer in range(self.n_polymers):
            main = list(self.start.keys())[list(self.start.values()).index(polymer)]
            if polymer not in cluster_choice:
                for i in range(self.length_of_polymer):
                    object = copy.deepcopy(self.space[main])
                    new_space[main] = object
                    main = str(copy.copy(object.next))
            else:
                for i in range(self.length_of_polymer):
                    object = copy.deepcopy(self.space[main])
                    original_coordinates = copy.copy(object.coordinates)
                    original_next = copy.copy(object.next)
                    object.coordinates[direction] += step
                    object.coordinates = self.periodic_coordinate(object.coordinates)
                    if i == 0:
                        if polymer in list(new_starts.values()):
                            del new_starts[str(original_coordinates)]
                        new_starts[str(object.coordinates)] = polymer
                    else:
                        object.previous[direction] += step
                        object.previous = self.periodic_coordinate(object.previous)
                    
                    if i == self.length_of_polymer - 1:
                        if polymer in list(new_lasts.values()):
                            del new_lasts[str(original_coordinates)]
                        new_lasts[str(object.coordinates)] = polymer
                    else:
                        object.next[direction] += step
                        object.next = self.periodic_coordinate(object.next)
                    new_space[str(object.coordinates)] = object
                    main = str(original_next)

        if len(new_starts) != len(new_lasts):
            print(direction, step)
            print(cluster_choice)
            print(new_starts)
            print(new_lasts)
            print(new_space)
            sys.exit('problem problem')

        individual_ncs, individual_nchcs, individual_energies, native_contacts, hydrophobic_contacts, energy = self.measure_system_energy(new_space, new_starts, new_lasts)

        deltaE = energy - self.energy
        if self.move_success(deltaE):
            if self.record_moves:
                self.acceptance_records.append(1)
            self.space = new_space
            self.start = new_starts
            self.last = new_lasts
            self.individual_ncs = individual_ncs
            self.individual_nchcs = individual_nchcs
            self.individual_energy = individual_energies
            self.native_contacts = native_contacts
            self.non_covalent_hydrophobic_contacts = hydrophobic_contacts
            self.energy = energy
            return True
        else:
            if self.record_moves:
                self.acceptance_records.append(0)
            return False
        
    def simulate(self, n_mcmc=10000, interval=100, record_intervals=False, anneal=True, beta_lower_bound=0, beta_upper_bound=1, beta_interval=0.05):
        """Run Metropolis Monte Carlo for ``n_mcmc`` steps.

        Each step applies one move chosen uniformly from: end, corner, corner-anywhere, corner flip,
        crankshaft, reptation and rotation, plus translation (and cluster moves, if
        ``allow_cluster_move``) when there is more than one chain.

        Parameters
        ----------
        n_mcmc : int, optional
            Number of steps.
        interval : int, optional
            Steps between recordings when ``record_intervals`` is True.
        record_intervals : bool, optional
            If True, every ``interval`` steps: append the energy, beta, contact counts, per-chain
            stats, number of clusters (``cluster_records``) and the clusters themselves
            (``cluster_size_records``) to the record lists, update the plateau check (variance of the last 29
            recorded energies below 2 sets ``plateaued`` and ``plateau_time``), and run
            ``validate_chain``.
        anneal : bool, optional
            If True, beta increases by ``beta_interval`` every
            ``n_mcmc * beta_interval / (beta_upper_bound - beta_lower_bound)`` steps. Raises ValueError
            if the two bounds are equal.
        beta_lower_bound : float, optional
            Beta at the start of the run. Used whether or not ``anneal`` is True.
        beta_upper_bound : float, optional
            Beta reached at the end of an annealed run.
        beta_interval : float, optional
            Size of each beta increment when annealing.
        """
        if anneal:
            if beta_upper_bound == beta_lower_bound:
                raise ValueError('beta_lower_bound and beta_upper_bound must be different when anneal=True.')
            substep = round(n_mcmc*beta_interval/(beta_upper_bound - beta_lower_bound), 0)
        else:
            substep = None
        self.beta = beta_lower_bound
        self.energy_records.append(copy.copy(round(self.energy,2)))
        self.beta_records.append(copy.copy(self.beta))
        self.native_contacts_records.append(copy.copy(self.native_contacts))
        self.nchc_records.append(copy.copy(self.non_covalent_hydrophobic_contacts))
        self.individual_ncs_records.append(copy.copy(self.individual_ncs))
        self.individual_nchcs_records.append(copy.copy(self.individual_nchcs))
        self.cluster_records.append(self.calculate_number_of_clusters())
        self.cluster_size_records.append(self.collect_networks())

        for step in range(1, n_mcmc+1):
            if anneal:
                if step%substep == 0:
                    self.beta += beta_interval
                    self.beta = round(self.beta, 2)
            all_functions = [self.end_move, self.corner_move, self.corner_move_anywhere, self.corner_flip, self.crankshaft_move, self.reptation_move, 
                            self.rotation_move]
            if self.n_polymers > 1:
                all_functions = all_functions + [self.transform_move]
                if self.allow_cluster_move:
                    all_functions = all_functions + [self.cluster_move]
            system_random.choice(all_functions)()
            self.n_mcmc += 1
            if record_intervals and step%interval == 0 and step > 0:
                self.energy_records.append(copy.copy(round(self.energy,2)))
                self.beta_records.append(copy.copy(self.beta))
                self.native_contacts_records.append(copy.copy(self.native_contacts))
                self.nchc_records.append(copy.copy(self.non_covalent_hydrophobic_contacts))
                self.individual_energy_records.append(copy.copy(self.individual_energy))
                self.individual_ncs_records.append(copy.copy(self.individual_ncs))
                self.individual_nchcs_records.append(copy.copy(self.individual_nchcs))
                self.cluster_records.append(self.calculate_number_of_clusters())
                self.cluster_size_records.append(self.collect_networks())
                if not self.plateaued and len(self.energy_records) > 30 and np.var(self.energy_records[-30:-1]) < 2:
                    self.plateau_time = copy.copy(self.n_mcmc)
                    self.plateaued = True
                self.validate_chain()
            if step%(n_mcmc/10) == 0:
                print(datetime.now())
                print('Completion: {}%'.format(step*100/n_mcmc))
                sys.stdout.flush()
            
            if self.record_moves and len(self.acceptance_records) != self.n_mcmc:
                raise ValueError('Recorded acceptance records are not equal to the number of MCMC steps')
        print('Fully Completed and Validated.')

    def energy_variation_graph(self, savefigure=False, figure_name='LatticePy_energy.png'):
        """Plot the recorded total energy against MC steps.

        Needs ``simulate(record_intervals=True)``. If ``savefigure`` is True the figure is saved to
        ``figure_name`` as a 300-dpi PNG; otherwise it is shown.
        """
        n_total = len(self.energy_records)
        if n_total == 0:
            raise RuntimeError('No energy values have been recorded so far. Please rerun the simulation with record_intervals = True in the simulate() function.')
        
        interval = self.n_mcmc/n_total
        x = [interval*i for i in range(n_total)]
        plt.figure(figsize=(12,6))
        plt.plot(x, self.energy_records, linestyle='solid', linewidth=2, markersize=0)
        plt.title('Variation of system energy over all MCMC steps')
        plt.xlabel('Number of MCMC steps')
        plt.ylabel('System Energy')
        if savefigure:
            plt.savefig(figure_name, dpi=300, format='png')
        else:
            plt.show()
    
    def native_contacts_over_time(self, savefigure=False, figure_name='LatticePy_native_contacts.png'):
        """Plot the recorded number of non-bonded contacts (``native_contacts``) against MC steps.

        If ``savefigure`` is True the figure is saved to ``figure_name`` as a 300-dpi PNG; otherwise
        it is shown.
        """
        n_total = len(self.native_contacts_records)
        interval = self.n_mcmc/n_total
        x = [interval*i for i in range(n_total)]
        plt.figure(figsize=(12,6))
        plt.plot(x, self.native_contacts_records, linestyle='solid', linewidth=2, markersize=0)
        plt.title('Variation of Native Contacts over all MCMC steps')
        plt.xlabel('Number of MCMC steps')
        plt.ylabel('Number of Native Contacts')
        if savefigure:
            plt.savefig(figure_name, dpi=300, format='png')
        else:
            plt.show()
    
    def nchc_over_time(self, savefigure=False, figure_name='LatticePy_nchc.png'):
        """Plot the recorded number of (-1, -1) contacts (NCHC) against MC steps.

        If ``savefigure`` is True the figure is saved to ``figure_name`` as a 300-dpi PNG; otherwise
        it is shown.
        """
        n_total = len(self.nchc_records)
        interval = self.n_mcmc/n_total
        x = [interval*i for i in range(n_total)]
        plt.figure(figsize=(12,6))
        plt.plot(x, self.nchc_records, linestyle='solid', linewidth=2, markersize=0)
        plt.title('Variation of NCHC Contacts over all MCMC steps')
        plt.xlabel('Number of MCMC steps')
        plt.ylabel('Number of Native Contacts')
        if savefigure:
            plt.savefig(figure_name, dpi=300, format='png')
        else:
            plt.show()

    def animate(self):
        """Show a 3-D plotly animation of recorded configurations.

        Expects ``self.records`` to hold one DataFrame per frame (as returned by
        ``visualize(simulating=True)``), with matching ``energy_records`` and ``beta_records``.
        ``simulate`` does not fill ``self.records``.
        """
        def make_figure(i):
            """Return the 3-D scatter trace for frame ``i`` of ``self.records``."""
            df = self.records[i]
            return go.Scatter3d(x=df.x, y=df.y, z=df.z, mode='markers+lines',
                    marker = dict(size = 7, color = df.c, opacity = 0.8 ), showlegend=False )

        all_data = [[make_figure(i)] for i in range(len(self.records))]
        frames = []

        for i in range(len(all_data)):
            df = self.records[i]
            x_min = min(df.x)-2
            y_min = min(df.y)-2
            z_max = max(df.z)+2
            energy = self.energy_records[i]
            beta = self.beta_records[i]
            annotations = [dict(x=x_min, y=y_min, z=z_max, text='Energy: {}'.format(energy), showarrow=False),  dict(x=x_min, y=y_min, z=z_max+2, text='Current Beta: {}'.format(beta), showarrow=False) ]  
            layout = go.layout(scene=dict(annotations=annotations))
            frame = go.Frame(data=all_data[i], layout=layout)
            frames.append(frame)

        first_frame = self.records[0]

        fig = go.Figure(data=[go.Scatter3d(x=first_frame.x, y=first_frame.y, z=first_frame.z, 
                            mode='markers+lines', marker = dict(size = 7, color = first_frame.c, opacity = 0.8 ), showlegend=False )],

            layout=go.Layout(updatemenus=[dict(type="buttons", buttons=[dict(label="Play", method="animate", args=[None])])]
            ), frames=frames)

        fig.show()
        
    def visualize(self, simulating=False, to_html=False, html_file_name='LatticePy_figure.html'):
        """Show all chains in 3-D, with beads coloured by polarity (+1 red, -1 blue).

        The plot is annotated with the energy, number of contacts, number of (-1, -1) contacts and
        beta.

        Parameters
        ----------
        simulating : bool, optional
            If True, return a DataFrame (x, y, z, colour) for the first chain in ``last`` instead of
            plotting.
        to_html : bool, optional
            If True, write the figure to ``html_file_name`` instead of showing it.
        html_file_name : str, optional
            Output file for ``to_html``.
        """
        data = []
        x_min = None
        y_min = None
        z_max = None
        for end_coordinates in list(self.last.keys()):
            x = []
            y = []
            z = [] 
            polarities = []
            current = self.space[str(end_coordinates)]
            
            while True:
                coordinates = current.coordinates
                x.append(coordinates[0])
                y.append(coordinates[1])
                z.append(coordinates[2])
                polarities.append(current.polarity)
                if current.previous is None:
                    break
                else:
                    current = self.space[str(current.previous)]
            
            colors = {}
            colors[1] = 'red'
            colors[-1] = 'blue'
            df = pd.DataFrame()
            df['x'] = x
            df['y'] = y
            df['z'] = z
            df['c'] = [colors[i] for i in polarities]

            if simulating:
                return df
            if x_min is None:
                x_min = min(df.x)
                y_min = min(df.y)
                z_max = min(df.z)

            if min(df.x) < x_min:
                x_min = min(df.x)
            if min(df.y) < y_min:
                y_min = min(df.y)
            if max(df.z) > z_max:
                z_max = max(df.z)
    
            fig1 = go.Scatter3d(x=df['x'], y=df['y'], z=df['z'], mode='markers+lines', marker = dict(
                                         size = 7,
                                         color = df['c'],
                                         opacity = 0.8 ), showlegend=False)
            
            data.append(fig1)

        annotations = [dict(x=x_min-2, y=y_min-2, z=z_max+1, text='Energy: {}'.format(round(self.energy,2)), showarrow=False), 
                       dict(x=x_min-2, y=y_min-2, z=z_max+3, text='Native Contacts: {}'.format(self.native_contacts), showarrow=False),
                       dict(x=x_min-2, y=y_min-2, z=z_max+5, text='Non Covalent Hydrophobic Contacts: {}'.format(self.non_covalent_hydrophobic_contacts), showarrow=False),
                       dict(x=x_min-2, y=y_min-2, z=z_max+7, text='Current Beta: {}'.format(round(self.beta,2)), showarrow=False)]          
        fig = go.Figure(data=data)
        fig.update_layout(scene=dict(annotations=annotations))
        if to_html:
            fig.write_html(html_file_name)
        else:
            plotly.offline.iplot(fig, filename='simple-3d-scatter')
