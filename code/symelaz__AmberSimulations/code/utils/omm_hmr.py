"""
omm_hmr.py

Hydrogen Mass Repartitioning (HMR) for OpenMM.

Works for AMBER (AmberPrmtopFile) and CHARMM (CharmmPsfFile) inputs: it only
uses top.topology, which both provide. (The original CHARMM-GUI version used
psf.bond_list and therefore crashed for AMBER topologies.)

Each hydrogen bonded to a heavy atom gets `hmass`; the added mass is taken from
the heavy atom, so the total mass is unchanged. Water is NOT repartitioned
(it is rigid, SETTLE handles it). Extra points / virtual sites (mass 0) are
never touched.

Usage: use a 4 fs time step (dt = 0.004) together with cons = HBonds.
"""

from __future__ import print_function

from simtk.unit import *
from simtk.openmm import *
from simtk.openmm.app import *

WATER_RESNAMES = {'HOH', 'WAT', 'TIP3', 'TIP3P', 'T3P', 'TP3', 'T4E', 'TP4', 'T4P', 'TIP4', 'OPC', 'SOL', 'SPC'}

def HydrogenMassRepartition(system, top, hmass=3.024):
    hydrogen = element.hydrogen
    n_h = 0
    for atom1, atom2 in top.topology.bonds():
        if atom1.residue.name.upper() in WATER_RESNAMES: continue
        if atom1.element == hydrogen:
            atom1, atom2 = atom2, atom1  # now atom2 is the hydrogen (if any)
        if atom2.element != hydrogen or atom1.element in (hydrogen, None): continue
        mass_h = system.getParticleMass(atom2.index)
        transfer = hmass*dalton - mass_h
        system.setParticleMass(atom2.index, hmass*dalton)
        system.setParticleMass(atom1.index, system.getParticleMass(atom1.index) - transfer)
        n_h += 1
    print("HMR: %d hydrogens set to %.3f Da (water excluded)" % (n_h, hmass))
    return system
