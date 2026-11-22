#!/usr/bin/env python
import numpy as np
from openmm.app import *
from openmm import *
#from openmmtools.integrators import NoseHooverChainVelocityVerletIntegrator  # Add this import
from simtk.unit import *
from sys import stdout
from argparse import ArgumentParser
import re


def build_parser():
    parser = ArgumentParser(description="Run openmm simulations")
    parser.add_argument('-ff', type=str, help="Force field (xml)", default="ff99SBildn.xml")
    parser.add_argument('-pdb', type=str, help="Input pdb file", default="protein.pdb")
    parser.add_argument('-xtc', type=str, help="Output xtc file", default="traj.xtc")
    # Integrator details.
    parser.add_argument('-temp', type=float, help="temperature in Kelvin", default=298.15)
    parser.add_argument('-friction', type=float, help="friction coefficient in langevin thermostating (/ps)", default=298.15)
    parser.add_argument('-tstep', type=float, help="timestep in ps (default=0.002)", default=0.002)
    # Simulation details.
    parser.add_argument('-Nsteps', type=int, help="Number of steps (default=20000)", default=20000)
    # Save intervals.
    parser.add_argument('-save_xtc', type=int, help="save frame every x step (default=500-->1ps)", default=500)
    parser.add_argument('-save_stdout', type=int, help="report every x step (default=500-->1ps)", default=500)
    # Correction map.
    parser.add_argument('-cmap', help="Input CMAP file (Optional). The units must be in kcal/mol in the file and are automatically converted in kJ/mol for openMM", default=None)
    # bAIes.
    parser.add_argument('-baies', help="Input BAIES file (Optional). Must contain the atom indexes and the parameters the same format as PLUMED input files", default=None)
    parser.add_argument('-model', help="bAIes histogram model (gauss or lognorm)", default="gauss")
    parser.add_argument('-prior', help="bAIes error prior (none, jeffreys (default) or cauchy)", default="jeffreys")
    return parser
    
    
def read_cmap_file(cmapfile):
    output_cmap = dict()
    with open(cmapfile, "r") as f:
        cmap_data = [el.rstrip("\n") for el in f.readlines()]
    type_ = None
    title = None
    output_matrix = []
    bins = []
    rama = []
    for el in cmap_data:
        if re.search("-\(", el) and re.search(", type \d+", el):
            if type_ is not None:
                rama = np.reshape(rama, (24, 24))
                output_cmap[title] = rama
            type_ = el.split()[-1]
            title = el.split()[1].rstrip(",")
            bins = []
            rama = []
        elif re.match("# ", el) and re.search("\d+", el):
            bins.append(float(el.split()[1]))
        elif type_ is not None and el != "":
            rama += [float(ell) for ell in el.split()]
            pass
    rama = np.reshape(rama, (24, 24))
    output_cmap[title] = rama
    return output_cmap


def shift_matrix_values(matrix, step_0=1, step_1=0):
    new_matrix = np.zeros(matrix.shape)
    for i in range(matrix.shape[0]):
        for j in range(matrix.shape[1]):
            new_matrix[i][j] = matrix[(i-step_0)%matrix.shape[0]][(j-step_1)%matrix.shape[1]]
    return new_matrix


if __name__ == "__main__":
    args = build_parser().parse_args()
    
    pdb = PDBFile(args.pdb)
    forcefield = ForceField(args.ff)
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=CutoffNonPeriodic, constraints=HBonds) #nonbondedCutoff=1*nanometer
    #system.setDefaultPeriodicBoxVectors((20, 0, 0), (0, 20, 0), (0, 0, 20))
        
    ### ADDING CMAPS
    if args.cmap:
        CMAP_force = CMAPTorsionForce()
        # Read CMAPs file
        CMAP_list = dict()
        CMAP_data = read_cmap_file(args.cmap) # format 'XXX-(XXX)'
        for k in CMAP_data.keys():
            print("adding CMAP force {}".format(k))
            matrix = CMAP_data[k] * 4.184 # (multiplication factor to convert from kcal to kJ) 
            # Adapt format for openMM: value range [0, 2pi] instead of [-pi, pi].
            formatted_matrix = shift_matrix_values(matrix, int(matrix.shape[0]/2), int(matrix.shape[1]/2)).transpose()
            index = CMAP_force.addMap(24, formatted_matrix.flatten()) # matrix[i+size*j] -> phi = i2pi/size ; psi = j2pi/size
            CMAP_list[k] = index
    
        # look for openmm index:
        residues_dict = dict()
        for atom_index, atom in enumerate(pdb.topology.atoms()):
            resid = int(atom.residue.id)
            resname = atom.residue.name
            name = atom.name
            chainID = atom.residue.chain.id
            #atom_info = (atom.name, int(atom.residue.id), atom.residue.chain.id)
            if (resid, chainID) not in residues_dict.keys():
                residues_dict[(resid, chainID)] = {"resname": resname}
            residues_dict[(resid, chainID)][name] = atom_index
    
        # Add torsion correction maps:   
        for i, (resid, chainID) in enumerate(residues_dict.keys()):
            if (resid-1, chainID) in residues_dict.keys() and "C" in residues_dict[(resid-1, chainID)].keys() and (resid+1, chainID) in residues_dict.keys() and "N" in residues_dict[(resid+1, chainID)].keys():
                print("Adding CMAP on residue {}".format(resid))
                # This residue can have a CMAP:
                a1 = residues_dict[(resid-1, chainID)]["C"]
                a2 = residues_dict[(resid, chainID)]["N"]
                a3 = residues_dict[(resid, chainID)]["CA"]
                a4 = residues_dict[(resid, chainID)]["C"]
                b1 = a2
                b2 = a3
                b3 = a4
                b4 = residues_dict[(resid+1, chainID)]["N"]
    
                resname_next = "PRO" if residues_dict[(resid+1, chainID)]['resname'] == "PRO" else "XXX"
                CMAP_force.addTorsion(CMAP_list['{}-({})'.format(residues_dict[(resid, chainID)]['resname'], resname_next)], a1, a2, a3, a4, b1, b2, b3, b4)
                # For CS-CMAP: the keys of CMAP_list will simply be the resids, which makes it simple to implement.
                
        system.addForce(CMAP_force)


    ### ADDING BAIES
    if args.baies:
        kT = args.temp * 8.314462618*1e-3 # in kJ/mol
        # Define the force
        if args.model == "gauss":
            if args.prior == "jeffreys":
                baies_force = CustomBondForce("- kT * log( erf((r - mu)/(sigma*sqrt(2)) ) / (r - mu) )")
            elif args.prior == "cauchy":
                baies_force = CustomBondForce("kT * log( 0.5 * (r - mu)^2 + sigma^2 )")
            else:
                baies_force = CustomBondForce("kT * ((r - mu)^2) / (2 * sigma^2) ")
        elif args.model == "lognorm":
            if args.prior == "jeffreys":
                baies_force = CustomBondForce("- kT * log( erf((log(r) - mu)/(sigma*sqrt(2)) ) / (log(r) - mu) )")
            elif args.prior == "cauchy":
                baies_force = CustomBondForce("kT * log( 0.5 * (log(r) - mu)^2 + sigma^2 )")
            else:
                baies_force = CustomBondForce("kT * ((log(r) - mu)^2) / (2 * sigma^2) ")
        baies_force.addGlobalParameter("kT", kT)
        baies_force.addPerBondParameter("mu")
        baies_force.addPerBondParameter("sigma")
        # Load the parameters
        baies_data = np.loadtxt(args.baies)
        if len(list(baies_data.shape)) == 1:
            atom_idx_1 = baies_data[1] - 1
            atom_idx_2 = baies_data[2] - 1
            mus = baies_data[3]
            sigmas = baies_data[4]
            baies_force.addBond(int(atom_idx_1), int(atom_idx_2), [mus, sigmas])
        else:
            atom_idx_1 = np.array([el - 1  for el in baies_data[:,1]])
            atom_idx_2 = np.array([el - 1  for el in baies_data[:,2]])
            mus = baies_data[:,3]
            sigmas = baies_data[:,4]
            for i in range(mus.shape[0]):
                baies_force.addBond(int(atom_idx_1[i]), int(atom_idx_2[i]), [mus[i], sigmas[i]])
        system.addForce(baies_force)

    
    ### THERMOSTAT

    ## Langevin thermostat
    integrator = LangevinMiddleIntegrator(args.temp*kelvin, args.friction/picosecond, args.tstep*picoseconds)

    # Platform (GPU (CUDA) or CPU (CPU))
    platform = Platform.getPlatformByName('CUDA')
    platform.setPropertyDefaultValue('Precision', 'double') # Important for bAIes
    #properties = {'CudaDeviceIndex': '0,1', 'CudaPrecision': 'double'}

    # Define simulation
    simulation = Simulation(pdb.topology, system, integrator, platform) #, properties)
    simulation.context.setPositions(pdb.positions)
    
    # Check force exclusions
    forces = system.getForces()
    for i, force in enumerate(forces):
        print(f"Force {i}: {force.__class__.__name__}")
        if hasattr(force, 'getNumExceptions'):
            print(f"  Exceptions: {force.getNumExceptions()}")
        if hasattr(force, 'getNumExclusions'):
            print(f"  Exclusions: {force.getNumExclusions()}")

    # Energy minimization
    simulation.minimizeEnergy()

    # Reporters
    simulation.reporters.append(XTCReporter(args.xtc, reportInterval=args.save_xtc)) # 5k = 10ps
    simulation.reporters.append(StateDataReporter(stdout, reportInterval=args.save_stdout, step=True, potentialEnergy=True, temperature=True))
    
    # Production
    simulation.step(args.Nsteps)











