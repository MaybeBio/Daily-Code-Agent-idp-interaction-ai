import time

import os
import argparse, sys
import numpy as np

import gsd.hoomd
import hoomd
import hoomd.md
#import h5py

import sys
from PythonFuncs import *
from hoomd import ashbaugh_plugin
import copy
#from hoomd import plugin_HPS_SS
#import ashbaugh_plugin as aplugin

PWD = os.getcwd()

def run(FolderPath, Restart=False, ExtendedSteps=0, prerun=False, resizeSteps=100000):
    ### Read Input Data
    ### All inputs are in lammps units, have to convert to 
    #new Parameter SimType Calvados3 default C2
    #new Parameter Domain
    Params = readParam(f"{FolderPath}/HOOMD_Setup/Params.txt")

    Seqs, NBeads, NChains, InputBonds, InputAngles, InputDihedrals = readSequences(f"{FolderPath}/HOOMD_Setup/Sequences.txt")

    InputPositions, InputTypes, InputCharges, InputMasses, _, Diameter, InputImage = readParticleData(f"{FolderPath}/HOOMD_Setup/Particles.txt", NBeads, Seqs)
    if Params["UseAngles"]:
        dihedral_eps, dihedral_dict, dihedral_list, dihedral_IDs, dihedral_AllIDs = readDihedrals(f"{FolderPath}/HOOMD_Setup/DihedralMap.txt", Seqs, InputTypes)
    else:
        dihedral_eps, dihedral_dict, dihedral_list, dihedral_IDs, dihedral_AllIDs=[],[],[],[],[]

    #Other LambdaDict for C2/C3
    IDS, IDToResName, IDToCharge, IDToMass, IDToSigma, IDToLambda = readDictionaries(f"{FolderPath}/HOOMD_Setup/Dictionaries.txt")

    ### constants 
    bondLength = 0.38

    kb = 0.00831446262
    kT = kb*Params["Temp"]
    forces = []

    ### Prepare HOOMD 
    tmp = []
    for i in IDS: 
        tmp.append(str(IDToResName[i]))
    Types = np.array(tmp)


    device = hoomd.device.GPU() if Params["Device"]=="GPU" else hoomd.device.CPU()
    sim = hoomd.Simulation(device=device, seed=Params["Seed"])

    if prerun:
        Params["NSteps"] = resizeSteps+8000
        Params["Create_Start_Config"] = False

    if Restart:
        TrajectoryNumber , NStepsOld = CountNumberOfTrajectoryFiles(FolderPath)
        NStepsOld *= Params["NOut"] 
        lastTrajectory = Params['Trajectory'] if TrajectoryNumber==1 else f"{Params['Trajectory'][:-4]}_{TrajectoryNumber-1}.gsd"

        RestartPath=f"{FolderPath}{Params['Simname']}_Restart_{TrajectoryNumber}.gsd"
        CopyLastFrameToRestartFile(FolderPath+lastTrajectory, RestartPath)

        Params["NSteps"] = Params["NSteps"]-NStepsOld if ExtendedSteps==0 else ExtendedSteps-NStepsOld ### Either extend to meet initial goal or extend simulations.
        #Params["NSteps"] = Params["NSteps"] if Params["NSteps"]>0 else 0

        print(f"NSteps {Params['NSteps']}")
        #Params["NSteps"] = int(NewGoal)-int(NStepsOld)  if int(NewGoal)-int(NStepsOld) >0 else 0  ### avoid negativ steps

        with open(f"{FolderPath}/HOOMD_Setup/{Params['Simname']} +_RestartLog.txt", 'a+') as f:
            f.write(f"Restart at timestep {NStepsOld} trying to extend to {Params['NSteps']+NStepsOld} by running {Params['NSteps']} additional steps.\n")
        sim.timestep = int(NStepsOld)
        sim.create_state_from_gsd(filename=RestartPath) ### takes the last frame
        Params["Trajectory"] = f"{Params['Trajectory'][:-4]}_{TrajectoryNumber}.gsd"
    else:
        if Params["Create_Start_Config"]:
            snapshot = gsd.hoomd.Frame()  
            #if Params["Use_Minimised_GSD"]:
                ### Specify particles
            snapshot.particles.N = NBeads
            snapshot.particles.position = InputPositions.astype(np.float32) ### convert to nm
            snapshot.particles.types = tuple(Types)
            snapshot.particles.typeid = InputTypes.astype(np.int32)
            snapshot.particles.image = InputImage

            snapshot.configuration.box = [Params["Lx"], Params["Ly"], Params["Lz"], 0, 0, 0]

            snapshot.particles.mass = InputMasses.astype(np.float32)
            snapshot.particles.charge = InputCharges.astype(np.float32)

            ### use undocumented feature to use "floppy bodies"
            snapshot.particles.body = np.genfromtxt(f"{FolderPath}/HOOMD_Setup/Groups.txt") 

            B_N = NBeads-NChains
            B_types = ['O-O']
            B_typeid = np.zeros(B_N, dtype=np.int32)
            B_group = InputBonds#.astype(np.uint32)

            if Params["SimulationType"]=="Calvados3":
                ## read the ENM_indice and backbone data
                B_N, B_types, B_typeid, B_group, ENMharmonic = read_ENM_HOOD_indices(f"{FolderPath}/HOOMD_Setup/ENM_indices.txt")


            snapshot.bonds.N = B_N
            snapshot.bonds.types = list(B_types)
            snapshot.bonds.typeid = B_typeid
            snapshot.bonds.group = B_group

            ## Create Angles
            snapshot.angles.N = NBeads-2*NChains
            snapshot.angles.types = ['O-O-O']
            snapshot.angles.typeid = np.zeros(snapshot.angles.N, dtype=int)
            snapshot.angles.group = InputAngles

            ## Create Dihedrals
            #snapshot.dihedrals.N =  NBeads-3*NChains 
            #snapshot.dihedrals.types = []#list(dihedral_list)
            #snapshot.dihedrals.typeid =  dihedral_AllIDs
            #snapshot.dihedrals.group = InputDihedrals
            if prerun:
                with gsd.hoomd.open(name=FolderPath + Params["Simname"] + "_" + "prerun.gsd", mode='w') as f:
                    f.append(snapshot)
                sim.create_state_from_gsd(filename=FolderPath + Params["Simname"] + "_" + "prerun.gsd")
            else:
                with gsd.hoomd.open(name=FolderPath + Params["Simname"] + ".gsd", mode='w') as f:
                    f.append(snapshot)
                sim.create_state_from_gsd(filename=FolderPath + Params["Simname"] + ".gsd")

        else:
            print(f"Creat start from: {FolderPath+Params['Simname']}.gsd")
            sim.create_state_from_gsd(filename=FolderPath+Params["Simname"]+".gsd")
    if Params["NSteps"] ==0:
        print("The number of necessary steps is zero. No simulation will be run.")
        return -1

    if Params["SimulationType"]!="Calvados3":
        ### Harmonic bonds
        harmonic = hoomd.md.bond.Harmonic()
        harmonic.params['O-O'] = dict(k=8033, r0=bondLength) ###calvados2: k=8033kJ/mol/nm^2 k=1000kJ/nm^2 = 10KJ/AA^2
        forces.append(harmonic)
    if Params["SimulationType"]=="Calvados3":
        ### Harmonic bonds
        harmonic = hoomd.md.bond.Harmonic()

        ## read the ENM_indice and backbone data
        B_N, B_types, B_typeid, B_group, ENMharmonic = read_ENM_HOOD_indices(f"{FolderPath}/HOOMD_Setup/ENM_indices.txt")
        
        for typ in B_types: 
            harmonic.params[typ] = dict(k=ENMharmonic[typ]["k"], r0=ENMharmonic[typ]["r"])

        forces.append(harmonic)

    if Params["UseAngles"]:
        cell2 = hoomd.md.nlist.Tree(buffer=0.4,default_r_cut=Params["AHCutoff"], exclusions=['bond', 'angle', 'dihedral', '1-3', '1-4'])
        cell  = hoomd.md.nlist.Tree(buffer=0.4,default_r_cut=0.0, exclusions=['bond', 'angle', 'dihedral', '1-3', '1-4'] )
    else:
        cell2 = hoomd.md.nlist.Tree(buffer=0.4,default_r_cut=Params["AHCutoff"], exclusions=['bond','body'])
        cell  = hoomd.md.nlist.Tree(buffer=0.4,default_r_cut=0.0, exclusions=['bond','body'])


    if Params["UseCharge"]:
        # # electrostatics forces
        yukawa = hoomd.md.pair.Yukawa(nlist=cell)
        for i in IDToResName.keys():
            res_i = IDToResName[i]
            for j in IDToResName.keys():
                res_j = IDToResName[j]
                yukawa.params[(res_i,res_j)] = dict(epsilon=IDToCharge[i]*IDToCharge[j]*Params["yk_prefactor"], kappa=Params["kappa"]) ### (KJ , 1/nm)d

                if ((IDToCharge[i]==0) or (IDToCharge[j]==0)):
                    yukawa.r_cut[(res_i,res_j)] = 0.0
                else:
                    yukawa.r_cut[(res_i,res_j)] = Params["YukawaCutoff"]
        forces.append(yukawa)
    
    # # nonbonded: ashbaugh-hatch potential
    ash = ashbaugh_plugin.pair.AshbaughPair(nlist=cell2, default_r_cut=2.0, default_r_on=0.0)
    for i in IDToResName.keys():
        res_i = IDToResName[i]
        for j in IDToResName.keys():
            res_j = IDToResName[j]
            ash.params[(res_i, res_j)] = {"epsilon":0.8368, "sigma":round((IDToSigma[i]+IDToSigma[j])/20.0,5), "lam":(IDToLambda[i]+IDToLambda[j])/2.0} ### convert to nm as well
            ash.r_cut[(res_i,res_j)] = Params["AHCutoff"] 
    forces.append(ash)


    if Params["UseAngles"]:
        # ANGLE POTENTIAL 
        anglePotential = getAnglePotential()
        forces.append(anglePotential)

        ### DIHEDRAL
        dihedralPotential = getDihedralPotential(dihedral_list, dihedral_IDs, dihedral_eps)
        forces.append(dihedralPotential)


    if Params["Alt_GSD_Start"]!="-":  ### not needed?
        sim = hoomd.Simulation(device=device, seed=Params["Seed"])
        sim.create_state_from_gsd(filename=FolderPath+Params["Alt_GSD_Start"])


    ### logging quantitites
    thermodynamic_properties = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermodynamic_properties)

    logger = hoomd.logging.Logger(categories=['scalar', 'string'])
    #logger.add(thermodynamic_properties)

    write_mode  ='wb'
    if prerun:
        filename = FolderPath+"/prerun_"+Params["Trajectory"]
        gsd_writer = hoomd.write.GSD(trigger=hoomd.trigger.Periodic(int(Params["NOut"])), filename=filename, filter=hoomd.filter.All(), mode=write_mode,dynamic=['particles/position', 'particles/image', 'configuration/box']) # -> trigger can be set to NSteps for saving only the last frame
        gsd_writer.log = logger

        #logger = hoomd.logging.Logger(categories=['scalar', 'string']) #'sequence'
        logger.add(sim, quantities=['timestep', 'walltime', 'tps'])
        table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(100000), logger=logger)
        sim.operations.writers.append(table)
    else:
        write_mode = 'ab' if Restart else 'wb'
        gsd_writer = hoomd.write.GSD(trigger=hoomd.trigger.Periodic(Params["NOut"]), filename=FolderPath+Params["Trajectory"], filter=hoomd.filter.All(),mode=write_mode,dynamic=['particles/position', 'particles/image'])
        #sim.operations.writers.append(gsd_writer)
        gsd_writer.log = logger

        #logger = hoomd.logging.Logger(categories=['scalar', 'string']) #'sequence'
        logger.add(sim, quantities=['timestep', 'walltime', 'tps'])
        table = hoomd.write.Table(trigger=hoomd.trigger.Periodic(100000), logger=logger)
        sim.operations.writers.append(table)

    hdf5_writer=[]
    if True:
    ### track presure
        logger_pres = hoomd.logging.Logger(categories=['scalar', 'sequence'])
        logger_pres.add(sim, quantities=["timestep"])
        logger_pres.add(sim, quantities=["initial_timestep"])
        logger_pres.add(thermodynamic_properties, quantities=['kinetic_temperature','kinetic_energy', 'potential_energy','pressure', 'pressure_tensor'])
        write_mode = 'a' if Restart else 'w'
        hdf5_writer = hoomd.write.HDF5Log(
            trigger=hoomd.trigger.Periodic(100), filename=FolderPath+'pressure.h5', mode=write_mode, logger=logger_pres
        )
    else:
        print("Thermodynamic Quantities are not tracked!")


    ### apply langevin
    integrator = hoomd.md.Integrator(dt=Params["dt"]) 

    if prerun:
        rho = 0.4 / (1.66/1000) # Mass/Volume = dalton/AA^3 = dalton/nm^3/1000
        A = Params["Lx"] * Params["Lz"] # nm^2
        TotalMass = sum(InputMasses) # dalton
        # dalton /((dalton / nm^3) * nm^2) = nm 
        L = TotalMass / (A*rho) # nm

        Ly_start = np.max(InputPositions[:,1])-np.min(InputPositions[:,1])*1.
        state_box = [Params["Lx"], Ly_start , Params["Lz"], 0, 0, 0]
        goalbox = [Params["Lx"], L , Params["Lz"], 0, 0, 0]

        box_resize = hoomd.update.BoxResize(trigger=hoomd.trigger.Periodic(10), box1=state_box, box2=goalbox, variant=hoomd.variant.Ramp(A=1.0, B=2.0, t_start=8000, t_ramp=resizeSteps))
        sim.operations.updaters.append(box_resize)


    nvt = hoomd.md.methods.Langevin(filter=hoomd.filter.All(), kT=kT) ### ps^-1
    integrator.methods = [nvt]
    integrator.dt = Params["dt"]

    ### apply increased friction to thermalized quicker in preruns as in HPS model
    for i in IDToResName.keys():
        name = IDToResName[i]
        nvt.gamma[name] = IDToMass[i]*10.0**-3
        nvt.gamma_r[name] = (0.0, 0.0, 0.0)
        
    sim.operations.integrator=integrator
    sim.operations.integrator.forces=forces

    print("Before simulation\n")

    sim.operations.writers.append(gsd_writer)
    sim.operations.writers.append(hdf5_writer)

    ### pre equilibrate the bonds by dissipating energy from the stretched bonds 
    if prerun:
        for fac in [50.0, 10.0,5.0,2.0, 1.5, 1.0]:
            for i in range(10):
                print(f"fac {fac} i {i}")
                integrator.dt = Params["dt"]/fac
                sim.run(100)
                sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=kT/fac)
    else:
        if not Restart:
            ### apply higher friction to thermalize quicker during the 0.1 mu s
            ### the excess potential energy created by unraveling the 
            ### proteins that cross the pbc from the prerun
            integrator.dt = Params["dt"]
            
            sim.state.thermalize_particle_momenta(filter=hoomd.filter.All(), kT=kT)
            sim.run(10_000_000)

    ### apply friction as in HPS model
    for i in IDToResName.keys():
        name = IDToResName[i]
        nvt.gamma[name] = IDToMass[i]*10.0**-5
        nvt.gamma_r[name] = (0.0, 0.0, 0.0)

    ### optimise cell list buffer
    now = sim.timestep
    if now is None:
        now = 0
    hoomd.md.tune.NeighborListBuffer(trigger=hoomd.trigger.Before(now+20000), nlist=cell , maximum_buffer=1.0, solver=hoomd.tune.GradientDescent())
    hoomd.md.tune.NeighborListBuffer(trigger=hoomd.trigger.Before(now+20000), nlist=cell2, maximum_buffer=1.0, solver=hoomd.tune.GradientDescent())

    ### start simulation

    integrator.dt = Params["dt"]
    print(f"Run simulation for {Params['NSteps']} steps")
    sim.run(Params["NSteps"], write_at_start=True)

    print(f"initial {sim.initial_timestep}")
    print(f"final   {sim.final_timestep}")

    
    print(f"TPS: {sim.tps:0.5g}")
    print(f"WallTime: {sim.walltime:0.5g}")

def str2bool(v):
    v = v.lower()
    if v in ('yes','true','t','y','1'):  return True
    if v in ('no','false','f','n','0'):  return False
    raise argparse.ArgumentTypeError('Boolean value expected.')

#def run(folder, flag, steps):
#    print(f'>>> run(folder={folder!r}, flag={flag}, steps={steps})')
#    # … your actual work …

def main():
    p = argparse.ArgumentParser()
    p.add_argument('folder')
    p.add_argument('flag', type=str2bool)
    p.add_argument('steps', type=int)
    a = p.parse_args()
    print(a.folder, a.flag, a.steps)   # echo for debugging
    run(a.folder, a.flag, a.steps)

if __name__ == '__main__':
    sys.exit(main())
