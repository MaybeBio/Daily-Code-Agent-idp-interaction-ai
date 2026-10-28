#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main function that is called to initialize and run DDFT
This code is the one that works well for non-salt case, this is derived from modifying v3.

"""

from __future__ import print_function
import fipy as fp
import os
import sys

#sys.path.insert(0,'/project/zerze/asilalah/software/RPA_ddft3/zero_flux/coalesce-20')
from fipy import * #finite volume numerical solver package
import fipy.solvers.scipy as sps
import numpy as np #math operators etc for Python
import matplotlib.pyplot as plt #plotting and graph functions for Python 
import argparse #unfamiliar with this, read more, parses input from files to text i think?
import Protein_RPA.utils.free_energy as f_en
from Protein_RPA.utils.input_parse import *
from Protein_RPA.utils.graphics import * #getting code from RPA calculation files from protein
import initial_density as density #obtaining initial condition files
from plot_all import * #plotting functions for results
import chempot_table as chempot #chemical potential as a function of density, learn more
from petsc4py import PETSc #numerical LA solver


#from boundary_condition import *mod
import timeit #python module to measure execution time
import subprocess #run child processes, error catching and handling, etc from this module
from Protein_RPA import seq_list as sl #sequence charge calculations from protein (?)
from BC import Boundary #boundary conditions, Robin BCs used
import pickle #serializing module for python
from initialize import Initialize #imports initialize class written for this code

def run_CDFT(args):
    """
    Function that calculate the volume density of monomer at different time steps
    """

    start = timeit.default_timer()

    input_parameters = input_parse(args.i);
  

    output_dir=args.o
    print(output_dir)
    state_dir=args.s
    elapsed = 0.0
    steps = 0
    petscwrapper = SerialPETScCommWrapper()

    #    sys.path.insert(0,'/project/zerze/asilalah/software/RPA_ddft3/zero_flux/coalesce-20')
    petscwrapper.Barrier()

    # Initialize the system and Boundary Condition (BC)
    System=Initialize(args)

    petscwrapper.Barrier()
    BC=Boundary(System)

    parallelComm.Barrier()
    alpha,beta,coeff_diff,coeff_impl,g,RobinCoeff=BC.parameters()

    start = timeit.default_timer()


    # Time discretization
    t = fp.Variable(0.0)
    dt = input_parameters['dt'];

    dt = input_parameters['dt'];
    dt_max = input_parameters['dt_max'];
    dt_min = input_parameters['dt_min'];
    tolerance = input_parameters['tolerance'];
    dimension=int(input_parameters['dimension']);

    steps=System.steps0
    t.value+=System.t0
    n_capture=input_parameters['n_capture']
    #file_state=os.getcwd()+'/'+state_dir+'state'+'.pickle'
    file_state=state_dir+'state'+'.pickle'
    dir_states=state_dir+'states/'

    n_save=input_parameters['state_save']
    xi_left=np.min(System.xi.globalValue)-0.05
    xi_right=np.max(System.xi.globalValue)+0.05


    #while (elapsed <= System.duration) and (steps <= System.total_steps) and (dt>System.dt_min):
        
    if parallelComm.procID == 0:
        sys.stdout.write("steps dt  "+ str(steps)+"  "+str(dt))    


    start2=timeit.default_timer()

    petscwrapper.Barrier()
    #First we calculate chemical potential, and chemical potential derivative
    # This is then followed with update/calculation of Boundary Condition parameters
    System.calc_chemical()
    BC.update_parameters(System)

    
    petscwrapper.Barrier()

    #Convection Equation with Convection BC:
    eqn0=fp.TransientTerm(coeff=1e-3,var=System.phi)==fp.ConvectionTerm(coeff=System.xi.faceGrad,var=System.phi) +alpha*(BC.RobinCoeff * g * coeff_diff).divergence + beta*fp.ImplicitSourceTerm(coeff=(BC.RobinCoeff * coeff_impl).divergence,var=System.phi)
    eqn1=fp.ImplicitSourceTerm(coeff=1.0,var=System.xi)== System.xi_now + fp.ImplicitSourceTerm(coeff=System.dxi_now,var=System.phi) -System.dxi_now*System.phi_now -fp.DiffusionTerm(coeff=System.kappa,var=System.phi) 


    #eqn is defined as both equations 
    eqn=eqn0 & eqn1


    start4=timeit.default_timer()
    #Every time both equations are solved, we need to make sure that the residue is smaller than tolerance
    #The impose on residue is done by sweep method untill residue (res0) satisfies the tolerance limit 
    res0=1 #Initial residual 
    sweep=0 # Counting number of sweeps

    while res0>tolerance:

        startsweep=timeit.default_timer()
        res0=eqn.sweep()

        if sweep>10:
            dt*=0.9
            System.phi[:]=System.phi.old
            System.xi[:]=System.xi.old
            sweep=0 # resetting sweep
        sweep+=1

    solver = sps.LinearLUSolver(tolerance = System.tolerance)
    #eqn.solve(dt=1.,solver=solver)
    #steps += 1
    #elapsed += dt
    #t.value = t.value +dt

    #dt *= 1.1
    #dt = min(dt, dt_max)
    System.phi.updateOld()
    System.xi.updateOld()
    
    if petscwrapper.procID == 0:
        sys.stdout.write('Wall time step c'+ str(steps)+'  '+str(round((timeit.default_timer()-start2),2))+' s'+'\n')

    petscwrapper.Barrier()
    start5=timeit.default_timer()
    System.global_phi()

    #Write output in png and movie formats
    if(dimension==2 and np.mod(steps,n_capture)==0):
        plot_density_2D(System,t,n_capture,steps,output_dir)
        plot_density_1D(System,t,n_capture,steps,output_dir)	
        plot_chempotential(System,t,n_capture,steps,output_dir,xi_left,xi_right)
    
    if(np.mod(steps+1,n_save)==0):
        with open(file_state, "wb") as f:
            pickle.dump([steps,t.value,System.input_parameters,System.phi_value,System.xi_value],f)
        f.close()

        with open(dir_states+'state-'+str(steps+1)+'.pickle', "wb") as f:
            pickle.dump([steps,t.value,System.input_parameters,System.phi_value,System.xi_value],f)
        f.close()

    #epcomm = Epetra.PyComm()
    n_proc=petscwrapper.Nproc
    if parallelComm.procID == 0:
        sys.stdout.write('number of processor '+str(n_proc)+'\n')
        sys.stdout.write('nx ny '+str(System.nx)+' '+str(System.ny)+'\n')
        sys.stdout.write('Wall time '+str(round((timeit.default_timer()-start),2))+' s'+'\n')
       

    with open(file_state, "wb") as f:
        pickle.dump([steps,t.value,System.input_parameters,System.phi_value,System.xi_value],f)
    f.close() 
    

if __name__ == "__main__":
     parser = argparse.ArgumentParser(description='Take output filename to run CH simulations')
     parser.add_argument('--i',help="Name of input params", required = True);
     parser.add_argument('--o',help="Name of output folder", required = True);
     parser.add_argument('--s',help="Name of state folder", required = True);
     args = parser.parse_args();

     run_CDFT(args);
