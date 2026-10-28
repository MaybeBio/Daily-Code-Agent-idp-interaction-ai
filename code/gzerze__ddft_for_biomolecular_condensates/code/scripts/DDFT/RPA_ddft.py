#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Main driver for initializing and running the Dynamic Density
Functional Theory (DDFT) simulation.
Evolves the conserved density field phi(r, t) of a phase-separating polymeric
or protein system using the FiPy finite-volume framework.


Numerical Method:
-----------------
- Spatial discretization: Cell-centered Finite Volume Method (FiPy).
- Boundary conditions: Generalized Robin-type boundary condition decomposition.
- Temporal integration: Implicit Euler with Picard / non-linear residual sweeps.
- Adaptive time-stepping: Backtracking with dt reduction on slow convergence;
  optimistic dt expansion upon convergence.
- Linear backend: PETSc iterative solvers with MPI synchronization.

This module evolves coupled spatial density fields over time using the FiPy 
finite-volume framework with PETSc-backed linear solvers. Time integration 
employs an adaptive time-stepping scheme with Picard/non-linear residual sweeps 
and coupled Robin-type boundary condition updates.
"""

from __future__ import print_function
import fipy as fp
import os
import sys

from fipy import * #finite volume numerical solver package
import fipy.tools.numerix as npx
import numpy as np #math operators etc for Python
import matplotlib.pyplot as plt #plotting and graph functions for Python 
import argparse
import Protein_RPA.utils.free_energy as f_en
from Protein_RPA.utils.input_parse import *
from Protein_RPA.utils.graphics import * #getting code from RPA calculation files from protein
import initial_density as density #obtaining initial condition files
from plot_all import * #plotting functions for results
import chempot_table as chempot #chemical potential as a function of density
from petsc4py import PETSc #numerical LA solver

import timeit #python module to measure execution time
import subprocess #run child processes, error catching and handling, etc from this module
from Protein_RPA import seq_list as sl #sequence charge calculations from protein (?)
from BC import Boundary #boundary conditions, Robin BCs used
import pickle #serializing module for python
from initialize import Initialize #imports initialize class written for this code

def run_DDFT(args):
    """
    Initialize and execute the DDFT simulation.

    The simulation evolves the density fields over time by solving
    a coupled system of numerical equations. The chemical potential,
    boundary-condition parameters, and density fields are updated
    at each time step.

    Parameters
    ----------
    args : argparse.Namespace
        Command-line interface arguments containing:
        - args.i : str, path to input configuration file.
        - args.o : str, path to directory for generated figures/plots.
        - args.s : str, path to directory for binary checkpoint files.
    """
    

    # 1. Initialize simulation timing and input parameters
    start = timeit.default_timer()

    # Read simulation parameters from the input configuration file
    input_parameters = input_parse(args.i);
  

    output_dir=args.o
    print(output_dir)
    state_dir=args.s
    elapsed = 0.0
    steps = 0
    petscwrapper = SerialPETScCommWrapper()

    petscwrapper.Barrier()


    #2. Initialize the simulation system and boundary conditions

    # Initialize the computational mesh, density fields, initial
    # conditions, and other physical and numerical parameters.
    System=Initialize(args)
    
    petscwrapper.Barrier()
    
    BC=Boundary(System)

    parallelComm.Barrier()
    alpha,beta,coeff_diff,coeff_impl,g,RobinCoeff=BC.parameters()


    #3. Configure time integration and output settings
    start = timeit.default_timer()


    # Define the simulation time variable used by FiPy
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

    #4. The main time-integration loop
    # Continue the simulation while all stopping criteria are satisfied:
    #     1. The elapsed simulation time has not exceeded the target duration.
    #     2. The total number of time steps has not been reached.
    #     3. The time step remains larger than the minimum allowed value.
    while (elapsed <= System.duration) and (steps <= System.total_steps) and (dt>System.dt_min):
        
        if parallelComm.procID == 0:
            sys.stdout.write("steps dt  "+ str(steps)+"  "+str(dt))    


        start2=timeit.default_timer()
    
        petscwrapper.Barrier()

        # Calculate the chemical potential and its relevant derivatives
        # This is then followed with update/calculation of Boundary Condition parameters
        System.calc_chemical()
        BC.update_parameters(System)

        
        petscwrapper.Barrier()

        # COUPLED PDE SYSTEM DISCRETIZATION
        #The implementation below is done so we can avoid solving the equations using forward euler causing the simulation to either blow up with numerical instability or take weeks to finish.
        # ---------------------------------------------------------------------
        # Equation 0 (Evolution of phi):
        # Solves d(phi)/dt =div(phi * grad(xi)) + Robin_BC_terms (boundary condition)
        # Implemented using a convective term with face-gradient velocity field grad(xi).

        # Equation 1 ( Equation for xi):
        # FiPy discretization & Linearization:
        # - LHS ImplicitSourceTerm(1.0, xi): Solves for xi as an implicit unknown
        # - RHS Linearized mu_local(phi): Taylor expansion around (phi_now, xi_now):
        #     mu_local(phi) ≈ xi_now + (dxi/dphi)_now * (phi - phi_now)
        #   Represented by an implicit coupling term for phi:
        #     + ImplicitSourceTerm(coeff=dxi_now, var=phi)
        #   and explicit load vector contributions:
        #     + xi_now - dxi_now * phi_now
        # - RHS -DiffusionTerm(kappa, phi): FiPy defines DiffusionTerm(G, u) as
        #   div(G * grad(u)). With the leading minus sign, this represents
        #   -div(kappa * grad(phi)) = -kappa * laplacian(phi).
        
        #Convection Equation with Convection BC:
        eqn0=fp.TransientTerm(coeff=1.0,var=System.phi)==fp.ConvectionTerm(coeff=System.xi.faceGrad,var=System.phi) +alpha*(BC.RobinCoeff * g * coeff_diff).divergence + beta*fp.ImplicitSourceTerm(coeff=(BC.RobinCoeff * coeff_impl).divergence,var=System.phi)
        eqn1=fp.ImplicitSourceTerm(coeff=1.0,var=System.xi)== System.xi_now + fp.ImplicitSourceTerm(coeff=System.dxi_now,var=System.phi) -System.dxi_now*System.phi_now -fp.DiffusionTerm(coeff=System.kappa,var=System.phi) 

        # Couple both PDEs into a single non-linear block system
        eqn=eqn0 & eqn1


        start4=timeit.default_timer()
        #Every time both equations are solved, we need to make sure that the residue is smaller than tolerance
        #The impose on residue is done by sweep method untill residue (res0) satisfies the tolerance limit 
        res0=1 #Initial residual 
        sweep=0 # Counting number of sweeps
        while res0>tolerance:

            startsweep=timeit.default_timer()
            res0=eqn.sweep(dt=dt)
 
            if sweep>10:
                dt*=0.9
                System.phi[:]=System.phi.old
                System.xi[:]=System.xi.old
                sweep=0 # resetting sweep
            sweep+=1

        steps += 1
        elapsed += dt
        t.value = t.value +dt

        dt *= 1.1
        dt = min(dt, dt_max)
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
        if(dimension==1 and np.mod(steps,n_capture)==0):
            #plot_density_2D(System,t,n_capture,steps,output_dir)
            plot_density_1Ds(System,t,n_capture,steps,output_dir)	
            #plot_chempotential(System,t,n_capture,steps,output_dir,xi_left,xi_right)
        
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
     parser = argparse.ArgumentParser(description='Run DDFT simulations using the specified parameters')
     parser.add_argument('--i',help="Path to the input parameter file", required = True);
     parser.add_argument('--o',help="Path to the output directory for plots and results", required = True);
     parser.add_argument('--s',help="Path to the directory for saving simulation states", required = True);
     args = parser.parse_args();

     run_DDFT(args);
