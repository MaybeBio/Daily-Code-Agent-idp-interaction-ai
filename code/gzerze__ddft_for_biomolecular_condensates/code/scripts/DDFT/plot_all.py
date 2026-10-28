import matplotlib.pyplot as plt
import numpy as np
import Protein_RPA.utils.global_vars as gv

"""
Plotting utilities for visualizing simulation results.

This module contains functions for generating 2D and 1D plots of the
simulated volume fraction and chemical potential. The functions use
the simulation mesh and field variables stored in the System object
and save the resulting figures to the output directory.
"""



def plot_density_2D(System,t,n_capture,steps,output_dir):
	"""
    Plot the 2D volume fraction of the system.
    Parameters:
    System : Simulation system containing the mesh and volume-fraction field.
    t : time interval
    n_capture : frequence of image capture
    steps : Simulation step number used in the output filename.
    output_dir : Directory where the figure will be saved.
    """
	
	fig, ax  =plt.subplots()
	cs = ax.tricontourf(System.mesh.x.globalValue,System.mesh.y.globalValue,gv.r_res*System.phi.globalValue,cmap=plt.cm.get_cmap("Blues"),levels=np.linspace(0,0.2,256))	
	fig.colorbar(cs)
	ax.legend(title=str(t.value))
	ax.set_xlabel('y')
	ax.set_ylabel('x')
	plt.title("Volume fraction")
	fig.savefig(fname=output_dir+'P_step_{step}.png'.format(step=steps),dpi=300,format='png')


def plot_chempotential(System,t,n_capture,steps,output_dir,xi_left,xi_right):
    """
    Plot the 2D chemical potential.
    Parameters:
    System : Simulation system containing the mesh and volume-fraction field.
    t : time interval
    n_capture : frequence of image capture
    steps : Simulation step number used in the output filename.
    output_dir : Directory where the figure will be saved.
    xi_left :  Minimum value of the chemical-potential color scale.
    xi_right : Maximum value of the chemical-potential color scale.
    """

	
	fig, ax  =plt.subplots()
	cs = ax.tricontourf(System.mesh.x.globalValue,System.mesh.y.globalValue,System.xi.globalValue,cmap=plt.cm.get_cmap("Blues"),levels=np.linspace(xi_left,xi_right,256))
	fig.colorbar(cs)
	ax.set_xlabel('y')
	ax.set_ylabel('x')
	plt.title("Chemical Potential")
	ax.legend(title=str(t.value))
	fig.savefig(fname=output_dir+'M_step_{step}.png'.format(step=steps),dpi=300,format='png')


def plot_density_1D(System,t,n_capture,steps,output_dir):
	"""
    Plot a 1D cross-section of the 2D volume fraction.
    The plot uses the center of the x-direction and plots the volume
    fraction along the corresponding y-direction.
    """

	
	fig, ax  =plt.subplots()
	xp=System.mesh.x.globalValue.reshape(System.nx,System.ny)
	ax.plot(xp[System.nx//2,:],gv.r_res*System.phi_value[System.ny//2,:])
	ax.set_xlabel('x')
	ax.set_ylim([0,0.16])
	plt.title("Volume fraction")
	ax.legend(title=str(t.value))
	fig.savefig(fname=output_dir+'D_step_{step}.png'.format(step=steps),dpi=300,format='png')

def plot_density_1Ds(System,t,n_capture,steps,output_dir):
	"""
    Plot the complete 1D volume-fraction data.
    This function assumes that System.mesh.x.globalValue and
    System.phi_value are already arranged consistently for plotting.
    """
	
	fig, ax  =plt.subplots()
	xp=System.mesh.x.globalValue.reshape(System.nx,System.ny)
	ax.plot(xp,gv.r_res*System.phi_value)
	ax.set_xlabel('x')
	ax.set_ylim([0,0.2])
	plt.title("Volume fraction")
	ax.legend(title=str(t.value))
	fig.savefig(fname=output_dir+'D_step_{step}.png'.format(step=steps),dpi=300,format='png')

def plot_densitylog_2D(System,t,n_capture,steps,output_dir):
	print()
