import matplotlib.pyplot as plt
import numpy as np
import Protein_RPA.utils.global_vars as gv
from matplotlib import ticker

"""
Plotting utilities for visualizing simulation results.

This module contains functions for generating 2D and 1D plots of the
simulated volume fraction and chemical potential. The functions use
the simulation mesh and field variables stored in the System object
and save the resulting figures to the output directory.
"""


def plot_densitylog_2D(mesh,phi,t,steps,output_dir):
	fig, ax  =plt.subplots()
	cs = ax.tricontourf(mesh.x.globalValue,mesh.y.globalValue,gv.r_res*phi.globalValue,cmap=plt.cm.get_cmap("Blues"),levels=np.linspace(0.001,0.2,256),locator=ticker.LogLocator())	
	fig.colorbar(cs,ticks=[0.0001,0.001,0.01,0.1])
	ax.legend(title=str(int(t)))
	ax.set_xlabel('y')
	ax.set_ylabel('x')
	plt.title("Volume fraction")
	fig.savefig(fname=output_dir+'Plog_step_{step}.png'.format(step=steps),dpi=300,format='png')
	plt.close(fig)

def plot_density_2D(mesh,phi,t,steps,output_dir,phi_range):
	fig, ax  =plt.subplots()
	cs = ax.tricontourf(mesh.x.value,mesh.y.value,gv.r_res*phi.value,cmap=plt.cm.get_cmap("Blues"),levels=np.linspace(phi_range[0],phi_range[1],256))	
	fig.colorbar(cs, cax=cbar_ax, label='$\phi$')
	ax.legend(title=str(round(float(t)*0.07,1))+" ns",title_fontsize=15)
	ax.set_xlabel('y')
	ax.set_ylabel('x')
	#plt.title("Volume fraction")
	fig.savefig(fname=output_dir+'P_step_{step}.png'.format(step=steps),dpi=300,format='png')
	plt.close(fig)


def plot_chempotential(mesh,xi,t,steps,output_dir,xi_range):
	fig, ax  =plt.subplots()
	xi_min,xi_max=np.min(xi.value),np.max(xi.value)
	fig, ax
	cs = ax.tricontourf(mesh.x.value,mesh.y.value,xi.value,cmap=plt.cm.get_cmap("Blues"),levels=np.linspace(xi_min,xi_max,256))
	fig.colorbar(cs, cax=cbar_ax, label='$\mu$')
	ax.set_xlabel('y')
	ax.set_ylabel('x')
	#plt.title("Chemical Potential")
	ax.legend(title=str(round(float(t)*0.07,1))+" ns",title_fontsize=15)
	fig.savefig(fname=output_dir+'M_step_{step}.png'.format(step=steps),dpi=300,format='png')
	plt.close(fig)


def plot_density_1D(mesh,phi_value,t,steps,output_dir,phi_range):
	fig, ax  =plt.subplots()
	xp=mesh.x.value.reshape(mesh.nx,mesh.ny)
	ax.plot(xp[mesh.nx//2,:],gv.r_res*phi_value[mesh.ny//2,:])
	ax.set_xlabel('x')
	ax.set_ylim([phi_range[0],phi_range[1]])
	plt.title("Volume fraction")
	ax.legend(title=str(t))
	fig.savefig(fname=output_dir+'D_step_{step}.png'.format(step=steps),dpi=300,format='png')
	plt.close(fig)

def plot_densitylog_1D(mesh,phi_value,t,steps,output_dir,phi_range):
	fig, ax  =plt.subplots()
	xp=mesh.x.value.reshape(mesh.nx,mesh.ny)
	ax.semilogy(xp[mesh.nx//2,:],gv.r_res*phi_value[mesh.ny//2,:])
	ax.set_xlabel('x')
	ax.set_ylim([phi_range[0],phi_range[1]])
	plt.title("Volume fraction")
	ax.legend(title=str(t))
	fig.savefig(fname=output_dir+'Dlog_step_{step}.png'.format(step=steps),dpi=300,format='png')
	plt.close(fig)
