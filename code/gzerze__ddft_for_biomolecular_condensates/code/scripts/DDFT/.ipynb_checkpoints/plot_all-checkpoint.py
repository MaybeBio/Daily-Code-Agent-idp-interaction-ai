import matplotlib.pyplot as plt
import numpy as np
import Protein_RPA.utils.global_vars as gv

def plot_density_2D(System,t,n_capture,steps,output_dir):
	fig, ax  =plt.subplots()
	cs = ax.tricontourf(System.mesh.x.globalValue,System.mesh.y.globalValue,gv.r_res*System.phi.globalValue,cmap=plt.cm.get_cmap("Blues"),levels=np.linspace(0,0.2,256))	
	fig.colorbar(cs)
	ax.legend(title=str(t.value))
	ax.set_xlabel('y')
	ax.set_ylabel('x')
	plt.title("Volume fraction")
	fig.savefig(fname=output_dir+'P_step_{step}.png'.format(step=steps),dpi=300,format='png')


def plot_chempotential(System,t,n_capture,steps,output_dir,xi_left,xi_right):
	fig, ax  =plt.subplots()
	cs = ax.tricontourf(System.mesh.x.globalValue,System.mesh.y.globalValue,System.xi.globalValue,cmap=plt.cm.get_cmap("Blues"),levels=np.linspace(xi_left,xi_right,256))
	fig.colorbar(cs)
	ax.set_xlabel('y')
	ax.set_ylabel('x')
	plt.title("Chemical Potential")
	ax.legend(title=str(t.value))
	fig.savefig(fname=output_dir+'M_step_{step}.png'.format(step=steps),dpi=300,format='png')


def plot_density_1D(System,t,n_capture,steps,output_dir):
	#phi_value=System.phi.globalValue.reshape(nx,ny)
	fig, ax  =plt.subplots()
	xp=System.mesh.x.globalValue.reshape(System.nx,System.ny)
	ax.plot(xp[System.nx//2,:],gv.r_res*System.phi_value[System.ny//2,:])
	ax.set_xlabel('x')
	ax.set_ylim([0,0.16])
	plt.title("Volume fraction")
	ax.legend(title=str(t.value))
	fig.savefig(fname=output_dir+'D_step_{step}.png'.format(step=steps),dpi=300,format='png')
