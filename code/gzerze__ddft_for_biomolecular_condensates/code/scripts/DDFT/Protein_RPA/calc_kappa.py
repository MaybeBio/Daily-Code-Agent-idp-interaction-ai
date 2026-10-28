# RPA+FH model for single charge sequence
# Short-range interaction contributes to only the k=0 FH term
# The FH term ehs parameters follow the definition in the PRL paper



import sys
import time
import multiprocessing as mp
import numpy as np

import f_min_solve_1p_1salt as fs
import thermal_1p_1salt as tt
import utils.global_vars as gv
import seq_list as sl

# a= radius of water (~2.4A)
# gv.eta= (b/a)^3=(3.8/2.4)^3=3.96~4
# gv.r_res=(rm/a)^3~1
# gv.r_sal=(rs/a)^3~1


# convert [Salt] in mM to [Salt]/[H2O]
phis_mM = float(sys.argv[5])
phis = phis_mM*0.001/(1000./18.)

# ehs must be a 2-element list: the first for entropy and the latter enthalpy
ehs = [float(x) for x in sys.argv[3:5]]

# use  phi-dependent permittivity or not 
ef = False

if sys.argv[2] == "cri_calc":
   cri_calc_end = 1
else:
   cri_calc_end = 0
   uref=float(sys.argv[2])
   umax = float(sys.argv[2])

du = 0.1

#=========================== Set up parameters ===========================
seq_name = sys.argv[1]  # Select a sequence in seq_list.py
sig, N, the_seq = sl.get_the_charge(seq_name)

HP = tt.RPAFH(sig, ehs=ehs, epsfun=ef)

#======================= Calculate critical point ========================

print('Seq:' , seq_name, '=' , the_seq ,'\nphi_s=', phis , \
      'r_res ='   , gv.r_res , ', r_con =' , gv.r_con , \
      ', r_sal =' , gv.r_sal , ', eta ='   , gv.eta, \
      '\nehs :' , ehs[0], ehs[1] )

t0 = time.time()

#critical_point
phi_cri, u_cri = fs.cri_calc(HP, phis)

print('Critical point found in', time.time() - t0 , 's')
print( 'u_cri =', '{:.4e}'.format(u_cri) , 
       'T*_cri = {:.4f}'.format(1/u_cri),
       ', phi_cri =','{:.8e}'.format(phi_cri) )

if(cri_calc_end):
    sys.exit()

#============================ Set up u range =============================
ddu = du/10;
umin = (np.floor(u_cri/ddu)+1)*ddu
uclose = (np.floor(u_cri/du)+2)*du
if umax < u_cri:
    umax = np.floor(u_cri*1.5) 
if uclose > umax:
    uclose = umax

umin=uref
umax=uref
uclose=uref

uall = np.append( np.arange(umin, uclose, ddu ), \
                  np.arange(uclose, umax+ddu, ddu ) )
print('sizzling')
print(uall, flush=True)


#==================== Parallel calculate multiple u's ====================


def bisp_parallel(u):
    sp1, sp2 = fs.ps_sp_solve( HP, phis, u, phi_cri  ) 
    print( u, sp1, sp2, 'sp done!', flush=True)
    bi1, bi2 = fs.ps_bi_solve( HP, phis, u, [sp1, sp2], phi_cri) 
    print( u, bi1, bi2, 'bi done!', flush=True)
    
    return sp1, sp2, bi1, bi2



if __name__ == '__main__':

	pool = mp.Pool(processes=4)
	sp1ss, sp2ss, bi1ss, bi2ss = zip(*pool.map(bisp_parallel, uall))

	ind_slc = np.where(np.array(bi1ss) > gv.phi_min_sys)[0]
	unew = uall[ind_slc]
	sp1s, sp2s = np.array(sp1ss)[ind_slc], np.array(sp2ss)[ind_slc]
	bi1s, bi2s = np.array(bi1ss)[ind_slc], np.array(bi2ss)[ind_slc]
	new_umax = np.max(unew)
	nnew = ind_slc.shape[0]

	sp_out, bi_out = np.zeros((2*nnew+1, 2)), np.zeros((2*nnew+1, 2))
	sp_out[:,0] = np.append(np.append(sp1s[::-1], phi_cri), sp2s)
	sp_out[:,1] = np.append(np.append(unew[::-1],u_cri), unew)
	bi_out[:,0] = np.append(np.append(bi1s[::-1],phi_cri), bi2s)
	bi_out[:,1] = sp_out[:,1] 

	print("sp_out bi_out")
	print(sp_out)
	print(bi_out)


	monosize = str(gv.r_res) + '_' + str(gv.r_con) + '_' + str(gv.r_sal)
	ehs_str = '_'.join(str(x) for x in ehs )

	calc_info = '_RPAFH_N{}_phis_{:.5f}_{}_eh{:.2f}_es{:.2f}_umax{:.2f}_du{:.2f}_ddu{:.2f}.txt'.format(
             N, phis, seq_name,ehs[0], ehs[1],new_umax,du,ddu)
	sp_file = './results/sp' + calc_info
	bi_file = './results/bi' + calc_info

	print(sp_file)
	print(bi_file)

	cri_info = "u_cri= " + str(u_cri) + " , phi_cri= " + str(phi_cri)

	#np.savetxt(sp_file, sp_out, fmt = '%.8e', header= cri_info )
	#np.savetxt(bi_file, bi_out, fmt = '%.8e', header= cri_info )


	# We want to calculate the influence parameter 'kappa' that is responsible for surface tension
	# After running 'python calc_kappa.py Ddx4_N1 4 0.15 -0.3 0' ( at u=1/T*=4), we have dilute and dense concentration
	dilute=bi_out[0,0]
	dense=bi_out[2,0]
	nsp=100
	#discretizing the phi space 
	ddense=(dense-dilute)/nsp
	phi_space=np.linspace(dilute,dense,nsp) 
	Iphi=0
	# Constructing the free energy of two-phase coexistence (Fbase)
	Fbase=tt.feng(HP,dilute,phis,uref)+ ( tt.feng(HP,dense,phis,uref)-tt.feng(HP,dilute,phis,uref)) *(phi_space-dilute)/(dense-dilute)
	# Free energy of mixed phase is always higher than Fbase
	Fmixed=[]
	for i in range(nsp):
		Fmixed.append(tt.feng(HP,phi_space[i],phis,uref))
		print(Fmixed[i]-Fbase[i])
		#Iphi+=np.sqrt(tt.feng(HP, phi_space[i], phis, uref)-Fbase[i])*ddense
		if (Fmixed[i]-Fbase[i]) >= 0:
			Iphi+=np.sqrt(Fmixed[i]-Fbase[i])*ddense
	print("total I",Iphi)
	gamma=0.104 #(in mJ/m^2)
	gamma=gamma*1e-3 #(in J/m^2)
	gv.a=3.8 # The size of lattice/smallest component(water molecule)
	gv.bm=3.8 # The size of monomer/smallest component(water molecule)
        # range of gv.bm from ? to ?

	kappa=(gamma)**2 *(gv.a**3)*(gv.bm**3)/(4*(Iphi**2))
	print("Integral of energy",Iphi)
	print("kappa",kappa)
	print(" dense dilute ",dense," ",dilute)
	print(" dense dilute ",uref,"  ",tt.dfeng(HP, dense, phis, uref))
	print(" dense dilute ",uref,"  ",tt.dfeng(HP, dilute, phis, uref))
	print(" dense dilute ",uref,"  ",(tt.feng(HP, dense, phis, uref)-tt.feng(HP,dilute,phis,uref)/(dense-dilute )) )
