# RPA+FH model for single charge sequence
# Short-range interaction contributes to only the k=0 FH term
# The FH term ehs parameters follow the definition in the PRL paper

# ver Git.1 Apr 14, 2020
# Upload to github

# ver0: May 27, 2018

import sys
import time
import multiprocessing as mp
import numpy as np

import f_min_solve_1p_1salt as fs
import thermal_1p_1salt as tt
import global_vars as gv
import seq_list as sl

gv.r_res = 1
gv.r_con = 1
gv.r_sal = 1
gv.eta = 1 


# convert [Salt] in mM to [Salt]/[H2O]
phis_mM = float(sys.argv[5])
phis = phis_mM*0.001/(1000./18.)

# ehs must be a 2-element list: the first for entropy and the latter enthalpy
ehs = [float(x) for x in sys.argv[3:5]]

# use  phi-dependent permittivity or not 
ef = False

u=float(sys.argv[2])


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

#======== Solve dilute and dense densities for one Temperature=======#
def bisp_parallel(u):
    sp1, sp2 = fs.ps_sp_solve( HP, phis, u, phi_cri  )
    print( u, sp1, sp2, 'sp done!', flush=True)
    bi1, bi2 = fs.ps_bi_solve( HP, phis, u, [sp1, sp2], phi_cri)
    print( u, bi1, bi2, 'bi done!', flush=True)

    return sp1, sp2, bi1, bi2



#==================== Free Energy Calculation ====================

utest=u

sp1,sp2,bi1,bi2=bisp_parallel(utest)
chem1=tt.chem_pot(HP,bi1,phis,utest)
chem2=tt.chem_pot(HP,bi2,phis,utest)



ener1=tt.feng(HP,bi1,phis,utest)
ener2=tt.feng(HP,bi2,phis,utest)


dener1=tt.dfeng(HP,bi1,phis,utest)
dener2=tt.dfeng(HP,bi2,phis,utest)

print("chem_pot ",chem1,chem2)
print("feng ",ener1+(1-bi1)*dener1,ener2+(1-bi2)*dener2 )



energy_file='energy.dat'

n_phi=100
phi_monomers1=np.linspace(0.8*bi1,1.2*bi1,50)
phi_monomers2=np.linspace(1.3*bi1,1.2*bi2,50)
phi_monomers=np.concatenate((phi_monomers1,phi_monomers2))

u=utest

f_en=[]
df_en=[]
chem_pot=[]
for i in range((n_phi)):
     phi_monomer=phi_monomers[i]
     free_energy=tt.feng(HP,phi_monomer,phis,u)
     dfree_energy=tt.dfeng(HP,phi_monomer,phis,u)
     chem_energy=tt.chem_pot(HP,phi_monomer,phis,u)

     f_en.append(free_energy)
     df_en.append(dfree_energy)
     chem_pot.append(chem_energy)
     print("energy",phi_monomer, free_energy,dfree_energy)

print(np.sum(HP['sig']))
print(tt.chi_calc(HP, u))
f_en=np.array(f_en)
df_en=np.array(df_en)

print(len(f_en),len(df_en))
np.savetxt(energy_file,np.transpose([phi_monomers,f_en,df_en,chem_pot]),fmt='%1.4e', newline='\n')
