#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Test function: test the utilization of free energy class and its functions
"""

from __future__ import print_function
import fipy as fp
import os
#from fipy.solvers.pysparse import LinearLUSolver as Solver
import numpy as np
import matplotlib.pyplot as plt

import utils.free_energy as f_en
import seq_list as sl

seq_name = "Ddx4_N1"
ehs=np.array([0.15, -0.3])
ef=False

sig, N, the_seq = sl.get_the_charge(seq_name)
HP=f_en.RPAFH(sig, ehs=ehs, epsfun=ef)
FE=f_en.free_energy(HP)

print(HP['ehs'])
#print(HP['pc'])
print(FE.Enp(0.1,0))
print(FE.feng(0.1,0,0.1))
print(FE.dfeng(0.1,0,0.1))
print(FE.ddfeng(0.1,0,0.1))

print(FE.chem_pot(0.001,0,4))





energy_file='energy2.dat'

bi1=0.000022
bi2=0.149
utest=4
phis=0

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
     free_energy=FE.feng(phi_monomer,phis,u)
     dfree_energy=FE.dfeng(phi_monomer,phis,u)
     chem_energy=FE.chem_pot(phi_monomer,phis,u)

     f_en.append(free_energy)
     df_en.append(dfree_energy)
     chem_pot.append(chem_energy)
     print("energy",phi_monomer, free_energy,dfree_energy)


f_en=np.array(f_en)
df_en=np.array(df_en)

print(len(f_en),len(df_en))
np.savetxt(energy_file,np.transpose([phi_monomers,f_en,df_en,chem_pot]),fmt='%1.4e', newline='\n')












