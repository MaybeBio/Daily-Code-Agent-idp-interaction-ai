import numpy as np
import Protein_RPA.utils.free_energy as f_en
from Protein_RPA.utils.input_parse import *
from initialize import Initialize
import Protein_RPA.utils.global_vars as gv
# Below is the function used to interpolate chemical potential at particular density (phi)
# In this version,salt density is set zero
def lookup_table(FE,u,phis):
    
    from scipy.interpolate import CubicSpline
    n1=200
    n2=200
    n3=200
    n4=200
    n5=200
    n6=200
    N=n1+n2+n3+n4+n5+n6
    phi_table1=np.linspace(1e-6,0.001,n1)/gv.r_res
    phi_table2=np.linspace(0.001+1e-4,0.01,n2)/gv.r_res
    phi_table3=np.linspace(0.01+1e-4,0.05,n3)/gv.r_res
    phi_table4=np.linspace(0.05+1e-4,0.2,n4)/gv.r_res
    phi_table5=np.linspace(0.2+1e-4,0.5,n5)/gv.r_res
    phi_table6=np.linspace(0.5+1e-4,0.85,n6)/gv.r_res
   
    args=(phi_table1,phi_table2,phi_table3,phi_table4,phi_table5,phi_table6)
    phi_table=np.concatenate(args)

    mu_table=[]
    for i in range(N):
        mu_table.append(FE.chem_pot(phi_table[i],phis,u))
    mu_table=np.array(mu_table)


    #for i in range(N):
    #    print(i,phi_table[i],mu_table[i])


    mu_spl = CubicSpline(phi_table,mu_table)

    return mu_spl




# Below is the function used to interpolate derivative of chemical potential w.r.t. phi at particular density (phi)
# This derivative is used in the implementation of backward/implicit euler
# In this version,salt density is set zero
def dmu_table(FE,u,phis):
    
    from scipy.interpolate import CubicSpline
    n1=200
    n2=200
    n3=200
    n4=200
    n5=200
    n6=200
    N=n1+n2+n3+n4+n5+n6
    phi_table1=np.linspace(1e-6,0.001,n1)/gv.r_res
    phi_table2=np.linspace(0.001+1e-4,0.01,n2)/gv.r_res
    phi_table3=np.linspace(0.01+1e-4,0.05,n3)/gv.r_res
    phi_table4=np.linspace(0.05+1e-4,0.2,n4)/gv.r_res
    phi_table5=np.linspace(0.2+1e-4,0.5,n5)/gv.r_res
    phi_table6=np.linspace(0.5+1e-4,0.85,n6)/gv.r_res
   
    args=(phi_table1,phi_table2,phi_table3,phi_table4,phi_table5,phi_table6)
    phi_table=np.concatenate(args)

    dmu_table=[]
    for i in range(N):
        dmu_table.append(FE.d_chem_pot(phi_table[i],phis,u))
    dmu_table=np.array(dmu_table)


    #for i in range(N):
    #    print(i,phi_table[i],dmu_table[i])


    dmu_spl = CubicSpline(phi_table,dmu_table)

    return dmu_spl











