import numpy as np
import Protein_RPA.utils.free_energy as f_en
import Protein_RPA.utils.free_energy_2comp as f_en2
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



    dmu_spl = CubicSpline(phi_table,dmu_table)

    return dmu_spl


#adding the spline function for chemical potential w.r.t. phi1 and phi2
def lookup_table2(FE,u,phis):
    
    from scipy.interpolate import CloughTocher2DInterpolator
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
   
    args=(phi_table1,phi1_table,phi1_table,phi1_table,phi1_table,phi1_table)
    phi1_table=np.concatenate(args)
    phi2_table=np.concatenate(args)


    mu1_table=np.empty([N,N])
    for i in range(N):
        for j in range(N):
            mu1_table[i,j]=FE.chem_pot1(phi1_table[i],phi2_table[j],phis,u)
    
    mu2_table=np.empty([N,N])
    for i in range(N):
        for j in range(N):
            mu2_table[i,j]=FE.chem_pot2(phi2_table[i],phi2_table[j],phis,u)



    mu1_spl = CloughTocher2DInterpolator(list(zip(phi1_table,phi2_table)),mu1_table)
    mu2_spl = CloughTocher2DInterpolator(list(zip(phi1_table,phi2_table)),mu2_table)

    return mu1_spl, mu2_spl

# Below is the function used to interpolate derivative of chemical potential w.r.t. phi1 and phi2 at particular density (phi)
# This derivative is used in the implementation of backward/implicit euler
# In this version,salt density is set zero
def dmu_table2(FE,u,phis):
    
    from scipy.interpolate import CubicSpline
    from scipy.interpolate import CloughTocher2DInterpolator
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
    phi1_table=np.concatenate(args)
    phi2_table=np.concatenate(args)

    dmu1_table=np.empty([N,N])
    for i in range(N):
        for j in range(N):
            dmu1_table[i,j]=FE.d_chem_pot(phi1_table[i],phis,u)



    dmu1_spl = CloughTocher2DInterpolator(list(zip(phi1_table,phi2_table)),dmu_table)




    return dmu1_spl






