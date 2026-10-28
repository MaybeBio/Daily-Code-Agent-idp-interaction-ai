#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Modules that are used to define the:


1.  Free Energy of protein in the presence of counterions, water, and salt
    This energy includes three terms:
    a. Entropy term
    b. Electrostatic energy
    c. Flory Huggins energy

2.  

"""
import numpy as np
#from . import global_vars as gv
import Protein_RPA.utils.global_vars as gv
import scipy.integrate as sci


#HP = RPAFH(sig, ehs=ehs, epsfun=ef)

def RPAFH(sig, ehs=[0,0], epsfun=False, eps_modify_ehs=False,
          eps_a = 18.931087269965023,
          eps_b = 84.51003476887941):
    #sequence parameters
    sig = np.array(sig) # charge pattern
    N   = sig.shape[0]  # sequence length
    pc  = np.abs(np.sum(sig))/N # prefactor for counterions
    Q   = np.sum(sig*sig)/N     # fraction of charged residues (sig=+/-1)

    # linear summation for S(k)
    mel = np.kron(sig, sig).reshape((N, N))
    Tel = np.array([ np.sum(mel.diagonal(n) + mel.diagonal(-n)) for n in range(N)])
    Tel[0] /= 2
    L = np.arange(N)

    HP = {  'sig': sig, \
            'N'  : N,   \
            'pc' : pc,  \
            'Q'  : Q,   \
            'L'  : L,   \
            'Tel': Tel, \
            'ehs': ehs  \
           }

    if epsfun:
        a, b = eps_a, eps_b
        HP['eps0'] = b
        flinear = lambda x: a*x + b*(1-x)
        HP['epsx']   = lambda x: b/flinear(x)
        HP['depsx']  = lambda x: -b*(a-b)/(flinear(x))**2
        HP['ddepsx'] = lambda x: 2*b*(a-b)*(a-b)/(flinear(x))**3
    else:
        HP['eps0']   = 1
        HP['epsx']   = lambda x: 1*(x==x)
        HP['depsx']  = lambda x: 0*x
        HP['ddepsx'] = lambda x: 0*x

    # using eps_r=eps0 to rescale ehs (assuming the input ehs is of eps_r=1)
    if eps_modify_ehs:
        ehs[0] *= HP['eps0']

    return HP



from .. import seq_list as sl
seq_name = "Ddx4_N1"
#Below is just default values in case the ehs parameters are not defined in input files
ehs=np.array([0.15, -0.3])
ef=False

sig, N, the_seq = sl.get_the_charge(seq_name)

HP = RPAFH(sig, ehs=ehs, epsfun=ef)

class free_energy():
    """
    
    Implementation of fre energy

    """

    def __init__(self,HP):
        self.c_smear=0
        self.HP=HP

        self.Gamma = 1    # Short-range cutoff factor
        self.intlim= 200
        self.NoSelfEnergy = True


#  Entropy part (needs to be rewritten for multiple components)

    def s_calc(self,x):
        return (x > gv.phi_min_sys )*x*np.log(x+(x<gv.phi_min_sys))

    def Enp(self,phi, phis):
        phic = self.HP['pc']*phi + phis
        return 1/self.HP['N']*self.s_calc(phi) + self.s_calc(phic) + self.s_calc(phis) \
            	+ self.s_calc(1-phi*gv.r_res-phic*gv.r_con-phis*gv.r_sal)

    def dEnp(self, phi, phis):
        phic = self.HP['pc']*phi + phis
        return ( 1 + np.log(phi) )/self.HP['N'] +self.HP['pc']*(1 + np.log(phic+(self.HP['pc']==0)) ) \
             - (gv.r_res + gv.r_con*self.HP['pc']) \
                *( 1 + np.log(1-phi*gv.r_res-phic*gv.r_con-phis*gv.r_sal) )

    def ddEnp(self, phi, phis):
        phic =self.HP['pc']*phi + phis
        return 1/self.HP['N']/phi + self.HP['pc']*self.HP['pc']/(phic + (self.HP['pc']==0))*(np.abs(self.HP['pc'])>0) \
           + (gv.r_res + gv.r_con*self.HP['pc'])*(gv.r_res + gv.r_con*self.HP['pc']) \
              /(1-phi*gv.r_res-phic*gv.r_con-phis*gv.r_sal)

#  Electrostatic Part

    def Uel(self,k,u):
        return 4*np.pi*u/(k*k*(1+self.Gamma*k*k))*np.exp(-self.c_smear*k*k)

    def Sk(self,k):
        return np.mean( self.HP['Tel']*np.exp(-k*k*self.HP['L']/6) )
    
    def fel(self,phi,phis,u):
        f1 = sci.quad(self.fel_toint1, 0, np.inf, args=(phi,phis,u), limit=self.intlim)[0]
        f2 = sci.quad(self.fel_toint2, 0, np.inf, args=(phi,phis,u), limit=self.intlim)[0]
        return f1+f2

    def fel_toint1(self,k, phi, phis, u):
        sk = self.Sk(k)
        lk  = self.Uel(k,u)
        epx = self.HP['epsx'](phi*gv.r_res)

        G1 = lk*epx*( 2*phis + phi*( self.HP['pc'] + sk )  )

        return 1/(4*np.pi*np.pi)*k*k*( 1/gv.eta*np.log(1+gv.eta*G1) - G1 + G1*G1/2 )

    def fel_toint2(self,k,phi, phis, u):
        sk = self.Sk(k)
        lk  = self.Uel(k,u)
        epx = self.HP['epsx'](phi*gv.r_res)

        G1 = lk*epx*( 2*phis + phi*( self.HP['pc'] + sk )  )
        G2 = lk*epx*( 2*phis + phi*( self.HP['pc'] + self. HP['Q'] ) )

        return 1/(4*np.pi*np.pi)*k*k*( G1-G1*G1/2-self.NoSelfEnergy*G2)


    def dfel(self,phi,phis,u):
        return sci.quad(self.dfel_toint, 0, np.inf, args=(phi,phis,u), limit=self.intlim )[0]

    def dfel_toint(self,k, phi, phis, u):
        sk   = self.Sk(k)
        lk   = self.Uel(k,u)
        epx  = self.HP['epsx'](phi*gv.r_res)
        depx = gv.r_res*self.HP['depsx'](phi*gv.r_res)

        SSk = self.HP['pc'] + sk
        G1 = lk*epx*( 2*phis + phi*SSk  )

        A1 = ( depx*(2*phis + phi*SSk ) + epx*SSk )/(1+gv.eta*G1)
        A2 = depx*( 2*phis + phi*( self.HP['pc'] + self.HP['Q'] ) ) + epx*( self.HP['pc'] + self.HP['Q'] )

        return 1/(4*np.pi*np.pi)*k*k*lk*( A1 - self.NoSelfEnergy*A2)

    def ddfel(self,phi,phis,u):
        return sci.quad(self.ddfel_toint, 0, np.inf, args=(phi,phis,u),limit=self.intlim )[0]

    def ddfel_toint(self,k, phi, phis, u):
        sk   = self.Sk(k)
        lk   = self.Uel(k,u)
        epx  = self.HP['epsx'](phi*gv.r_res)
        depx = gv.r_res*self.HP['depsx'](phi*gv.r_res)
        ddepx = gv.r_res*gv.r_res*self.HP['ddepsx'](phi*gv.r_res)

        SSk  = self.HP['pc'] + sk
        G1   = lk*epx*( 2*phis + phi*SSk  )
        dG1  = lk*( depx*(2*phis + phi*SSk ) + epx*SSk )
        ddG1 = lk*( ddepx*(2*phis+phi*SSk) + 2*depx*SSk)

        A11 = dG1/(1+gv.eta*G1)
        A12 = ddG1/(1+gv.eta*G1)

        A2  = ddepx*( 2*phis + phi*( self.HP['pc'] + self.HP['Q'] ) ) \
            + 2*depx*( self.HP['pc'] + self.HP['Q'] )

        return 1/(4*np.pi*np.pi)*k*k*( -gv.eta*A11*A11 + A12 - self.NoSelfEnergy*lk*A2 )

    #f_FH functions
    
    def chi_calc(self,u):
        ehs = self.HP['ehs']
        return gv.r_res*gv.r_res*( ehs[0]*u + ehs[1] )

    # free energy functions
    def feng(self, phi, phis, u):
        return self.Enp( phi, phis) + self.fel(phi,phis,u) +gv.r_res**(-1)*self.chi_calc(u)*phi*(1-gv.r_res*phi)

    def dfeng(self,phi, phis, u):
        return self.dEnp(phi, phis) + self.dfel( phi,phis,u) + gv.r_res**(-1)*self.chi_calc(u)*(1-2*gv.r_res*phi)

    def ddfeng(self, phi, phis, u):
        return self.ddEnp(phi, phis) + self.ddfel(phi,phis,u) - 2*self.chi_calc( u)

    def chem_pot(self,phi,phis,u):
        ri=gv.r_res
        ti=1.0
        return ri*self.feng(phi,phis,u) + ti*self.dfeng(phi,phis,u)-ri*phi*self.dfeng(phi,phis,u)


    def d_chem_pot(self,phi,phis,u):
        ri=gv.r_res
        ti=1.0
        return ri*self.dfeng(phi,phis,u) + ti*self.ddfeng(phi,phis,u)-ri*phi*self.ddfeng(phi,phis,u)-ri*self.dfeng(phi,phis,u)


   


