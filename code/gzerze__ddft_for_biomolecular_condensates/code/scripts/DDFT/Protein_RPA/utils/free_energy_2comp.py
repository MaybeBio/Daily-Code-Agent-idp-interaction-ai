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
          eps_b = 84.51003476887941  ):
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
def RPAFH_12(sig1,sig2):
    sig1 = np.array(sig1) # charge pattern
    N1   = sig1.shape[0]  # sequence length
    pc1  = np.abs(np.sum(sig1))/N1 # prefactor for counterions
    sig2 = np.array(sig2) # charge pattern
    N2   = sig2.shape[0]  # sequence length
    pc2  = np.abs(np.sum(sig2))/N2 # prefactor for counterions
    if (len(sig2) > len(sig1)):
        zeroes = np.zeros((len(sig2)-len(sig1)))
        sig1=np.concatenate((sig1,zeroes),axis=0)
    elif (len(sig2) < len(sig1)):
        zeroes = np.zeros((len(sig1)-len(sig2)))
        sig2=np.concatenate((sig2,zeroes),axis=0)   
    Q   = np.sum(sig1*sig2)/np.sqrt((N1+N2)/2)     # fraction of charged residues (sig=+/-1)
    Gamma = np.linalg.norm((sig1-sig2),1)/np.mean((np.linalg.norm(sig1,1),np.linalg.norm(sig2,1)))
    # linear summation for S(k)
    mel = np.kron(sig1, sig2).reshape((max(N1,N2), max(N1,N2)))
    Tel1 = np.array([ np.sum(mel.diagonal(n) + mel.diagonal(-n)) for n in range(max(N1,N2))])
    Tel1[0] /= 2
    Tel2 = np.array([ np.sum(mel.diagonal(n) + mel.diagonal(-n)) for n in range(max(N1,N2))])
    Tel2[0] /= 2
    L1 = np.arange(max(N1,N2))
    L2 = np.arange(max(N1,N2))
    L=(L1+L2)/2

    HP12 = {    'sig1': sig1, \
                'N1'  : N1,   \
                'pc1' : pc1,  \
                'sig2': sig2, \
                'N2'  : N2,   \
                'pc2' : pc2,  \
                'Q'  : Q,   \
                'L1'  : L1,   \
                'Tel1': Tel1, \
                'L2'  : L2,   \
                'Tel2': Tel2, \
                'L': L,  \
                'Gamma': Gamma \
           }
    return HP12

from .. import seq_list as sl
seq_name1 = "Ddx4_N1"
seq_name2 = "Ddx4_N1"
#Below is just default values in case the ehs parameters are not defined in input files
ehs1=np.array([0.15, -0.3])
ehs2=np.array([0.15, -0.3])
ef=False

sig1, N1, the_seq1 = sl.get_the_charge(seq_name1)
sig2, N1, the_seq2 = sl.get_the_charge(seq_name2)

HP1 = RPAFH(sig1, ehs=ehs1, epsfun=ef)
HP2 = RPAFH(sig2, ehs=ehs2, epsfun=ef)
HP = [HP1,HP2]

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
 

    def Enp(self, phi1, phi2, phis):
        phic = self.HP[0]['pc']*phi1 + phis+self.HP[1]['pc']*phi2
        return 1/self.HP[0]['N']*self.s_calc(phi1)+1/self.HP[1]['N']*self.s_calc(phi2) + self.s_calc(phic) + self.s_calc(phis) \
            	+ self.s_calc(1-phi1*gv.r_res-phic*gv.r_con-phis*gv.r_sal-phi2*gv.r_res)

    def d1Enp(self, phi1, phi2, phis):
        phic = self.HP[0]['pc']*phi1 + phis + self.HP[1]['pc']*phi2
        return ( 1 + np.log(phi1) )/self.HP[0]['N'] +self.HP[0]['pc']*(1 + np.log(phic+(self.HP[0]['pc']==0)) ) \
             - (gv.r_res + gv.r_con*self.HP[0]['pc']+ gv.r_con*self.HP[1]['pc']) \
                *( 1 + np.log(1-phi1*gv.r_res-phic*gv.r_con-phis*gv.r_sal-phi2*gv.r_res) )
    
    def d2Enp(self, phi1, phi2, phis):
        phic = self.HP[0]['pc']*phi1 + phis + self.HP[1]['pc']*phi2
        return ( 1 + np.log(phi2) )/self.HP[1]['N'] +self.HP[1]['pc']*(1 + np.log(phic+(self.HP[1]['pc']==0)) ) \
             - (gv.r_res + gv.r_con*self.HP[0]['pc']+ gv.r_con*self.HP[1]['pc']) \
                *( 1 + np.log(1-phi1*gv.r_res-phic*gv.r_con-phis*gv.r_sal-phi2*gv.r_res) )

    def dd1Enp(self, phi1, phi2, phis):
        phic =self.HP[0]['pc']*phi1 + phis + self.HP[1]['pc']*phi2
        return 1/self.HP[0]['N']/phi1 + self.HP[0]['pc']*self.HP[0]['pc']/(phic + (self.HP[0]['pc']==0))*(np.abs(self.HP[0]['pc'])>0) \
           + (gv.r_res + gv.r_con*self.HP[0]['pc'])*(gv.r_res + gv.r_con*self.HP[0]['pc']) \
              /(1-phi1*gv.r_res-phic*gv.r_con-phis*gv.r_sal-phi2*gv.r_res)
    
    def dd2Enp(self, phi1, phi2, phis):
        phic =self.HP[0]['pc']*phi2 + phis + self.HP[1]['pc']*phi2
        return 1/self.HP[1]['N']/phi2 + self.HP[1]['pc']*self.HP[0]['pc']/(phic + (self.HP[1]['pc']==0))*(np.abs(self.HP[1]['pc'])>0) \
           + (gv.r_res + gv.r_con*self.HP[1]['pc'])*(gv.r_res + gv.r_con*self.HP[1]['pc']) \
              /(1-phi1*gv.r_res-phic*gv.r_con-phis*gv.r_sal-phi2*gv.r_res)    

#  Electrostatic Part for single polymer component 1

    def U1el(self,k,u):
        return 4*np.pi*u/(k*k*(1+self.Gamma*k*k))*np.exp(-self.c_smear*k*k)

    def S1k(self,k):
        return np.mean(self.HP[0]['Tel']*np.exp(-k*k*self.HP[0]['L']/6) )
    
    def fel1(self,phi1,phis,u):
        f1 = sci.quad(self.fel1_toint1, 0, np.inf, args=(phi1,phis,u), limit=self.intlim)[0]
        f2 = sci.quad(self.fel1_toint2, 0, np.inf, args=(phi1,phis,u), limit=self.intlim)[0]
        return f1+f2

    def fel1_toint1(self,k, phi1, phis, u):
        sk = self.S1k(k)
        lk  = self.U1el(k,u)
        epx = self.HP[0]['epsx'](phi1*gv.r_res)

        G1 = lk*epx*( 2*phis + phi1*( self.HP[0]['pc'] + sk )  )

        return 1/(4*np.pi*np.pi)*k*k*( 1/gv.eta*np.log(1+gv.eta*G1) - G1 + G1*G1/2 )

    def fel1_toint2(self,k,phi1, phis, u):
        sk = self.S1k(k)
        lk  = self.U1el(k,u)
        epx = self.HP[0]['epsx'](phi1*gv.r_res)

        G1 = lk*epx*( 2*phis + phi1*( self.HP[0]['pc'] + sk ))
        G2 = lk*epx*( 2*phis + phi1*( self.HP[0]['pc'] +self. HP[0]['Q'] ) )

        return 1/(4*np.pi*np.pi)*k*k*( G1-G1*G1/2-self.NoSelfEnergy*G2)


    def d1fel(self,phi1,phis,u):
        return sci.quad(self.d1fel_toint, 0, np.inf, args=(phi1,phis,u), limit=self.intlim )[0]

    def d1fel_toint(self,k, phi1, phis, u):
        sk   = self.S1k(k)
        lk   = self.U1el(k,u)
        epx  = self.HP[0]['epsx'](phi1*gv.r_res)
        depx = gv.r_res*self.HP[0]['depsx'](phi1*gv.r_res)

        SSk = self.HP[0]['pc'] + sk
        G1 = lk*epx*( 2*phis + phi1*SSk  )

        A1 = ( depx*(2*phis + phi1*SSk ) + epx*SSk )/(1+gv.eta*G1)
        A2 = depx*( 2*phis + phi1*( self.HP[0]['pc'] + self.HP[0]['Q'] ) ) + epx*( self.HP[0]['pc'] + self.HP[0]['Q'] )

        return 1/(4*np.pi*np.pi)*k*k*lk*( A1- self.NoSelfEnergy*A2)

    def dd1fel(self,phi1,phis,u):
        return sci.quad(self.dd1fel_toint, 0, np.inf, args=(phi1,phis,u),limit=self.intlim )[0]

    def dd1fel_toint(self,k, phi1, phis, u):
        sk   = self.S1k(k)
        lk   = self.U1el(k,u)
        epx  = self.HP[0]['epsx'](phi1*gv.r_res)
        depx = gv.r_res*self.HP[0]['depsx'](phi1*gv.r_res)
        ddepx = gv.r_res*gv.r_res*self.HP[0]['ddepsx'](phi1*gv.r_res)

        SSk  = self.HP[0]['pc'] + sk
        G1   = lk*epx*( 2*phis + phi1*SSk  )
        dG1  = lk*( depx*(2*phis + phi1*SSk ) + epx*SSk )
        ddG1 = lk*( ddepx*(2*phis+phi1*SSk) + 2*depx*SSk)

        A11 = dG1/(1+gv.eta*G1)
        A12 = ddG1/(1+gv.eta*G1)

        A2  = ddepx*( 2*phis + phi1*( self.HP[0]['pc'] + self.HP[0]['Q'] ) ) \
            + 2*depx*( self.HP[0]['pc'] + self.HP[0]['Q'] )

        return 1/(4*np.pi*np.pi)*k*k*( -gv.eta*A11*A11 + A12 - self.NoSelfEnergy*lk*A2 )
    
    #  Electrostatic Part for single polymer component 2


    def U2el(self,k,u):
        return 4*np.pi*u/(k*k*(1+self.Gamma*k*k))*np.exp(-self.c_smear*k*k)

    def S2k(self,k):
        return np.mean( self.HP[1]['Tel']*np.exp(-k*k*self.HP[1]['L']/6) )
    
    def fel2(self,phi2,phis,u):
        f1 = sci.quad(self.fel2_toint1, 0, np.inf, args=(phi2,phis,u), limit=self.intlim)[0]
        f2 = sci.quad(self.fel2_toint2, 0, np.inf, args=(phi2,phis,u), limit=self.intlim)[0]
        return f1+f2

    def fel2_toint1(self,k, phi2, phis, u):
        sk = self.S2k(k)
        lk  = self.U2el(k,u)
        epx = self.HP[1]['epsx'](phi2*gv.r_res)

        G1 = lk*epx*( 2*phis + phi2*( self.HP[1]['pc'] + sk )  )

        return 1/(4*np.pi*np.pi)*k*k*( 1/gv.eta*np.log(1+gv.eta*G1) - G1 + G1*G1/2 )

    def fel2_toint2(self,k,phi2, phis, u):
        sk = self.S2k(k)
        lk  = self.U2el(k,u)
        epx = self.HP[1]['epsx'](phi2*gv.r_res)

        G1 = lk*epx*( 2*phis + phi2*( self.HP[1]['pc'] + sk ))
        G2 = lk*epx*( 2*phis + phi2*( self.HP[1]['pc'] +self. HP[1]['Q'] ) )

        return 1/(4*np.pi*np.pi)*k*k*( G1-G1*G1/2-self.NoSelfEnergy*G2)


    def d2fel(self,phi2,phis,u):
        return sci.quad(self.d2fel_toint, 0, np.inf, args=(phi2,phis,u), limit=self.intlim )[0]

    def d2fel_toint(self,k, phi2, phis, u):
        sk   = self.Sk(k)
        lk   = self.Uel(k,u)
        epx  = self.HP[1]['epsx'](phi2*gv.r_res)
        depx = gv.r_res*self.HP[1]['depsx'](phi2*gv.r_res)

        SSk = self.HP[1]['pc'] + sk
        G1 = lk*epx*( 2*phis + phi2*SSk  )

        A1 = ( depx*(2*phis + phi2*SSk ) + epx*SSk )/(1+gv.eta*G1)
        A2 = depx*( 2*phis + phi2*( self.HP[1]['pc'] + self.HP[1]['Q'] ) ) + epx*( self.HP[1]['pc'] + self.HP[1]['Q'] )

        return 1/(4*np.pi*np.pi)*k*k*lk*( A1- self.NoSelfEnergy*A2)

    def dd2fel(self,phi2,phis,u):
        return sci.quad(self.dd2fel_toint, 0, np.inf, args=(phi2,phis,u),limit=self.intlim )[0]

    def dd2fel_toint(self,k, phi2, phis, u):
        sk   = self.S2k(k)
        lk   = self.U2el(k,u)
        epx  = self.HP[1]['epsx'](phi2*gv.r_res)
        depx = gv.r_res*self.HP[1]['depsx'](phi2*gv.r_res)
        ddepx = gv.r_res*gv.r_res*self.HP[0]['ddepsx'](phi2*gv.r_res)

        SSk  = self.HP[0]['pc'] + sk
        G1   = lk*epx*( 2*phis + phi2*SSk  )
        dG1  = lk*( depx*(2*phis + phi2*SSk ) + epx*SSk )
        ddG1 = lk*( ddepx*(2*phis+phi2*SSk) + 2*depx*SSk)

        A11 = dG1/(1+gv.eta*G1)
        A12 = ddG1/(1+gv.eta*G1)

        A2  = ddepx*( 2*phis + phi2*( self.HP[1]['pc'] + self.HP[1]['Q'] ) ) \
            + 2*depx*( self.HP[1]['pc'] + self.HP[1]['Q'] )

        return 1/(4*np.pi*np.pi)*k*k*( -gv.eta*A11*A11 + A12 - self.NoSelfEnergy*lk*A2 )

    
    #1-2 electrostatic interactions, modeled approximately as the same energy contribution as intramolecular interactions

    def U12el(self,k,u):
        return 4*np.pi*u/(k*k*(1+self.Gamma*k*k))*np.exp(-self.c_smear*k*k)

    def S12k(self,k):
        return np.mean( [self.HP[2]['Tel1']*np.exp(-k*k*self.HP[2]['L1']/6),self.HP[2]['Tel2']*np.exp(-k*k*self.HP[2]['L2']/6)] )

    def fel12(self,phi1,phi2,u):
        f1 = sci.quad(self.fel12_toint1, 0, np.inf, args=(phi1,phi2,u), limit=self.intlim)[0]
        f2 = sci.quad(self.fel12_toint2, 0, np.inf, args=(phi1,phi2,u), limit=self.intlim)[0]
        return f1+f2

    def fel12_toint1(self,k, phi1, phi2, u):
        sk = self.S12k(k)
        lk  = self.U12el(k,u)
        epx = ((phi1+phi2)*gv.r_res)

        #Gamma = np.linalg.norm((HP[2]['sig1']-HP[2]['sig2']),1)/np.mean((np.linalg.norm(HP[2]['sig1'],1),np.linalg.norm(HP[2]['sig2'],1)))

        G1 = self.HP[2]['Gamma']*lk*epx*(phi1*( self.HP[0]['pc']) + phi2*( self.HP[1]['pc']))  

        return 1/(4*np.pi*np.pi)*k*k*( 1/gv.eta*np.log(1+gv.eta*G1) - G1 + G1*G1/2 )

    def fel12_toint2(self,k,phi1, phi2, u):
        sk = self.S12k(k)
        lk  = self.U12el(k,u)
        epx = (phi1*gv.r_res)*phi2
        
        #Gamma = np.linalg.norm((HP[2]['sig1']-HP[2]['sig2']),1)/np.mean((np.linalg.norm(HP[2]['sig1'],1),np.linalg.norm(HP[2]['sig2'],1)))

        G1 = self.HP[2]['Gamma']*lk*epx*(phi1*self.HP[0]['pc'] + phi2*self.HP[1]['pc'] + sk )
        G2 = lk*epx*(phi1*self.HP[0]['pc'] +phi2*self.HP[1]['pc']+ self.HP[2]['Q'] )

        return 1/(4*np.pi*np.pi)*k*k*( G1-G1*G1/2-self.NoSelfEnergy*G2)

    #f_FH functions
    
    def chi_calc1(self,u):
        ehs = self.HP[0]['ehs']
        return gv.r_res*gv.r_res*( ehs1[0]*u + ehs1[1] )
    
    def chi_calc2(self,u):
        ehs = self.HP[1]['ehs']
        return gv.r_res*gv.r_res*( ehs2[0]*u + ehs2[1] )  

    # free energy functions
    def feng(self, phi1, phi2, phis, u):
        return self.Enp( phi1, phi2,phis) + self.fel1(phi1,phis,u)+ self.fel2(phi2,phis,u) -self.fel12(phi1,phi2,u)+gv.r_res**(-1)*self.chi_calc1(u)*phi1*(1-gv.r_res*phi1-gv.r_res*phi2)+gv.r_res**(-1)*self.chi_calc2(u)*phi2*(1-gv.r_res*phi1-gv.r_res*phi2)

    def d1feng(self,phi1, phi2, phis, u):
        return self.d1Enp(phi1,phi2,phis) + self.dfel(phi2,phis,u) + gv.r_res**(-1)*self.chi_calc1(u)*(1-2*gv.r_res*phi1)- gv.r_res**(-1)*self.chi_calc2(u)*(gv.r_res*phi1*phi2)
    
    def d2feng(self,phi1, phi2, phis, u):
        return self.d2Enp(phi1,phi2,phis) + self.dfel(phi2,phis,u) + gv.r_res**(-1)*self.chi_calc2(u)*(1-2*gv.r_res*phi2)- gv.r_res**(-1)*self.chi_calc1(u)*(gv.r_res*phi1*phi2)

    def dd1feng(self, phi1, phi2, phis, u):
        return self.dd1Enp(phi1,phi2,phis) + self.dd1fel(phi1,phis,u) - 2*self.chi_calc1(u) - gv.r_res**(-1)*self.chi_calc2(u)*(gv.r_res*phi2)
    
    def dd2feng(self, phi1, phi2, phis, u):
        return self.dd1Enp(phi1,phi2,phis) + self.dd1fel(phi2,phis,u) - 2*self.chi_calc2(u) - gv.r_res**(-1)*self.chi_calc1(u)*(gv.r_res*phi1)

    def chem_pot1(self,phi1,phi2,phis,u):
        ri=gv.r_res
        ti=1.0
        return ri*self.feng(phi1,phi2,phis,u) + ti*self.d1feng(phi1,phi2,phis,u)-ri*phi1*self.d1feng(phi1,phi2,phis,u)
    
    def chem_pot2(self,phi1,phi2,phis,u):
        ri=gv.r_res
        ti=1.0
        return ri*self.feng(phi1,phi2,phis,u) + ti*self.d2feng(phi1,phi2,phis,u)-ri*phi2*self.d2feng(phi1,phi2,phis,u)

    def d_chem_pot2(self,phi1,phi2,phis,u):
        ri=gv.r_res
        ti=1.0
        return ri*self.d2feng(phi1,phi2,phis,u) + ti*self.dd2feng(phi1,phi2,phis,u)-ri*phi2*self.dd2feng(phi1,phi2,phis,u)-ri*self.d2feng(phi1,phi2,phis,u)
    
    def d_chem_pot1(self,phi1,phi2,phis,u):
        ri=gv.r_res
        ti=1.0
        return ri*self.d1feng(phi1,phi2,phis,u) + ti*self.dd1feng(phi1,phi2,phis,u)-ri*phi1*self.dd1feng(phi1,phi2,phis,u)-ri*self.d1feng(phi1,phi2,phis,u)


   


