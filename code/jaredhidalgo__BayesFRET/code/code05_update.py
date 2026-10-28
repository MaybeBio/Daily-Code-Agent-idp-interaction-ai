"""
BayesFRET: File 6
===
An Experiment-Adjusted HDP-HMM to Analyze Surface-Immobilized smFRET Data

About
---
`code05_update.py`: The file to handle each MCMC iteration.

Function Hierarchy
---
* `sampler_update()`
    * `UPDATE_STATES()`
    * `UPDATE_PHOTOPHYSICS()`
    * `UPDATE_TRANSITION_PROBS()`
    * `UPDATE_EMISSIONS()`
        * `HMC()`
            * `get_u_grad()`
            * `get_u()`
        * `MG()`
        * `FLIPPER()`
    * `UPDATE_AUX()` is part of the `Sample` class in `code01_classes.py` since all parameters are already there.

Author
---
Code by Jared Hidalgo. 

Inspired by MATLAB code from Ioannis Sgouralis, Shreya Madaan, Franky Djutanta, Rachael Kha, Rizal F. Hariadi, and Steve Pressé for "A Bayesian Nonparametric Approach to Single Molecule Förster Resonance Energy Transfer".
"""
# Import internal packages for use.
from copy import deepcopy
# Import external packages for use.
import numpy as np
from numpy.matlib import repmat
# Import program files for documentation.
from code01_classes import Universal, Sample, Params






def sampler_update(sample: Sample, params: Params, U: Universal, v: np.int64):
    """
    The parent function to update the sample data.
    
    Function Hierarchy
    ---
    * `UPDATE_STATES()`
    * `UPDATE_PHOTOPHYSICS()`
    * `UPDATE_TRANSITION_PROBS()`
    * `UPDATE_EMISSIONS()`
        * `HMC()`
            * `get_u_grad()`
            * `get_u()`
        * `MG()`
        * `FLIPPER()`
    * `UPDATE_AUX()` is part of the `Sample` class in `code01_classes.py` since all parameters are already there.
    """
    sample.n += 1
    UPDATE_STATES( sample, params, U, v )
    UPDATE_PHOTOPHYSICS( sample, params, U, v )
    UPDATE_TRANSITION_PROBS( sample, params, U, v )
    UPDATE_EMISSIONS( sample, params,  U, v )
    sample.UPDATE_AUX()






def UPDATE_STATES(sample: Sample, params: Params, U: Universal, v: np.int64):
    """
    Update #1
    ===
    Updates photo-states of each dye (`ft_D`, `ft_A`) and calculated conformational states (`st`) throughout the experiment.

    Args:
        v (int): Index of the current simulation.
    """
    # FORM TRANSITION PROBS FOR DYES ONLY
    p = np.empty((4,4))
    p[np.ix_([0, 2], [0, 2])] =                              1 - sample.wi_D[0]      # D: 0 -> 0
    p[np.ix_([0, 2], [1, 3])] =                                  sample.wi_D[0]      # D: 0 -> 1
    p[np.ix_([1, 3], [0, 2])] =                              1 - sample.wi_D[1]      # D: 1 -> 0
    p[np.ix_([1, 3], [1, 3])] =                                  sample.wi_D[1]      # D: 1 -> 1
    p[np.ix_([0, 1], [0, 1])] = p[np.ix_([0, 1], [0, 1])] * (1 - sample.wi_A[0])     # A: 0 -> 0
    p[np.ix_([0, 1], [2, 3])] = p[np.ix_([0, 1], [2, 3])] *      sample.wi_A[0]      # A: 0 -> 1
    p[np.ix_([2, 3], [0, 1])] = p[np.ix_([2, 3], [0, 1])] * (1 - sample.wi_A[1])     # A: 1 -> 0
    p[np.ix_([2, 3], [2, 3])] = p[np.ix_([2, 3], [2, 3])] *      sample.wi_A[1]      # A: 1 -> 1
    
    o = np.empty(4)
    o[[0, 2]] =                                              1 - sample.wi_D[2]      # D: * -> 0
    o[[1, 3]] =                                                  sample.wi_D[2]      # D: * -> 1
    o[[0, 1]] =                                 o[[0, 1]] * (1 - sample.wi_A[2])     # A: * -> 0
    o[[2, 3]] =                                 o[[2, 3]] *      sample.wi_A[2]      # A: * -> 1
    
    # FORM TRANSITION PROBS IN COMBINED STATE-SPACE
    P = np.empty((4*params.K_lim, 4*params.K_lim))
    O = np.empty( 4*params.K_lim )
    for k in range(params.K_lim):
        for m in range(params.K_lim):
            P[ (4*m):(4*m+4) , (4*k):(4*k+4) ] = p * sample.pm[m,k]
        O[ (4*k):(4*k+4) ] = o * sample.ps[k]
    
    # COMPUTE LIKELIHOOD IN COMBINED STATE-SPACE
    A_forw = np.empty( (4*params.K_lim, params.T) )
    f_D    = np.array( [[0],[1],[0],[1]] )
    f_A    = np.array( [[0],[0],[1],[1]] )
    for k in range(params.K_lim):
        tmp1 = f_D * (f_A * sample.lam_D[k] + ( 1-f_A ) * sample.lam_Z)
        tmp2 = f_D *  f_A * sample.lam_A[k]
        mu_D = params.dD * params.qD * ( sample.ksi_D + params.cDD * tmp1 + params.cAD * tmp2 )
        mu_A = params.dD * params.qA * ( sample.ksi_A + params.cDA * tmp1 + params.cAA * tmp2 )
        A_forw[ (4*k):(4*k+4) ] = ( params.It_D * np.log(mu_D) + params.It_A * np.log(mu_A) ) - ( mu_D + mu_A )

    A_forw = np.maximum( np.exp( A_forw - A_forw.max(0) ), np.finfo(float).tiny )
    
    # FILTER FORWARD IN COMBINED STATE-SPACE
    tmp         = A_forw[:,0] * O
    A_forw[:,0] = tmp / np.sum( tmp )
    for t in range(1, params.T):
        tmp         = A_forw[:,t] * ( np.transpose(P) @ A_forw[:,t-1] )
        A_forw[:,t] = tmp / np.sum( tmp )
    
    # SAMPLE BACKWARD IN COMBINED STATE-SPACE
    ct = np.empty( params.T, dtype = np.int64 )
    ct[-1] = U.get_samples( A_forw[:, -1], v )
    for t in range(params.T-2, -1, -1):
        ct[t] = U.get_samples( A_forw[:, t] * P[:, ct[t+1]], v )
    ct += 1 # NOTE: Adjust to mimic original code
    
    # RECOVER ORIGINAL STATE-SPACE
    sample.ft_D = np.mod(            ct+1, 2 )
    sample.ft_A = np.mod(   np.ceil( ct/2-1 ), 2 )
    sample.st   = np.int64( np.ceil( ct/4 ) ) - 1






def UPDATE_PHOTOPHYSICS(sample: Sample, params: Params, U: Universal, v: np.int64):
    """
    Update #2
    ===
    Updates transition probabilities of the photo-states for each dye (`wi_D`, `wi_A`).

    Args:
        v (int): Index of the current simulation.
    """
    # UPDATE DONOR
    idx_D_1 = 1 + np.where(               sample.ft_D[:-1] )[0]
    idx_D_0 = 1 + np.where(np.logical_not(sample.ft_D[:-1]))[0]
    
    Da =  params.wi_D_prior_eta  + [ np.sum( sample.ft_D[idx_D_0] ), np.sum( sample.ft_D[idx_D_1] ), sample.ft_D[0]  ]
    Db = (params.wi_D_prior_zeta + [ np.sum( 1*np.logical_not( sample.ft_D[idx_D_0] ) ),
                                     np.sum( 1*np.logical_not( sample.ft_D[idx_D_1] ) ), 
                                             1*np.logical_not( sample.ft_D[0] ) ])
    sample.wi_D = U.RNGs[v].betarnd(Da, Db)
    
    # UPDATE ACCEPTOR
    idx_A_1 = 1 + np.where(               sample.ft_A[:-1] )[0]
    idx_A_0 = 1 + np.where(np.logical_not(sample.ft_A[:-1]))[0]
    
    Aa =  params.wi_A_prior_eta  + [ np.sum( sample.ft_A[idx_A_0] ), np.sum( sample.ft_A[idx_A_1] ), sample.ft_A[0]  ]
    Ab = (params.wi_A_prior_zeta + [ np.sum( 1*np.logical_not( sample.ft_A[idx_A_0] ) ),
                                     np.sum( 1*np.logical_not( sample.ft_A[idx_A_1] ) ), 
                                             1*np.logical_not( sample.ft_A[0] ) ] )
    sample.wi_A = U.RNGs[v].betarnd(Aa, Ab)






def UPDATE_TRANSITION_PROBS(sample: Sample, params: Params, U: Universal, v: np.int64):
    """
    Update #3
    ===
    Updates the nonhierarchical prior (`bm`) for transition probabilities of conformational states (`pm`, `ps`).

    Args:
        v (int): Index of the current simulation.
    """
    # FIRST GATHER COUNTS
    pm = np.zeros( (params.K_lim, params.K_lim), dtype = np.int64 )
    for t in range(1, params.T): 
        pm[ sample.st[t-1], sample.st[t] ] += 1
    pm = np.moveaxis(pm, 0, 1) # Adjusted.
    
    ps = np.zeros( params.K_lim, dtype = np.int64 )
    jk = sample.st[0]
    ps[jk] = 1

    x = np.where( pm>0 )
    x = zip(x[0], x[1])
    
    # THEN UPDATE BASE
    for rep in range(U.rep_bm): # Default: range(15)
        
        # RESAMPLE COUNTS M
        M = np.zeros( (params.K_lim+1, params.K_lim), dtype = np.float64 )
        Arm: np.ndarray[np.float64] = repmat( params.alpha*sample.bm, params.K_lim+1, 1 )
        for a, b in x:
            M[a,b] = 1 + ( np.sum( (U.RNGs[v].rand(pm[a,b]-1)*(   np.arange(1, pm[a,b]) + Arm[a,b]     )) < Arm[a,b],     axis = 0 ) )
        M.flat[jk] = 1 + ( np.sum( (U.RNGs[v].rand(ps.flat[jk])*( np.arange(1, ps[jk])  + Arm.flat[jk] )) < Arm.flat[jk], axis = 0 ) )
        
        # RESAMPLE BASE
        M = np.moveaxis( np.reshape( M, (params.K_lim, params.K_lim+1) ), 0, 1 ) # Adjusted.
        sample.bm = U.RNGs[v].dirrnd( params.gamma/params.K_lim + np.sum( M, axis = 0 ) )
        
    # LAST SAMPLE TRANSITION PROBS
    sample.pm = U.RNGs[v].dirrnd( pm + params.alpha * np.tile(sample.bm, (params.K_lim, 1)) )
    sample.ps = U.RNGs[v].dirrnd( ps + params.alpha * sample.bm )






def UPDATE_EMISSIONS(sample: Sample, params: Params, U: Universal, v: np.int64):
    """
    Update #4
    ===
    Updates the photoemission multiplier (`tht`) and photoemission rate priors (`tht`, `rho_D`, `rho_A`, `kap_D`, `kap_A`, `kap_Z`) based on `st`.
    
    Uses Hamiltonian Monte Carlo (`HMC()`), Metropolis-within-Gibbs (`MG()`), and Flipper (`FLIPPER()`) functions.

    Args:
        v (int): Index of the current simulation.
    """
    tht_phi_post = params.tht_prior_phi + np.sum( params.It_D + params.It_A )
    sample.K_set, sample.K_loc, sample.K_cnt = np.unique( sample.st, return_inverse=True, return_counts=True )
    sample.K_sz  = sample.K_set.size
    
    for rep in range(U.rep_tht): # Default: range(5)
        
        # UPDATE tht (direct)
        sample.tht = ( params.tht_prior_phi / params.tht_prior_psi 
                       + params.dD * ( params.T*( params.qD*sample.rho_D + params.qA*sample.rho_A )
                                       + ( params.qD*params.cDD + params.qA*params.cDA ) * np.sum( sample.ft_D * (sample.ft_A * sample.kap_D[sample.st] + (1-sample.ft_A) * sample.kap_Z) )
                                       + ( params.qD*params.cAD + params.qA*params.cAA ) * np.sum( sample.ft_D *  sample.ft_A * sample.kap_A[sample.st] ) ) )
        sample.tht = U.RNGs[v].gamrnd1( tht_phi_post )/sample.tht
        
        # UPDATE BACKGROUND + ACTIVE MULTIPLIERS (HMC) --> Default: range(125)
        if U.HMC_L != 0: HMC( sample, params, U, v )
        
        # UPDATE BACKGROUND + ACTIVE MULTIPLIERS (MG) --> Default: range(5)
        if U.MG_L != 0: MG( sample, params, U, v )
        
        # FLIP DONOR ACTIVE MULTIPLIERS (MG) --> once
        FLIPPER( sample, params, U, v )
    
    # UPDATE INACTIVE MULTIPLIERS (direct)
    Sm = np.setdiff1d( np.arange(params.K_lim), sample.K_set )
    sample.kap_D[Sm] = U.RNGs[v].gamrnd( params.kap_D_prior_phi, params.kap_D_prior_psi/params.kap_D_prior_phi, Sm.size )
    sample.kap_A[Sm] = U.RNGs[v].gamrnd( params.kap_A_prior_phi, params.kap_A_prior_psi/params.kap_A_prior_phi, Sm.size )






def HMC(sample: Sample, params: Params, U: Universal, v: np.int64):
    """
    Update Emissions > Hamiltonian Monte Carlo
    ===

    Updates the photoemission multiplier (`tht`) and photoemission rate priors (`tht`, `rho_D`, `rho_A`, `kap_D`, `kap_A`, `kap_Z`).

    Args:
        v (int): Index of the current simulation.
    """
    HMC_eps = U.HMC_eps * U.RNGs[v].rand1()

    # MASSES & MOMENTA
    m_RD = 1;                     sample_RD = m_RD * U.RNGs[v].randn1()
    m_RA = 1;                     sample_RA = m_RA * U.RNGs[v].randn1()
    m_KD = np.ones(sample.K_sz);  sample_KD = m_KD * U.RNGs[v].randn(sample.K_sz)
    m_KA = np.ones(sample.K_sz);  sample_KA = m_KA * U.RNGs[v].randn(sample.K_sz)
    m_KZ = 1;                     sample_KZ = m_KZ * U.RNGs[v].randn1()
    
    # PRE-ALLOCATION
    # 1: upsilon / Potential energy / main parameters. 2: P / Kinetic energy / rate of change.
    propos_rD = np.empty((U.HMC_L, 1));            propos_RD = np.empty((U.HMC_L, 1))
    propos_rA = np.empty((U.HMC_L, 1));            propos_RA = np.empty((U.HMC_L, 1))
    propos_kD = np.empty((U.HMC_L, sample.K_sz));  propos_KD = np.empty((U.HMC_L, sample.K_sz))
    propos_kA = np.empty((U.HMC_L, sample.K_sz));  propos_KA = np.empty((U.HMC_L, sample.K_sz))
    propos_kZ = np.empty((U.HMC_L, 1));            propos_KZ = np.empty((U.HMC_L, 1))

    # INIT STEP
    l = 0
    propos_rD[l] = sample.rho_D
    propos_rA[l] = sample.rho_A
    propos_kD[l] = sample.kap_D[sample.K_set]
    propos_kA[l] = sample.kap_A[sample.K_set]
    propos_kZ[l] = sample.kap_Z

    [U_rD, U_rA, U_kD, U_kA, U_kZ] = get_u_grad( sample, params, propos_rD[l], propos_rA[l], propos_kD[l], propos_kA[l], propos_kZ[l] )
    
    propos_RD[l] = sample_RD - 0.5 * HMC_eps * U_rD
    propos_RA[l] = sample_RA - 0.5 * HMC_eps * U_rA
    propos_KD[l] = sample_KD - 0.5 * HMC_eps * U_kD 
    propos_KA[l] = sample_KA - 0.5 * HMC_eps * U_kA 
    propos_KZ[l] = sample_KZ - 0.5 * HMC_eps * U_kZ

    # LEAP-FROG FORWARD
    for l in range(1, U.HMC_L - 1): # HMC_L (default) = 125
        propos_rD[l] = propos_rD[l-1] + HMC_eps * propos_RD[l-1] / m_RD
        propos_rA[l] = propos_rA[l-1] + HMC_eps * propos_RA[l-1] / m_RA
        propos_kD[l] = propos_kD[l-1] + HMC_eps * propos_KD[l-1] / m_KD 
        propos_kA[l] = propos_kA[l-1] + HMC_eps * propos_KA[l-1] / m_KA 
        propos_kZ[l] = propos_kZ[l-1] + HMC_eps * propos_KZ[l-1] / m_KZ
        
        [U_rD, U_rA, U_kD, U_kA, U_kZ] = get_u_grad( sample, params, propos_rD[l], propos_rA[l], propos_kD[l], propos_kA[l], propos_kZ[l] )
        
        propos_RD[l] = propos_RD[l-1] - HMC_eps * U_rD
        propos_RA[l] = propos_RA[l-1] - HMC_eps * U_rA
        propos_KD[l] = propos_KD[l-1] - HMC_eps * U_kD 
        propos_KA[l] = propos_KA[l-1] - HMC_eps * U_kA 
        propos_KZ[l] = propos_KZ[l-1] - HMC_eps * U_kZ

    # TERM STEP
    l = U.HMC_L - 1
    propos_rD[l] = propos_rD[l-1] + HMC_eps * propos_RD[l-1] / m_RD
    propos_rA[l] = propos_rA[l-1] + HMC_eps * propos_RA[l-1] / m_RA
    propos_kD[l] = propos_kD[l-1] + HMC_eps * propos_KD[l-1] / m_KD 
    propos_kA[l] = propos_kA[l-1] + HMC_eps * propos_KA[l-1] / m_KA 
    propos_kZ[l] = propos_kZ[l-1] + HMC_eps * propos_KZ[l-1] / m_KZ

    [U_rD, U_rA, U_kD, U_kA, U_kZ] = get_u_grad( sample, params, propos_rD[l], propos_rA[l], propos_kD[l], propos_kA[l], propos_kZ[l] )
    
    propos_RD[l] = propos_RD[l-1] - 0.5 * HMC_eps * U_rD
    propos_RA[l] = propos_RA[l-1] - 0.5 * HMC_eps * U_rA
    propos_KD[l] = propos_KD[l-1] - 0.5 * HMC_eps * U_kD
    propos_KA[l] = propos_KA[l-1] - 0.5 * HMC_eps * U_kA
    propos_KZ[l] = propos_KZ[l-1] - 0.5 * HMC_eps * U_kZ
    
    # ACCEPTANCE
    log_a = ( get_u( sample, params, sample.rho_D, sample.rho_A, sample.kap_D[sample.K_set], sample.kap_A[sample.K_set], sample.kap_Z )
            - get_u( sample, params, propos_rD[-1], propos_rA[-1], propos_kD[-1], propos_kA[-1], propos_kZ[-1]  )
            + 0.5 * ( (sample_RD**2         - propos_RD[-1]**2)         / m_RD
                    + (sample_RA**2         - propos_RA[-1]**2)         / m_RA
                    + (np.square(sample_KD) - np.square(propos_KD[-1])) / m_KD
                    + (np.square(sample_KA) - np.square(propos_KA[-1])) / m_KA
                    + (sample_KZ**2         - propos_KZ[-1]**2)         / m_KZ ) )
    
    if np.any( np.log(U.RNGs[v].rand1()) < log_a ):
        sample.rho_D               = propos_rD[-1][0]
        sample.rho_A               = propos_rA[-1][0]
        sample.kap_D[sample.K_set] = propos_kD[-1]
        sample.kap_A[sample.K_set] = propos_kA[-1]
        sample.kap_Z               = propos_kZ[-1][0]
        sample.rec[0][0]           = sample.rec[0][0] + 1
    sample.rec[0][1] = sample.rec[0][1] + 1






def get_u_grad(sample: Sample, params: Params, rD: np.float64, rA: np.float64, kD: np.ndarray[np.float64], kA: np.ndarray[np.float64], kZ: np.float64):
    """
    Update Emissions > Hamiltonian Monte Carlo > Gradient of the Log Posterior Density
    ===

    ∇_emission log f(emission)

    Args:
        rD (float): The donor dye's background emission.
        rA (float): The acceptor dye's background emission.
        kD (array[float]): The donor dye's emissions with FRET.
        kA (array[float]): The acceptor dye's emissions with FRET.
        kZ (float): The donor dye's emission without FRET.
    
    Returns:
        :U_rD: (float) Adjusted donor dye's background emission.
        :U_rA: (float) Adjusted acceptor dye's background emission.
        :U_kD: (array[float]) Adjusted donor dye's emissions with FRET.
        :U_kA: (array[float]) Adjusted acceptor dye's emissions with FRET.
        :U_kZ: (float) Adjusted donor dye's emission without FRET.
    """
    tmp1 = sample.ft_D * ( sample.ft_A * kD[sample.K_loc] + (1-sample.ft_A) * kZ )
    tmp2 = sample.ft_D *   sample.ft_A * kA[sample.K_loc]
    ID_hat_D = params.It_D / ( rD + params.cDD*tmp1 + params.cAD*tmp2 )
    IA_hat_A = params.It_A / ( rA + params.cDA*tmp1 + params.cAA*tmp2 )
    
    tmp1 = sample.tht * params.dD * params.T
    U_rD = ( - np.sum(ID_hat_D) 
             + params.qD * tmp1
             + (1 - params.rho_D_prior_phi)/rD
             +      params.rho_D_prior_phi /params.rho_D_prior_psi )
    
    U_rA = ( - np.sum(IA_hat_A) 
             + params.qA * tmp1
             + (1 - params.rho_A_prior_phi)/rA
             +      params.rho_A_prior_phi /params.rho_A_prior_psi )
    
    U_kD = np.empty(sample.K_sz)
    U_kA = np.empty(sample.K_sz)
    
    for k in range(sample.K_sz):
        idx = np.where((sample.ft_D == 1) & (sample.ft_A == 1) & (sample.K_loc == k))[0]
        tmp1 = sample.tht * params.dD * idx.size
        tmp2 = np.sum(ID_hat_D[idx])
        tmp3 = np.sum(IA_hat_A[idx])
        
        U_kD[k] = ( - params.cDD * tmp2 
                    - params.cDA * tmp3
                    + ( params.qD*params.cDD + params.qA*params.cDA ) * tmp1
                    + ( 1 - params.kap_D_prior_phi )/kD[k]
                    +       params.kap_D_prior_phi / params.kap_D_prior_psi )
        
        U_kA[k] = ( - params.cAD * tmp2 
                    - params.cAA * tmp3
                    + ( params.qD*params.cAD + params.qA*params.cAA ) * tmp1
                    + ( 1 - params.kap_A_prior_phi )/kA[k]
                    +       params.kap_A_prior_phi / params.kap_A_prior_psi )


    idx = np.where((sample.ft_D == 1) & (sample.ft_A == 0))[0]

    U_kZ = ( - params.cDD * np.sum(ID_hat_D[idx]) 
             - params.cDA * np.sum(IA_hat_A[idx])
             + ( params.qD*params.cDD + params.qA*params.cDA ) * sample.tht * params.dD * idx.size
             + ( 1 - params.kap_Z_prior_phi )/kZ
             +       params.kap_Z_prior_phi / params.kap_Z_prior_psi )
    
    return U_rD, U_rA, U_kD, U_kA, U_kZ






def get_u(sample: Sample, params: Params, rD: np.float64, rA: np.float64, kD: np.ndarray[np.float64], kA: np.ndarray[np.float64], kZ: np.float64):
    """
    Update Emissions > Hamiltonian Monte Carlo > Log Posterior Density
    ===

    f(emission)

    Args:
        rD (float): The donor dye's background emission.
        rA (float): The acceptor dye's background emission.
        kD (array[float]): The donor dye's emissions with FRET.
        kA (array[float]): The acceptor dye's emissions with FRET.
        kZ (float): The donor dye's emission without FRET.
    
    Returns:
        :float:
    """
    if ((rD<0) or (rA<0) or np.any(kD<0) or np.any(kA<0) or np.any(kZ<0)):
        return np.inf
    else:
        tmp1 = sample.ft_D * ( sample.ft_A * kD[sample.K_loc] + (1-sample.ft_A) * kZ )
        tmp2 = sample.ft_D *   sample.ft_A * kA[sample.K_loc]
        hat_D = rD + params.cDD*tmp1 + params.cAD*tmp2 # mu_D / tht
        hat_A = rA + params.cDA*tmp1 + params.cAA*tmp2 # mu_A / tht
        return ( sample.tht*params.dD*( params.qD * np.sum(hat_D) + params.qA * np.sum(hat_A) )
                  - np.sum( params.It_D * np.log(hat_D) + params.It_A * np.log(hat_A) )
                     + (1 - params.rho_D_prior_phi)*       np.log(rD)                + params.rho_D_prior_phi/params.rho_D_prior_psi*       rD
                     + (1 - params.rho_A_prior_phi)*       np.log(rA)                + params.rho_A_prior_phi/params.rho_A_prior_psi*       rA
                     + (1 - params.kap_D_prior_phi)*np.sum(np.log(kD))               + params.kap_D_prior_phi/params.kap_D_prior_psi*np.sum(kD)
                     + (1 - params.kap_A_prior_phi)*np.sum(np.log(kA))               + params.kap_A_prior_phi/params.kap_A_prior_psi*np.sum(kA)
                     + (1 - params.kap_Z_prior_phi)*       np.log(kZ, where=kZ != 0) + params.kap_Z_prior_phi/params.kap_Z_prior_psi*       kZ  )






def MG(sample: Sample, params: Params, U: Universal, v: np.int64):
    """
    Update Emissions > Metropolis-within-Gibbs
    ===

    Updates `rho_D`, `rho_A`, `kap_D`, `kap_A`, `kap_Z`, `rec`.

    Args:
        v (int): Index of the current simulation.
    """
    tmp1 = sample.tht * params.dD
    A = np.array([[params.cDD, params.cAD], [params.cDA, params.cAA]])
    
    tmp2 = sample.ft_D * ( sample.ft_A * sample.kap_D[sample.K_loc] + (1-sample.ft_A) * sample.kap_Z )
    tmp3 = sample.ft_D *   sample.ft_A * sample.kap_A[sample.K_loc]
    sample_mu_D = params.qD * tmp1 * ( sample.rho_D + params.cDD*tmp2 + params.cAD*tmp3 ) # Current photodetection rate of the donor dye.
    sample_mu_A = params.qA * tmp1 * ( sample.rho_A + params.cDA*tmp2 + params.cAA*tmp3 ) # Current photodetection rate of the acceptor dye.
    
    for l in range(U.MG_L): # Default: range(5)
        
        ri = U.RNGs[v].randi1(1, 3, True)
        propos_rho_D = sample.rho_D * ( 1 if ri == 2 else U.RNGs[v].gamrnd1( params.MG_a_D )/params.MG_a_D ) # Proposed scaling factor for background emission rate of the donor dye.
        propos_rho_A = sample.rho_A * ( 1 if ri == 3 else U.RNGs[v].gamrnd1( params.MG_a_A )/params.MG_a_A ) # Proposed scaling factor for background emission rate of the acceptor dye.
        
        b = np.array([sample.rho_D - propos_rho_D, sample.rho_A - propos_rho_A])
        x = np.linalg.solve(A, b)
        propos_kap_D: np.ndarray[np.float64] = sample.kap_D + x[0] # Proposed scaling factor for donor dye emission rate with FRET.
        propos_kap_A: np.ndarray[np.float64] = sample.kap_A + x[1] # Proposed scaling factor for acceptor dye emission rate with FRET.
        propos_kap_Z: np.float64             = sample.kap_Z + x[0] # Proposed scaling factor for donor dye emission rate without FRET.
        
        tmp2 = sample.ft_D * ( sample.ft_A * propos_kap_D[sample.K_loc] + (1-sample.ft_A) * propos_kap_Z )
        tmp3 = sample.ft_D *   sample.ft_A * propos_kap_A[sample.K_loc]
        propos_mu_D = params.qD * tmp1 * ( propos_rho_D + params.cDD*tmp2 + params.cAD*tmp3 ) # Proposed photodetection rate of the donor dye.
        propos_mu_A = params.qA * tmp1 * ( propos_rho_A + params.cDA*tmp2 + params.cAA*tmp3 ) # Proposed photodetection rate of the acceptor dye.
        
        # ACCEPTANCE
        if (np.all(propos_kap_D > 0) and np.all(propos_kap_A > 0) and (propos_kap_Z > 0)):
            
            tmp4A = propos_rho_D/sample.rho_D
            tmp4B = sample.rho_D/propos_rho_D
            tmp5A = propos_rho_A/sample.rho_A
            tmp5B = sample.rho_A/propos_rho_A
            log_a = ( np.sum((sample_mu_D-propos_mu_D) + (sample_mu_A-propos_mu_A)
                            + params.It_D * np.log( propos_mu_D / sample_mu_D )
                            + params.It_A * np.log( propos_mu_A / sample_mu_A ))
                     + ( params.rho_D_prior_phi-1 )*       np.log( tmp4A                      ) + params.rho_D_prior_phi*      ( sample.rho_D-propos_rho_D )/params.rho_D_prior_psi
                     + ( params.rho_A_prior_phi-1 )*       np.log( tmp5A                      ) + params.rho_A_prior_phi*      ( sample.rho_A-propos_rho_A )/params.rho_A_prior_psi
                     + ( params.kap_D_prior_phi-1 )*np.sum(np.log( propos_kap_D/sample.kap_D )) + params.kap_D_prior_phi*np.sum( sample.kap_D-propos_kap_D )/params.kap_D_prior_psi
                     + ( params.kap_A_prior_phi-1 )*np.sum(np.log( propos_kap_A/sample.kap_A )) + params.kap_A_prior_phi*np.sum( sample.kap_A-propos_kap_A )/params.kap_A_prior_psi
                     + ( params.kap_Z_prior_phi-1 )*       np.log( propos_kap_Z/sample.kap_Z )  + params.kap_Z_prior_phi*      ( sample.kap_Z-propos_kap_Z )/params.kap_Z_prior_psi
                     + ( 2*params.MG_a_D-1 )*np.log( tmp4B ) + params.MG_a_D*( tmp4A - tmp4B )
                     + ( 2*params.MG_a_A-1 )*np.log( tmp5B ) + params.MG_a_A*( tmp5A - tmp5B ) )
            will_accept = np.any( np.log(U.RNGs[v].rand1()) < log_a )
            
        else: will_accept = False
        
        if will_accept:
            sample.rho_D = propos_rho_D
            sample.rho_A = propos_rho_A
            sample.kap_D = propos_kap_D
            sample.kap_A = propos_kap_A
            sample.kap_Z = propos_kap_Z
            sample.rec[1][0] = sample.rec[1][0] + 1
            sample_mu_D  = propos_mu_D
            sample_mu_A  = propos_mu_A
        sample.rec[1][1] = sample.rec[1][1] + 1






def FLIPPER(sample: Sample, params: Params, U: Universal, v: np.int64):
    """
    Update Emissions > Metropolis-within-Gibbs > Flipper
    ===

    Swaps the donor emission without FRET `kap_Z` with a random donor emission with FRET `kap_D[k]` and tests.

    Updates `kap_D`, `kap_Z`, `rec`.

    Args:
        v (int): Index of the current simulation.
    """
    tmp1 = sample.tht * params.dD
    tmp3 = sample.ft_D * sample.ft_A * sample.kap_A[sample.K_loc]
    
    tmp2 = sample.ft_D * ( sample.ft_A * sample.kap_D[sample.K_loc] + (1-sample.ft_A) * sample.kap_Z )
    sample_mu_D = params.qD * tmp1 * ( sample.rho_D + params.cDD*tmp2 + params.cAD*tmp3 ) # Current photodetection rate of the donor dye.
    sample_mu_A = params.qA * tmp1 * ( sample.rho_A + params.cDA*tmp2 + params.cAA*tmp3 ) # Current photodetection rate of the acceptor dye.
        
    # FLIP HERE
    k = U.RNGs[v].randi1(0, params.K_lim)
    propos_kap_D: np.ndarray[np.float64] = deepcopy( sample.kap_D )     # Proposed scaling factor for donor dye emission rate with FRET.
    propos_kap_D[k]                      =           sample.kap_Z
    propos_kap_Z: np.float64             = deepcopy( sample.kap_D[k] )  # Proposed scaling factor for donor dye emission rate without FRET.
    
    tmp2 = sample.ft_D * ( sample.ft_A * propos_kap_D[sample.K_loc] + (1-sample.ft_A) * propos_kap_Z )
    propos_mu_D = params.qD * tmp1 * ( sample.rho_D + params.cDD*tmp2 + params.cAD*tmp3 ) # Proposed photodetection rate of the donor dye.
    propos_mu_A = params.qA * tmp1 * ( sample.rho_A + params.cDA*tmp2 + params.cAA*tmp3 ) # Proposed photodetection rate of the acceptor dye.
    
    # ACCEPTANCE
    log_a = ( np.sum( (sample_mu_D-propos_mu_D) + (sample_mu_A-propos_mu_A)
                     + params.It_D * np.log( propos_mu_D / sample_mu_D )
                     + params.It_A * np.log( propos_mu_A / sample_mu_A ) )
             + ( params.kap_D_prior_phi-1 )*np.log( propos_kap_D[k]/sample.kap_D[k] ) + params.kap_D_prior_phi*( sample.kap_D[k]-propos_kap_D[k] )/params.kap_D_prior_psi
             + ( params.kap_Z_prior_phi-1 )*np.log( propos_kap_Z   /sample.kap_Z    ) + params.kap_Z_prior_phi*( sample.kap_Z   -propos_kap_Z    )/params.kap_Z_prior_psi )
    
    if np.any( np.log(U.RNGs[v].rand1()) < log_a ):
        sample.kap_D = propos_kap_D
        sample.kap_Z = propos_kap_Z
        sample.rec[2][0] = sample.rec[2][0] + 1
    sample.rec[2][1] = sample.rec[2][1] + 1


