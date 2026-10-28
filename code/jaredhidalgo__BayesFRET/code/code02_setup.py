"""
BayesFRET: File 3
===
An Experiment-Adjusted HDP-HMM to Analyze Surface-Immobilized smFRET Data

About
---
code02_setup.py: The file to generate & save Figure 1 + synthetic data.

Author
---
Code by Jared Hidalgo. 

Inspired by MATLAB code from Ioannis Sgouralis, Shreya Madaan, Franky Djutanta, Rachael Kha, Rizal F. Hariadi, and Steve Pressé for "A Bayesian Nonparametric Approach to Single Molecule Förster Resonance Energy Transfer".
"""
# Import dependencies.
import os, sys
spl = __file__.split( os.sep )
DIR_PROGRAM = f"{os.sep}".join( spl[:-1] )
sys.path.append( DIR_PROGRAM )

# Import internal packages for use.
import warnings
# Import external packages for use.
import numpy as np
import matplotlib.ticker as ticker
from scipy.stats import norm, gamma
from scipy.special import beta as betaS
from scipy.special import gamma as gammaS
from scipy.special import kn as besselk
# Import packages for documentation.
from matplotlib.axes import Axes
from code01_classes import Params, True_Samples, Universal

mum = 6
red_a = 0.7





def SHOW_EXPERIMENTAL_DATA(units_t: str, units_I: str, It_D: np.ndarray[np.float64], It_A: np.ndarray[np.float64], T: np.int64, dt: np.float64, U: Universal):
    """
    Generate & save Figure 1. 
    
    Graphs
    ---
    * Photon intensities `It_A`, `It_D` (graph + density)
    * Apparent FRET efficiency `It_A/(It_D+It_A)` (graph + density)

    Args:
        units_t: Unit of time for figures.
        units_I: Unit of intensity for figures.
        It_D: Original data: Set of photon intensities of donor dye (size `T`).
        It_A: Original data: Set of photon intensities of acceptor dye (size `T`).
        T: Size of the dataset (`size(It_D) = size(It_A)`).
        dt: δt: Frame rate (a.k.a. measurement acquisition period, or the frequency of data) in seconds, regardless of the unit of time.

    Returns:
        :U: `Universal` object to handle Figure objects.
    """
    fig1 = U.make_figure( "1 | Preview: Experimental Data", True )
    tn_bnd = np.arange(T) * dt
    xlim = (-1, T*dt+1)
    bns = np.int64(T/10)

    # Graph photon intensities.
    ax1 = fig1.add_subplot( 3, mum, (1, 2*mum-1),
                            xlim        = xlim, 
                            ylabel      = "Intensities (" + units_I + ")" )
    ax1.grid( True, alpha=0.2 )
    ax1.xaxis.set_major_locator(ticker.MultipleLocator( np.int64(T*dt/10) ))
    ax1 = plot_stairs( ax1, tn_bnd, 
                            [It_D, It_A], 
                            [U.col_D, U.col_A], 
                            ['Donor channel', 'Acceptor channel'], U )
    ax1.legend( ncol=2, loc='upper right', shadow=False, framealpha=U.fa, fontsize=U.fs_l )

    # Graph densities of photon intensities.
    ax1D = fig1.add_subplot( 3, mum, (mum, 2*mum),
                             xticks      = [], 
                             xticklabels = [],
                             sharey      = ax1 )
    ax1D.grid( True, alpha=0.2 )
    ax1D.yaxis.tick_right()
    ax1D.yaxis.set_label_position("right")
    ax1D.hist( It_D, color=U.col_D, bins=bns, density=True, orientation="horizontal" )
    ax1D.hist( It_A, color=U.col_A, bins=bns, density=True, orientation="horizontal", alpha=red_a )
    
    # Graph apparent FRET efficiency.
    ax2 = fig1.add_subplot( 3, mum, (2*mum+1, 3*mum-1),
                            sharex = ax1, 
                            xlabel = f"Time in smFRET experiment ({units_t})",
                            ylim   = (-10, 110),
                            ylabel = "Apparent FRET\nefficiency (%)" )
    ax2.grid( True, alpha=0.2 )
    ax2.xaxis.set_major_locator(ticker.MultipleLocator( np.int64(T*dt/10) ))
    ax2 = plot_stairs( ax2, tn_bnd, 
                            100*It_A/(It_D+It_A),
                            U.col_m,
                            'Apparent', U )

    # Graph density of apparent FRET efficiency.
    ax2D = fig1.add_subplot( 3, mum, 3*mum,
                             xticks      = [], 
                             xticklabels = [], 
                             xlabel      = "Density", 
                             sharey      = ax2 )
    ax2D.grid( True, alpha=0.2 )
    ax2D.yaxis.tick_right()
    ax2D.yaxis.set_label_position("right")
    ax2D.hist( 100*It_A/(It_D+It_A), color=U.col_m, bins=bns, density=True, orientation="horizontal" )
    
    if U.show:
        fig1.show()
        fig1.canvas.flush_events()
    fig1.savefig(U.func_getActivePath( 'BayesFRET_fig01 preview_data.png' ))
    return U






def GENERATE_SYNTHETIC_DATA(U: Universal):
    """
    Generate & store data + Figure 1.
    
    Graphs
    ---
    * States of observed molecule `st`
    * States of donor + acceptor dyes `ft_D`, `ft_A`
    * Photon intensities `It_A`, `It_D` (graph + density)
    * Apparent FRET efficiency `It_A/(It_D+It_A)` (graph + density)

    Returns:
        :It_D: Set of photon intensities of donor dye.
        :It_A: Set of photon intensities of acceptor dye.
        :units_t: Unit of time for figures.
        :units_I: Unit of intensity for figures.
        :dt: δt: Frame rate (a.k.a. measurement acquisition period, or the frequency of data) in seconds, regardless of the unit of time.
        :dD: δτ = δt (1-d): Exposure period, or how much time within the frame rate (δt/`dt`) is dedicated to capturing the intensity using the dead time (d).
        :cDD: Cross-talk proportion from donor dye to donor channel. `cDD + cDA = 1`
        :cAA: Cross-talk proportion from acceptor dye to acceptor channel. `cAA + cAD = 1`
        :qD: Detector quantum efficiency (photodetection percentage) for the donor dye.
        :qA: Detector quantum efficiency (photodetection percentage) for the acceptor dye.
        :T: Size of the dataset (`size(It_D) = size(It_A)`).
        :U: `Universal` object to handle Figure objects.
    """
    units_I = "photons"
    units_t = "s"
    f     = 10/1                            # Frame rate 1/[t]
    d     = 0.01                            # Fraction of dead time per measurement acquisition period (dt)
    T     = 1000                            # Total number of steps
    c_D   = 0.9                             # Cross-talk coefficients
    c_A   = 0.75
    q_D   = 0.85                            # Quantum efficiency
    q_A   = 0.75
    tht   = 1e3                             # Overall emission rate [I]/[t]
    rho_D = 0.05                            # Background emission multipliers
    rho_A = 0.10
    kap_Z = 6.0                             # Dye emission multipliers
    kap_D = np.array( [3.5, 2.5, 1.5], dtype=np.float64 )
    kap_A = np.array( [1.0, 2.0, 3.0], dtype=np.float64 )
    pm    = np.array([[0.96, 0.04, 0.00],   # Transtion probs
                      [0.04, 0.92, 0.04], 
                      [0.00, 0.04, 0.96]], dtype=np.float64 )
    ps    = np.array( [0.50, 0.50, 0.00], dtype=np.float64 )
    w0_D  = 0.25                            # 0 -> 1
    w1_D  = 0.98                            # 1 -> 1
    ws_D  = 0.80                            # * -> 1
    w0_A  = 0.15                            # 0 -> 1
    w1_A  = 0.96                            # 1 -> 1
    ws_A  = 0.90                            # * -> 1
    
    D     = T/f                             # Total duration [t]
    ksi_D = rho_D * tht                     # Background emission rates [I]/[t]
    ksi_A = rho_A * tht
    lam_D = kap_D * tht                     # Dye emission rates [I]/[t]
    lam_A = kap_A * tht
    lam_Z = kap_Z * tht
    
    # GENERATE STATE TRACES
    st      = np.zeros(T, dtype = np.int64)
    ft_D    = np.zeros(T)
    ft_A    = np.zeros(T)
    st[0]   = U.get_samples(ps, -1)
    ft_D[0] = U.RNGs[-1].rand1() <= ws_D
    ft_A[0] = U.RNGs[-1].rand1() <= ws_A
    for t in range(1, T):
        st[t]   = U.get_samples(pm[ st[t-1] ], -1)
        ft_D[t] = U.RNGs[-1].rand1() <= (w1_D if ft_D[t-1] else w0_D)
        ft_A[t] = U.RNGs[-1].rand1() <= (w1_A if ft_A[t-1] else w0_A)
    
    # GENERATE INTENSITIES
    dt   = D/T                              # [t]
    dD   = (1-d) * dt                       # [t]
    It_D = ksi_D + c_D * ft_D * ( ft_A*lam_D[st] + (1-ft_A)*lam_Z ) + ( 1-c_A ) * ft_D * ft_A * lam_A[st]
    It_D = U.RNGs[-1].poissrnd( q_D*dD*It_D )
    It_D = np.maximum( It_D, np.full_like(It_D, np.finfo(float).tiny) )
    It_A = ksi_A + ( 1-c_D ) * ft_D * ( ft_A*lam_D[st] + (1-ft_A)*lam_Z ) + c_A * ft_D * ft_A * lam_A[st]
    It_A = U.RNGs[-1].poissrnd( q_A*dD*It_A )
    It_A = np.maximum( It_A, np.full_like(It_A, np.finfo(float).tiny) )
    
    # Store true values for comparison.
    U.TS.true_line = 100*lam_A[st]/(lam_D[st]+lam_A[st])
    U.TS.true_st = st
    U.TS.true_K_set, U.TS.true_K_loc, U.TS.true_K_cnt = np.unique( st, return_inverse=True, return_counts=True )
    U.TS.true_K_sz  = U.TS.true_K_set.size
    U.TS.true_tht = tht
    U.TS.true_rho_D = rho_D
    U.TS.true_rho_A = rho_A
    U.TS.true_kap_D = kap_D
    U.TS.true_kap_A = kap_A
    U.TS.true_kap_Z = kap_Z
    U.TS.true_ksi_D = ksi_D
    U.TS.true_ksi_A = ksi_A
    U.TS.true_lam_D = lam_D
    U.TS.true_lam_A = lam_A
    U.TS.true_lam_Z = lam_Z
    U.TS.true_ft_D = ft_D
    U.TS.true_ft_A = ft_A
    U.TS.true_wi_D = [w0_D, w1_D, ws_D]
    U.TS.true_wi_A = [w0_A, w1_A, ws_A]
    U.TS.true_pm = pm
    U.TS.true_ps = ps
    
    # GRAPH
    U = plot_synthetic_data(kap_D, T, dt, st, ft_D, ft_A, It_D, It_A, units_t, units_I, lam_D, lam_A, U)

    return It_D, It_A, units_t, units_I, dt, dD, c_D, c_A, q_D, q_A, T, U






def REUSE_SYNTHETIC_DATA(params: Params, true: True_Samples, U: Universal):
    """
    Setup existing synthetic data.
    """
    U.TS = true
    U = plot_synthetic_data( true.true_kap_D, params.T, params.dt, 
                             true.true_st, true.true_ft_D, true.true_ft_A, params.It_D, params.It_A,
                             params.units_t, params.units_I, true.true_lam_D, true.true_lam_A, U )
    return U






def plot_synthetic_data(kap_D, T, dt, st, ft_D, ft_A, It_D, It_A, units_t, units_I, lam_D, lam_A, U: Universal):
    """
    Create Figure 1.
    """
    # GRAPH
    M = kap_D.size
    tn_bnd = np.arange(T) * dt              # [t]
    bns = np.int64(T/10)
    xlim = (-1, T*dt+1)
    fig1 = U.make_figure( "1 | Preview: Synthetic Data", True )
    
    # Graph molecular states.
    ax1 = fig1.add_subplot( 5, mum, (1, mum-1),
                            xlim        = xlim, 
                            ylim        = (0, M+1),
                            ylabel      = "Conformational\nstate $k$",
                            yticks      = range(1, M+1),
                            yticklabels = [i for i in range(1, M+1)] )
    ax1.grid( True, alpha=0.2 )
    ax1.xaxis.set_major_locator(ticker.MultipleLocator( np.int64(T*dt/10) ))
    ax1 = plot_stairs( ax1, tn_bnd, st+1, U.col_M, [], U )
    
    # Graph states of fluorophores.
    ax2 = fig1.add_subplot( 5, mum, (mum+1, 2*mum-1),
                            sharex      = ax1,
                            ylim        = (-1, 2),
                            yticks      = (0, 1),
                            yticklabels = ['Dark', 'Bright'],
                            ylabel      = "Photo-State" )
    ax2.grid( True, alpha=0.2 )
    ax2 = plot_stairs( ax2, tn_bnd, 
                            [ft_D + 0.02, ft_A - 0.02],
                            [U.col_D, U.col_A], 
                            ['Donor', 'Acceptor'], U )

    # Put legend here.
    ax2D = fig1.add_subplot( 5, mum, 2*mum,
                             alpha = 0, frameon = False )
    ax2D.plot( [], [], color=U.col_D, label="Donor" )
    ax2D.plot( [], [], color=U.col_A, label="Acceptor" )
    ax2D.legend( loc="center left", shadow=False, framealpha=U.fa, fontsize=U.fs_l )
    ax2D.set_xticks([])
    ax2D.set_yticks([])
    
    # Graph photon intensities.
    ax3 = fig1.add_subplot( 5, mum, (2*mum+1, 4*mum-1),
                            sharex      = ax1, 
                            ylabel      = "Intensities (" + units_I + ")" )
    ax3.grid( True, alpha=0.2 )
    ax3 = plot_stairs( ax3, tn_bnd, 
                            [It_D, It_A], 
                            [U.col_D, U.col_A], 
                            [None, None], U ) # ['Donor channel', 'Acceptor channel']
    
    # Graph densities of photon intensities.
    ax3D = fig1.add_subplot( 5, mum, (3*mum, 4*mum),
                             xticks      = [], 
                             xticklabels = [], 
                             sharey      = ax3 )
    ax3D.grid( True, alpha=0.2 )
    ax3D.yaxis.tick_right()
    ax3D.yaxis.set_label_position("right")
    ax3D.hist( It_D, color=U.col_D, bins=bns, density=True, orientation="horizontal" )
    ax3D.hist( It_A, color=U.col_A, bins=bns, density=True, orientation="horizontal", alpha=red_a )
    
    # Graph apparent FRET efficiency.
    ax4 = fig1.add_subplot( 5, mum, (4*mum+1, 5*mum-1),
                            sharex = ax1, 
                            xlabel = f"Time in smFRET experiment ({units_t})",
                            ylim   = (-10, 140),
                            ylabel = "FRET\nefficiency (%)" )
    ax4.grid( True, alpha=0.2 )
    ax4.xaxis.set_major_locator(ticker.MultipleLocator( np.int64(T*dt/10) ))
    ax4 = plot_stairs( ax4, tn_bnd, 
                            [100*It_A/(It_D+It_A), 100*lam_A[st]/(lam_D[st]+lam_A[st])],
                            [U.col_m, U.col_M], 
                            ['Apparent', 'True'], U )
    ax4.legend( ncol=2, loc='upper right', shadow=False, framealpha=U.fa, fontsize=U.fs_l )

    # Graph density of apparent FRET efficiency.
    ax4D = fig1.add_subplot( 5, mum, 5*mum, 
                             xlabel      = "Density", 
                             xticks      = [], 
                             xticklabels = [], 
                             sharey      = ax4 )
    ax4D.grid( True, alpha=0.2 )
    ax4D.yaxis.tick_right()
    ax4D.yaxis.set_label_position("right")
    ax4D.hist( 100*It_A/(It_D+It_A), color=U.col_m, bins=bns, density=True, orientation="horizontal" )
    
    if U.show:
        fig1.show()
        fig1.canvas.flush_events()
    fig1.savefig(U.func_getActivePath( 'BayesFRET_fig01 preview_data.png' ))
    return U






def plot_stairs(ax: Axes, r_edges: np.ndarray[np.float64], w_centers: list[np.ndarray[np.float64]], colors: list[tuple], labels: list[str], U: Universal) -> Axes:
    """
    Template to make step graph.
    """
    if len(r_edges) == 0 and len(w_centers) == 0:
        r_edges   = np.array([[0], [1], [2.5], [3], [4]])
        w_centers = U.RNGs[-1].rand((r_edges.size - 1, 5))
    
    if len(colors) == 3:
        ax.step( r_edges, w_centers,    linewidth=0.5, color=colors,    label=labels )
    elif len(w_centers) == 2:
        ax.step( r_edges, w_centers[0], linewidth=0.5, color=colors[0], label=labels[0] )
        ax.step( r_edges, w_centers[1], linewidth=0.5, color=colors[1], label=labels[1] )
    
    return ax





def GRAPH_PRIORS(params: Params, U: Universal):
    """
    Generates & saves Figure 2a.
    """
    fig2 = U.make_figure( "2 | Priors", True )
    app        = np.array([params.It_D, params.It_A]) / params.dD
    lambda_ref = np.mean( app )
    w_lim      = np.array([0, 1.2*np.max( app )])
    num        = 3
    mum        = 4
    
    s11 = fig2.add_subplot( num, mum, (0*mum+1, (num-2)*mum+1),
                            ylim   = w_lim,
                            ylabel = "Photoemission (PE) rate (" + params.units_I + "/" + params.units_t + ")",
                            title  = "Apparent Dye PE" )
    s11.hist( params.It_D/params.dD, bins='auto', density=True, orientation='horizontal', facecolor=U.col_D, alpha=0.9, ec='black', label="$\\lambda_*^D$" )
    s11.hist( params.It_A/params.dD, bins='auto', density=True, orientation='horizontal', facecolor=U.col_A, alpha=0.5, ec='black', label="$\\lambda_*^A$" )
    s11.plot( s11.get_xlim(), [lambda_ref]*2, ':', color='k', label="$\\lambda_*$ mean" )
    s11.legend( loc='upper right', framealpha=0, fontsize=U.fs_l )
    s11.grid( True, alpha=0.3 )
    s11.ticklabel_format( style='sci', axis='x', scilimits=(0,0) )
    
    s12 = fig2.add_subplot( num, mum, (0*mum+2, (num-2)*mum+2),
                            xlim  = (0, 1),
                            ylim  = w_lim,
                            title = "Multiplier Prior" )
    s12 = show_prior( s12, 'gamma', '-', 'b', params.tht_prior_phi, params.tht_prior_psi/params.tht_prior_phi, U.eps, 1, 0, 0, "$\\theta$" )
    s12.plot( s12.get_xlim(), [lambda_ref]*2, ':', color='k' )
    s12.legend( loc='upper right', framealpha=0, fontsize=U.fs_l )
    s12.grid( True, alpha=0.3 )
    
    s13 = fig2.add_subplot( num, mum, (0*mum+3, (num-2)*mum+3),
                            xlim  = (0, 1),
                            ylim  = w_lim,
                            title = "Background PE Priors" )
    s13 = show_prior( s13, 'besselK', '--', U.col_D, [params.tht_prior_phi, params.rho_D_prior_phi], [params.tht_prior_psi, params.rho_D_prior_psi], U.eps, 1, 0, 0,   "$\\xi^D$" )
    s13 = show_prior( s13, 'besselK', '--', U.col_A, [params.tht_prior_phi, params.rho_A_prior_phi], [params.tht_prior_psi, params.rho_A_prior_psi], U.eps, 0.7, 0, 0, "$\\xi^A$" )
    s13.plot( s13.get_xlim(), [lambda_ref]*2, ':', color='k' )
    s13.legend( loc='upper right', framealpha=0, fontsize=U.fs_l )
    s13.grid( True, alpha=0.3 )
    
    lb = (0,1,1)
    s14 = fig2.add_subplot( num, mum, (0*mum+4, (num-2)*mum+4),
                            xlim   = (0, 1),
                            ylim   = w_lim,
                            title  = "Dye PE Priors" )
    s14 = show_prior( s14, 'besselK', '--', U.col_D, [params.tht_prior_phi, params.kap_D_prior_phi], [params.tht_prior_psi, params.kap_D_prior_psi], U.eps, 1, 0, 0,   "$\\lambda^D$" )
    s14 = show_prior( s14, 'besselK', '--', U.col_A, [params.tht_prior_phi, params.kap_A_prior_phi], [params.tht_prior_psi, params.kap_A_prior_psi], U.eps, 0.7, 0, 0, "$\\lambda^A$" )
    s14 = show_prior( s14, 'besselK', '--', lb,      [params.tht_prior_phi, params.kap_Z_prior_phi], [params.tht_prior_psi, params.kap_Z_prior_psi], U.eps, 0.4, 0, 0, "$\\lambda^Z$" )
    s14.plot( s14.get_xlim(), [lambda_ref]*2, ':', color='k' )
    s14.legend( loc='upper right', framealpha=0, fontsize=U.fs_l )
    s14.grid( True, alpha=0.3 )
    
    s21 = fig2.add_subplot( num, mum, (num-1)*mum+1,
                            xlim   = (0, 1),
                            ylim   = (-10, 110),
                            yticks = range(0, 120, 20),
                            ylabel = "Efficiency (%)",
                            xlabel = "Probability density",
                            title  = "FRET Efficiency Prior" )
    s21 = show_prior( s21, 'eff_custom', '--', 'b', [params.kap_D_prior_phi, params.kap_D_prior_psi], [params.kap_A_prior_phi, params.kap_A_prior_psi], U.eps, 1, 0.05, 0, "$E$" )
    xlim = s21.get_xlim()
    s21.plot( xlim, [0]*2,   ':', color='k', label="Limits" )
    s21.plot( xlim, [100]*2, ':', color='k' )
    s21.legend( loc='center left', framealpha=0, fontsize=U.fs_l )
    s21.grid( True, alpha=0.3 )
    
    s22 = fig2.add_subplot( num, mum, (num-1)*mum+2,
                              xlim   = (0, 1),
                              ylim   = (0, 5),
                              yticks = range(0, 6),
                              xlabel = "Probability density",
                              ylabel = "Value (1)",
                              title  = "Signal-to-Noise Ratio Priors" )
    BS_temp = np.linspace( s22.get_ylim()[0], s22.get_ylim()[1], num=1000, endpoint=True )
    fD_temp = params.rho_D_prior_psi / params.kap_D_prior_psi * params.kap_D_prior_phi / params.rho_D_prior_phi
    fA_temp = params.rho_A_prior_psi / params.kap_A_prior_psi * params.kap_A_prior_phi / params.rho_A_prior_phi
    pD_temp = fD_temp * np.power( fD_temp*BS_temp, params.kap_D_prior_phi-1 ) * np.power( 1+(fD_temp*BS_temp), -params.kap_D_prior_phi-params.rho_D_prior_phi ) / betaS( params.kap_D_prior_phi, params.rho_D_prior_phi )
    pA_temp = fA_temp * np.power( fA_temp*BS_temp, params.kap_A_prior_phi-1 ) * np.power( 1+(fA_temp*BS_temp), -params.kap_A_prior_phi-params.rho_A_prior_phi ) / betaS( params.kap_A_prior_phi, params.rho_A_prior_phi )
    xlim = s22.get_xlim()
    pD_temp = np.amin(xlim) + 1*(np.amax(xlim)-np.amin(xlim)) * pD_temp / np.amax(pD_temp)
    pA_temp = np.amin(xlim) + 1*(np.amax(xlim)-np.amin(xlim)) * pA_temp / np.amax(pA_temp)
    s22.plot( pD_temp, BS_temp, '--', color=U.col_D, label="$SNR^D$", alpha=1 )
    s22.plot( pA_temp, BS_temp, '--', color=U.col_A, label="$SNR^A$", alpha=0.7 )
    s22.plot( xlim, [1]*2, ':', color='k', label="" )
    s22.legend( loc='upper right', framealpha=0, fontsize=U.fs_l )
    s22.grid( True, alpha=0.3 )
    
    s23 = fig2.add_subplot( num, mum, (num-1)*mum+3,
                              xlim   = (0, 1),
                              ylim   = (0, 5),
                              yticks = range(0, 6),
                              xlabel = "Probability density",
                              title  = "Background Scaling Priors" )
    s23 = show_prior( s23, 'gamma', '-', U.col_D, params.rho_D_prior_phi, params.rho_D_prior_psi/params.rho_D_prior_phi, U.eps, 1, 0, 0,   "$\\rho^D$" )
    s23 = show_prior( s23, 'gamma', '-', U.col_A, params.rho_A_prior_phi, params.rho_A_prior_psi/params.rho_A_prior_phi, U.eps, 0.7, 0, 0, "$\\rho^A$" )
    s23.plot( s23.get_xlim(), [1]*2, ':', color='k' )
    s23.legend( loc='upper right', framealpha=0, fontsize=U.fs_l )
    s23.grid( True, alpha=0.3 )
    
    s24 = fig2.add_subplot( num, mum, (num-1)*mum+4,
                              xlim   = (0, 1),
                              ylim   = (0, 5),
                              yticks = range(0, 6),
                              xlabel = "Probability density",
                              title  = "Dye Scale Priors" )
    s24 = show_prior( s24, 'gamma', '-', U.col_D, params.kap_D_prior_phi, params.kap_D_prior_psi/params.kap_D_prior_phi, U.eps, 1, 0, 0,   "$\\kappa^D$" )
    s24 = show_prior( s24, 'gamma', '-', U.col_A, params.kap_A_prior_phi, params.kap_A_prior_psi/params.kap_A_prior_phi, U.eps, 0.7, 0, 0, "$\\kappa^A$" )
    s24 = show_prior( s24, 'gamma', '-', lb,      params.kap_Z_prior_phi, params.kap_Z_prior_psi/params.kap_Z_prior_phi, U.eps, 0.4, 0, 0, "$\\kappa^Z$" )
    s24.plot( s24.get_xlim(), [1]*2, ':', color='k' )
    s24.legend( loc='upper right', framealpha=0, fontsize=U.fs_l )
    s24.grid( True, alpha=0.3 )
    
    if U.show:
        fig2.show()
        fig2.canvas.flush_events()
    fig2.savefig(U.func_getActivePath( f"BayesFRET_fig02 priors.png" ))






def show_prior(ax: Axes, tag: str, line: str, color, p1, p2, eps, a: np.float64, exX: np.float64, exY: np.float64, label: str):
    
    y_range = np.linspace( ax.get_ylim()[0], ax.get_ylim()[1], num = 1000, endpoint = True )
    
    if tag == 'norm':            # p1 = mu, p2 = sg
        x_range = norm.pdf(y_range, p1, p2) + eps
        
    elif tag == 'symnorm':       # p1 = mu, p2 = sg
        x_range = 0.5 * norm.pdf(y_range, p1, p2) + 0.5 * norm.pdf(y_range, -p1, p2) + eps
        
    elif tag == 'gamma':         # p1 = a, p2 = b
        x_range = gamma.pdf(y_range, a=p1, scale=p2)
        
    elif tag == 'eff_custom':    # p1 = [phi_D,psi_D], p2 = [phi_A,psi_A]
        y_range = y_range / 100
        tmp = 1 / y_range-1
        x_range = ( np.power( p2[1] / p1[1] * p1[0] / p2[0], p1[0] ) 
                    / betaS( p1[0], p2[0] ) 
                    * np.power( tmp, p1[0]-1 ) 
                    * np.power( 1 + p2[1] / p1[1] * p2[0] / p1[0] * tmp, -p1[0]-p2[0] ) 
                    / np.square(y_range) )
        x_range[y_range < 0] = 0
        x_range[y_range > 1] = 0
        y_range = y_range * 100
        
    elif tag == 'besselK':       # p1[0:1] = phi[0:1], p2 = psi[0:1]
        warnings.filterwarnings('ignore')
        tmp = p1[0]*p1[1] / ( p2[0]*p2[1] ) * y_range
        x_range = ( eps + 2/( gammaS(p1[0]) * gammaS(p1[1]) ) 
                    * np.power( tmp, 0.5*(p1[0] + p1[1]) ) 
                    * besselk( p1[0] - p1[1], 2*np.sqrt(tmp) ) / y_range )
        x_range[y_range < 0] = 0
        warnings.resetwarnings()
        
    else: raise Exception( "INTERNAL ERROR: Unknown tag for code02_setup.py > show_prior()! Try again!" )
    
    x_range = ax.get_xlim()[0] + 1*( ax.get_xlim()[1] - ax.get_xlim()[0] ) * x_range / np.amax( x_range, where=~np.isnan(x_range), initial=-1 )
    ax.plot( x_range, y_range, line, color=color, alpha=a, label=label )
    if exX != 0: ax.set_xlim(( ax.get_xlim()[0] + exX, ax.get_xlim()[1] + exX ))
    if exY != 0: ax.set_ylim(( ax.get_ylim()[0] + exY, ax.get_ylim()[1] + exY ))
    return ax