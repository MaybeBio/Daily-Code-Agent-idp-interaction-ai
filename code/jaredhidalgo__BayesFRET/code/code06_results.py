"""
BayesFRET: File 7
===
An Experiment-Adjusted HDP-HMM to Analyze Surface-Immobilized smFRET Data

About
---
`code06_results.py`: The file to generate & save the rest of the figures after the simulations.

Function Hierarchy
---
* EFFICIENCY_AND_EMISSIONS()
* PHOTOPHYSICS()
* CONVERGENCE()
* FINAL_GRAPHS()
    * plot_best_and_quantiles()
    * make_colquad()

Author
---
Code by Jared Hidalgo. 

Inspired by MATLAB code from Ioannis Sgouralis, Shreya Madaan, Franky Djutanta, Rachael Kha, Rizal F. Hariadi, and Steve Pressé for "A Bayesian Nonparametric Approach to Single Molecule Förster Resonance Energy Transfer".
"""
# Import internal packages for use.
from collections import Counter
from pickle import load
# Import external packages for use.
import numpy as np
import matplotlib.ticker as ticker
from scipy.stats.mstats import mquantiles
# Import program files for documentation.
from code01_classes import Chain_History, Universal
from code03_mcmc import Params
from matplotlib.axes import Axes 






def GRAPH_ALL(params: Params, U: Universal):
    """
    Parent class to graph all results.

    Function Hierarchy
    ---
    * EFFICIENCY_AND_EMISSIONS()
    * PHOTOPHYSICS()
    * CONVERGENCE()
    * FINAL_GRAPHS()
        * plot_best_and_quantiles()
        * make_colquad()
    """
    S: list[Chain_History] = []
    for v in U.range_seeds:
        with open(U.filenames[v], "rb") as f:
            S.append( load(f) )
    idxA = np.arange( 0, params.N+1, 1, dtype=np.int64 ) # Range of all samples.
    burnin_idx = np.ceil( U.burn_in*params.N ) # End of the burn-in period.
    idxB = np.arange( burnin_idx, params.N+1, 1, dtype=np.int64 ) # Range of samples after the burn-in period.
    eff_star = 100 * params.It_A / (params.It_A + params.It_D)
    
    EFFICIENCY_AND_EMISSIONS(S, params, U, idxB, eff_star)
    PHOTOPHYSICS(S, params, U, idxB)
    CONVERGENCE(S, params, U, idxA)
    FINAL_GRAPHS(S, params, U, idxA, burnin_idx, idxB, eff_star)





def EFFICIENCY_AND_EMISSIONS(S: list[Chain_History], params: Params, U: Universal, idxB: np.ndarray, eff_star: np.ndarray[float]):
    """
    Generates & saves Figures 4, 5A, and 5B.

    Args:
        idxB (array): Range of samples after the burn-in period (default: [301, 1000]).
        eff_star (array[float]): The apparent FRET efficiency.
    """
    lam_D_star = params.It_D / params.dD
    lam_A_star = params.It_A / params.dD
    if U.is_syn:
        true_eff_mcmc = 100 * U.TS.true_kap_A / (U.TS.true_kap_A + U.TS.true_kap_D)
    

    #
    # FIGURE 4
    #
    U._updateStatus2( "Graphing Figure 4", "" )
    fig = U.make_figure( "4 | Densities: FRET Efficiencies", True )
    bins = 75
    a = 0.5
    at = 0.5
    lwt = 0.8

    ylim = (0, 100)
    yticks = range(0, 110, 10)
    for i in U.range_seeds:
        ax = fig.add_subplot( 2, 2, i+1, 
                               xlabel = f"Probability density of RNG Seed {U.RNGs[i].seed_str}",
                               ylim   = ylim,
                               yticks = yticks )
        eff_mcmc = 100 * S[i].eff[idxB]
        ax.hist( eff_mcmc, bins=bins, density=True, stacked=True, orientation='horizontal', color=[U.col_R[i]]*eff_mcmc.shape[1], label="Estimated" )
        ax.hist( eff_star, bins=bins, density=True, stacked=True, orientation='horizontal', color=U.col_m, alpha=a, label=("Apparent" if i == 0 else None) )
        if U.is_syn:
            ax.plot( ax.get_xlim(), [true_eff_mcmc[0]]*2,  'k--', lw=lwt, alpha=at, label=("True" if i == 0 else None) )
            ax.plot( ax.get_xlim(), [true_eff_mcmc[1:]]*2, 'k--', lw=lwt, alpha=at )
        if i % 2 == 0: ax.set_ylabel( "FRET efficiency (%)" )

    if U.show:
        fig.show()
        fig.canvas.flush_events()
    fig.savefig(U.func_getActivePath( f"BayesFRET_fig04 FRET_efficiency_densities.png" ))
    U._updateStatus2( "Done!", "\n" )
    

    #
    # FIGURE 5a
    #
    U._updateStatus2( "Graphing Figure 5a", "" )
    fig = U.make_figure( "5a | Densities: Donor Emission Rates", True )
    bins1 = 125
    bins2 = 45
    bins3 = 25
    a = 0.5

    # Bottom graphs
    for i in U.range_seeds: 
        ax = fig.add_subplot( 3, 4, i+9, xlabel = "Probability density" )
        # Set y-label only for the first subplot.
        if i == 0: ax.set_ylabel( "$\\xi^D$: Background PE rate ({}/{})".format(params.units_I, params.units_t) )
        # Main histogram.
        ax.hist( S[i].ksi_D[idxB], bins=bins3, density=True, orientation='horizontal', color=U.col_R[i] )
        # Plot true_ksi_D if synthetic data is available.
        if U.is_syn: ax.plot( ax.get_xlim(), [U.TS.true_ksi_D]*2, 'k--', lw=lwt, alpha=at )
        # Set scientific notation for x-axis.
        ax.ticklabel_format( style='sci', axis='x', scilimits=(0,0) )
    
    # Upper graphs
    lam_D_star = params.It_D / params.dD
    #lam_D_mcmcs = [np.moveaxis( S[i].lam_D[idxB], 0, 1 ) for i in U.range_seeds]
    max_lam_D = np.amax([ 3*np.mean(l, dtype = l.dtype) for l in [ np.moveaxis( S[i].lam_D[idxB], 0, 1 ) for i in U.range_seeds ] ])
    for i in U.range_seeds: 
        if_first = i == 0
        ax = fig.add_subplot( 3, 4, (i+1, i+5),
                                 ylim  = (0, max_lam_D),
                                 title = f"RNG Seed {U.RNGs[i].seed_str}" )
        # Set y-label only for the first subplot.
        if if_first: ax.set_ylabel( "$\\lambda^D$: Donor PE rate ({}/{})".format(params.units_I, params.units_t) )
        # Main histograms.
        lam_D_mcmc = np.moveaxis( S[i].lam_D[idxB], 0, 1 )
        ax.hist( lam_D_mcmc, bins=bins1, density=True, stacked=True, orientation='horizontal', color=[U.col_R[i]]*lam_D_mcmc.shape[1], label="Estimated" )
        ax.hist( lam_D_star, bins=bins2, density=True, stacked=True, orientation='horizontal', color=U.col_m, alpha=a, label=("Apparent" if if_first else None) )
        # Handle true values and their labels.
        if U.is_syn:
            if if_first:
                ax.plot( ax.get_xlim(), [U.TS.true_lam_D[0]]*2,  'k--', lw=lwt, alpha=at, label="True" )
                ax.plot( ax.get_xlim(), [U.TS.true_lam_D[1:]]*2, 'k--', lw=lwt, alpha=at )
            else:
                ax.plot( ax.get_xlim(), [U.TS.true_lam_D]*2,  'k--', lw=lwt, alpha=at )
        # Set scientific notation for x-axis.
        ax.ticklabel_format( style='sci', axis='x', scilimits=(0,0) )
        ax.legend( loc='upper right', framealpha=0, fontsize=U.fs_l )

    if U.show:
        fig.show()
        fig.canvas.flush_events()
    fig.savefig(U.func_getActivePath( f"BayesFRET_fig05a donor_emission_rates.png" ))
    U._updateStatus2( "Done!", "\n" )
    

    #
    # FIGURE 5b
    #
    U._updateStatus2( "Graphing Figure 5b", "" )
    fig = U.make_figure( "5b | Densities: Acceptor Emission Rates", True )

    # Bottom graphs
    for i in U.range_seeds: 
        ax = fig.add_subplot( 3, 4, i+9, xlabel = "Probability density" )
        # Set y-label only for the first subplot.
        if i == 0: ax.set_ylabel( "$\\xi^A$: Background PE rate ({}/{})".format(params.units_I, params.units_t) )
        # Main histogram.
        ax.hist( S[i].ksi_A[idxB], bins=bins3, density=True, orientation='horizontal', color=U.col_R[i] )
        # Plot true_ksi_A if synthetic data is available.
        if U.is_syn: ax.plot( ax.get_xlim(), [U.TS.true_ksi_A]*2, 'k--', lw=lwt, alpha=at )
        # Set scientific notation for x-axis.
        ax.ticklabel_format( style='sci', axis='x', scilimits=(0,0) )
    
    # Upper graphs
    lam_A_star = params.It_A / params.dD
    #lam_A_mcmcs = [ np.moveaxis( S[i].lam_A[idxB], 0, 1 ) for i in U.range_seeds ]
    max_lam_A = np.amax([ 3*np.mean(l, dtype = l.dtype) for l in [ np.moveaxis( S[i].lam_A[idxB], 0, 1 ) for i in U.range_seeds ] ])
    for i in U.range_seeds: 
        if_first = i == 0
        ax = fig.add_subplot( 3, 4, (i+1, i+5),
                                 ylim  = (0, max_lam_A),
                                 title = f"RNG Seed {U.RNGs[i].seed_str}" )
        # Set y-label only for the first subplot.
        if if_first: ax.set_ylabel( "$\\lambda^A$: Acceptor PE rate ({}/{})".format(params.units_I, params.units_t) )
        # Main histograms.
        lam_A_mcmc = np.moveaxis( S[i].lam_A[idxB], 0, 1 )
        ax.hist( lam_A_mcmc, bins=bins1, density=True, stacked=True, orientation='horizontal', color=[U.col_R[i]]*lam_A_mcmc.shape[1], label="Estimated" )
        ax.hist( lam_A_star,     bins=bins2, density=True, stacked=True, orientation='horizontal', color=U.col_m, alpha=a, label=("Apparent" if if_first else None) )
        # Handle true values and their labels.
        if U.is_syn:
            if if_first:
                ax.plot( ax.get_xlim(), [U.TS.true_lam_A[0]]*2,  'k--', lw=lwt, alpha=at, label="True" )
                ax.plot( ax.get_xlim(), [U.TS.true_lam_A[1:]]*2, 'k--', lw=lwt, alpha=at )
            else:
                ax.plot( ax.get_xlim(), [U.TS.true_lam_A]*2,  'k--', lw=lwt, alpha=at )
        # Set scientific notation for x-axis.
        ax.ticklabel_format( style='sci', axis='x', scilimits=(0,0) )
        ax.legend( loc='upper right', framealpha=0, fontsize=U.fs_l )

    if U.show:
        fig.show()
        fig.canvas.flush_events()
    fig.savefig(U.func_getActivePath( f"BayesFRET_fig05b acceptor_emission_rates.png" ))
    U._updateStatus2( "Done!", "\n" )






def PHOTOPHYSICS(S: list[Chain_History], params: Params, U: Universal, idxB: np.ndarray):
    """
    Generates & saves Figures 6a - 6d.

    Args:
        idxB (array): Range of samples after the burn-in period (default: [301, 1000]).
    """
    U._updateStatus2( "Gathering data for Figure 6a - 6d", "\r" )
    if U.is_syn:
        true_w0_D_mcmc   = 1 - U.TS.true_wi_D[0]
        true_w1_D_mcmc   =     U.TS.true_wi_D[1]
        true_w0_A_mcmc   = 1 - U.TS.true_wi_A[0]
        true_w1_A_mcmc   =     U.TS.true_wi_A[1]
        true_tau0_D_mean = params.dt/(1 - true_w0_D_mcmc)
        true_tau1_D_mean = params.dt/(1 - true_w1_D_mcmc)
        true_tau0_A_mean = params.dt/(1 - true_w0_A_mcmc)
        true_tau1_A_mean = params.dt/(1 - true_w1_A_mcmc)
    
    w_bnd = np.linspace(0, 1, 175)
    o     = (1, 0.5, 0)
    g     = (0, 1, 0)
    m     = (1, 0, 1)
    r     = (1, 0, 0)
    al    = 1
    at    = 1
    lwt   = 0.8
    fig_labels = ["6a", "6b", "6c", "6d"]

    for i in U.range_seeds:
        U._updateStatus2( f"Graphing Figure {fig_labels[i]}", "" )
        w0_D_mcmc = 1 - S[i].wi_D[idxB, 0]
        w1_D_mcmc =     S[i].wi_D[idxB, 1]
        w0_A_mcmc = 1 - S[i].wi_A[idxB, 0]
        w1_A_mcmc =     S[i].wi_A[idxB, 1]
        fig = U.make_figure( f"{fig_labels[i]} | Densities: Photophysics of RNG Seed {U.RNGs[i].seed_str}", True )
        ax1 = fig.add_subplot( 2, 3, 1,
                               xlim = (0, 1),
                               ylabel = "Donor probability density" )
        ax1.hist( w0_D_mcmc, w_bnd, density=True, facecolor=o, alpha=al, label="$\\omega^D_0$" )
        ax1.hist( w1_D_mcmc, w_bnd, density=True, facecolor=g, alpha=al, label="$\\omega^D_1$" )
        ylim = ax1.get_ylim()
        if U.is_syn:
            ax1.plot( [true_w0_D_mcmc]*2, ylim, 'k--', lw=lwt+0.1, alpha=at )
            ax1.plot( [true_w0_D_mcmc]*2, ylim, '--', color=o, lw=lwt, alpha=at, label="True $\\omega^D_0$" )
            ax1.plot( [true_w1_D_mcmc]*2, ylim, 'k--', lw=lwt+0.1, alpha=at )
            ax1.plot( [true_w1_D_mcmc]*2, ylim, '--', color=g, lw=lwt, alpha=at, label="True $\\omega^D_1$" )
        ax1.set_ylim(( ylim[0], 1.3*ylim[1] )) # Resize y-axis for legend.
        ax1.legend( loc='upper center', ncol=2, framealpha=0, fontsize=U.fs_l )
    
        ax4 = fig.add_subplot( 2, 3, 4,
                                sharex = ax1,
                                xlabel = "$\\omega$: Self-transition photoswitching prob.",
                                ylabel = "Acceptor probability density" )
        ax4.hist( w0_A_mcmc, w_bnd, density=True, facecolor=m, alpha=al, label="$\\omega^A_0$" )
        ax4.hist( w1_A_mcmc, w_bnd, density=True, facecolor=r, alpha=al, label="$\\omega^A_1$" )
        ylim = ax4.get_ylim()
        if U.is_syn:
            ax4.plot( [true_w0_A_mcmc]*2, ylim, 'k--', lw=lwt+0.1, alpha=at )
            ax4.plot( [true_w0_A_mcmc]*2, ylim, '--', color=m, lw=lwt, alpha=at, label="True $\\omega^A_0$" )
            ax4.plot( [true_w1_A_mcmc]*2, ylim, 'k--', lw=lwt+0.1, alpha=at )
            ax4.plot( [true_w1_A_mcmc]*2, ylim, '--', color=r, lw=lwt, alpha=at, label="True $\\omega^A_1$" )
        ax4.set_ylim(( ylim[0], 1.3*ylim[1] )) # Resize y-axis for legend.
        ax4.legend( loc='upper center', ncol=2, framealpha=0, fontsize=U.fs_l )
        
        t_bnd = np.logspace( np.log10(0.5*params.dt), np.log10(2*params.N*params.dt), 75 )
        w1 = np.zeros_like( params.dt/(1-w0_D_mcmc) ) + 1. / w0_D_mcmc.size
        w2 = np.zeros_like( params.dt/(1-w1_D_mcmc) ) + 1. / w1_D_mcmc.size
        w3 = np.zeros_like( params.dt/(1-w0_A_mcmc) ) + 1. / w0_A_mcmc.size
        w4 = np.zeros_like( params.dt/(1-w1_A_mcmc) ) + 1. / w1_A_mcmc.size
        
        ax2 = fig.add_subplot( 2, 3, (2, 3),
                                xlim = (params.dt * 0.5, params.dt * 2 * params.N),
                                xscale = 'log',
                                ylabel = "Donor probability (norm.)" )
        ax2.hist( params.dt/(1 - w0_D_mcmc), t_bnd, weights=w1, density=True, stacked=True, facecolor=o, label="$\\tau^D_0$" )
        ax2.hist( params.dt/(1 - w1_D_mcmc), t_bnd, weights=w2, density=True, stacked=True, facecolor=g, label="$\\tau^D_1$" )
        ylim = ax2.get_ylim()
        ax2.plot( [params.dt]*2,          ylim, 'k--', label="$\\delta t$, $N$ x $\\delta t$" )
        ax2.plot( [params.N*params.dt]*2, ylim, 'k--' )
        if U.is_syn:
            ax2.plot( [true_tau0_D_mean]*2, ylim, 'k--', lw=lwt+0.1, alpha=at )
            ax2.plot( [true_tau0_D_mean]*2, ylim, '--', color=o, lw=lwt, alpha=at, label="True $\\tau^D_0$" )
            ax2.plot( [true_tau1_D_mean]*2, ylim, 'k--', lw=lwt+0.1, alpha=at )
            ax2.plot( [true_tau1_D_mean]*2, ylim, '--', color=g, lw=lwt, alpha=at, label="True $\\tau^D_1$" )
        ax2.set_ylim(( ylim[0], 1.15*ylim[1] )) # Resize y-axis for legend.
        ax2.legend( loc='upper center', ncol=5, framealpha=0, fontsize=U.fs_l )
        
        ax5 = fig.add_subplot( 2, 3, (5, 6),
                                sharex = ax2,
                                xlabel = "$\\tau$: Mean dwell time ({}) FOR SEED {}".format( params.units_t, U.RNGs[i].seed_str ),
                                ylabel = "Acceptor probability (norm.)" )
        ax5.hist( params.dt/(1 - w0_A_mcmc), t_bnd, weights=w3, density=True, stacked=True, facecolor=m, label="$\\tau^A_0$" )
        ax5.hist( params.dt/(1 - w1_A_mcmc), t_bnd, weights=w4, density=True, stacked=True, facecolor=r, label="$\\tau^A_1$" )
        ylim = ax5.get_ylim()
        ax5.plot( [params.dt]*2,          ylim, 'k--', label="$\\delta t$, $N$ x $\\delta t$" )
        ax5.plot( [params.N*params.dt]*2, ylim, 'k--' )
        if U.is_syn:
            ax5.plot( [true_tau0_A_mean]*2, ylim, 'k--', lw=lwt+0.1, alpha=at )
            ax5.plot( [true_tau0_A_mean]*2, ylim, '--', color=m, lw=lwt, alpha=at, label="True $\\tau^A_0$" )
            ax5.plot( [true_tau1_A_mean]*2, ylim, 'k--', lw=lwt+0.1, alpha=at )
            ax5.plot( [true_tau1_A_mean]*2, ylim, '--', color=r, lw=lwt, alpha=at, label="True $\\tau^A_1$" )
        ax5.set_ylim(( ylim[0], 1.15*ylim[1] )) # Resize y-axis for legend.
        ax5.legend( loc='upper center', ncol=5, framealpha=0, fontsize=U.fs_l )

        if U.show:
            fig.show()
            fig.canvas.flush_events()
        fig.savefig(U.func_getActivePath( f"BayesFRET_fig0{fig_labels[i]} photophysics_of_seed_{U.RNGs[i].seed_str}.png" ))
        U._updateStatus2( "Done!", "\n" )






def CONVERGENCE(S: list[Chain_History], params: Params, U: Universal, idxA: np.ndarray):
    """
    Generates & saves Figure 7. Visualizes the convergence of key variables.

    Args:
        idxA (array): All samples (default: [0, 1000]).
    """
    U._updateStatus2( "Graphing Figure 7", "" )
    fig7 = U.make_figure( "7 | Results: Sample Convergence", True )
    iln = len(idxA)
    lw = 0.2
    psz = 0.3

    # Initialize plots.
    ax1 = fig7.add_subplot( 2, 4, 1, 
                            ylabel = "$\\omega^D_0$: Donor dark to bright prob. (%)" )
    ax1.xaxis.set_major_locator(ticker.MultipleLocator( np.int64(params.N/4) ))
    ax2 = fig7.add_subplot( 2, 4, 2, 
                            sharex = ax1, 
                            ylabel = "$\\omega^D_1$: Donor bright to bright prob. (%)" )
    ax3 = fig7.add_subplot( 2, 4, 3, 
                            sharex = ax1, 
                            ylabel = "$\\omega^A_0$: Acceptor dark to bright prob. (%)" )
    ax4 = fig7.add_subplot( 2, 4, 4, 
                            sharex = ax1, 
                            ylabel = "$\\omega^A_1$: Acceptor bright to bright prob. (%)" )
    ax5 = fig7.add_subplot( 2, 4, 5, 
                            sharex = ax1,
                            xlabel = "$n$: MCMC iteration", 
                            ylabel = "$\\xi^D$: Donor background rate ("+params.units_I+"/"+params.units_t+")" )
    ax6 = fig7.add_subplot( 2, 4, 6, 
                            sharex = ax1,
                            xlabel = "$n$: MCMC iteration", 
                            ylabel = "$\\xi^A$: Acceptor background rate ("+params.units_I+"/"+params.units_t+")" )
    ax7 = fig7.add_subplot( 2, 4, 7, 
                            sharex = ax1,
                            xlabel = "$n$: MCMC iteration", 
                            ylabel = "$\\lambda^Z$: Donor no FRET rate ("+params.units_I+"/"+params.units_t+")" )
    ax8 = fig7.add_subplot( 2, 4, 8,
                            sharex = ax1,
                            xlabel = "$n$: MCMC iteration",
                            ylabel = "$K$: Cumulative mode of states" )

    # Graph variable for each simulation.
    alphas = np.linspace(0.6, 0.3, 4)
    for i, sample in enumerate(S):
        ax1.plot( idxA, sample.wi_D[idxA,0].flat.copy()*100, '-o', markersize=psz, lw=lw, color=U.col_R[i], alpha=alphas[i], 
                  label=f"{U.RNGs[i].seed_str}" )
        ax2.plot( idxA, sample.wi_D[idxA,1].flat.copy()*100, '-o', markersize=psz, lw=lw, color=U.col_R[i], alpha=alphas[i] )
        ax3.plot( idxA, sample.wi_A[idxA,0].flat.copy()*100, '-o', markersize=psz, lw=lw, color=U.col_R[i], alpha=alphas[i] )
        ax4.plot( idxA, sample.wi_A[idxA,1].flat.copy()*100, '-o', markersize=psz, lw=lw, color=U.col_R[i], alpha=alphas[i] )
        ax5.plot( idxA, sample.ksi_D[idxA],                  '-o', markersize=psz, lw=lw, color=U.col_R[i], alpha=alphas[i] )
        ax6.plot( idxA, sample.ksi_A[idxA],                  '-o', markersize=psz, lw=lw, color=U.col_R[i], alpha=alphas[i] )
        ax7.plot( idxA, sample.lam_Z[idxA],                  '-o', markersize=psz, lw=lw, color=U.col_R[i], alpha=alphas[i] )
        ax8.plot( idxA, [Counter(sample.K_sz[:n]).most_common(1)[0][0] for n in range(1, iln+1)], 
                                                           '-o', markersize=psz, lw=lw, color=U.col_R[i], alpha=alphas[i] )
    
    # If synthetic, add the true values used in generating the synthetic data.
    if U.is_syn:
        lww = 0.8
        lwb = lww-0.1
        a = 0.8
        ax1.plot( idxA, np.repeat(U.TS.true_wi_D[0]*100, iln), 'w--', lw=lww, alpha=a )
        ax1.plot( idxA, np.repeat(U.TS.true_wi_D[0]*100, iln), 'k--', lw=lwb, alpha=a, label="True" )
        ax2.plot( idxA, np.repeat(U.TS.true_wi_D[1]*100, iln), 'w--', lw=lww, alpha=a )
        ax2.plot( idxA, np.repeat(U.TS.true_wi_D[1]*100, iln), 'k--', lw=lwb, alpha=a )
        ax3.plot( idxA, np.repeat(U.TS.true_wi_A[0]*100, iln), 'w--', lw=lww, alpha=a )
        ax3.plot( idxA, np.repeat(U.TS.true_wi_A[0]*100, iln), 'k--', lw=lwb, alpha=a )
        ax4.plot( idxA, np.repeat(U.TS.true_wi_A[1]*100, iln), 'w--', lw=lww, alpha=a )
        ax4.plot( idxA, np.repeat(U.TS.true_wi_A[1]*100, iln), 'k--', lw=lwb, alpha=a )
        ax5.plot( idxA, np.repeat(U.TS.true_ksi_D, iln), 'w--', lw=lww, alpha=a )
        ax5.plot( idxA, np.repeat(U.TS.true_ksi_D, iln), 'k--', lw=lwb, alpha=a )
        ax6.plot( idxA, np.repeat(U.TS.true_ksi_A, iln), 'w--', lw=lww, alpha=a )
        ax6.plot( idxA, np.repeat(U.TS.true_ksi_A, iln), 'k--', lw=lwb, alpha=a )
        ax7.plot( idxA, np.repeat(U.TS.true_lam_Z, iln), 'w--', lw=lww, alpha=a )
        ax7.plot( idxA, np.repeat(U.TS.true_lam_Z, iln), 'k--', lw=lwb, alpha=a )
        ax8.plot( idxA, np.repeat(U.TS.true_K_sz, iln), 'w--', lw=lww, alpha=a )
        ax8.plot( idxA, np.repeat(U.TS.true_K_sz, iln), 'k--', lw=lwb, alpha=a )
    
    ax1.legend( title="RNG Seeds", loc='upper right', framealpha=0, fontsize=U.fs_l )
    
    if U.show:
        fig7.show()
        fig7.canvas.flush_events()
    fig7.savefig(U.func_getActivePath( f"BayesFRET_fig07 sample_convergence.png" ))
    U._updateStatus2( "Done!", "\n" )
    






def FINAL_GRAPHS(S: list[Chain_History], params: Params, U: Universal, idxA: np.ndarray, burnin_idx: int, idxB: np.ndarray, eff_star: np.ndarray[float]):
    """
    Generates & saves Figures 8a + 8b.

    Args:
        idxA (array): All samples (default: [0, 1000]).
        burnin_idx (int): End of the burn-in period (default: 300).
        idxB (array): Range of samples after the burn-in period (default: [301, 1000]).
        eff_star (array[float]): The apparent FRET efficiency.
    """
    #
    # SETUP DATA
    #
    U._updateStatus2( "Gathering + writing data for Figures 8a + 8b", "\r" )
    ticker_mult = np.int64( params.T*params.dt/10 )
    t = np.arange( 0.5, params.T ) * params.dt
    eff_mcmc_all = [[ 100 * c.eff[ idxB[j], c.st[idxB[j]] ] 
                      for j in range(idxB.size) ]
                    for c in S]

    # Best FRET efficiency traces.
    eff_medn_all = [ np.median( mc, axis=0 ) for mc in eff_mcmc_all ] # Median of all FRET traces for each simulation.
    ix_all = [ np.argmin( np.sum( (med - mc) ** 2, 1 ) ) for med, mc in zip(eff_medn_all, eff_mcmc_all) ] # Indices of best traces for sets ignoring burn-in period (`mc`).
    eff_medn_best_all = [ mc[ix] for mc, ix in zip(eff_mcmc_all, ix_all) ] # Best traces.
    ix_adj = [int(i + burnin_idx) for i in ix_all] # Adjusted indices of best traces for sets with all samples.
    
    # Extract data from best traces.
    for v in U.range_seeds:
        with open( U.func_getActivePath( f"BayesFRET_data_sim_{U.RNGs[v].seed_str}_preview.txt" ), "a" ) as file:
            a = f"Best Sample (n = {ix_adj[v]})"
            b = "".join( ["-"]*len(a) )
            file.write( f"\n\n{a}\n{b}\n" )
            for obj in [S[v]]:
                for attr in vars(obj):
                    if attr not in ["seed_int", "runtime", "rec_ord"]:
                        file.write( f"{attr}: {getattr(obj, attr, 'NaN')[ix_adj[v]]}\n\n\n" )
    
    # Posterior quantiles.
    col_start = [(1, 190/255, 190/255), (1, 1, 190/255), (190/255, 1, 190/255), (190/255, 190/255, 255/255)]
    col_end = [ ( np.abs(U.col_R[i][0]-10/255), np.abs(U.col_R[i][1]-10/255), np.abs(U.col_R[i][2]-10/255) ) for i in U.range_seeds ]
    q = np.delete( np.linspace( 0, 1, 2*U.Q+3, endpoint=True ), [0, U.Q+1, -1], None )
    tflat = np.array([t, np.flip(t)]).flat
    colq_all = [ make_colquad( U.Q, col_start[i], col_end[i] ) for i in U.range_seeds ]
    eff_quan_all = [ mquantiles( eff_mcmc_all[i], q, alphap=0.5, betap=0.5, axis=0 ) for i in U.range_seeds ]
    flip_eq_all = [ np.flip(eff_quan) for eff_quan in eff_quan_all ]
    qi = np.around( 100*q, 1 ) # Labels
    lbs = [ "{}% - {}%".format(qi[i], qi[qi.size-i-1]) for i in range(U.Q) ]
    

    #
    # FIGURE 8a
    #
    U._updateStatus2( "Graphing Figure 8a", "" )
    fig8a = U.make_figure( "8a | Results: State Convergence", True )
    bins = np.int64( params.N / (10 if params.N >= 50 else 2) )
    num = U.len_seeds
    lw1 = 0.5
    lw2a = 0.4
    lw2b = lw2a - 0.1
    mum = 4
    a = 0.5
    
    # Red seed: HMC, MH, Flipper
    i = 0
    ax11 = fig8a.add_subplot( num, mum, 1,
                              ylim = (-10, 140),
                              title = "Algorithm Cum. Accep. (%)" )
    ax11.xaxis.set_major_locator(ticker.MultipleLocator( np.int64(params.N/4) ))
    ax11.grid( True, alpha=0.3 )
    ax11.plot( idxA, S[i].rec_ord[:, 1], '-', lw=lw1, color=U.col_HM, label="HMC" )
    ax11.plot( idxA, S[i].rec_ord[:, 2], '-', lw=lw1, color=U.col_MG, label="MG" )
    ax11.plot( idxA, S[i].rec_ord[:, 3], '-', lw=lw1, color=U.col_FL, label="FL" )
    ax11.legend( ncol=3, loc="upper right", framealpha=0, fontsize=U.fs_l-0.1 )

    # Red seed: K density
    ax21 = fig8a.add_subplot( num, mum, 2,
                              title = "Density by RNG Seed" )
    ax21.grid( False )
    ax21.hist( S[i].K_sz[idxB], bins=bins, density=True, stacked=True, color=U.col_R[i], label=U.RNGs[i].seed_str )
    ylim = ax21.get_ylim()
    if U.is_syn:
        ax21.plot( [U.TS.true_K_sz]*2, ylim, 'w--', lw=lw2a )
        ax21.plot( [U.TS.true_K_sz]*2, ylim, 'k--', lw=lw2b, label="True" )
    ax21.set_ylim(( ylim[0], 1.2*ylim[1] ))
    ax21.legend( ncol=2, loc="upper right", framealpha=0, fontsize=U.fs_l-0.1 )

    # Red seed: FRET efficiency trace
    ax31 = fig8a.add_subplot( num, mum, (3, 4),
                              title = "Best FRET Efficiency Trace with Posterior Quantiles (%)" )
    ax31.xaxis.set_major_locator(ticker.MultipleLocator( ticker_mult ))
    ax31.plot( t, eff_star, linewidth=0.3, color=U.col_m, alpha=a, label="Apparent" )
    plot_best_and_quantiles( U, ax31, i, t, tflat, eff_quan_all, flip_eq_all, colq_all, eff_medn_best_all, True )
    
    # Yellow seed: HMC, MH, Flipper
    i = 1
    ax12 = fig8a.add_subplot( num, mum, 5,
                              sharex = ax11,
                              sharey = ax11 )
    ax12.grid( True, alpha=0.3 )
    ax12.plot( idxA, S[i].rec_ord[:, 1], '-', lw=lw1, color=U.col_HM )
    ax12.plot( idxA, S[i].rec_ord[:, 2], '-', lw=lw1, color=U.col_MG )
    ax12.plot( idxA, S[i].rec_ord[:, 3], '-', lw=lw1, color=U.col_FL )

    # Yellow seed: K density
    ax22 = fig8a.add_subplot( num, mum, 6 )
    ax22.grid( False )
    ax22.hist( S[i].K_sz[idxB], bins=bins, density=True, stacked=True, color=U.col_R[i], label=U.RNGs[i].seed_str )
    ylim = ax22.get_ylim()
    if U.is_syn:
        ax22.plot( [U.TS.true_K_sz]*2, ylim, 'w--', lw=lw2a )
        ax22.plot( [U.TS.true_K_sz]*2, ylim, 'k--', lw=lw2b )
    ax22.set_ylim(( ylim[0], 1.2*ylim[1] ))
    ax22.legend( ncol=1, loc="upper right", framealpha=0, fontsize=U.fs_l-0.1 )

    # Yellow seed: FRET efficiency trace
    ax32 = fig8a.add_subplot( num, mum, (7, 8),
                              sharex = ax31 )
    ax32.plot( t, eff_star, linewidth=0.3, color=U.col_m, alpha=a )
    plot_best_and_quantiles( U, ax32, i, t, tflat, eff_quan_all, flip_eq_all, colq_all, eff_medn_best_all, False )
    
    # Green seed: HMC, MH, Flipper
    i = 2
    ax13 = fig8a.add_subplot( num, mum, 9,
                              sharex = ax11,
                              sharey = ax11)
    ax13.grid( True, alpha=0.3 )
    ax13.plot( idxA, S[i].rec_ord[:, 1], '-', lw=lw1, color=U.col_HM )
    ax13.plot( idxA, S[i].rec_ord[:, 2], '-', lw=lw1, color=U.col_MG )
    ax13.plot( idxA, S[i].rec_ord[:, 3], '-', lw=lw1, color=U.col_FL )

    # Green seed: K density
    ax23 = fig8a.add_subplot( num, mum, 10 )
    ax23.grid( False )
    ax23.hist( S[i].K_sz[idxB], bins=bins, density=True, stacked=True, color=U.col_R[i], label=U.RNGs[i].seed_str )
    ylim = ax23.get_ylim()
    if U.is_syn:
        ax23.plot( [U.TS.true_K_sz]*2, ylim, 'w--', lw=lw2a )
        ax23.plot( [U.TS.true_K_sz]*2, ylim, 'k--', lw=lw2b )
    ax23.set_ylim(( ylim[0], 1.2*ylim[1] ))
    ax23.legend( ncol=1, loc="upper right", framealpha=0, fontsize=U.fs_l-0.1 )

    # Green seed: FRET efficiency trace
    ax33 = fig8a.add_subplot( num, mum, (11, 12),
                              sharex = ax31 )
    ax33.plot( t, eff_star, linewidth=0.3, color=U.col_m, alpha=a )
    plot_best_and_quantiles( U, ax33, i, t, tflat, eff_quan_all, flip_eq_all, colq_all, eff_medn_best_all, False )
    
    # Blue seed: HMC, MH, Flipper
    i = 3
    ax14 = fig8a.add_subplot( num, mum, 13,
                              sharex = ax11,
                              sharey = ax11,
                              xlabel = "$n$: MCMC iteration" )
    ax14.grid( True, alpha=0.3 )
    ax14.plot( idxA, S[i].rec_ord[:, 1], '-', lw=lw1, color=U.col_HM )
    ax14.plot( idxA, S[i].rec_ord[:, 2], '-', lw=lw1, color=U.col_MG )
    ax14.plot( idxA, S[i].rec_ord[:, 3], '-', lw=lw1, color=U.col_FL )

    # Blue seed: K density
    ax24 = fig8a.add_subplot( num, mum, 14,
                              xlabel = "$K$: Number of states" )
    ax24.grid( False )
    ax24.hist( S[i].K_sz[idxB], bins=bins, density=True, stacked=True, color=U.col_R[i], label=U.RNGs[i].seed_str )
    ylim = ax24.get_ylim()
    if U.is_syn:
        ax24.plot( [U.TS.true_K_sz]*2, ylim, 'w--', lw=lw2a )
        ax24.plot( [U.TS.true_K_sz]*2, ylim, 'k--', lw=lw2b )
    ax24.set_ylim(( ylim[0], 1.2*ylim[1] ))
    ax24.legend( ncol=1, loc="upper right", framealpha=0, fontsize=U.fs_l-0.1 )

    # Blue seed: FRET efficiency trace
    ax34 = fig8a.add_subplot( num, mum, (15, 16),
                              sharex = ax31,
                              xlabel = f"Time in smFRET experiment ({params.units_t})" )
    ax34.plot( t, eff_star, linewidth=0.3, color=U.col_m, alpha=a )
    plot_best_and_quantiles( U, ax34, i, t, tflat, eff_quan_all, flip_eq_all, colq_all, eff_medn_best_all, False )

    if U.show:
        fig8a.show()
        fig8a.canvas.flush_events()
    fig8a.savefig(U.func_getActivePath( f"BayesFRET_fig08a state_convergence.png" ))
    U._updateStatus2( "Done!", "\n" )


    #
    # FIGURE 8b
    #
    U._updateStatus2( "Graphing Figure 8b", "" )
    fig8b = U.make_figure( "8b | Best FRET Efficiency Traces", False, 1 )
    ax = fig8b.add_subplot( 1, 5, (1, 4),
                            anchor = 'W',
                            ylim   = (0, 100),
                            yticks = range(0, 110, 10),
                            xlabel = f"Time in smFRET experiment ({params.units_t})",
                            ylabel = "FRET efficiency (%)" )
    ax.xaxis.set_major_locator(ticker.MultipleLocator( ticker_mult ))
    ax.grid(False)

    # Add the apparent FRET efficiency.
    P0, = ax.plot( t, eff_star, linewidth=0.3, color=U.col_m, alpha=a )
    allP = [P0]
    allL = ["Apparent"]

    # Add the posterior quantiles for each seed.
    P1 = []; P2 = []; P3 = []; P4 = []
    P = [P1, P2, P3, P4]
    for i in U.range_seeds:
        for j in range(U.Q):
            P[i].append( ax.fill( tflat, np.array([eff_quan_all[i][j], flip_eq_all[i][j]]).flat,
                                  color=colq_all[i][j], edgecolor=colq_all[i][j], alpha=0.5 )[0] )

    # Add best FRET efficiency traces.
    for i in U.range_seeds:
        p, = ax.plot( t, eff_medn_best_all[i], '-', linewidth=0.75, color=U.col_R[i] )
        allP.append(p)
        allL.append( f"Seed {U.RNGs[i].seed_str}, $n$: {ix_adj[i]}, $K$: {S[i].K_sz[ix_adj[i]]}" )

    # If synthetic, add true line on top of the traces.
    if U.is_syn:
        ax.plot( t, U.TS.true_line, 'w--', linewidth=lw2a ) # white background
        h3, = ax.plot( t, U.TS.true_line, 'k--', linewidth=lw2b ) # black dotted line
        allP.append(h3)
        allL.append("True")
    
    # Organize legends.
    ax.add_artist( ax.legend( allP, allL,   bbox_to_anchor=(1.02, 1.0),  loc="upper left", title="Best Traces", title_fontsize=7, fancybox=True ) )
    ax.add_artist( ax.legend( [], [],       bbox_to_anchor=(1.10, 0.65), loc="upper left", framealpha=0, title="Posterior Quantiles") )
    ax.add_artist( ax.legend( P1, [""]*U.Q, bbox_to_anchor=(1.02, 0.6),  loc="upper left", framealpha=0 ) )
    ax.add_artist( ax.legend( P2, [""]*U.Q, bbox_to_anchor=(1.07, 0.6),  loc="upper left", framealpha=0 ) )
    ax.add_artist( ax.legend( P3, [""]*U.Q, bbox_to_anchor=(1.12, 0.6),  loc="upper left", framealpha=0 ) )
    ax.add_artist( ax.legend( P4, lbs,      bbox_to_anchor=(1.17, 0.6),  loc="upper left", framealpha=0 ) )

    if U.show:
        fig8b.show()
        fig8b.canvas.flush_events()
    fig8b.savefig(U.func_getActivePath( f"BayesFRET_fig08b FRET_efficiency_best_traces.png" ))
    U._updateStatus2( "Done!", "\n" )
    



    


def plot_best_and_quantiles(U: Universal, ax: Axes, i: int, t, tflat, eff_quan_all, flip_eq_all, colq_all, eff_medn_best_all, add_label: bool ):
    """
    Plot the best FRET efficiency trace and the posterior quantiles of each simulation.

    Args:
        i (int): Simulation index in [0, 3].
    """
    ax.grid( False )
    for j in range(U.Q):
       ax.fill( tflat, np.array([eff_quan_all[i][j], flip_eq_all[i][j]]).flat,
                color=colq_all[i][j], edgecolor=colq_all[i][j], alpha=0.5 )
    ax.plot( t, eff_medn_best_all[i], '-o', markersize=0.3, lw=0.2, color=U.col_R[i] )
    if U.is_syn:
        ax.plot( t, U.TS.true_line, 'w--', linewidth=0.4 ) # white background
        ax.plot( t, U.TS.true_line, 'k--', linewidth=0.3 ) # black dotted line
    if add_label:
        ax.legend( loc="upper right", framealpha=0, fontsize=U.fs_l-0.1 )






def make_colquad(Q: np.int64, col_1: tuple, col_2: tuple):
    """
    Creates a set of colors (a gradient) between RGB tuples `col_1` and `col_2` by `Q` colors.

    ex. `col_1 = (0.5, 1, 188/255)`
    """
    return np.moveaxis(np.array([ np.linspace(col_1[0], col_2[0], Q, endpoint=True),
                                  np.linspace(col_1[1], col_2[1], Q, endpoint=True),
                                  np.linspace(col_1[2], col_2[2], Q, endpoint=True) ]), 0, 1)


