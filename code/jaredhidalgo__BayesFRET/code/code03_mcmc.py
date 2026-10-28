"""
BayesFRET: File 4
===
An Experiment-Adjusted HDP-HMM to Analyze Surface-Immobilized smFRET Data

About
---
code03_mcmc.py: The file to handle the chain simulations (`Chain_Main`). Uses code04_chart.py and code05_update.py.

Author
---
Code by Jared Hidalgo. 

Inspired by MATLAB code from Ioannis Sgouralis, Shreya Madaan, Franky Djutanta, Rachael Kha, Rizal F. Hariadi, and Steve Pressé for "A Bayesian Nonparametric Approach to Single Molecule Förster Resonance Energy Transfer".
"""
# Import internal packages for use.
from copy import deepcopy
import datetime, pickle, time
# Import external packages for use.
import numpy as np
# Import program files for use.
from code04_chart import Chain_Visualize
from code05_update import sampler_update
# Import program files for documentation.
from code01_classes import Universal, Sample, Params, Chain_History






class Chain_Main:
    """
    During simulations: Updates Figure 2, Chain_History objects.

    After simulations: Stores all simulation data + Figure 2.
    """
    runtimes = []

    def __init__(self, params: Params, samples: list[Sample], U: Universal):
        """
        Args:
            samples (list[Sample]): The current samples of each seed with respect to the current MCMC iteration `n`.
        """
        lstICO = [ f" - {i} of 6.ico" for i in range(1, 6) ]
        if U.show: G = Chain_Visualize( params, samples, U, True, True )
        history = [ Chain_History( params, samples[i], U.RNGs[i].seed_int ) for i in U.range_seeds ]
        
        digits_n = np.int64(np.log10(params.N)) + 1
        digits_s = np.int64(np.max([ np.int64(np.log10(r.seed_int))+1 for r in U.RNGs[:-1] ]))
        digits_sValue = "4" if digits_s <= 4 else str(digits_s)
        digits_K = np.int64(np.log10(params.K_lim)) + 1

        print(f"Running simulations............................... Stop @ n = {params.N}.")
        self.prt1 = "".join([" "]*(digits_s-4)) + "Seed  |  " + "".join([" "]*(digits_n-1)) + "n  |       HMC  |        MG  |  Flipper  |  " + "".join([" "]*(digits_K-1)) + "K  |  Runtime"
        print( self.prt1 )
        self.prt1 = self.prt1[:-12]
        U.menu_t3.setText( self.prt1 )
        self.prt2_f = str("{:" + digits_sValue + "d}  |  {:" + str(digits_n) + "d}  |  {:3.2f}% {}  |  {:3.2f}% {}  |  {:3.2f}% {}  |  {:" + str(digits_K) + "d}  ")

        # For each distinct MCMC simulation...
        for v, sample in enumerate(samples):
            U.tray_state = lstICO[v]
            U.tray.setIcon( U.func_getIcon( f"{U.THEME}{U.tray_state}" ) )
            tme = time.time()
            first_samp = deepcopy( sample )

            for n in range(1, params.N+1):
                # Sample for each step in the MCMC simulation.
                sampler_update( sample, params, U, v )
                # Option to record samples at certain intervals with `params.i_skip`.
                if ( n % params.i_skip == 0 ):
                    history[v].st   [n] = sample.st
                    history[v].ft_D [n] = sample.ft_D
                    history[v].ft_A [n] = sample.ft_A
                    history[v].eff  [n] = sample.eff
                    history[v].eff_select.append(sample.eff[sample.K_set])
                    history[v].K_set.append(sample.K_set)
                    history[v].K_cnt.append(sample.K_cnt)
                    history[v].K_loc.append(sample.K_loc)
                    history[v].K_sz [n] = sample.K_sz
                    history[v].wi_D [n] = sample.wi_D
                    history[v].wi_A [n] = sample.wi_A
                    history[v].bm   [n] = sample.bm
                    history[v].pm   [n] = sample.pm.flatten('F')
                    history[v].ps   [n] = sample.ps
                    history[v].ksi_D[n] = sample.ksi_D
                    history[v].ksi_A[n] = sample.ksi_A
                    history[v].lam_Z[n] = sample.lam_Z
                    history[v].lam_D[n] = sample.lam_D
                    history[v].lam_A[n] = sample.lam_A
                    history[v].tht  [n] = sample.tht
                    history[v].rho_D[n] = sample.rho_D
                    history[v].rho_A[n] = sample.rho_A
                    history[v].kap_D[n] = sample.kap_D
                    history[v].kap_A[n] = sample.kap_A
                    history[v].kap_Z[n] = sample.kap_Z
                    
                    history[v].rec_ord[n] = [ n, 
                                                sample.rec[0, 0]/sample.rec[0, 1] * 100, 
                                                sample.rec[1, 0]/sample.rec[1, 1] * 100, 
                                                sample.rec[2, 0]/sample.rec[2, 1] * 100 ]
                    if U.show: G.update( sample, U, v )
                    self.func_updateStatus3( history[v].rec_ord, n, U.RNGs[v].seed_int, history[v].K_sz[n], U )
                    print("\r", end="")

            # If showing graphs, save background.
            if U.show: G.bg = G.fig3.canvas.copy_from_bbox( G.fig3.bbox )

            # Save runtime and finish simulation.
            runtime = str(datetime.timedelta(seconds = time.time() - tme))
            self.func_updateStatus3( history[v].rec_ord, n-1, U.RNGs[v].seed_int, history[v].K_sz[-1], U )
            print(f"|  {runtime}")
            history[v].runtime = runtime
            self.runtimes.append(runtime)

            # Export data (pickle).
            with open( U.filenames[v], "wb" ) as file:
                pickle.dump( history[v], file )
            
            # Export data (txt).
            with open( U.func_getActivePath( f"BayesFRET_data_sim_{U.RNGs[v].seed_str}_preview.txt" ), "w+" ) as file:
                file.write( f"{U.title}\nSAMPLE DATA FOR SIMULATION WITH RNG SEED {U.RNGs[v].seed_str}\n" +  
                            f"====================================================================================\n\n\n" + 
                            f"Runtime: {runtime}\n\n\nFirst Sample\n------------\n" )
                for obj in [first_samp]:
                    for attr in vars(obj):
                        file.write( f"{attr}: {getattr(obj, attr, 'NaN')}\n\n\n" )
                a = f"Final Sample (n = {params.N})"
                b = "".join( ["-"]*len(a) )
                file.write( f"\n\n{a}\n{b}\n" )
                for obj in [sample]:
                    for attr in vars(obj):
                        file.write( f"{attr}: {getattr(obj, attr, 'NaN')}\n\n\n" )
            del first_samp
        
        print()
        if U.show:
            # Close temporary window for Figure 3.
            U.plt.close( G.fig3 )
        U.menu_t3.setText( "-- Status --" )
        # (Re)make Figure 3.
        G_done = Chain_Visualize( params, samples, U, U.show, False )
        G_done.restore( samples, U )
        # Notify.
        U._updateStatus2( "Exporting Figure 3 (5/5)", "" )
        G_done.fig3.savefig( U.func_getActivePath( f"BayesFRET_fig03 FRET_efficiency_final_traces.png" ) )
        U._updateStatus2( "Done!", "\n" )
        # Update tray icon.
        U.tray_state = lstICO[-1]
        U.tray.setIcon( U.func_getIcon( f"{U.THEME}{U.tray_state}" ) )
    


    def func_updateStatus3(self, rec_ord: np.ndarray[np.float64], n: np.int64, s: np.int64, K_sz: np.int64, U: Universal):
        """
        Standard for printing the progress of each iteration on the command line, context menu, and tray tooltip.

        Args:
            rec_ord (array[float]): Copy of cumulative ratios to accept/run MG (1) + HMC (2) + Flipper (3) algorithms by time t (0).
            n (int): Current iteration of chain.
            s (int): The seed of the simulation
            K_sz_med (int): The median of the number of conformational states up to `n`.
        """
        p  = rec_ord[n]
        pr = rec_ord[n-1]
        p1 = u"\u2191" if p[1] > pr[1] else u"\u2193" if p[1] < pr[1] else u"\u003d"
        p2 = u"\u2191" if p[2] > pr[2] else u"\u2193" if p[2] < pr[2] else u"\u003d"
        p3 = u"\u2191" if p[3] > pr[3] else u"\u2193" if p[3] < pr[3] else u"\u003d"
        prt2 = self.prt2_f.format( s, np.int64(p[0]), p[1], p1, p[2], p2, p[3], p3, K_sz )
        print( prt2, end="" )
        U.menu_run.setText( prt2 )
        U.tray.setToolTip( self.prt1 + "\n" + prt2 )


