"""
BayesFRET: File 5
===
An Experiment-Adjusted HDP-HMM to Analyze Surface-Immobilized smFRET Data

About
---
`code04_chart.py`: The file to visualize progress of all simulations (`Chain_Visualize`).

Authors
---
Code by Jared Hidalgo. 

Inspired by MATLAB code from Ioannis Sgouralis, Shreya Madaan, Franky Djutanta, Rachael Kha, Rizal F. Hariadi, and Steve Pressé for "A Bayesian Nonparametric Approach to Single Molecule Förster Resonance Energy Transfer".
"""
# Import external packages for use.
import numpy as np
import matplotlib.ticker as ticker
# Import packages for documentation.
from matplotlib.axes import Axes
from matplotlib.lines import Line2D
from code01_classes import Universal, Sample, Params

df = 0.02
bins_active = 25
bins_still = 50

psz = 0.4
l_thin = 0.2
l_conv = 0.6
l_stat = 1






class Chain_Visualize:
    """
    Simulations for Figure 3
    ---
    * Before: Creates background.
    * During: Shows most recent sample.
    * After: Shows final sample + exports Figure 3.
    """
    def __init__(self, params: Params, samples: list[Sample], U: Universal, show_fig: bool, curr_samp: bool):
        """
        Creates the background of Figure 3.

        Args:
            show_fig (bool): If this Figure shows on a window or not.
            curr_samp (bool): If this class shows the current sample or not.
        """
        self.params = params
        num = U.len_seeds
        mum = 6
        typ = "Current" if curr_samp else "Final"
        self.fig3 = U.make_figure( f"3 | FRET Efficiency: {typ} Samples", True )
        if curr_samp:
            win = self.fig3.canvas.window()
            win.setFixedSize( win.size() )
        """Figure object for Figure 3."""
        t       = np.arange(0.5, params.T) * params.dt
        xlim_E  = (-1, params.T * params.dt + 1)
        xlim_C  = (0, 0.4)
        app_eff = 100 * params.It_A / (params.It_D + params.It_A + U.eps)


        # 1st RNG - FRET
        ax_rngE_1 = self.fig3.add_subplot( num, mum, (1, mum-1), # 4, 6, (1, 5)
                                           xlim        = xlim_E,
                                           ylim        = (-10, 130),
                                           ylabel      = "FRET eff. (%)" )
        ax_rngE_1.grid( True, alpha=0.3 )
        ax_rngE_1.xaxis.set_major_locator(ticker.MultipleLocator( np.int64(params.T*params.dt/10 )))
        ax_rngE_1.yaxis.set_major_locator(ticker.MultipleLocator( 50 ))
        #
        ax_rngE_1.plot( t, app_eff, linewidth=l_thin, color=U.col_m, label="Apparent" )
        if U.is_syn:
            ax_rngE_1.plot( t, U.TS.true_line, linewidth=l_thin, color='k', label="True" )
        label = f"RNG Seed {U.RNGs[0].seed_str}" + ("" if curr_samp else f", $n$: {params.N}, $K$: {samples[0].K_sz}")
        ax_rngE_1.plot( t, [-np.inf]*t, '-o', markersize=psz, linewidth=l_thin, color=U.col_R[0], label=label )
        ax_rngE_1.legend( ncol=3, loc='upper right', framealpha=U.fa, fontsize=U.fs_l )


        # 1st RNG - Additional
        ax_rngC_1 = self.fig3.add_subplot( num, mum, mum,  # 4, 6, 6
                                           xlim    = xlim_C,
                                           sharey  = ax_rngE_1 )
        ax_rngC_1.grid( True, alpha=0.3 )
        ax_rngC_1.xaxis.set_major_locator(ticker.MultipleLocator( 0.1 ))
        ax_rngC_1.yaxis.tick_right()
        ax_rngC_1.yaxis.set_label_position("right")
        # Results - Density
        ax_rngC_1.hist( app_eff, bins=params.T, density=True, orientation='horizontal', color=U.col_m )
        if U.is_syn:
            ax_rngC_1.hist( U.TS.true_line, bins=params.K_lim, density=True, orientation='horizontal', color='black', alpha=0.8 )


        # 2nd RNG - FRET
        ax_rngE_2 = self.fig3.add_subplot( num, mum, (mum+1, 2*mum-1), # 4, 6, (7, 11)
                                           sharex = ax_rngE_1,
                                           sharey = ax_rngE_1,
                                           ylabel = "FRET eff. (%)" )
        ax_rngE_2.grid( True, alpha=0.3 )
        #
        ax_rngE_2.plot( t, app_eff, linewidth=l_thin, color=U.col_m )
        if U.is_syn:
            ax_rngE_2.plot( t, U.TS.true_line, linewidth=l_thin, color='k' )
        label = f"RNG Seed {U.RNGs[1].seed_str}" + ("" if curr_samp else f", $n$: {params.N}, $K$: {samples[1].K_sz}")
        ax_rngE_2.plot( t, [-np.inf]*t, '-o', markersize=psz, linewidth=l_thin, color=U.col_R[1], label=label )
        ax_rngE_2.legend( ncol=1, loc='upper right', framealpha=U.fa, fontsize=U.fs_l )


        # 2nd RNG - Additional
        ax_rngC_2 = self.fig3.add_subplot( num, mum, 2*mum,  # 4, 6, 12
                                           sharex = ax_rngC_1,
                                           sharey = ax_rngC_1 )
        ax_rngC_2.grid( True, alpha=0.3 )
        ax_rngC_2.yaxis.tick_right()
        ax_rngC_2.yaxis.set_label_position("right")
        # Results - Density
        ax_rngC_2.hist( app_eff, bins=params.T, density=True, orientation='horizontal', color=U.col_m, label="Apparent" )
        if U.is_syn:
            ax_rngC_2.hist( U.TS.true_line, bins=params.K_lim, density=True, orientation='horizontal', color='black', alpha=0.8 )


        # 3rd RNG - FRET
        ax_rngE_3 = self.fig3.add_subplot( num, mum, (2*mum+1, 3*mum-1), # 4, 6, (13, 17)
                                           sharex = ax_rngE_1,
                                           sharey = ax_rngE_1,
                                           ylabel = "FRET eff. (%)" )
        ax_rngE_3.grid( True, alpha=0.3 )
        #
        ax_rngE_3.plot( t, app_eff, linewidth=l_thin, color=U.col_m )
        if U.is_syn:
            ax_rngE_3.plot( t, U.TS.true_line, linewidth=l_thin, color='k' )
        label = f"RNG Seed {U.RNGs[2].seed_str}" + ("" if curr_samp else f", $n$: {params.N}, $K$: {samples[2].K_sz}")
        ax_rngE_3.plot( t, [-np.inf]*t, '-o', markersize=psz, linewidth=l_thin, color=U.col_R[2], label=label )
        ax_rngE_3.legend( ncol=1, loc='upper right', framealpha=U.fa, fontsize=U.fs_l )


        # 3rd RNG - Additional
        ax_rngC_3 = self.fig3.add_subplot( num, mum, 3*mum,  # 4, 6, 18
                                           sharex = ax_rngC_1,
                                           sharey = ax_rngC_1 )
        ax_rngC_3.grid( True, alpha=0.3 )
        ax_rngC_3.yaxis.tick_right()
        ax_rngC_3.yaxis.set_label_position("right")
        # Results - Density
        ax_rngC_3.hist( app_eff, bins=params.T, density=True, orientation='horizontal', color=U.col_m, label="Apparent" )
        if U.is_syn:
            ax_rngC_3.hist( U.TS.true_line, bins=params.K_lim, density=True, orientation='horizontal', color='black', alpha=0.8 )


        # 4th RNG - FRET
        ax_rngE_4 = self.fig3.add_subplot( num, mum, (3*mum+1, 4*mum-1), # 4, 6, (19, 23)
                                           sharex = ax_rngE_1,
                                           sharey = ax_rngE_1,
                                           xlabel = f"Time in smFRET experiment ({params.units_t})",
                                           ylabel = "FRET eff. (%)" )
        ax_rngE_4.grid( True, alpha=0.3 )
        # 
        ax_rngE_4.plot( t, app_eff, linewidth=l_thin, color=U.col_m )
        if U.is_syn:
            ax_rngE_4.plot( t, U.TS.true_line, linewidth=l_thin, color='k' )
        label = f"RNG Seed {U.RNGs[3].seed_str}" + ("" if curr_samp else f", $n$: {params.N}, $K$: {samples[3].K_sz}")
        ax_rngE_4.plot( t, [-np.inf]*t, '-o', markersize=psz, linewidth=l_thin, color=U.col_R[3], label=label )
        ax_rngE_4.legend( ncol=1, loc='upper right', framealpha=U.fa, fontsize=U.fs_l )


        # 4th RNG - Additional
        ax_rngC_4 = self.fig3.add_subplot( num, mum, 4*mum,  # 4, 6, 24
                                           sharex = ax_rngC_1,
                                           sharey = ax_rngC_1,
                                           xlabel = "Density" )
        ax_rngC_4.grid( True, alpha=0.3 )
        ax_rngC_4.yaxis.tick_right()
        ax_rngC_4.yaxis.set_label_position("right")
        # Results - Density
        ax_rngC_4.hist( app_eff, bins=params.T, density=True, orientation='horizontal', color=U.col_m )
        if U.is_syn:
            ax_rngC_4.hist( U.TS.true_line, bins=params.K_lim, density=True, orientation='horizontal', color='black', alpha=0.8 )

        # Combine by seeds
        self.axs_eff = [ax_rngE_1, ax_rngE_2, ax_rngE_3, ax_rngE_4]
        """List of Axes for the FRET efficiencies."""
        self.lns_eff = [ax.get_lines() for ax in self.axs_eff]
        """For each Axis, a list of Lines."""
        self.axs_den = [ax_rngC_1, ax_rngC_2, ax_rngC_3, ax_rngC_4]
        """List of Axes for the densities of FRET efficiencies."""
        
        if show_fig:
            U.plt.show(block=False)
            U.plt.pause(0.1)
        
        self.bg = self.fig3.canvas.copy_from_bbox( self.fig3.bbox )
        """Copy of background of Figure 3 for updating."""
        self.fig3.canvas.blit( self.fig3.bbox )
        self.fig3.canvas.flush_events()
    

    
    def update(self, sample: Sample, U: Universal, v: np.int64 ):
        """
        Updates Figure 3 with blitting.

        Args:
            v (int): Simulation version in `axs_eff` and `lns_eff`.
        """
        self.fig3.canvas.restore_region( self.bg )

        data = 100 * sample.eff[sample.st]
        self.set_data( self.axs_eff[v], self.lns_eff[v][U.line_indx], ydata=data ) # Left
        self.draw_hist( self.axs_den[v], data, U.col_R[v] ) # Right

        self.fig3.canvas.blit( self.fig3.bbox )
        self.fig3.canvas.flush_events()
    


    def restore(self, samples: list[Sample], U: Universal):
        """
        Full update for Figure 3 via blitting.
        """
        self.fig3.canvas.restore_region( self.bg )
        status = "Exporting Figure 3 ({:d}/5)"

        for v, sample in enumerate(samples):
            U._updateStatus2( status.format(v+1), "\r" )
            data = 100 * sample.eff[sample.st]
            self.set_data( self.axs_eff[v], self.lns_eff[v][U.line_indx], ydata=data ) # Left
            self.draw_hist( self.axs_den[v], data, U.col_R[v] ) # Right
            self.axs_eff[v].figure.canvas.draw()
            self.axs_den[v].figure.canvas.draw()
        
        self.fig3.canvas.blit( self.fig3.bbox )
        self.fig3.canvas.flush_events()
    


    def set_data(self, ax: Axes, line: Line2D, xdata = [], ydata = []):
        xempty = len(xdata) == 0
        yempty = len(ydata) == 0
        if xempty and yempty:
            raise Exception("INTERNAL ERROR! `set_data(..., xdata, ydata)` in `code04_chart.py` must have one or both inputs of data.")
        else:
            if not xempty: line.set_xdata( xdata )
            if not yempty: line.set_ydata( ydata )
            ax.draw_artist( line )
    


    def draw_hist(self, ax: Axes, data: np.ndarray, color: tuple):
        n, bins, rects = ax.hist( data, bins=self.params.K_lim, density=True, orientation='horizontal', color=color, ec=color, alpha=0.5 )
        for r in rects: ax.draw_artist( r )