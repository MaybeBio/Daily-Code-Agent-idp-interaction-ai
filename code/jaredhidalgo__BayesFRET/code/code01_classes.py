"""
BayesFRET: File 2
===
An Experiment-Adjusted HDP-HMM to Analyze Surface-Immobilized smFRET Data

About
---
code01_classes.py: The file to host custom classes used across the rest of the files.

Classes: `RNG`, `True Params`, `Universal`, `Params`, `Sample`, `Chain_History`

Author
---
Code by Jared Hidalgo. 

Inspired by MATLAB code from Ioannis Sgouralis, Shreya Madaan, Franky Djutanta, Rachael Kha, Rizal F. Hariadi, and Steve Pressé for "A Bayesian Nonparametric Approach to Single Molecule Förster Resonance Energy Transfer".
"""
# Import internal packages for use.
import os, platform
# Import external packages for use.
import numpy as np
from numpy.matlib import repmat
import matplotlib.pyplot as plt
# Import external packages for documentation.
from PySide6.QtGui import (QAction, QIcon)
from PySide6.QtWidgets import QSystemTrayIcon
from matplotlib.figure import Figure


DPI = 200
plt.rcParams.update({
    "font.size": 5,
    "figure.dpi": DPI,
    "xtick.labelsize": 4,
    "ytick.labelsize": 4
})





class RNG:
    """
    Manage random number generators in MATLAB style.
    """


    
    def __init__(self, seed: np.int64 = None):
        self.seed_int = seed
        """The specific RNG seed (int) for the Generator class."""
        self.seed_str = str(self.seed_int) if seed != None else "[None]"
        """The specific RNG seed (str) for the Generator class."""
        self.rng = np.random.default_rng( seed )
        """The Generator class hosting the RNG functions."""
        self.matx_randg = np.vectorize( self.rng.gamma )
        """Vectorization of the Gamma distribution."""
        self.matx_poiss = np.vectorize( self.rng.poisson )
        """Vectorization of the Poisson distribution."""
    


    def rand1(self):
        """Draws a sample from the “continuous uniform” distribution over [0.0, 1.0)."""
        return self.rng.random()
    


    def rand(self, sz: np.int64|tuple[np.int64]):
        """
        Draws samples from the “continuous uniform” distribution over [0.0, 1.0).
        
        Args:
            sz: The length (`np.int64`) or shape (`tuple[np.int64]`) of the resulting array. Value(s) must be `> 1`.
        """
        return self.rng.random( size = sz )
    


    def randi1(self, low: np.int64, high: np.int64, endpoint = False) -> np.int64:
        """
        Draws a random integer from the "discrete uniform" distribution.
        
        Args:
            low: Minimum integer.
            high: Maximum integer.
            endpoint: Include `high`, optional.
        """
        return self.rng.integers( low, high, endpoint = endpoint )
    


    def randi(self, low: np.int64, high: np.int64, sz: np.int64|tuple[np.int64], endpoint = False) -> np.ndarray[np.int64]:
        """
        Draws random integers from the "discrete uniform" distribution.

        Args:
            low: Minimum integer.
            high: Maximum integer.
            sz: Length (np.int64) or shape (tuple[np.int64]) of the resulting array.
            endpoint: Include `high`, optional.
        """
        return self.rng.integers( low, high, size = sz, endpoint = endpoint )
    


    def randn1(self):
        """Draws a sample from the standard Normal distribution over [0.0, 1.0)."""
        return self.rng.standard_normal()
    


    def randn(self, sz: np.int64|tuple[np.int64]) -> np.ndarray[np.float64]:
        """
        Draws samples from the standard Normal distribution (mean=0, stdev=1).

        Args:
            sz: The length (`np.int64`) or shape (`tuple[np.int64]`) of the resulting array. Value(s) must be `> 1`.
        """
        return self.rng.standard_normal( size = sz )
    


    def poissrnd(self, lambd: np.ndarray[np.float64]) -> np.ndarray[np.float64]:
        """
        Draws samples from Pois(`lambd`).
        
        Args:
            lambd: Parameter λ of the Poisson distribution. `lambd > 0`.
        """
        return self.matx_poiss( lambd )
    


    def betarnd(self, a: np.ndarray[np.float64|np.bool_], b: np.ndarray[np.float64|np.bool_]) -> np.ndarray[np.float64]:
        """
        Draws samples from Beta(`a`, `b`). `a`, `b` can be arrays of the same size or simply floats.
        
        Args:
            a: Shapes of the Beta distribution α (`all(a > 0)`)
            b: Scales of the Beta distribution β (`all(b > 0)`)
        """
        return self.rng.beta( a, b )
    


    def gamrnd1(self, a: np.float64, b = 1.0):
        """
        Draws a sample from Gamma(`a`, `b`). If scale = 1, then leave `b` blank.
        
        Args:
            a: Shape of the Gamma distribution α (`a > 0`)
            b: Scale of the Gamma distribution β (`b > 0`), optional
        """
        return self.rng.gamma( a, b )
    


    def gamrnd(self, a: np.float64, b: np.float64, sz: np.int64|tuple[np.int64]) -> np.ndarray[np.float64]:
        """
        Draws samples from Gamma(`a`, `b`).
        
        Args:
            a: Shape of the Gamma distribution α (`a > 0`)
            b: Scale of the Gamma distribution β (`b > 0`)
            sz: The length (`np.int64`) or shape (`tuple[np.int64]`) of the resulting array. Value(s) must be `>= 2`.
        """
        return self.rng.gamma( a, b, size = sz )
    


    def dirrnd(self, a: np.ndarray[np.float64]) -> np.ndarray[np.float64]:
        """
        Custom function that draws samples from the Dirichlet distribution.
        """
        p = self.matx_randg(a)
        return p / (np.sum(p) if p.shape[0] == p.size else np.sum(p, 1)[:, None])






class True_Samples:
    """
    Save and categorize true values of the samples. For synthetic data only.
    """
    true_line: np.ndarray[np.float64]
    """True FRET efficiencies for graphing."""
    true_st: np.ndarray[np.float64]
    """s_n True: Conformational states of the observed molecule."""
    true_K_set: np.ndarray[np.int64]
    """set(K) True: The set of unique conformational states from `unique(sample.st)`."""
    true_K_loc: np.ndarray[np.int64]
    """location(K) True: The locations of the unique conformational states from `unique(sample.st, return_inverse=True)`."""
    true_K_cnt: np.ndarray[np.int64]
    """count(K) True: The counts of the unique conformational states from `unique(sample.st, return_counts=True)."""
    true_K_sz: np.int64
    """K True: Number of true conformational states."""
    true_ft_D: np.ndarray[np.bool_]
    """f_t^D True: Photo-states of the donor dye throughout the experiment."""
    true_ft_A: np.ndarray[np.bool_]
    """f_t^A True: Photo-states of the acceptor dye throughout the experiment."""
    true_tht: np.float64
    """θ True: Photoemission multiplier."""
    true_rho_D: np.float64
    """ρ^D True: Scaling factor of background emission rates for the donor dye."""
    true_rho_A: np.float64
    """ρ^A True: Scaling factor of background emission rates for the acceptor dye."""
    true_kap_D: np.ndarray[np.float64]
    """κ^D True: Scaling factor of photoemission rates for the donor dye with FRET."""
    true_kap_A: np.ndarray[np.float64]
    """κ^A True: Scaling factor of photoemission rates for the acceptor dye with FRET."""
    true_kap_Z: np.float64
    """κ^Z True: Scaling factor of photoemission rate for the donor dye without FRET."""
    true_ksi_D: np.float64
    """ξ^D True: Prior of background emission rates for the donor dye."""
    true_ksi_A: np.float64
    """ξ^A True: Prior of background emission rates for the acceptor dye."""
    true_lam_D: np.ndarray[np.float64]
    """λ^D True: Prior of photoemission rates for the donor dye with FRET."""
    true_lam_A: np.ndarray[np.float64]
    """λ^A True: Prior of photoemission rates for the acceptor dye with FRET."""
    true_lam_Z: np.float64
    """λ^Z True: Prior of photoemission rate for the donor dye without FRET."""
    true_wi_D: list[np.float64]
    """ω^D True: Transition probabilities of the photo-states for the donor dye."""
    true_wi_A: list[np.float64]
    """ω^A True: Transition probabilities of the photo-states for the acceptor dye."""
    true_pm: np.ndarray[np.float64]
    """π True: The transition matrix between K states.\n\nSize: `K*^2`"""
    true_ps: np.ndarray[np.float64]
    """π* Current: The probability array that state K will occur in the initial state.\n\nSize: `K`"""



    def __init__(self):
        pass






class Universal:
    """
    Object to carry variables + functions across files + classes.
    """
    plt = plt
    OS = platform.system()
    """The current operating system."""
    DIR_PROGRAM = os.getcwd()
    """Current directory of the BayesFRET program.\n\nCHANGE THIS IF FOLDER MOVED AFTER SIMULATION!"""
    dpi = DPI
    """Dots per inch of the figure."""
    eps = np.spacing(1) # = np.finfo(float).eps

    col_M  = (0, 0, 0.5)
    """Color for main molecule - Blue"""
    col_m  = (0, 1, 1)
    """Color for auxiliary molecule - Light blue"""
    col_D  = (0.3, 1.0, 0.3)
    """Color for main donor - Light green"""
    col_A  = (1.0, 0.3, 0.3)
    """Color for main acceptor - Light red"""
    col_HM = (0, 0, 0)
    """Color for Hamiltonian Monte Carlo - Black"""
    col_MG = (1, 0, 1)
    """Color for Metropolis-Hastings - Pink"""
    col_FL = (1, 0.5, 0)
    """Color for Flipper - Orange"""
    col_R  = [(1, 0, 0),
              (204/255, 204/255, 0.0),
              (0, 204/255, 0),
              (0, 0, 1)]
    """Colors for random number generators - Red, Yellow, Green, Blue"""
    fs_l   = 5
    """Font size for legend text."""
    fa     = 0
    """Alpha for legend frame."""
    figs: list[Figure] = []



    def __init__( self, title: str, show: bool, large: bool, is_syn: bool, RNGs: list[RNG], ad: str, rd: str, theme: str, 
                  Q: np.int64, bi: np.float64, rep_tht: np.int64, rep_bm: np.int64, MG_L: np.int64, HMC_L: np.int64, HMC_eps: np.float64, 
                  mt3: QAction, mr: QAction, tray: QSystemTrayIcon, ts: str ):
        self.title = title
        """Title of the BayesFRET program."""
        self.show = show
        """Boolean to open windows for figures (`True`) or not (`False`)."""
        self.large = large
        """Boolean if figures are large (`True`) or not (`False`)."""
        self.w = 3200 if large else 1280
        """Width of figures."""
        self.h = 1800 if large else  720
        """Height of figures."""
        self.is_syn = is_syn
        """Boolean if data is synthetic (`True`) or experimental (`False`)."""
        self.TS = True_Samples() if is_syn else None
        """Collection of true parameters (`True_Samples`). For synthetic data + internal use only. `None` when exported."""
        self.line_indx = 2 if is_syn else 1
        """Index of line for sampled data. 1 if experimental. 2 if synthetic."""
        self.dt = "syn" if is_syn else "exp"
        """Data type for image filename."""
        self.RNGs = RNGs
        """List of RNG classes. `RNGs[-1]` is strictly to generate synthetic data (no seed). The rest are for each simulation."""
        self.len_seeds = len(RNGs) - 1
        """Number of RNG seeds, one for each simulation, not including `RNGs[-1]`."""
        self.range_seeds = range(self.len_seeds)
        """Range for seed values (0-3)."""
        self.list_seeds = [RNGs[i].seed_str for i in self.range_seeds]
        """List of seed values in string format."""
        self.DIR_ACTIVE = ad
        """The full path to the active directory to store output. ex. C:\\...\\parentDir\\BayesFRET_{curr_date}\n\nCHANGE THIS IF FOLDER MOVED AFTER SIMULATION!"""
        self.DIR_RESOURCES = rd
        """Current directory of the BayesFRET program.\n\nCHANGE THIS IF THE BAYESFRET FOLDER IS MOVED!"""
        self.filenames = [ self.func_getActivePath( f"BayesFRET_data_sim_{self.RNGs[v].seed_str}_history.p" ) for v in self.range_seeds ]
        """The filenames of the RNG history pickle files."""
        self.THEME = theme
        """The dark/light theme for icons: 'BayesFRET dark icon' or 'BayesFRET light icon'"""
        self.Q = Q
        """The number of posterior quantiles to show in the last graph. Posterior quantilies demonstrate the sampling frequency of the simulation's FRET efficiency trace."""
        self.burn_in = bi
        """Ignores the first percentage of results to hide outliers."""
        self.rep_tht = rep_tht
        """Number of times to calculate `tht`: θ ~ Gamma( φ_θ, ψ_θ ) and other photoemission priors."""
        self.rep_bm = rep_bm
        """Number of times to calculate `bm`: β ~ Dir( γ/K, ..., γ/K )."""
        self.MG_L = MG_L
        """Length to repeat the Metropolis-Hastings algorithm."""
        self.HMC_L = HMC_L
        """Length to run the Hamiltonian Monte Carlo algorithm."""
        self.HMC_eps = HMC_eps
        """Leap size for the Hamiltonian Monte Carlo algorithm."""
        self.menu_t3 = mt3
        """Topmost context menu button ("Status") to update status. `None` when exported."""
        self.menu_run = mr
        """Second topmost context menu button ("Run") to update status. `None` when exported."""
        self.tray = tray
        """Icon for system tray. `None` when exported."""
        self.tray_state = ts
        """Second half of QIcon filename (Universal)."""
    


    def get_samples(self, p: np.ndarray, v: np.int64):
        """
        Categorical distribution for conformational states based on probabilities.
        
        Args:
            p: The probability array.
            v: The version of the simulation. `v=-1` is strictly for synthetic data.
        """
        ind = np.argsort(-1*p)
        P = np.cumsum( np.sort(p)[::-1] )
        return ind[ np.argmax( P[-1] * self.RNGs[v].rand1() <= P ) ]



    def make_figure(self, title: str, tl: bool, pad=0):
        """
        Start to make window for Figure.

        Args:
            title: Only for `self.show=True`.
            tl: Tight layout or not.
            pad: Padding of tight layout (only for Figure 5).
        """
        if self.show:
            fig = self.plt.figure( title, dpi=self.dpi, tight_layout=tl ) # NOTE: Do we need self.dpi??
            if self.large:  self.plt.get_current_fig_manager().window.showMaximized()
            else:           self.plt.get_current_fig_manager().resize( self.w, self.h )
        else:
            fig = self.plt.figure( title, dpi=self.dpi, tight_layout=tl, figsize=( self.w/self.dpi, self.h/self.dpi ) )
        self.plt.get_current_fig_manager().window.setWindowIcon( self.func_getIcon( f"{self.THEME}.ico" ) )
        if pad != 0:
            self.plt.tight_layout( pad=pad )
        self.figs.append(fig)
        return fig



    def func_getActivePath(self, fname: str):
        """Returns the full filepath of the active directory merged with the file `fname`."""
        return os.path.join( self.DIR_ACTIVE, fname )
    


    def func_getIcon(self, file: str):
        """
        Shortcut to return formatted icon from the resources folder.
        """
        return QIcon( os.path.join( self.DIR_RESOURCES, file ) )



    def _updateStatus2(self, status: str, end: str):
        """
        Updates the one-time status on the command line, context menu, and tray tooltip. DO NOT USE OUTSIDE OF SIMULATION.
        """
        if "Done!" not in status:
            dots = "."*(50 - len(status))
            print(f"{status}{dots} ", end=end)
            self.menu_run.setText( status )
            self.tray.setToolTip( status )
        else:
            print(f"{status}", end=end)
            self.menu_run.setText( self.menu_run.text() + ": Done!" )
            self.tray.setToolTip( self.tray.toolTip() + ": Done!" )



    def func_delForPickle(self):
        """
        Remove objects that can't be pickled.
        """
        del(self.TS)
        del(self.menu_t3)
        del(self.menu_run)
        del(self.tray)
        del(self.tray_state)






class Params:
    """
    Stores and labels all the static parameters.
    """
    i_skip  = 1
    """Accepts sample on `n`th iteration when `n mod i_skip = 0`. `i_skip = 1` means no skipping."""
    MG_a_D = 250
    """α^D for the proposed background emission multiplier for the donor dye in the MG algorithm: ρ^D_* = ρ^D_{n-1} [ Gamma( α^D, 1 ) / α^D ]."""
    MG_a_A = 250
    """α^A for the proposed background emission multiplier for the acceptor dye in the MG algorithm: ρ^A_* = ρ^A_{n-1} [ Gamma( α^A, 1 ) / α^A ]."""



    def __init__( self, It_D: np.ndarray[np.float64], It_A: np.ndarray[np.float64], T: np.int64, units_t: str, units_I: str, 
                  dt: np.float64, dD: np.float64, cDD: np.float64, cAA: np.float64, qD: np.float64, qA: np.float64, 
                  N: np.int64, K_lim: np.int64, alpha: np.float64, gamma: np.float64,
                  wi_D_prior_eta: np.ndarray[np.float64], wi_D_prior_zeta: np.ndarray[np.float64], wi_A_prior_eta: np.ndarray[np.float64], wi_A_prior_zeta: np.ndarray[np.float64],
                  rho_D_prior_phi: np.float64, rho_D_prior_psi: np.float64, rho_A_prior_phi: np.float64, rho_A_prior_psi: np.float64, tht_prior_phi: np.float64, tht_prior_psi: np.float64,
                  kap_D_prior_phi: np.float64, kap_D_prior_psi: np.float64, kap_A_prior_phi: np.float64, kap_A_prior_psi: np.float64, kap_Z_prior_phi: np.float64, kap_Z_prior_psi: np.float64 ):
        self.It_D = It_D
        """Original data: Set of photon intensities of donor dye.\n\nSize: `T`"""
        self.It_A = It_A
        """Original data: Set of photon intensities of acceptor dye.\n\nSize: `T`"""
        self.T = T
        """Size of the dataset (`size(It_D) = size(It_A)`)."""
        self.units_t = units_t
        """Unit of time for figures."""
        self.units_I = units_I
        """Unit of intensity for figures."""
        self.dt = dt
        """δt: Frame rate (a.k.a. measurement acquisition period, or the frequency of data) in seconds, regardless of the unit of time."""
        self.dD = dD
        """δτ = δt (1-d): Exposure period, or how much time within the frame rate (δt/`dt`) is dedicated to capturing the intensity using the dead time (d)."""
        self.cDD = cDD
        """Cross-talk proportion from donor dye to donor channel. `cDD + cDA = 1`"""
        self.cDA = 1 - cDD
        """Cross-talk proportion from donor dye to acceptor channel. `cDD + cDA = 1`"""
        self.cAA = cAA
        """Cross-talk proportion from acceptor dye to acceptor channel. `cAA + cAD = 1`"""
        self.cAD = 1 - cAA
        """Cross-talk proportion from acceptor dye to donor channel. `cAA + cAD = 1`"""
        self.qD = qD
        """Detector quantum efficiency (photodetection percentage) for the donor dye."""
        self.qA = qA
        """Detector quantum efficiency (photodetection percentage) for the acceptor dye."""
        self.N = N
        """Length (a.k.a. number of iterations) for each simulation."""
        self.K_lim = K_lim
        """K: The weak limit to estimating the number of conformational states. Also contributes to the nonparametric hierarchical prior `bm`: β ~ Dir( γ/K, ..., γ/K )."""
        self.alpha = alpha
        """Concentration parameter α for the conformational state probability prior `ps`: ~π ~ DP( α, β )."""
        self.gamma = gamma
        """γ for the nonparametric hierarchical prior `bm`: β ~ Dir( γ/K, ..., γ/K )."""
        self.wi_D_prior_eta = wi_D_prior_eta
        """Shape η^D for photo-state transition priors of the donor dye: ω^D ~ Beta( η^D, ζ^D ).\n\nSize 3"""
        self.wi_D_prior_zeta = wi_D_prior_zeta
        """Scale ζ^D for photo-state transition priors of the donor dye: ω^D ~ Beta( η^D, ζ^D ).\n\nSize 3"""
        self.wi_A_prior_eta = wi_A_prior_eta
        """Shape η^A for photo-state transition priors of the acceptor dye: ω^A ~ Beta( η^A, ζ^A ).\n\nSize 3"""
        self.wi_A_prior_zeta = wi_A_prior_zeta
        """Scale ζ^A for photo-state transition priors of the acceptor dye: ω^A ~ Beta( η^A, ζ^A ).\n\nSize 3"""
        self.rho_D_prior_phi = rho_D_prior_phi
        """φ^D_ρ in scaling factor of background emission rates for the donor dye: ρ^D ~ Gamma( φ^D_ρ, ψ^D_ρ / φ^D_ρ )."""
        self.rho_D_prior_psi = rho_D_prior_psi
        """ψ^D_ρ in scaling factor of background emission rates for the donor dye: ρ^D ~ Gamma( φ^D_ρ, ψ^D_ρ / φ^D_ρ )."""
        self.rho_A_prior_phi = rho_A_prior_phi
        """φ^A_ρ in scaling factor of background emission rates for the acceptor dye: ρ^A ~ Gamma( φ^A_ρ, ψ^A_ρ / φ^A_ρ )."""
        self.rho_A_prior_psi = rho_A_prior_psi
        """ψ^A_ρ in scaling factor of background emission rates for the acceptor dye: ρ^A ~ Gamma( φ^A_ρ, ψ^A_ρ / φ^A_ρ )."""
        self.tht_prior_phi = tht_prior_phi
        """Shape φ_θ for the photoemission multiplier: θ ~ Gamma( φ_θ, ψ_θ )."""
        self.tht_prior_psi = tht_prior_psi
        """Scale ψ_θ for the photoemission multiplier: θ ~ Gamma( φ_θ, ψ_θ )."""
        self.kap_D_prior_phi = kap_D_prior_phi
        """φ^D_κ in scaling factor of photoemission rates for the donor dye with FRET: κ^D ~ Gamma( φ^D_κ, ψ^D_κ / φ^D_κ )."""
        self.kap_D_prior_psi = kap_D_prior_psi
        """ψ^D_κ in scaling factor of photoemission rates for the donor dye with FRET: κ^D ~ Gamma( φ^D_κ, ψ^D_κ / φ^D_κ )."""
        self.kap_A_prior_phi = kap_A_prior_phi
        """φ^A_κ in scaling factor of photoemission rates for the acceptor dye with FRET: κ^A ~ Gamma( φ^A_κ, ψ^A_κ / φ^A_κ )."""
        self.kap_A_prior_psi = kap_A_prior_psi
        """ψ^A_κ in scaling factor of photoemission rates for the acceptor dye with FRET: κ^A ~ Gamma( φ^A_κ, ψ^A_κ / φ^A_κ )."""
        self.kap_Z_prior_phi = kap_Z_prior_phi
        """φ^Z_κ in scaling factor of photoemission rates for the donor dye without FRET: κ^Z ~ Gamma( φ^Z_κ, ψ^Z_κ / φ^Z_κ )."""
        self.kap_Z_prior_psi = kap_Z_prior_psi
        """ψ^Z_κ in scaling factor of photoemission rates for the donor dye without FRET: κ^Z ~ Gamma( φ^Z_κ, ψ^Z_κ / φ^Z_κ )."""






class Sample:
    """
    Calculates and stores the initial priors and states.
    """
    n = 0
    """Current iteration in chain."""
    rec: np.ndarray[np.float64] = repmat([0, np.finfo(float).tiny], 3, 1)
    """Current ratios to accept [i][0] / run [i][1] the HMC [0] + MG [1] + Flipper [2] algorithms."""
    ksi_D: np.float64
    """ξ^D Current: Background photoemission rates for the donor dye."""
    ksi_A: np.float64
    """ξ^A Current: Background photoemission rates for the acceptor dye."""
    lam_D: np.ndarray[np.float64]
    """λ^D Current: Donor dye photoemission rates with FRET.\n\nSize: `K`"""
    lam_A: np.ndarray[np.float64]
    """λ^A Current: Acceptor dye photoemission rates with FRET.\n\nSize: `K`"""
    lam_Z: np.float64
    """λ^Z Current: Donor dye photoemission rate without FRET."""
    eff: np.ndarray[np.float64]
    """E Current: Apparent FRET efficiencies of the `K` states.\n\nSize: `K`"""
    eff_select: np.ndarray[np.float64]
    """E Selected: Apparent FRET efficiency of the unique states.\n\nSize: `K_set`"""
    K_set: np.ndarray[np.int64]
    """set(K): The set of unique conformational states from `unique(sample.st)`. Filters `kap_A` and `kap_D` for updating."""
    K_cnt: np.ndarray[np.intp]
    """count(K): The counts of the unique conformational states from `unique(sample.st, return_counts=True)`."""
    K_loc: np.ndarray[np.intp]
    """location(K): The locations of the unique conformational states from `unique(sample.st, return_inverse=True)`. Filters `kap_A` and `kap_D` for `mu`."""


    
    def __init__(self, params: Params, U: Universal, v: np.int64):
        """
        Calculates and stores the initial priors and states.

        Args:
            v: Index of the current simulation.
        """

        # states
        self.st   = np.ones( params.T, dtype = np.int64 )
        """s_n Current: Conformational states of the observed molecule throughout the experiment.\n\nSize: `T`"""
        self.K_set, self.K_loc, self.K_cnt = np.unique( self.st, return_inverse = True, return_counts = True )
        self.K_sz = self.K_set.size
        """size(K) Current: Number of unique conformational states."""
        self.ft_D = np.ones( params.T, dtype = np.bool_ )
        """f_0^D Current: The photo-states (bright/dark) of the donor dye throughout the experiment.\n\nSize: `T`"""
        self.ft_A = np.ones( params.T, dtype = np.bool_ )
        """f_0^A Current: The photo-states (bright/dark) of the acceptor dye throughout the experiment.\n\nSize: `T`"""
        
        # transition probabilities
        self.bm    = params.gamma * np.ones(params.K_lim) / params.K_lim
        """Beta dist. Current: The nonparametric hierarchical prior.\n\nSize: `K`"""
        self.pm    = U.RNGs[v].dirrnd( params.alpha * np.tile( self.bm[0], (params.K_lim, params.K_lim) ) )
        """π Current: The transition matrix between K states.\n\nSize: `K*^2`"""
        self.ps    = U.RNGs[v].dirrnd( params.alpha * self.bm )
        """π* Current: The probability array that state K will occur in the initial state.\n\nSize: `K`"""

        # photophysics
        self.wi_D  = U.RNGs[v].betarnd( params.wi_D_prior_eta, params.wi_D_prior_zeta )
        """ω^D Current: The transition probabilities of the photo-states for the donor dye where ω^D = {ω^D_0, ω^D_1, ω^D_*}."""
        self.wi_A  = U.RNGs[v].betarnd( params.wi_A_prior_eta, params.wi_A_prior_zeta )
        """ω^A Current: The transition probabilities of the photo-states for the acceptor dye where ω^A = {ω^A_0, ω^A_1, ω^A_*}."""
        
        # emission rates + multipliers
        self.tht   = U.RNGs[v].gamrnd1( params.tht_prior_phi,   params.tht_prior_psi   / params.tht_prior_phi   )
        """θ Current: The photoemission multiplier."""
        self.rho_D = U.RNGs[v].gamrnd1( params.rho_D_prior_phi, params.rho_D_prior_psi / params.rho_D_prior_phi )
        """ρ^D Current: Scaling factor for background photoemission rates of the donor dye."""
        self.rho_A = U.RNGs[v].gamrnd1( params.rho_A_prior_phi, params.rho_A_prior_psi / params.rho_A_prior_phi )
        """ρ^A Current: Scaling factor for background photoemission rates of the acceptor dye."""
        self.kap_D = U.RNGs[v].gamrnd(  params.kap_D_prior_phi, params.kap_D_prior_psi / params.kap_D_prior_phi, params.K_lim )
        """κ^D Current: Scaling factors for photoemission rates of the donor dye with FRET.\n\nSize: `K`"""
        self.kap_A = U.RNGs[v].gamrnd(  params.kap_A_prior_phi, params.kap_A_prior_psi / params.kap_A_prior_phi, params.K_lim )
        """κ^A Current: Scaling factors for photoemission rates of the acceptor dye with FRET.\n\nSize: `K`"""
        self.kap_Z = U.RNGs[v].gamrnd1( params.kap_Z_prior_phi, params.kap_Z_prior_psi / params.kap_Z_prior_phi )
        """κ^Z Current: Scaling factor for photoemission rates of the donor dye without FRET."""

        self.UPDATE_AUX()
    


    def UPDATE_AUX(self):
        """
        Update Auxiliary Variables
        ===
        Calculate photoemission rates from priors.

        Updates `ksi_D`, `ksi_A`, `lam_D`, `lam_A`, `lam_Z`, `eff`.
        """
        self.ksi_D = self.tht * self.rho_D
        self.ksi_A = self.tht * self.rho_A
        self.lam_D = self.tht * self.kap_D
        self.lam_A = self.tht * self.kap_A
        self.lam_Z = self.tht * self.kap_Z
        self.eff   = self.kap_A / (self.kap_A + self.kap_D)
        self.eff_select = self.eff[self.K_set]






class Chain_History:
    """
    Keeps a numerical record of this simulation.
    """
    runtime: str
    """String of simulation runtime."""

    def __init__(self, params: Params, sample: Sample, seed_int: np.int64):
        self.rec_ord = np.zeros((params.N+1, 4), dtype = np.float64)
        """Copy of cumulative ratios to accept/run HMC (1) + MG (2) + Flipper (3) algorithms by time $ (0)."""
        self.seed_int = seed_int
        """The specific RNG seed (int) for the Generator class."""
        self.K_sz  = np.concatenate( ([sample.K_sz],  np.zeros((params.N, ),                dtype = np.int64)),   dtype = np.int64,   axis = 0 )
        """K Array: The number of conformational states throughout the simulation at iteration `n`.\n\nSize: `N+1`"""

        # states
        self.st    = np.concatenate( ([sample.st],    np.zeros((params.N, params.T),        dtype = np.int64)),   dtype = np.int64,   axis = 0 )
        """s_t Array: Conformational states of the observed molecule throughout the experiment at iteration `n`.\n\nShape: `N+1, T`"""
        K_set, K_loc, K_cnt = np.unique( self.st, return_inverse=True, return_counts=True )
        self.K_set = [K_set]
        """
        set(K) Array: The set of unique conformational states from `unique(st)` at iteration `n`. Filters `kap_A` and `kap_D` for updating.

        Shape: List of `N` samples, each of size `self.K_sz[n]`.
        """
        self.K_loc = [K_loc]
        """
        location(K) Array: The locations of the unique conformational states from `unique(st, return_inverse=True)` at iteration `n`. Filters `kap_A` and `kap_D` for `mu`.

        Shape: List of `N` samples, each of size `self.K_sz[n]`.
        """
        self.K_cnt = [K_cnt]
        """
        count(K) Array: The counts of the unique conformational states from `unique(st, return_counts=True)` at iteration `n`.

        Shape: List of `N` samples, each of size `self.K_sz[n]`.
        """
        self.ft_D  = np.concatenate( ([sample.ft_D],  np.zeros((params.N, params.T),        dtype = np.bool_)),   dtype = np.bool_,   axis = 0 )
        """f_t^D Array: The photo-states (bright/dark) of the donor dye throughout the experiment at iteration `n`.\n\nShape: `N+1, T`"""
        self.ft_A  = np.concatenate( ([sample.ft_A],  np.zeros((params.N, params.T),        dtype = np.bool_)),   dtype = np.bool_,   axis = 0 )
        """f_t^A Array: The photo-states (bright/dark) of the acceptor dye throughout the experiment at iteration `n`.\n\nShape: `N+1, T`"""
        self.eff   = np.concatenate( ([sample.eff],   np.empty((params.N, params.K_lim),    dtype = np.float64)), dtype = np.float64, axis = 0 )
        """E Array: The apparent FRET efficiency at iteration `n`.\n\nShape: `N+1, K_lim`"""
        self.eff_select = [self.eff[K_set]]
        """E Selected: Apparent FRET efficiency of the unique states at iteration `n`.\n\nShape: `N+1, K_set`"""
        
        # photophysics
        self.wi_D  = np.concatenate( ([sample.wi_D],  np.empty((params.N, 3),               dtype = np.float64)), dtype = np.float64, axis = 0 )
        """ω^D = {ω^D_0, ω^D_1, ω^D_*} Array: The transition probabilities of the photo-states for the donor dye at iteration `n`.\n\nShape: `N+1, 3`"""
        self.wi_A  = np.concatenate( ([sample.wi_A],  np.empty((params.N, 3),               dtype = np.float64)), dtype = np.float64, axis = 0 )
        """ω^A = {ω^A_0, ω^A_1, ω^A_*} Array: The transition probabilities of the photo-states for the acceptor dye at iteration `n`.\n\nShape: `N+1, 3`"""
        
        # transition probabilities
        self.bm    = np.concatenate( ([sample.bm],    np.empty((params.N, params.K_lim),    dtype = np.float64)), dtype = np.float64, axis = 0 )
        """Beta distribution Array: The nonparametric hierarchical prior at iteration `n`.\n\nShape: `N+1, K_lim`"""
        self.pm    = np.concatenate( ([sample.pm.flatten('F')], np.empty((params.N, params.K_lim**2), dtype = np.float64)), dtype = np.float64, axis = 0 )
        """π Array: The transition matrix between K states at iteration `n`.\n\nSize: `K^2`"""
        self.ps    = np.concatenate( ([sample.ps],    np.empty((params.N, params.K_lim),    dtype = np.float64)), dtype = np.float64, axis = 0 )
        """~π Array: The probability array that state K will occur in the initial state at iteration `n`.\n\nShape: `N+1, K_lim`"""
        
        # emission rates + multipliers
        self.tht   = np.concatenate( ([sample.tht],   np.empty((params.N, ),                dtype = np.float64)), dtype = np.float64, axis = 0 )
        """θ Array: The photoemission multiplier at iteration `n`.\n\nSize: `N+1`"""
        self.rho_D = np.concatenate( ([sample.rho_D], np.empty((params.N, ),                dtype = np.float64)), dtype = np.float64, axis = 0 )
        """ρ^D Array: The scaling factor of background emission rates for the donor dye at iteration `n`.\n\nSize: `N+1`"""
        self.rho_A = np.concatenate( ([sample.rho_A], np.empty((params.N, ),                dtype = np.float64)), dtype = np.float64, axis = 0 )
        """ρ^A Array: The scaling factor of background emission rates for the acceptor dye at iteration `n`.\n\nSize: `N+1`"""
        self.kap_D = np.concatenate( ([sample.kap_D], np.empty((params.N, params.K_lim),    dtype = np.float64)), dtype = np.float64, axis = 0 )
        """κ^D Array: The scaling factor of photoemission rates for the donor dye with FRET at iteration `n`.\n\nShape: `N+1, K_lim`"""
        self.kap_A = np.concatenate( ([sample.kap_A], np.empty((params.N, params.K_lim),    dtype = np.float64)), dtype = np.float64, axis = 0 )
        """κ^A Array: The scaling factor of photoemission rates for the acceptor dye with FRET at iteration `n`.\n\nShape: `N+1, K_lim`"""
        self.kap_Z = np.concatenate( ([sample.kap_Z], np.empty((params.N, ),                dtype = np.float64)), dtype = np.float64, axis = 0 )
        """κ^Z Array: The scaling factor of photoemission rates for the donor dye without FRET at iteration `n`.\n\nSize: `N+1`"""
        
        # aux
        self.ksi_D = np.concatenate( ([sample.ksi_D], np.empty((params.N, ),                dtype = np.float64)), dtype = np.float64, axis = 0 )
        """ξ^D Array: The background photoemission rates for the donor dye at iteration `n`.\n\nSize: `N+1`"""
        self.ksi_A = np.concatenate( ([sample.ksi_A], np.empty((params.N, ),                dtype = np.float64)), dtype = np.float64, axis = 0 )
        """ξ^A Array: The background photoemission rates for the acceptor dye at iteration `n`.\n\nSize: `N+1`"""
        self.lam_D = np.concatenate( ([sample.lam_D], np.empty((params.N, params.K_lim),    dtype = np.float64)), dtype = np.float64, axis = 0 )
        """λ^D Array: The donor dye photoemission rates with FRET at iteration `n`.\n\nShape: `N+1, K_lim`"""
        self.lam_A = np.concatenate( ([sample.lam_A], np.empty((params.N, params.K_lim),    dtype = np.float64)), dtype = np.float64, axis = 0 )
        """λ^A Array: The acceptor dye photoemission rates with FRET at iteration `n`.\n\nShape: `N+1, K_lim`"""
        self.lam_Z = np.concatenate( ([sample.lam_Z], np.empty((params.N, ),                dtype = np.float64)), dtype = np.float64, axis = 0 )
        """λ^Z Array: The donor dye photoemission rate without FRET at iteration `n`.\n\nSize: `N+1`"""


