"""
BayesFRET: File 1
===
An Experiment-Adjusted HDP-HMM to Analyze Surface-Immobilized smFRET Data

About
---
begin.py: Generates GUI and binds commands.

`k in [1, K]`: Conformational states.
`t in [1, T]`: Data point of photon intensities.
`n in [0, N]`: MCMC steps.

Author
---
Code by Jared Hidalgo. 

Inspired by MATLAB code from Ioannis Sgouralis, Shreya Madaan, Franky Djutanta, Rachael Kha, Rizal F. Hariadi, and Steve Pressé for "A Bayesian Nonparametric Approach to Single Molecule Förster Resonance Energy Transfer".
"""
title = "BayesFRET: An Experiment-Adjusted HDP-HMM to Analyze Surface-Immobilized smFRET Data"
print("\nChecking requirements.............................\n")

# Import internal packages.
import os, pickle, sys, traceback
from datetime import datetime
from functools import partial
from importlib.metadata import distributions
from platform import system
from shutil import rmtree
from subprocess import call, Popen
from time import sleep

# Get directories.
spl = __file__.split( os.sep )
DIR_PROGRAM = f"{os.sep}".join( spl[:-1] )
DIR_RESOURCES = os.path.join( DIR_PROGRAM, "resources" )
DIR_RESUTLS = f"{os.sep}".join( spl[:-2] )

# Set command to open a directory.
od_dict = {
    "Windows": "explorer", 
    "Darwin": "open"
}
open_directory = od_dict[system()] if system() in od_dict else "xdg-open"

# File dependency check.
req_files = ["code01_classes.py", "code02_setup.py", "code03_mcmc.py", "code04_chart.py", "code05_update.py", "code06_results.py", "requirements.txt"]
lis_files = [f for f in os.listdir(DIR_PROGRAM) if f in req_files]
if len(lis_files) != len(req_files):
    sys.exit( "\n\nERROR: Missing Python files. --> Can't run program." )

# Package check: Online attempt.
try:
    call( f"{sys.executable} -m pip install -U pip", shell=True )
    call( f"{sys.executable} -m pip install -r \"{DIR_PROGRAM}{os.sep}requirements.txt\"", shell=True )
except:
    print( "WARNING: Either no internet or not using independent Python compiler." )

# Package check: Offline check.
with open(f"{DIR_PROGRAM}{os.sep}requirements.txt", "r") as r:
    req_pkgs = [i.split("==" if "==" in i else "\n")[0] for i in r.readlines()]
curr_pkgs    = [dist.metadata['Name'] for dist in distributions()]
missing_pkgs = [x for x in req_pkgs if x not in curr_pkgs]
if len(missing_pkgs) > 0:
    sys.exit( "\n\nERROR: Can't install missing Python packages. --> Can't run program." )

print("\n\nDone!\nBooting BayesFRET................................. ", end="")


# Import external packages.
from PySide6.QtCore import (Qt)
from PySide6.QtGui import (QAction, QIcon, QImage, QPixmap, QScreen)
from PySide6.QtWidgets import (QApplication, QFileDialog, QMainWindow, QMenu, QMessageBox, QSplashScreen, QSystemTrayIcon)

# Get dark/light theme.
import darkdetect
theme_switch = {
    "Dark": "BayesFRET light logo",
    "Light": "BayesFRET dark logo",
}
"""Type of logo to return based on theme switch."""

# Temporarily import program directory to PATH for importing program files from any directory.
sys.path.append( DIR_PROGRAM )

# Start splash screen.
app = QApplication( sys.argv )
main_icon = QImage(os.path.join( DIR_RESOURCES, f"{theme_switch[darkdetect.theme()]}.ico" ))
main_pixmap = QPixmap( main_icon )
splash = QSplashScreen( main_pixmap )
splash.show()
splash.setWindowIcon( main_pixmap )

# Continue importing external packages.
import chime, matplotlib
chime.theme('material')
matplotlib.use('QtAgg')
import numpy as np

# Import program files.
from code01_classes import RNG, Universal, Params, Sample, True_Samples
import code02_setup as setup
from code03_mcmc import Chain_Main
import code06_results as res
from gui01_mainwindow import Ui_MainWindow






class BayesFRET(Ui_MainWindow):
    """
    GUI wrapper with interactive functions and variables.
    """
    title = title
    U: Universal = None
    DIR_PROGRAM = DIR_PROGRAM
    """The directory from the latest input file."""
    DIR_SOURCE = DIR_PROGRAM
    """After importing a file, this var saves its directory."""
    DIR_ACTIVE = ""
    """The full directory to put results while running. ex. `C:\\...\\parentDir\\BayesFRET_{curr_date}`"""
    DIR_ACTIVE_TITLE = ""
    """The directory name to put results while running. ex. `BayesFRET_{curr_date}`"""
    It_D_exp = []
    """Dataset of photon intensities of the donor dye from experimental data."""
    It_A_exp = []
    """Dataset of photon intensities of the acceptor dye from experimental data."""
    syn_params: Params
    """Parameters reused for synthetic data."""
    syn_true: True_Samples
    """True samples reused for synthetic data."""
    tray_state = ".ico"
    """Second half of QIcon filename (GUI)."""
    THEME = theme_switch[darkdetect.theme()]
    
    b_show_graphs = True
    """If the program will show figures (`True`) or hide figures (`False`).

    `b_large` | `b_show_graphs` | Program Action
    :---------|:----------------|:----------------------------------------------------
    `False`   | `True`          | Show + save figures at 1280 x 720.
    `False`   | `False`         | Save figures at 1280 x 720.
    `True`    | `True`          | Show + save figures at maximized screen resolution.
    `True`    | `False`         | Save figures at 3200 x 1800."""
    b_large = False
    """If the figures will be large (`True`) or not (`False`).
    
    `b_large` | `b_show_graphs` | Program Action
    :---------|:----------------|:----------------------------------------------------
    `False`   | `True`          | Show + save figures at 1280 x 720.
    `False`   | `False`         | Save figures at 1280 x 720.
    `True`    | `True`          | Show + save figures at maximized screen resolution.
    `True`    | `False`         | Save figures at 3200 x 1800."""
    b_is_syn = False
    """If the program will USE synthetic (`True`) or experimental (`False`) data."""
    b_reuse_syn = False
    """If the program will REUSE synthetic data (`True`) or not (`False`)."""
    b_already_ran = False
    """If running simulations again with windows open, close windows."""
    b_valid_options = True
    """If options have valid inputs (`True`) or not (`False`)"""
    b_reset_was_enabled = False
    """Keeps track if the Reset button for Settings was enabled (`True`) or not (`False`)"""
    b_exp_donor_ready = False
    """If the donor file for experimental data is compatible (`True`) or not (`False`)"""
    b_exp_accep_ready = False
    """If the acceptor file for experimental data is compatible (`True`) or not (`False`)"""
    b_syn_ready = False
    """If the pickle file for synthetic data is ready (`True`) or not (`False`)"""
    b_is_running = False
    """If the program is running or not."""
    b_has_error = False
    """If the program has an error or not."""
    b_quitting = False
    """If the application is quitting or not. For `func_closeEvent()`."""

    options_grid                = np.ones( 46, dtype=np.int64 )
    """Array to track states of Settings input. `0` = Valid, `1` = Default, `2` = Invalid"""
    rng1_default               = 5
    rng2_default               = 44
    rng3_default               = 356
    rng4_default               = 9918
    all_unique_before          = True
    """If all RNG seeds were unique (`True`) or not (`False`) before processing the current seed change."""
    rng_defaults               = [rng1_default, rng2_default, rng3_default, rng4_default]
    """Internal array of default seeds."""
    rng_max_opts               = [1e2, 1e3, 1e4, 1e5, 1e6]
    """Array of maximum values for RNG seeds."""
    units_t_default            = "s"
    units_I_default            = "photons"
    dt_default                 = 0.1
    dD_default                 = 0.099
    cDD_default                = 0.9
    cAA_default                = 0.75
    qD_default                 = 0.85
    qA_default                 = 0.75
    N_default                  = 1000
    K_lim_default              = 25
    rep_tht_default            = 5
    alpha_default              = 1
    gamma_default              = 1
    rep_bm_default             = 15
    MG_L_default               = 5
    HMC_L_default              = 125
    HMC_eps_default            = 0.0001
    Q_default                  = 9
    burn_in_default            = 0.3
    rho_D_prior_phi_default    = 1
    rho_D_prior_psi_default    = 1
    rho_A_prior_phi_default    = 1
    rho_A_prior_psi_default    = 1
    tht_prior_phi_default      = 2
    tht_prior_psi_default      = -1
    """The actual value is `0.5 * np.sum(It_D + It_A)/(T * dD)` after data import."""
    kap_D_prior_phi_default    = 1
    kap_D_prior_psi_default    = 1
    kap_A_prior_phi_default    = 1
    kap_A_prior_psi_default    = 1
    kap_Z_prior_phi_default    = 1
    kap_Z_prior_psi_default    = 1
    wi_D_prior_eta_default     = np.array([1, 1, 1])
    wi_D_prior_zeta_default    = np.array([1, 1, 1])
    wi_A_prior_eta_default     = np.array([1, 1, 1])
    wi_A_prior_zeta_default    = np.array([1, 1, 1])

    windows_actions: list[QAction] = []
    """List of QActions for every active monitor."""



    def __init__(self, window: QMainWindow):
        super().setupUi(window)
        self.window = window
        self.window.closeEvent = self.func_closeEvent
        
        # Set the logo icon.
        self.logo = QIcon( main_pixmap )
        self.window.setWindowIcon( self.logo )

        # Create tray icon + menu.
        self.tray = QSystemTrayIcon( self.logo )
        self.menu = QMenu()

        # Continue setting up menu + tray.
        self.menu_t3 = self.menu.addAction( "-- Status --" )
        self.menu_t3.setDisabled( True )
        self.menu_run = QAction( self.func_getIcon( "play_music_icon_231499.ico" ), "Run BayesFRET" )
        self.menu_run.triggered.connect( self.func_SETUP_run )
        self.menu.addAction( self.menu_run )
        self.menu_run.setEnabled( False )
        self.menu_quit = QAction( self.func_getIcon( "vcsconflicting_93497.ico" ), "Exit BayesFRET" )
        self.menu_quit.triggered.connect( self.func_closeEvent )
        self.menu.addAction( self.menu_quit )
        self.menu.addSeparator()
        
        menu_t2 = self.menu.addAction( "-- Settings --" )
        menu_t2.setDisabled( True )
        self.menu1_donor = QAction( self.func_getIcon( "foldergreen_93329.ico" ), "Experimental Data: Import Donor File" )
        self.menu1_donor.triggered.connect( self.func_CONTEXT_import_exp_donor )
        self.menu.addAction( self.menu1_donor )
        self.menu2_acceptor = QAction( self.func_getIcon( "folderred_93207.ico" ), "Experimental Data: Import Acceptor File" )
        self.menu2_acceptor.triggered.connect( self.func_CONTEXT_import_exp_acceptor )
        self.menu.addAction( self.menu2_acceptor )
        self.menu3_syn = QAction( self.func_getIcon( "folderblue_92960.ico" ), "Synthetic Data: Import Pickle File" )
        self.menu3_syn.triggered.connect( self.func_CONTEXT_import_syn_data )
        self.menu.addAction( self.menu3_syn )
        self.menu4_show = QAction( self.func_getIcon( "disable_eye_hidden_hide_internet_security_view_icon_127055.ico" ), "Hide Graphs During Run" )
        self.menu4_show.triggered.connect( self.checkBox.click )
        self.menu.addAction( self.menu4_show )
        self.menu5_size = QAction( self.func_getIcon( "arrow_expand_full_fullscreen_internet_screen_security_icon_127065.ico" ), "Maximize Graph Size" )
        self.menu5_size.triggered.connect( self.func_CONTEXT_toggle_size )
        self.menu.addAction( self.menu5_size )
        self.menu6_random = QAction( self.func_getIcon( "multimedia_option_change_exchange_random_arrows_shuffle_icon_258796.ico" ), "Randomize RNG Seeds" )
        self.menu6_random.triggered.connect( self.func_SEEDS_randomize )
        self.menu.addAction( self.menu6_random )
        self.menu7_reset = QAction( self.func_getIcon( "arrow_back_previous_left_return_undo_icon_258802.ico" ), "Reset Settings + RNG Seeds" )
        self.menu7_reset.triggered.connect( self.func_OPTIONS_reset )
        self.menu.addAction( self.menu7_reset )
        self.menu7_reset.setEnabled( False )
        
        self.func_updateMonitors(None, False)
        app.screenRemoved.connect( partial( self.func_updateMonitors, True) )
        app.screenAdded.connect( partial( self.func_updateMonitors, True) )

        # Add menu to tray.
        self.tray.setContextMenu( self.menu )
        self.tray.setVisible(True)
        self.tray.setToolTip( "BayesFRET" )
        self.tray.activated.connect( self.func_raise_window )
        
        # Set button stylesheets.
        self.style_norm_r = str(self.pushButton_2.styleSheet)
        """Stylesheet to decorate the normal button style."""
        self.style_reset_on = ("QPushButton { background-color: red;" +
                                             "color: white;" + 
                                             "font-weight: bold;" +
                                             "border: 1px solid black;" +
                                             "border-style: outset;" +
                                             "border-radius: 6px;}" +
                               "QPushButton:hover { border: 2px solid rgb(0, 120, 212);" +
                                                   "background-color: rgb(245, 0, 100) }" +
                               "QPushButton:pressed { background-color: rgb(224, 0, 0);" +
                                                     "border-style: inset;}")
        """Stylesheet to decorate the Reset button when activated."""
        self.style_accepted = ("QPushButton { background-color: rgb(80, 245, 100);" +
                                             "color: black;" + 
                                             "border: 0px;" +
                                             "border-radius: 6px;}" +
                               "QPushButton:hover { border: 2px solid rgb(0, 120, 212);}" +
                               "QPushButton:pressed { border-style: inset;" +
                                                     "background-color: rgb(0, 224, 0);}")
        """Stylesheet to decorate the Donor and Acceptor buttons when input is accepted."""
        self.style_run = ("QPushButton { background-color: lime;" +
                                        "color: black;" + 
                                        "font-weight: bold;" +
                                        "border: 1px solid black;" +
                                        "border-style: outset;" +
                                        "border-radius: 6px;}" +
                          "QPushButton:hover { border: 2px solid rgb(0, 120, 212);" +
                                              "background-color: rgb(80, 245, 100);}" +
                          "QPushButton:pressed { border-style: inset;" +
                                                "background-color: rgb(0, 224, 0);}")
        """Stylesheet to decorate the Run button when activated."""
        
        # Disable buttons that aren't ready.
        self.pushButton_3.setEnabled(False)
        self.pushButton_2.setEnabled(False)
        self.pushButton_9.setEnabled(False)
        # Set all default values
        self.func_OPTIONS_reset()

        # Bind commands: Main buttons
        self.pushButton_2.clicked.connect( self.func_SETUP_run )
        self.pushButton_3.clicked.connect( self.func_OPTIONS_reset )
        # Bind commands: Data Source
        self.radioButton.released.connect( self.func_DATASOURCE_experimental )
        self.radioButton_2.released.connect( self.func_DATASOURCE_synthetic_new )
        self.radioButton_4.released.connect( self.func_DATASOURCE_synthetic_reuse )
        self.pushButton_5.clicked.connect( self.func_DATASOURCE_import_exp_donor )
        self.pushButton_4.clicked.connect( self.func_DATASOURCE_import_exp_acceptor )
        self.pushButton_9.clicked.connect( self.func_DATASOURCE_import_syn_data )
        # Bind commands: Graph Options
        self.checkBox.stateChanged.connect( self.func_GRAPH_show )
        self.radioButton_3.toggled.connect( self.func_GRAPH_large )
        self.radioButton_5.toggled.connect( self.func_GRAPH_small )
        # Bind commands: RNG Seeds
        self.lineEdit.textChanged.connect( self.check_RNG1 )
        self.lineEdit_2.textChanged.connect( self.check_RNG2 )
        self.lineEdit_3.textChanged.connect( self.check_RNG3 )
        self.lineEdit_4.textChanged.connect( self.check_RNG4 )
        self.pushButton_6.clicked.connect( self.func_SEEDS_randomize )
        # Bind commands and values: smFRET Parameters
        self.lineEdit_9.textChanged.connect(  self.check_units_t )
        self.lineEdit_10.textChanged.connect( self.check_units_I )
        self.lineEdit_11.textChanged.connect( self.check_dt )
        self.lineEdit_12.textChanged.connect( self.check_dD )
        self.lineEdit_14.textChanged.connect( self.check_cDD )
        self.lineEdit_16.textChanged.connect( self.check_cAA )
        self.lineEdit_15.textChanged.connect( self.check_qD )
        self.lineEdit_13.textChanged.connect( self.check_qA )
        # Bind commands and values: Algorithm Settings
        self.lineEdit_25.textChanged.connect( self.check_N ) # MCMC
        self.lineEdit_27.textChanged.connect( self.check_K_lim )
        self.lineEdit_71.textChanged.connect( self.check_rep_tht )
        self.lineEdit_26.textChanged.connect( self.check_alpha ) # Dirichlet
        self.lineEdit_28.textChanged.connect( self.check_gamma )
        self.lineEdit_70.textChanged.connect( self.check_rep_bm )
        self.lineEdit_55.textChanged.connect( self.check_MG_L ) # Methods
        self.lineEdit_53.textChanged.connect( self.check_HMC_L )
        self.lineEdit_54.textChanged.connect( self.check_HMC_eps )
        self.lineEdit_56.textChanged.connect( self.check_Q ) # Results
        self.lineEdit_69.textChanged.connect( self.check_burn_in )
        # Bind commands and values: Photoemission Priors
        self.lineEdit_43.textChanged.connect( self.check_rho_D_prior_phi ) # Background (rho)
        self.lineEdit_45.textChanged.connect( self.check_rho_D_prior_psi )
        self.lineEdit_44.textChanged.connect( self.check_rho_A_prior_phi )
        self.lineEdit_46.textChanged.connect( self.check_rho_A_prior_psi )
        self.lineEdit_30.textChanged.connect( self.check_tht_prior_phi ) # Multiplier (theta)
        self.lineEdit_31.textChanged.connect( self.check_tht_prior_psi )
        self.lineEdit_50.textChanged.connect( self.check_kap_D_prior_phi ) # Dyes (kappa)
        self.lineEdit_48.textChanged.connect( self.check_kap_D_prior_psi )
        self.lineEdit_47.textChanged.connect( self.check_kap_A_prior_phi )
        self.lineEdit_49.textChanged.connect( self.check_kap_A_prior_psi )
        self.lineEdit_51.textChanged.connect( self.check_kap_Z_prior_phi )
        self.lineEdit_52.textChanged.connect( self.check_kap_Z_prior_psi )
        # Bind commands and values: Photophysics Priors
        self.lineEdit_65.textChanged.connect( self.check_d_eta_0 )
        self.lineEdit_58.textChanged.connect( self.check_d_eta_1 )
        self.lineEdit_63.textChanged.connect( self.check_d_eta_z )
        self.lineEdit_67.textChanged.connect( self.check_d_zeta_0 )
        self.lineEdit_57.textChanged.connect( self.check_d_zeta_1 )
        self.lineEdit_62.textChanged.connect( self.check_d_zeta_z )
        self.lineEdit_68.textChanged.connect( self.check_a_eta_0 )
        self.lineEdit_61.textChanged.connect( self.check_a_eta_1 )
        self.lineEdit_60.textChanged.connect( self.check_a_eta_z )
        self.lineEdit_66.textChanged.connect( self.check_a_zeta_0 )
        self.lineEdit_64.textChanged.connect( self.check_a_zeta_1 )
        self.lineEdit_59.textChanged.connect( self.check_a_zeta_z )
        
        print(f"Done!\n\n\n\n                  BayesFRET is loaded!                  \n--------------------------------------------------------\n")
    
    
    
    def func_getIcon(self, filename: str):
        """
        Shortcut to return formatted icon from the resources folder.
        """
        return QIcon( os.path.join( DIR_RESOURCES, filename ) )
    
    
    
    def func_closeEvent(self, event):
        """
        Confirm closing the program while not running. Replaces closeEvent from `QMainWindow`.
        """
        chime.warning()
        self.func_raise_window()
        has_active_dir = self.b_is_running and self.DIR_ACTIVE != "" and os.path.exists( self.DIR_ACTIVE )
        # Display a message box to confirm the close event.
        txt = ( f"Are you sure you want to quit the simulations and close BayesFRET? This will also delete the current results folder:\n{self.DIR_ACTIVE_TITLE}" 
                if has_active_dir else "Are you sure you want to close BayesFRET?" )
        reply = QMessageBox( QMessageBox.Icon.Question,
                                "Quit Simulations and Close App?" if self.b_is_running else "Close App?",
                                txt,
                                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                                self.window,
                                Qt.WindowType.WindowStaysOnTopHint )
        r = reply.exec()

        if r == QMessageBox.StandardButton.Yes:
            self.b_quitting = True
            print("\n\nClosing...")
            # Delete active directory.
            if self.b_is_running:
                if self.DIR_ACTIVE != "" and os.path.exists( self.DIR_ACTIVE ):
                    rmtree( self.DIR_ACTIVE )
                    print(f"Active directory deleted: {self.DIR_ACTIVE}\n")
                    print("---< Run cancelled @ " + datetime.strftime( datetime.now(), "%a, %d %b %Y | %I:%M:%S %p" ) + " >---")
            # Close all figures that exist.
            try: self.U.plt.close("all")
            except: pass
            # NOTE on Mac: "The cached device pixel ratio value was stale on window expose.  Please file a QTBUG which explains how to reproduce."
            # Closing message.
            display_text = "Shutting down BayesFRET." # Length must be even number to align.
            lines = "".join( ["-"]*56 )
            spaces_num = (56 - len(display_text))//2
            spaces = "".join( [" "]*spaces_num )
            print(f"{lines}\n{spaces}{display_text}{spaces}\n\n\n")
            exit()
        else:
            try: event.ignore()
            except: pass



    def func_updateMonitors(self, changed_screen: QScreen, is_not_first: bool):
        """
        Reconstruct context menu to adapt for screens.
        """
        if is_not_first:
            for x in self.windows_actions:
                self.menu.removeAction(x)
        self.windows_actions = []
        for i, s in enumerate(app.screens()):
            self.windows_actions.append( self.menu.addSeparator() )
            self.windows_actions.append( self.menu.addAction( f"-- Move Panel to Monitor {s.name()} --" ) )
            self.windows_actions[-1].setDisabled(True)
            self.windows_actions.append( QAction( self.func_getIcon( "workspaceswitchertopleft_94786.ico" ), "Top-Left Corner" ) )
            self.windows_actions[-1].triggered.connect( partial(self.func_moveWindow_topLeft, i) )
            self.menu.addAction( self.windows_actions[-1] )
            self.windows_actions.append( QAction( self.func_getIcon( "workspaceswitcherrighttop_93848.ico" ), "Top-Right Corner" ) )
            self.windows_actions[-1].triggered.connect( partial(self.func_moveWindow_topRight, i) )
            self.menu.addAction( self.windows_actions[-1] )
            self.windows_actions.append( QAction( self.func_getIcon( "workspaceswitcherleftbottom_93857.ico" ), "Bottom-Left Corner" ) )
            self.windows_actions[-1].triggered.connect( partial(self.func_moveWindow_bottomLeft, i) )
            self.menu.addAction( self.windows_actions[-1] )
            self.windows_actions.append( QAction( self.func_getIcon( "workspaceswitcherrightbottom_93959.ico" ), "Bottom-Right Corner" ) )
            self.windows_actions[-1].triggered.connect( partial(self.func_moveWindow_bottomRight, i) )
            self.menu.addAction( self.windows_actions[-1] )



    def func_moveWindow_topLeft(self, monitor: int):
        """
        Move the control panel to the top-left corner of the chosen monitor. Useful while program is running.
        """
        self.window.move( app.screens()[monitor].availableGeometry().topLeft() )



    def func_moveWindow_topRight(self, monitor: int):
        """
        Move the control panel to the top-right corner of the chosen monitor. Useful while program is running.
        """
        move = app.screens()[monitor].availableGeometry().topRight()
        (ww, wh) = self.window.frameGeometry().size().toTuple()
        self.window.move( move.x()-ww, move.y() )



    def func_moveWindow_bottomLeft(self, monitor: int):
        """
        Move the control panel to the bottom-left corner of the chosen monitor. Useful while program is running.
        """
        move = app.screens()[monitor].availableGeometry().bottomLeft()
        (ww, wh) = self.window.frameGeometry().size().toTuple()
        self.window.move( move.x(), move.y()-wh )



    def func_moveWindow_bottomRight(self, monitor: int):
        """
        Move the control panel to the bottom-right corner of the chosen monitor. Useful while program is running.
        """
        move = app.screens()[monitor].availableGeometry().bottomRight()
        (ww, wh) = self.window.frameGeometry().size().toTuple()
        self.window.move( move.x()-ww, move.y()-wh )


    
    def func_raise_window(self):
        self.window.activateWindow()
        self.window.raise_()

    
    
    def _set_tooltip(self, object, txt: str):
        """
        Set the tooltip of a QObject.
        """
        object.setToolTip( "<html><head/><body><p><span style='font-size:12pt;'>" + txt + "</span></p></body></html>" )

    
    
    def _handle_run_button(self, is_ready: bool):
        """
        Handle and decorate a Run button.
        """
        self.menu_run.setEnabled( is_ready )
        self.pushButton_2.setEnabled( is_ready )
        self.pushButton_2.setStyleSheet( self.style_run if is_ready else self.style_norm_r )
        if is_ready: chime.info()
    

    
    def _invalid_opts_first(self):
        """
        Handle invalid options when it first happens. For `self._check_all_entries()`.
        """
        self.b_valid_options = False
        # Handle "Reset" button.
        self.menu7_reset.setEnabled( True )
        self.pushButton_3.setEnabled( True )
        self.pushButton_3.setStyleSheet( self.style_reset_on )
        self.b_reset_was_enabled = True
        # Handle "Run" button.
        self._handle_run_button( False )
        chime.error()
    
    
    
    def _complete_opts_restored(self, options_choice: str):
        """
        Handle options when valid. For `self._check_all_entries()`.
        """
        self.b_valid_options = True
        # Handle "Reset" button.
        self.b_reset_was_enabled = options_choice == "Valid"
        self.menu7_reset.setEnabled( self.b_reset_was_enabled )
        self.pushButton_3.setEnabled( self.b_reset_was_enabled )
        self.pushButton_3.setStyleSheet( self.style_norm_r )
        # Handle "Run" button.
        self._handle_run_button( 
            self.radioButton_2.isChecked() or (self.radioButton.isChecked() and self.b_exp_donor_ready and self.b_exp_accep_ready) )
    
    

    def func_DATASOURCE_experimental(self):
        """
        Select experimental data. Check if donor + acceptor data is already imported.
        """
        self.b_is_syn = False
        self.b_reuse_syn = False
        # Enable "smFRET" group.
        self.groupBox_4.setEnabled(True)
        # Enable "Experimental" buttons.
        self.pushButton_4.setEnabled(True)
        self.pushButton_5.setEnabled(True)
        # Disable "Synthetic - Reuse" buttons.
        self.pushButton_9.setEnabled(False)
        # Handle "Run" button.
        self._handle_run_button( self.b_exp_donor_ready and self.b_exp_accep_ready and self.b_valid_options )
    
    

    def func_DATASOURCE_synthetic_new(self):
        """
        Select synthetic mode: generates and evaluates synthetic data.
        """
        self.b_is_syn = True
        self.b_reuse_syn = False
        # Disable "smFRET" group.
        self.groupBox_4.setEnabled(False)
        # Disable "Experimental" buttons.
        self.pushButton_4.setEnabled(False)
        self.pushButton_5.setEnabled(False)
        # Disable "Synthetic - Reuse" buttons.
        self.pushButton_9.setEnabled(False)
        # Handle "Run" button.
        self._handle_run_button( self.b_valid_options )
    
    

    def func_DATASOURCE_synthetic_reuse(self):
        """
        Select synthetic mode: reuses existing synthetic data.
        """
        self.b_is_syn = True
        self.b_reuse_syn = True
        # Disable "smFRET" group.
        self.groupBox_4.setEnabled(False)
        # Disable "Experimental" buttons.
        self.pushButton_4.setEnabled(False)
        self.pushButton_5.setEnabled(False)
        # Disable "Synthetic - Reuse" buttons.
        self.pushButton_9.setEnabled(True)
        # Handle "Run" button.
        self._handle_run_button( self.b_valid_options and self.b_syn_ready )
    


    def func_CONTEXT_import_exp_donor(self):
        self.radioButton.click()
        self.func_DATASOURCE_import_exp_donor()
    


    def func_DATASOURCE_import_exp_donor(self):
        """
        Imports donor data, saves input directory.
        """
        donor_file = QFileDialog.getOpenFileName( 
            self.window, "Select the donor file", self.DIR_SOURCE, "Text files (*.txt *.csv)" )
        if donor_file[0] != "":
            self.DIR_SOURCE = os.path.dirname( str(donor_file[0].split(f"{os.sep}")[:-1]) )
            try:
                with open( donor_file[0], "r" ) as d:
                    self.It_D_exp = np.array( d.read().split( None if donor_file[0].split(".")[-1] == "txt" else "," ), dtype = np.float64 )
                b1 = not self.b_exp_accep_ready                                     # True when acceptor file is not imported.
                b2 = len(self.It_D_exp) == len(self.It_A_exp)                       # True when acceptor file is imported and file sizes are same.
                b3 = not np.all( self.It_D_exp == self.It_A_exp ) if b2 else False  # True when acceptor file is imported, file sizes are the same, and not all elements are the same.
                if b1 or (b2 and b3):
                    self._import_exp_donor( donor_file )
                elif self.b_exp_accep_ready:
                    if b2 and not b3:
                        self._import_exp_donor( donor_file )
                        chime.warning()
                        QMessageBox.warning( self.window, "WARNING: Duplicate Files!", 
                                             f"Your donor dataset is an exact duplicate of the acceptor dataset." + 
                                             "\nThis is not an impossibility, but check your files anyway!", 
                                             QMessageBox.StandardButton.Ok )
                    elif not b2:
                        self.It_D_exp = []
                        chime.error()
                        QMessageBox.critical( self.window, "ERROR: File Size Difference!", 
                                              f"Your donor dataset size must be the same as your acceptor dataset size." + 
                                              f"\nDonor size: {str(len(self.It_D_exp))} / Acceptor size: {str(len(self.It_A_exp))}\nThis file will not be imported.", 
                                              QMessageBox.StandardButton.Ok )
            except:
                chime.error()
                QMessageBox.critical( self.window, "ERROR: Unable to Import File!",
                                      "Sorry, but this file isn't compatible. Try a different file.",
                                      QMessageBox.StandardButton.Ok )
    
    

    def _import_exp_donor(self, donor_file: tuple[str, str]):
        """
        Successful import of the donor file.
        """
        self.pushButton_5.setStyleSheet( self.style_accepted )
        self.b_exp_donor_ready = True
        self._handle_run_button( self.b_valid_options and self.b_exp_accep_ready )
        self._set_tooltip( self.pushButton_5, donor_file[0] )
        self.menu1_donor.setText( "*Experimental Data: Replace Donor File" )
    


    def func_CONTEXT_import_exp_acceptor(self):
        self.radioButton.click()
        self.func_DATASOURCE_import_exp_acceptor()
    
    
    
    def func_DATASOURCE_import_exp_acceptor(self):
        """
        Imports acceptor data, saves input directory.
        """
        accep_file = QFileDialog.getOpenFileName( self.window, "Select the acceptor file", self.DIR_SOURCE, "Text files (*.txt *.csv)" )
        if accep_file[0] != "":
            self.DIR_SOURCE = os.path.dirname( str(accep_file[0].split(f"{os.sep}")[:-1]) )
            try:
                with open( accep_file[0], "r" ) as d:
                    self.It_A_exp = np.array( d.read().split( None if accep_file[0].split(".")[-1] == "txt" else "," ), dtype = np.float64 )
                b1 = not self.b_exp_donor_ready                                     # True when acceptor file is not imported.
                b2 = len(self.It_D_exp) == len(self.It_A_exp)                       # True when acceptor file is imported and file sizes are same.
                b3 = not np.all( self.It_D_exp == self.It_A_exp ) if b2 else False  # True when acceptor file is imported, file sizes are the same, and not all elements are the same.
                if b1 or (b2 and b3):
                    self._import_exp_acceptor( accep_file )
                elif self.b_exp_donor_ready:
                    if b2 and not b3:
                        self._import_exp_acceptor( accep_file )
                        chime.warning()
                        QMessageBox.warning( self.window, "WARNING: Duplicate Files!", 
                                             f"Your acceptor dataset is an exact duplicate of the donor dataset." + 
                                             "\nThis is not an impossibility, but check your files anyway!", 
                                             QMessageBox.StandardButton.Ok )
                    elif not b2:
                        self.It_A_exp = []
                        chime.error()
                        QMessageBox.critical( self.window, "ERROR: File Size Difference!", 
                                              f"Your acceptor dataset size must be the same as your donor dataset size." + 
                                              f"\nAcceptor size: {str(len(self.It_A_exp))} / Donor size: {str(len(self.It_D_exp))}", 
                                              QMessageBox.StandardButton.Ok )
            except:
                chime.error()
                QMessageBox.critical( self.window, "ERROR: Unable to Import File!", 
                                      "Sorry, but this file isn't compatible. Try a different file.", 
                                      QMessageBox.StandardButton.Ok )
    
    

    def _import_exp_acceptor(self, accep_file: tuple[str, str]):
        """
        Successful import of the acceptor file.
        """
        self.pushButton_4.setStyleSheet( self.style_accepted )
        self.b_exp_accep_ready = True
        self._handle_run_button( self.b_valid_options and self.b_exp_donor_ready )
        self._set_tooltip( self.pushButton_4, accep_file[0] )
        self.menu2_acceptor.setText( "*Experimental Data: Replace Acceptor File" )
    


    def func_CONTEXT_import_syn_data(self):
        self.radioButton_4.click()
        self.func_DATASOURCE_import_syn_data()
    
    
    
    def func_DATASOURCE_import_syn_data(self):
        """
        Imports all required synthetic data.
        """
        syn_file = QFileDialog.getOpenFileName( self.window, "Select \"BayesFRET_data_params_and_true.p\"", self.DIR_SOURCE, "Pickle files (*.p)" )
        if syn_file[0] != "":
            dir_split = os.path.split(syn_file[0])
            if dir_split[-1] == "BayesFRET_data_params_and_true.p":
                self.DIR_SOURCE = os.path.dirname( f"{os.sep}".join(dir_split[:-1]) )
                try:
                    with open( syn_file[0], "rb" ) as s:
                        x = pickle.load( s )
                        self.syn_params: Params = x[0]
                        self.syn_true: True_Samples = x[1]
                    self.pushButton_9.setStyleSheet( self.style_accepted )
                    self.b_syn_ready = True
                    self._handle_run_button( self.b_valid_options )
                    self._set_tooltip( self.pushButton_9, syn_file[0] )
                    self.menu3_syn.setText( "*Synthetic Data: Replace Pickle File" )
                except:
                    chime.error()
                    QMessageBox.critical( self.window, "ERROR: Unreadable File!", 
                                          "Something went wrong with importing this file. It could be corrupted.", 
                                          QMessageBox.StandardButton.Ok )
                    # Erase any partially saved data.
                    self.syn_params = None
                    self.syn_true = None
            else:
                chime.error()
                QMessageBox.critical( self.window, "ERROR: Unable to Import File!", 
                                      "Sorry, but the only acceptable filename is \"BayesFRET_data_params_and_true.p\". If you changed the filename, please change it back.", 
                                      QMessageBox.StandardButton.Ok )
    

    
    def func_GRAPH_show(self):
        """
        Show figures during the run (or not).
        """
        self.b_show_graphs = self.checkBox.isChecked()
        self.radioButton_3.setText( "Maximized" if self.b_show_graphs else "3200 x 1800" )
        self._set_tooltip( self.radioButton_3, "Based on active monitor." if self.b_show_graphs else "Resolution (px)" )
        strt = "Show" if self.b_show_graphs else "Hide"
        self._set_tooltip( self.checkBox, f"{strt} graphs while running BayesFRET. Regardless, graphs will be saved as PNG files." )
        iconC = "disable_eye_hidden_hide_internet_security_view_icon_127055.ico" if self.b_show_graphs else "business_eye_focus_internet_security_view_vision_icon_127037.ico"
        self.menu4_show.setIcon( self.func_getIcon( iconC ) )
        self.menu4_show.setText( ("Hide" if self.b_show_graphs else "Show") + " Graphs During Run" )
        if self.b_large:
            self.menu5_size.setText( "Decrease Graph Size to 1280 x 720" )
        else:
            self.menu5_size.setText( "Maximize Graph Size" if self.b_show_graphs else "Increase Graph Size to 3200 x 1800" )
    
    
    
    def func_GRAPH_large(self):
        """
        Graphs will be outputted and saved as 3200 x 1800 (`b_show_graphs=False`) or maximized (`b_show_graphs=True`).
        """
        self.b_large = True
        self.menu5_size.setIcon(self.func_getIcon( "arrow_exit_internet_minimize_reduce_screen_security_icon_127081.ico" ))
        self.menu5_size.setText( "Decrease Graph Size to 1280 x 720" )
    
    
    
    def func_GRAPH_small(self):
        """
        Graphs will be outputted and saved as 1280 x 720.
        """
        self.b_large = False
        self.menu5_size.setIcon(self.func_getIcon( "arrow_expand_full_fullscreen_internet_screen_security_icon_127065.ico" ))
        self.menu5_size.setText( "Maximize Graph Size" if self.b_show_graphs else "Increase Graph Size to 3200 x 1800" )
    
    
    
    def func_CONTEXT_toggle_size(self):
        """
        Toggle the size of the figures from the context menu.
        """
        if self.b_large:    self.radioButton_5.click()
        else:               self.radioButton_3.click()
    

    
    def func_SEEDS_randomize(self):
        """
        For each randomized seed, select a random maximum value from the list `self.rng_max_opts`.
        """
        all_diff = False
        while not all_diff: # Loop until all values are unique.
            ints = [ RNG().randi1(0, self.rng_max_opts[x]) for x in RNG().randi(0, len(self.rng_max_opts), 4) ]
            if len(np.unique(ints)) == len(ints):
                all_diff = True
        self.lineEdit.setText( str(ints[0]) )
        self.lineEdit_2.setText( str(ints[1]) )
        self.lineEdit_3.setText( str(ints[2]) )
        self.lineEdit_4.setText( str(ints[3]) )
    
    
    
    def func_OPTIONS_reset(self):
        """
        Reset the Settings panel and RNG seed values.
        """
        # RNG Seeds
        self.lineEdit.setText( str(self.rng1_default) )
        self.lineEdit_2.setText( str(self.rng2_default) )
        self.lineEdit_3.setText( str(self.rng3_default) )
        self.lineEdit_4.setText( str(self.rng4_default) )
        # smFRET Parameters
        self.lineEdit_9.setText( str(self.units_t_default) ) # Top
        self.lineEdit_10.setText( str(self.units_I_default) )
        self.lineEdit_11.setText( str(self.dt_default) )
        self.lineEdit_12.setText( str(self.dD_default) )
        self.lineEdit_14.setText( str(self.cDD_default) ) # Bottom
        self.lineEdit_16.setText( str(self.cAA_default) )
        self.lineEdit_15.setText( str(self.qD_default) )
        self.lineEdit_13.setText( str(self.qA_default) )
        # Algorithm Settings
        self.lineEdit_25.setText( str(self.N_default) ) # MCMC
        self.lineEdit_27.setText( str(self.K_lim_default) )
        self.lineEdit_71.setText( str(self.rep_tht_default) )
        self.lineEdit_26.setText( str(self.alpha_default) ) # Dirichlet
        self.lineEdit_28.setText( str(self.gamma_default) )
        self.lineEdit_70.setText( str(self.rep_bm_default) )
        self.lineEdit_55.setText( str(self.MG_L_default) ) # Methods
        self.lineEdit_53.setText( str(self.HMC_L_default) ) 
        self.lineEdit_54.setText( str(self.HMC_eps_default) )
        self.lineEdit_56.setText( str(self.Q_default) ) # Results
        self.lineEdit_69.setText( str(self.burn_in_default) )
        # Photoemission Priors
        self.lineEdit_43.setText( str(self.rho_D_prior_phi_default) ) # Background (rho)
        self.lineEdit_45.setText( str(self.rho_D_prior_psi_default) )
        self.lineEdit_44.setText( str(self.rho_A_prior_phi_default) )
        self.lineEdit_46.setText( str(self.rho_A_prior_psi_default) )
        self.lineEdit_30.setText( str(self.tht_prior_phi_default) ) # Multiplier (theta)
        self.lineEdit_31.setText( str(self.tht_prior_psi_default) )
        self.lineEdit_50.setText( str(self.kap_D_prior_phi_default) ) # Dyes (kappa)
        self.lineEdit_48.setText( str(self.kap_D_prior_psi_default) )
        self.lineEdit_47.setText( str(self.kap_A_prior_phi_default) )
        self.lineEdit_49.setText( str(self.kap_A_prior_psi_default) )
        self.lineEdit_51.setText( str(self.kap_Z_prior_phi_default) )
        self.lineEdit_52.setText( str(self.kap_Z_prior_psi_default) )
        # Photophysics Priors
        self.lineEdit_65.setText( str(self.wi_D_prior_eta_default[0]) ) # η^D column
        self.lineEdit_58.setText( str(self.wi_D_prior_eta_default[1]) )
        self.lineEdit_63.setText( str(self.wi_D_prior_eta_default[2]) )
        self.lineEdit_67.setText( str(self.wi_D_prior_zeta_default[0]) ) # ζ^D column
        self.lineEdit_57.setText( str(self.wi_D_prior_zeta_default[1]) )
        self.lineEdit_62.setText( str(self.wi_D_prior_zeta_default[2]) )
        self.lineEdit_68.setText( str(self.wi_A_prior_eta_default[0]) ) # η^A column
        self.lineEdit_61.setText( str(self.wi_A_prior_eta_default[1]) )
        self.lineEdit_60.setText( str(self.wi_A_prior_eta_default[2]) )
        self.lineEdit_66.setText( str(self.wi_A_prior_zeta_default[0]) ) # ζ^A column
        self.lineEdit_64.setText( str(self.wi_A_prior_zeta_default[1]) )
        self.lineEdit_59.setText( str(self.wi_A_prior_zeta_default[2]) )
        self.b_valid_options = True
        # Reset "Reset" button
        self.menu7_reset.setEnabled( False )
        self.pushButton_3.setEnabled( False )
        self.pushButton_3.setStyleSheet( self.style_norm_r )
        self.b_reset_was_enabled = False
    
    
    
    def func_SETUP_run(self):
        """
        Runs the simulations, produces graphs, exports data.
        """
        try:
            # Get + apply time stamp.
            curr_reset_state = self.pushButton_3.isEnabled()
            first_tme = datetime.now()
            start_datetime = datetime.strftime( first_tme, "%a, %d %b %Y | %I:%M:%S %p" )
            print( f"\n\n----< Run began at {start_datetime}  >----\n" )
            time_folder = datetime.strftime( first_tme, "%Y %m %d, %H %M %S %f" )
            time_folder = ("syn" if self.b_is_syn else "exp") + " - " + time_folder
            # Create active folder.
            self.func_updateStatus1( "Creating directory" )
            self.DIR_ACTIVE = os.path.join( DIR_RESUTLS, f"BayesFRET_{time_folder}" )
            self.DIR_ACTIVE_TITLE = f"BayesFRET_{time_folder}"
            os.makedirs(self.DIR_ACTIVE)
            self.func_updateStatus1( "Done!" )
            self.b_is_running = True

            # Handle plots.
            if self.b_already_ran:
                self.U.plt.close( "all" )
                self.tray_state = ".ico"
                self.tray.setIcon(self.func_getIcon( f"{self.THEME}{self.tray_state}" ))
            # Close GUI.
            self.centralwidget.setEnabled( False )
            # Close "Run" button.
            self.pushButton_2.setStyleSheet( self.style_reset_on )
            self.pushButton_2.setText( "Running" )
            # Handle context menu.
            self.menu1_donor.setEnabled( False )
            self.menu2_acceptor.setEnabled( False )
            self.menu3_syn.setEnabled( False )
            self.menu6_random.setEnabled( False )
            self.menu4_show.setEnabled( False )
            self.menu5_size.setEnabled( False )
            self.menu7_reset.setEnabled( False )
            self.menu_run.setEnabled( False )
            self.menu_run.setIcon(self.func_getIcon( "clock_time_watch_icon_181567.ico" ))
            txt = "Synthetic" if self.b_is_syn else "Experimental"
            self.menu_quit.setText( "Exit BayesFRET + " + txt + " Simulation" )
            sleep(1)

            # Setup Universal object.
            self.func_updateStatus1( "Setting up variables" )
            seeds = [ np.int64(self.lineEdit.text()),
                      np.int64(self.lineEdit_2.text()),
                      np.int64(self.lineEdit_3.text()),
                      np.int64(self.lineEdit_4.text()) ]
            RNGs = [ RNG(seeds[0]), RNG(seeds[1]), RNG(seeds[2]), RNG(seeds[3]), RNG() ]
            Q                  = np.int64(self.lineEdit_56.text())
            burn_in            = np.float64(self.lineEdit_69.text())

            # Save variables.
            N                  = np.int64(self.lineEdit_25.text())
            K_lim              = np.int64(self.lineEdit_27.text())
            rep_tht            = np.int64(self.lineEdit_71.text())
            alpha              = np.float64(self.lineEdit_26.text())
            gamma              = np.float64(self.lineEdit_28.text())
            rep_bm             = np.int64(self.lineEdit_70.text())
            MG_L               = np.int64(self.lineEdit_55.text())
            HMC_L              = np.int64(self.lineEdit_53.text())
            HMC_eps            = np.float64(self.lineEdit_54.text())
            rho_D_prior_phi    = np.float64(self.lineEdit_43.text())
            rho_D_prior_psi    = np.float64(self.lineEdit_45.text())
            rho_A_prior_phi    = np.float64(self.lineEdit_44.text())
            rho_A_prior_psi    = np.float64(self.lineEdit_46.text())
            tht_prior_phi      = np.float64(self.lineEdit_30.text())
            kap_D_prior_phi    = np.float64(self.lineEdit_50.text())
            kap_D_prior_psi    = np.float64(self.lineEdit_48.text())
            kap_A_prior_phi    = np.float64(self.lineEdit_47.text())
            kap_A_prior_psi    = np.float64(self.lineEdit_49.text())
            kap_Z_prior_phi    = np.float64(self.lineEdit_51.text())
            kap_Z_prior_psi    = np.float64(self.lineEdit_52.text())
            wi_D_prior_eta     = np.array([ self.lineEdit_65.text(), self.lineEdit_58.text(), self.lineEdit_63.text() ], dtype=np.float64)
            wi_D_prior_zeta    = np.array([ self.lineEdit_67.text(), self.lineEdit_57.text(), self.lineEdit_62.text() ], dtype=np.float64)
            wi_A_prior_eta     = np.array([ self.lineEdit_68.text(), self.lineEdit_61.text(), self.lineEdit_60.text() ], dtype=np.float64)
            wi_A_prior_zeta    = np.array([ self.lineEdit_66.text(), self.lineEdit_64.text(), self.lineEdit_59.text() ], dtype=np.float64)
            
            self.U = Universal( title, self.b_show_graphs, self.b_large, self.b_is_syn, RNGs, self.DIR_ACTIVE, DIR_RESOURCES,  
                                self.THEME, Q, burn_in, rep_tht, rep_bm, MG_L, HMC_L, HMC_eps, self.menu_t3, self.menu_run, self.tray, self.tray_state )
            self.func_updateStatus1( "Done!" )
            
            # Junction between experimental and synthetic data.
            if self.b_is_syn and not self.b_reuse_syn:
                self.func_updateStatus1( "Generating, graphing, exporting synthetic data" )
                It_D, It_A, units_t, units_I, dt, dD, cDD, cAA, qD, qA, T, self.U = setup.GENERATE_SYNTHETIC_DATA( self.U )
            elif self.b_is_syn and self.b_reuse_syn:
                self.func_updateStatus1( "Reimplementing synthetic data" )
                # Only copy "smFRET Parameters".
                It_D           = self.syn_params.It_D 
                It_A           = self.syn_params.It_A
                units_t        = self.syn_params.units_t
                units_I        = self.syn_params.units_I
                dt             = self.syn_params.dt
                dD             = self.syn_params.dD
                cDD            = self.syn_params.cDD
                cAA            = self.syn_params.cAA
                qD             = self.syn_params.qD
                qA             = self.syn_params.qA
                T              = np.size(self.syn_params.It_D )
                self.U         = setup.REUSE_SYNTHETIC_DATA( self.syn_params, self.syn_true, self.U )
            else:
                self.func_updateStatus1( "Graphing experimental data" )
                It_D           = self.It_D_exp
                It_A           = self.It_A_exp
                units_t        = self.lineEdit_9.text()
                units_I        = self.lineEdit_10.text()
                dt             = np.float64(self.lineEdit_11.text())
                dD             = np.float64(self.lineEdit_12.text())
                cDD            = np.float64(self.lineEdit_14.text())
                cAA            = np.float64(self.lineEdit_16.text())
                qD             = np.float64(self.lineEdit_15.text())
                qA             = np.float64(self.lineEdit_13.text())
                T              = np.size(self.It_D_exp)
                self.U         = setup.SHOW_EXPERIMENTAL_DATA( units_t, units_I, It_D, It_A, T, dt, self.U )
            
            tps                = np.float64(self.lineEdit_31.text())
            tht_prior_psi      = np.sum(It_D + It_A)/(2*T*dD) if tps == self.tht_prior_psi_default else tps
            self.func_updateStatus1( "Done!" )
            
            # Setup and export parameters.
            self.func_updateStatus1( "Exporting parameters and graphing priors" )
            params = Params( It_D, It_A, T, units_t, units_I, dt, dD, cDD, cAA, qD, qA, N, K_lim, alpha, gamma, 
                            wi_D_prior_eta, wi_D_prior_zeta, wi_A_prior_eta, wi_A_prior_zeta, 
                            rho_D_prior_phi, rho_D_prior_psi, rho_A_prior_phi, rho_A_prior_psi, tht_prior_phi, tht_prior_psi, 
                            kap_D_prior_phi, kap_D_prior_psi, kap_A_prior_phi, kap_A_prior_psi, kap_Z_prior_phi, kap_Z_prior_psi )
            
            # Export data (pickle). Bundled for reuse.
            fname_params = "BayesFRET_data_params" + ("_and_true" if self.b_is_syn else "")
            fp_params = self.U.func_getActivePath( f"{fname_params}.p" )
            with open(fp_params, "wb") as file:
                if self.U.is_syn: pickle.dump( [params, self.U.TS], file )
                else:             pickle.dump( params, file )
            
            # Export data (txt).
            DATA = { "Parameters": {attr: getattr(obj, attr, 'NaN') for obj in [params] for attr in vars(obj)} }
            if self.b_is_syn:
                DATA["True Values"] = {attr: getattr(obj, attr, 'NaN') for obj in [self.U.TS] for attr in vars(obj)}
            with open( self.U.func_getActivePath( f"{fname_params}.txt" ), "w+" ) as file:
                txt1 = "PARAMETERS AND TRUE VALUES" if self.b_is_syn else "PARAMETERS"
                txt2 = "="*len(self.title)
                file.write( f"{self.title}\n{txt1}\n{txt2}\n" )
                for key, val in DATA.items():
                    file.write(f'\n\n{key}\n----------\n')
                    for attr, obj in val.items():
                        file.write( f"{attr}: {obj}\n\n" )
            
            # Graph priors.
            setup.GRAPH_PRIORS( params, self.U )
            self.func_updateStatus1( "Done!\n" )

            # Setup samples.
            samples = [ Sample( params, self.U, i ) for i in self.U.range_seeds ]
            
            # MCMC work and graph results.
            chain = Chain_Main( params, samples, self.U )
            res.GRAPH_ALL( params, self.U )

            # Export Universal file.
            self.func_updateStatus1( "Exporting Universal file" )
            self.U.func_delForPickle()
            with open(self.U.func_getActivePath( "BayesFRET_Universal_class.p" ), "wb") as file:
                pickle.dump( self.U, file )
            self.func_updateStatus1( "Done!" )

            # Open file directory of saved files.
            self.func_updateStatus1( "Opening directory" )
            x = Popen( [open_directory, self.DIR_ACTIVE] ) 
            sleep(1)
            x.kill()
            # NOTE: Catch resource warning.
            
            # Handle GUI quickly.
            if not self.b_has_error:
                self.pushButton_2.setStyleSheet( self.style_run )
                self.pushButton_2.setText( "Run" )
                self.centralwidget.setEnabled( True )
                self.menu1_donor.setEnabled( True )
                self.menu2_acceptor.setEnabled( True )
                self.menu3_syn.setEnabled( True )
                self.menu6_random.setEnabled( True )
                self.menu4_show.setEnabled( True )
                self.menu5_size.setEnabled( True )
                self.menu7_reset.setEnabled( curr_reset_state )
            self.menu_quit.setText( "Exit BayesFRET" )
            self.b_already_ran = True
            self.DIR_ACTIVE = ""
            self.DIR_ACTIVE_TITLE = ""
            self.b_is_running = False

            # Wrap up program.
            chime.success()
            self.func_updateStatus1( "Done!" )
            print("FINISHED! All functions executed & all data exported.")
            self.menu_run.setEnabled( True )
            self.menu_run.setText( "Run BayesFRET" )
            self.menu_run.setIcon(self.func_getIcon( "play_music_icon_231499.ico" ))
            self.tray.setToolTip( "BayesFRET" )
            self.tray_state = " - 6 of 6.ico"
            self.tray.setIcon(self.func_getIcon( f"{self.THEME}{self.tray_state}" ))

            last_tme = datetime.now()
            runtime = str(last_tme - first_tme)
            last_datetime = datetime.strftime( last_tme, "%a, %d %b %Y | %I:%M:%S %p" )
            print(f"Total runtime: {runtime}\n\n----< Run cleared @ {last_datetime} >----\n\n\n> Ready for use.")
            with open( self.U.func_getActivePath("BayesFRET_performance.txt"), "w+" ) as f:
                f.write( f"Start: {start_datetime}\nFinish: {last_datetime}\nOverall Runtime: {runtime}\n\n" + 
                         "".join([f"Seed {seeds[i]} Runtime: {chain.runtimes[i]}\n" for i in range(4)]) )

        # Capture any errors and restart app.
        except BaseException:
            if not self.b_quitting:
                chime.error()
                self.b_has_error = True
                ex_type, ex_value, ex_traceback = sys.exc_info()
                err_title = f"{ex_type.__name__}: {ex_value}"
                stack_trace = [ "File: %s\nLine: %d\nFunction Name: %s\nMessage: %s" % (t[0], t[1], t[2], t[3]) 
                                for t in traceback.extract_tb(ex_traceback) ]
                dashes = "-"*25
                print(f"\n\n\n{dashes}ERROR!{dashes}")
                traceback.print_exc()
                print(f"{dashes}ERROR!{dashes}")
                try:
                    self.U.plt.close( "all" )
                    self.U = None
                except: pass
                self.tray.setToolTip( err_title )
                self.tray_state = " - error.ico"
                self.tray.setIcon(self.func_getIcon( f"{self.THEME}{self.tray_state}" ))
                self.menu_t3.setText( "-- Status --" )
                self.menu_run.setText( err_title )
                self.menu_run.setIcon(self.func_getIcon( "vcsconflicting_93497.ico" ))
                self.menu_quit.setText( "Exit BayesFRET" )
                self.b_already_ran = False
                self.b_is_running = False
                print("\n\n---< Run cancelled @ " + datetime.strftime( datetime.now(), "%a, %d %b %Y | %I:%M:%S %p" )
                + " >---")
                add = f" and delete the active folder {self.DIR_ACTIVE}." if self.DIR_ACTIVE != "" else "."
                reply = QMessageBox( QMessageBox.Icon.Critical,
                                    err_title,
                                    f"{stack_trace[-1]}\n\nBayesFRET is reset.\nCheck the command line for more details.\nContinue to open the app{add}",
                                    QMessageBox.StandardButton.Ok,
                                    self.window,
                                    Qt.WindowType.WindowStaysOnTopHint )
                r = reply.exec()
                if r in [QMessageBox.StandardButton.Ok, QMessageBox.StandardButton.Close]:
                    if self.DIR_ACTIVE != "" and os.path.exists( self.DIR_ACTIVE ):
                        rmtree( self.DIR_ACTIVE )
                        self.DIR_ACTIVE = ""
                    self.pushButton_2.setStyleSheet( self.style_run )
                    self.pushButton_2.setText( "Run" )
                    self.tray.setToolTip( "BayesFRET" )
                    self.tray.setIcon( main_pixmap )
                    self.centralwidget.setEnabled( True )
                    self.menu1_donor.setEnabled( True )
                    self.menu2_acceptor.setEnabled( True )
                    self.menu3_syn.setEnabled( True )
                    self.menu4_show.setEnabled( True )
                    self.menu5_size.setEnabled( True )
                    self.menu6_random.setEnabled( True )
                    self.menu7_reset.setEnabled( curr_reset_state )
                    self.menu_run.setEnabled( True )
                    self.menu_run.setText( "Run BayesFRET" )
                    self.menu_run.setIcon(self.func_getIcon( "play_music_icon_231499.ico" ))
                    self.b_has_error = False
                    print("\n\n> Ready for use.")
    


    def func_updateStatus1(self, status: str):
        """
        Updates the one-time status in the CMD, context menu, and tray tooltip.
        """
        if "Done!" not in status:
            dots = "."*(50 - len(status))
            print(f"{status}{dots} ", end="")
            self.menu_run.setText( status )
            self.tray.setToolTip( status )
        else:
            print(f"{status}")
            self.menu_run.setText( self.menu_run.text() + ": Done!" )
            self.tray.setToolTip( self.tray.toolTip() + ": Done!" )
    

    #
    # Internal: Auto-check any and every option.
    #
    def check_RNG1(self, s: str):            self._validate_rng_seed( s, self.rng1_default, 0 )
    def check_RNG2(self, s: str):            self._validate_rng_seed( s, self.rng2_default, 1 )
    def check_RNG3(self, s: str):            self._validate_rng_seed( s, self.rng3_default, 2 )
    def check_RNG4(self, s: str):            self._validate_rng_seed( s, self.rng4_default, 3 )
    def check_units_t(self, s: str):         self.options_grid[4]  = self._validate_str( s, self.units_t_default ); self._check_all_entries()
    def check_units_I(self, s: str):         self.options_grid[5]  = self._validate_str( s, self.units_I_default ); self._check_all_entries()
    def check_dt(self, s: str):              self.options_grid[6]  = self._validate_flt_0i( s, self.dt_default, ">" ); self._check_all_entries()
    def check_dD(self, s: str):              self.options_grid[7]  = self._validate_flt_0i( s, self.dD_default, ">" ); self._check_all_entries()
    def check_cDD(self, s: str):             self.options_grid[8]  = self._validate_flt_01( s, self.cDD_default, "[]" ); self._check_all_entries()
    def check_cAA(self, s: str):             self.options_grid[9]  = self._validate_flt_01( s, self.cAA_default, "[]" ); self._check_all_entries()
    def check_qD(self, s: str):              self.options_grid[10] = self._validate_flt_01( s, self.qD_default, "[]" ); self._check_all_entries()
    def check_qA(self, s: str):              self.options_grid[11] = self._validate_flt_01( s, self.qA_default, "[]" ); self._check_all_entries()
    def check_N(self, s: str):               self.options_grid[12] = self._validate_int_ge( s, self.N_default, 10 ); self._check_all_entries()
    def check_K_lim(self, s: str):           self.options_grid[13] = self._validate_int_ge( s, self.K_lim_default, 10 ); self._check_all_entries()
    def check_rep_tht(self, s: str):         self.options_grid[14] = self._validate_int_ge( s, self.rep_tht_default, 1 ); self._check_all_entries()
    def check_alpha(self, s: str):           self.options_grid[15] = self._validate_flt_0i( s, self.alpha_default, ">" ); self._check_all_entries()
    def check_gamma(self, s: str):           self.options_grid[16] = self._validate_flt_0i( s, self.gamma_default, ">" ); self._check_all_entries()
    def check_rep_bm(self, s: str):          self.options_grid[17] = self._validate_int_ge( s, self.rep_bm_default, 1 ); self._check_all_entries()
    def check_MG_L(self, s: str):            self._validate_algorithm_lengths(); # 18
    def check_HMC_L(self, s: str):           self._validate_algorithm_lengths(); # 19
    def check_HMC_eps(self, s: str):         self.options_grid[20] = self._validate_flt_0i( s, self.HMC_eps_default, ">" ); self._check_all_entries()
    def check_Q(self, s: str):               self.options_grid[21] = self._validate_int_ra( s, self.Q_default, (2, 15) ); self._check_all_entries()
    def check_burn_in(self, s: str):         self.options_grid[22] = self._validate_flt_01( s, self.burn_in_default, "[)" ); self._check_all_entries()
    def check_rho_D_prior_phi(self, s: str): self.options_grid[23] = self._validate_flt_0i( s, self.rho_D_prior_phi_default, ">" ); self._check_all_entries()
    def check_rho_D_prior_psi(self, s: str): self.options_grid[24] = self._validate_flt_0i( s, self.rho_D_prior_psi_default, ">" ); self._check_all_entries()
    def check_rho_A_prior_phi(self, s: str): self.options_grid[25] = self._validate_flt_0i( s, self.rho_A_prior_phi_default, ">" ); self._check_all_entries()
    def check_rho_A_prior_psi(self, s: str): self.options_grid[26] = self._validate_flt_0i( s, self.rho_A_prior_psi_default, ">" ); self._check_all_entries()
    def check_tht_prior_phi(self, s: str):   self.options_grid[27] = self._validate_flt_0i( s, self.tht_prior_phi_default, ">" ); self._check_all_entries()
    def check_tht_prior_psi(self, s: str):   self.options_grid[28] = self._validate_tht_prior_psi( s ); self._check_all_entries()
    def check_kap_D_prior_phi(self, s: str): self.options_grid[29] = self._validate_flt_0i( s, self.kap_D_prior_phi_default, ">" ); self._check_all_entries()
    def check_kap_D_prior_psi(self, s: str): self.options_grid[30] = self._validate_flt_0i( s, self.kap_D_prior_psi_default, ">" ); self._check_all_entries()
    def check_kap_A_prior_phi(self, s: str): self.options_grid[31] = self._validate_flt_0i( s, self.kap_A_prior_phi_default, ">" ); self._check_all_entries()
    def check_kap_A_prior_psi(self, s: str): self.options_grid[32] = self._validate_flt_0i( s, self.kap_A_prior_psi_default, ">" ); self._check_all_entries()
    def check_kap_Z_prior_phi(self, s: str): self.options_grid[33] = self._validate_flt_0i( s, self.kap_Z_prior_phi_default, ">" ); self._check_all_entries()
    def check_kap_Z_prior_psi(self, s: str): self.options_grid[34] = self._validate_flt_0i( s, self.kap_Z_prior_psi_default, ">" ); self._check_all_entries()
    def check_d_eta_0(self, s: str):         self.options_grid[35] = self._validate_flt_0i( s, self.wi_D_prior_eta_default[0], ">" ); self._check_all_entries()
    def check_d_eta_1(self, s: str):         self.options_grid[36] = self._validate_flt_0i( s, self.wi_D_prior_eta_default[1], ">" ); self._check_all_entries()
    def check_d_eta_z(self, s: str):         self.options_grid[37] = self._validate_flt_0i( s, self.wi_D_prior_eta_default[2], ">" ); self._check_all_entries()
    def check_d_zeta_0(self, s: str):        self.options_grid[38] = self._validate_flt_0i( s, self.wi_D_prior_zeta_default[0], ">" ); self._check_all_entries()
    def check_d_zeta_1(self, s: str):        self.options_grid[39] = self._validate_flt_0i( s, self.wi_D_prior_zeta_default[1], ">" ); self._check_all_entries()
    def check_d_zeta_z(self, s: str):        self.options_grid[40] = self._validate_flt_0i( s, self.wi_D_prior_zeta_default[2], ">" ); self._check_all_entries()
    def check_a_eta_0(self, s: str):         self.options_grid[41] = self._validate_flt_0i( s, self.wi_A_prior_eta_default[0], ">" ); self._check_all_entries()
    def check_a_eta_1(self, s: str):         self.options_grid[42] = self._validate_flt_0i( s, self.wi_A_prior_eta_default[1], ">" ); self._check_all_entries()
    def check_a_eta_z(self, s: str):         self.options_grid[43] = self._validate_flt_0i( s, self.wi_A_prior_eta_default[2], ">" ); self._check_all_entries()
    def check_a_zeta_0(self, s: str):        self.options_grid[44] = self._validate_flt_0i( s, self.wi_A_prior_zeta_default[0], ">" ); self._check_all_entries()
    def check_a_zeta_1(self, s: str):        self.options_grid[45] = self._validate_flt_0i( s, self.wi_A_prior_zeta_default[1], ">" ); self._check_all_entries()
    def check_a_zeta_z(self, s: str):        self.options_grid[46] = self._validate_flt_0i( s, self.wi_A_prior_zeta_default[2], ">" ); self._check_all_entries()
    
    
    
    def _validate_rng_seed(self, s: str, default: np.int64, indx: np.int64):
        """
        Checks input change of an RNG seed in Settings. Internally uses `_validate_int_ge()` and `_check_all_entries()`.
        
        All RNG seeds must be unique + non-negative.

        Args:
            s: The string from lineEdit to analyze.
            default: The default value to check against.
            indx: Index of the status of the current RNG seed in `self.options_grid`.

        Returns:
            :self.options_grid:
                0: Valid: `np.int64(s)` must be non-negative and not the same as the other seeds.
                1: Default value.
                2: Invalid: `np.int64(s)` is not an integer or is equal to another seed.
        """
        curr_seeds = [lE.text() for lE in [self.lineEdit, self.lineEdit_2, self.lineEdit_3, self.lineEdit_4]]
        all_unique_now = len(np.unique(curr_seeds)) == len(curr_seeds)
        if all_unique_now and self.all_unique_before:           # If the RNG seeds were and still are normal...
            self.options_grid[indx] = self._validate_int_ge(s, default, 0) # Check only the current RNG seed.
            self._check_all_entries()
        elif all_unique_now and not self.all_unique_before:     # If the RNG seeds were previously not unique, but it’s all good now...
            for i in range(4): # Check all RNG seeds.
                self.options_grid[i] = self._validate_int_ge(curr_seeds[i], self.rng_defaults[i], 0)
            self.all_unique_before = True
            self._check_all_entries()
        elif not all_unique_now and self.all_unique_before:     # If the RNG seeds were previously unique, but now they’re not...
            self.options_grid[:4] = 2 # Invalidate all RNG seeds.
            self.all_unique_before = False
            self._invalid_opts_first()
        else:                                                   # If the RNG seeds were previously not unique (`self.all_unique_before = False`) and still aren’t (`all_unique_now = False`)...
            pass
    
    
    
    def _validate_str(self, s: str, default: str):
        """
        Checks input change of a string variable in Settings.

        Parameters
        ---
        s: The string from lineEdit to analyze.
        default: The default value to check against.

        Returns:
            :self.options_grid:
                0: Valid: Any string.
                1: Default value.
                2: Invalid: Empty string.
        """
        return 2 if s == "" else np.int64( s == default )
    
    
    
    def _validate_int_ra(self, s: str, default: np.int64, min_max: tuple[np.int64]):
        """
        Checks input change of an integer variable in Settings. Must be within a range where `min` and `max` are included.

        Parameters
        ---
        s: The string from lineEdit to analyze.
        default: The default value to check against.
        range: A tuple of the minimum and maximum values to accept.

        Returns:
            :self.options_grid:
                0: Valid: `min <= np.int64(s) <= max`.
                1: Default value.
                2: Invalid: Out of range or not an integer.
        """
        try:
            inp = np.int64(s)
            if inp == default: return 1
            elif min_max[0] <= inp and inp <= min_max[1]: return 0
            else: return 2
        except: return 2
    
    
    
    def _validate_int_ge(self, s: str, default: np.int64, min: np.int64):
        """
        Checks input change of an integer variable in Settings. Must be greater than or equal to `min`.

        Parameters
        ---
        s: The string from lineEdit to analyze.
        default: The default value to check against.
        min: The minimum value to accept or reject input.

        Returns:
            :self.options_grid:
                0: Valid: `np.int64(s) >= min`.
                1: Default value.
                2: Invalid: `np.int64(s) < min` or not an integer.
        """
        try:    return np.int64( np.int64(s) == default ) if np.int64(s) >= min else 2
        except: return 2
    
    
    
    def _validate_algorithm_lengths(self):
        """
        Checks input change of an algorithm (MG/HMC) length in Settings. Internally uses `_validate_int_ge()`.
        
        One length can be 0, but not both.

        Returns:
            :self.options_grid:
                0: Valid: `np.int64(s)` must be non-negative or 0 without the other string equal to 0.
                1: Default value.
                2: Invalid: `np.int64(s)` is not non-negative or 0 with the other string equal to 0.
        """
        curr = [lE.text() for lE in [self.lineEdit_53, self.lineEdit_55]]
        try:
            HMC_L = np.int64(curr[0])
            MG_L  = np.int64(curr[1])
            both_invalid = HMC_L <= 0 and MG_L <= 0
            self.options_grid[18] = 2 if both_invalid else self._validate_int_ge( curr[1], self.MG_L_default, 0 )
            self.options_grid[19] = 2 if both_invalid else 0 if HMC_L == 0 else self._validate_int_ge( curr[0], self.HMC_L_default, 3 )
        except:
            self.options_grid[18] = self._validate_int_ge( curr[1], self.MG_L_default, 0 )
            try:
                HMC_L = np.int64(curr[0])
                self.options_grid[19] = 0 if HMC_L == 0 else self._validate_int_ge( curr[0], self.HMC_L_default, 3 )
            except: self.options_grid[19] = 2
        self._check_all_entries()
    

    
    def _validate_flt_01(self, s: str, default: np.float64, cond: str):
        """
        Check input change of np.float64 variable in Settings.
        
        Must be between 0 and 1, but `cond` determines if 0 and/or 1 are included.

        Args:
            s: The string from lineEdit to analyze.
            default: The default value to check against.
            cond: The condition to determine if 0 and/or 1 are included.
        
        Returns:
            :self.options_grid:
                0: Valid
                    `cond = "[]"`: Range is `0 <= np.float64(s) <= 1`
                    `cond = "(]"`: Range is `0 <  np.float64(s) <= 1`
                    `cond = "[)"`: Range is `0 <= np.float64(s) <  1`
                    `cond = "()"`: Range is `0 <  np.float64(s) <  1`
                1: Default value.
                2: Invalid: `s` is not a np.float64 or is out of range.
        """
        try: 
            inp = np.float64(s)
            if cond == "[]":    return np.int64( inp == default ) if inp >= 0 and inp <= 1 else 2
            elif cond == "(]":  return np.int64( inp == default ) if inp >  0 and inp <= 1 else 2
            elif cond == "[)":  return np.int64( inp == default ) if inp >= 0 and inp <  1 else 2
            elif cond == "()":  return np.int64( inp == default ) if inp >  0 and inp <  1 else 2
            else: raise Exception("INTERNAL ERROR! Variable `cond = " + cond + "` for function `_validate_flt_01` is invalid!")
        except: return 2
    
    
    
    def _validate_flt_0i(self, s: str, default: np.float64, cond: str):
        """
        Check input change of np.float64 variable in Settings. 
        
        Must be a non-negative float/int, and `cond` determines if 0 is included or not.

        Args:
            s: The string from lineEdit to analyze.
            default: The default value to check against.
            cond: The condition to determine if 0 is included or not.

        Returns:
            :self.options_grid:
                0: Valid
                    `cond = ">"`:  `np.float64(s) > 0`
                    `cond = ">="`: `np.float64(s) >= 0`
                1: Default value.
                2: Invalid: `s` is not a float/int or is out of range.
        """
        try: 
            inp = np.float64(s)
            if cond == ">":     return np.int64( inp == default ) if inp > 0  else 2
            elif cond == ">=":  return np.int64( inp == default ) if inp >= 0 else 2
            else: raise Exception("INTERNAL ERROR! Variable `cond = " + cond + "` for function `_validate_flt_0i` is invalid!")
        except: return 2
    
    
    
    def _validate_tht_prior_psi(self, s: str):  
        """
        Check input change of `tht_prior_psi` in Settings. 

        Args:
            s (str): The string from lineEdit to analyze.

        Returns:
            :self.options_grid:
                0: Valid: Any positive float/int.
                1: Default value. The actual value is `0.5 * np.sum(It_D + It_A)/(T * dD)`
                2: Invalid: Not a positive float/int.
        """
        if s == "-1": return 1
        else:
            try: return 0 if np.float64(s) > 0 else 2
            except: return 2

    
    
    def _check_all_entries(self):
        """
        For every change in Settings, check if Settings can be reset (valid + not default) or if change is invalid.
        """
        is_all_default = all(self.options_grid == 1)
        self.b_reset_was_enabled = not is_all_default
        self.pushButton_3.setEnabled( self.b_reset_was_enabled )
        self.menu7_reset.setEnabled( self.b_reset_was_enabled )
        if is_all_default:                  # If all values are default...
            if self.b_valid_options:        # Default: If all options were valid and are still valid.
                pass
            else:                           # Default: If options were invalid, but not anymore.
                self._complete_opts_restored( "Default" )
        else:
            if any(self.options_grid == 2): # If any values are not default...
                if self.b_valid_options:    # Valid: If options were valid, but are now invalid for the first time.
                    self._invalid_opts_first()
                else:                       # Valid: If options were already invalid and are still invalid.
                    pass
            else:
                if self.b_valid_options:    # Valid: If all options are valid.
                    pass
                else:                       # Valid: If options were invalid, but not anymore.
                    self._complete_opts_restored( "Valid" )





#
# RUN PROGRAM
#
if __name__ == '__main__': 
    win = QMainWindow()
    ui = BayesFRET(win)
    win.show()
    sleep(1)
    splash.finish(win)
    sys.exit( app.exec() )