<p align="center"><img src="https://github.com/jaredhidalgo/BayesFRET/raw/HEAD/resources/BayesFRET%20-%20banner.png"/></p>

Have you ever finished [an smFRET experiment](https://en.wikipedia.org/wiki/Single-molecule_FRET) on [an intrinsically-disordered protein](https://en.wikipedia.org/wiki/Intrinsically_disordered_proteins) and wondered "What do I do with this noisy data?" Maybe you'd like to know how many configurations (conformatoinal states) your enzyme made or how those states behaved during the experiment. Then this Python app is for you!

**BayesFRET** is a Python GUI program of a nonparametric Bayesian method (HDP-HMM) to estimate the number and behavior of distinct conformational states of an intrinsically-disordered molecule (ex. enzymes, Holliday junctions) from a surface-immobilized single-molecule Forster Resonance Energy Transfer (smFRET) experiment. This was originally a Ph.D.-level NSF project that became my [M.S. thesis](https://www.researchgate.net/publication/410979214_An_Experiment-Informed_Hierarchical_Dirichlet_Process-Hidden_Markov_Model_to_Analyze_Surface-Immobilized_smFRET_Data).

* **Data**: Real (experimental) or mimicked (synthetic) time-binned photon intensities and experiment parameters of a molecule. Prior distributions optional.
* **Method**: 4 Markov Chain Monte Carlo (MCMC) simulations on an experiment-informed Hierarchical Dirichlet Process-Hidden Markov Model (HDP-HMM).
* **Results**: A directory of graphs (PNG) and numerical data (Pickle, TXT files) of simulation history and convergence.

For example, an enzyme usually has at least 2 conformational states (open for business, closed for processing), but we don't know the number of states by default.​

Original work and MATLAB app were made by [Sgouralis et al](https://pubs.acs.org/doi/10.1021/acs.jpcb.8b09752). This GUI enhances and expands capabilities from the original work, taking initiative to increase the number of MCMC simulations from 1 to 4 and to increase app efficiency.


# Installation

1. **Compatible Operating Systems**: Successful testing on Windows 10/11 and macOS. Untested for Linux.
1. **Python Interpreter & Environment**
    * Compatible versions: 3.11, 3.12, 3.13
    * This documentation will assume the Python CMD call function is `python`.
1. **Download Files**: Assign a directory specifically for the files of this program (ex. `C:\Users\Cyvu37\Documents\BayesFRET`). 
    * Make sure the parent directory of the program directory (ex. `C:\Users\Cyvu37\Documents`) doesn't require administrator access. BayesFRET must be able to produce subdirectories and files!
    * Processing power may vary based on the directory's location.
1. **Install Packages**: [chime](https://github.com/MaxHalford/chime), [darkdetect](https://github.com/albertosottile/darkdetect), [matplotlib](https://matplotlib.org/), [numpy](https://numpy.org/), [PySide6](https://doc.qt.io/qtforpython-6/), [scipy](https://scipy.org/)
    * If internet is available when BayesFRET opens, the app will auto-run the `python -m pip install -r requirements.txt` command to verify the required packages. Make sure you have [pip](https://pip.pypa.io/en/stable/installation/) and [setuptools](https://pypi.org/project/setuptools/) in your Python environment.
    * The app will also auto-check your Python system for the required packages and the current directory for file dependencies (ex. `code01_classes.py`). Anything missing will cause the program to close.

***

**NOTE: The following content will be moved to the Wiki section.**


# Before Run: The Control Panel

<p align="center"><img src="https://github.com/jaredhidalgo/BayesFRET/raw/HEAD/resources/BayesFRET%20-%20control%20panel.jpg"/></p>

The program can either use data from real experiments ("Experimental") or it can run by itself using randomly generated data simulating data from an smFRET experiment ("Synthetic"). Either way, the program will run the four MCMC simulations one at a time and visualize the results. Graphs and numerical data are saved as soon as they are produced.

## Data Source Requirements
* **Synthetic Mode**: None. The program generates its own data and compares it with the MCMC simulations.
* **Experimental Mode**
    * 2 text/CSV files, each with a single array of numerical values separated by whitespace/commas.
        * 1 for the intensities of the *acceptor* dye throughout the experiment.
        * 1 for the intensities of the *donor* dye throughout the experiment.
    * Parameters of the experiment for the Options section of the GUI.

## GUI Features
TBA

## GUI: The System Tray Icon

<p align="center"><img src="https://github.com/jaredhidalgo/BayesFRET/raw/HEAD/resources/BayesFRET%20-%20system%20icon.jpg"/></p>

When opened, BayesFRET produces a system tray icon in Windows and macOS. The _Status_ and _Settings_ sections of the context menu are shortcuts to commands in the GUI. 

Context Menu CMD | GUI Equivalent
:----------------|:----------------
Run BayesFRET | "Run" button in Settings section (when enabled)
Exit BayesFRET | Closing the GUI
Experimental Data: Import Donor File | Select "Experimental", then "Donor" button
Experimental Data: Import Acceptor File | Select "Experimental", then "Acceptor" button
Synthetic Data: Import Pickle File | Select "Synthetic - Reuse", then "Parameters + True Values" button
Hide/Show Graphs During Run | Boolean to "Show Windows" checkbox
Graph Size*** | Changing radio buttons in Graph Options group
Randomize RNG Seeds | "Randomize!" button in RNG Seeds group
Reset Settings + RNG Seeds | "Reset" button in Options section (when enabled)

**Move GUI**: The commands below the _--Status--_ and _--Settings--_ sections of the context menu deal with positioning the GUI panel to one of the four corners of each active monitor in the local computer, accounting for the task bar in Windows and the menu bar in macOS. Each section corresponds to each active monitor. If you're showing graphs during a run, then these commands can move the GUI away from the open windows. 

**Changing Monitors**: When adding or removing an active monitor on the local computer, the GUI will automatically refresh the icon's context menu with commands for each active monitor.

**Status**: During a run, hovering the mouse over the system tray icon will produce a tooltip showing the status of the program. If the icon's context menu appears during the run, the _--Status--_ header and "Run BayesFRET" command will be replaced with status messages. Moreover, the commands in the _--Settings--_ section will be disabled.

### ***Handling the Graph Size Button

The 5th context menu button under _--Settings--_ may be confusing, so here are two perspectives.

1. **What It Means**: In the _Graph Options_ section of the GUI, when the radio button [A] is selected and the "Show Windows" checkbox is [B], then the text for the 5th context menu button under _--Settings--_ is [C]:

A | B | C
:------------|:--------------|:--------------
"1280 x 720" | enabled | "Maximize Graph Size"
"1280 x 720" | disabled | "Increase Graph Size to 3200 x 1800"
"Maximized" | enabled or disabled | "Decrease Graph Size to 1280 x 720"

2. **What It Does**: When the text for the 5th context menu button under _--Settings--_ is [A], selecting that button will change the radio button to [B] while keeping the "Show Windows" checkbox as [C]:

A | B | C
:------------|:--------------|:--------------
"Maximize Graph Size" | "Maximized" | enabled
"Increase Graph Size to 3200 x 1800" | "3200 x 1800" | disabled
"Decrease Graph Size to 1280 x 720" | "1280 x 720" | enabled or disabled


# During Run: Running the Program

When BayesFRET runs the simulations, the program stores the information inputted through the GUI into a directory with the exact timestamp of when the run started, or **the active directory**. 

BayesFRET uses a Universal class to store various miscellaneous data, including absolute URIs to their original directories. After the program finishes its run, data from the Universal class is saved to replicate the program's output. See the [Custom Classes](#custom-classes) section for more information.

## Avoid These Things While Running

* Don't modify the BayesFRET or active directories, including their URIs and internal files, while BayesFRET runs the simulations!
* Don't interact with the GUI, its windows, or the system tray icon while the simulations are running! This version doesn't support multithreading yet.


# After Run: Interpreting the Results

Each pickle file (`.p`) produced by this program contains data that can be instantly used in a Python environment. However, they contain custom Python classes from `code01_classes.py`. Copy `code01_classes.py` to your current directory to access the documentation of the Python classes. An advanced code editor like Visual Studio Code with extensions can provide the user with previews of documentation regarding those custom classes.

## Importing the Data

The following commands import all of the data produced by a run from the program:

```
import os, pickle

# Full directory path.
DIR = "...{os.sep}BayesFRET_syn - 2024 04 23, 13 12 01 614687"

# Import documentation.
# For this line to work, you must have a copy of `code01_classes.py` in your current directory.
from code01_classes import Chain_History, Params, True_Samples, Universal

# Import the Universal class (`U`).
with open( os.join( DIR, "BayesFRET_Universal_class.p" ), "rb" ) as file:
    U: Universal = pickle.load(file)

# Import the parameters class (`params`). If synthetic, also import the true samples class (`TS`).
with open( os.join( DIR, "BayesFRET_data_params_and_true.p" ), "rb" ) as file:
    if U.is_syn:
        X = pickle.load(file)
        params: Params = X[0]
        TS: True_Samples = X[1]
    else:
        params: Params = pickle.load(file)

# Import the sampling history for each simulation (`S[0]` - `S[3]`).
S: list[Chain_History] = []
FILENAMES = [os.join( DIR, "BayesFRET_data_sim_{i}_history.p" ) for i in U.list_seeds]
for v in U.range_seeds:
    with open(FILENAMES[v], "rb") as f:
        S.append( pickle.load(f) )
```

After that, you can use the variables `U`, `params`, `TS`, and `S` to your liking.

## Custom Classes

* `Chain_History`: History of samples regarding key variables throughout the MCMC simulation.
    * Experimental + synthetic data: Found in `BayesFRET_data_sim_{i}_history.p` where `{i}` is the RNG seed value.
* `Params`: Fixed parameters such as the options from the Algorithm Settings, Photoemission Priors, and Photophysics Priors sections of the GUI. 
    * Experimental data: Found in `BayesFRET_data_params.p`.
    * Synthetic data: Found in `BayesFRET_data_params_and_true.p`.
* `True_Samples`: Variables used in making the synthetic data from `code02_setup.py`.
    * Synthetic data only: Found in `BayesFRET_data_params_and_true.p`.
* `Universal`: Holds all other miscellaneous variables, including...
    * Filenames (not filepaths) of the Chain_History pickle files (`U.filenames`).\*
    * List of the unique RNG seed values (`U.list_seeds`).
    * List of the RNG classes corresponding to each RNG seed value (`U.RNGs`).
    * Number of Posterior Quantiles (`U.Q`).
    * The full path of the active directory (`U.DIR_ACTIVE`)\* and a function to merge the path with a filename (`U.func_getActivePath(...)`).\*
    * The full path of the program's original directory (`U.DIR_PROGRAM`).\*

\*Changes in a filename or a directory path affect these variables and functions. The following code demonstrates how to save those changes in the Universal file (`BayesFRET_Universal_class.p`).

```
# >> How to Handle Name Changes <<
# directory name change: DIR --> "...\BayesFRET_syn FINAL"
# Filename change for 3rd seed: 'BayesFRET_data_sim_627_history.p' --> 'BayesFRET_627.p'

DIR2 = "...{os.sep}BayesFRET_syn FINAL"
with open( os.join( DIR2, "BayesFRET_Universal_class.p" ), "rb" ) as file:
    U: Universal = pickle.load(file)
U.DIR_ACTIVE = DIR2                 # This updates the `FILEPATH` function.
U.filenames[2] = 'BayesFRET_627.p'  # This allows you to use the previous for loop to iterate all filenames.
...
with open(U.FILEPATH( os.join( DIR2, "BayesFRET_Universal_class.p" ) ), "wb") as file:
    pickle.dump( U, file )
```


## Usage Examples

1. Using the default RNG seeds (5, 44, 356, 9918), use the following commands to extract the FRET efficiencies `eff` of each unique state at the last MCMC iteration `n=1000` for the MCMC simulation with RNG seed 5 (`S[0]`).

```
n = 1000

# The set of unique conformational states at the MCMC iteration n=1000.
ks = S[0].K_set[n]

# Only extract the FRET efficiencies (`eff`) of the unique conformational states (`ks`) at n=1000. Results between 0 and 1.
return S[0].eff[n, ks]
```

2. Find the evaluated conformational states at time t (`st`) of the acceptor dye's photoemission rates $\lambda^A_{s_t}$ at the MCMC iteration `n=649` for the MCMC simulation with RNG seed 9918 (`S[3]`).

```
return S[3].lam_A[ S[3].st[649] ]
```
