# bAIes-IDP-OpenMM

You can run bAIes with openMM in this directory.

You will find here the following directories:
* `scripts`, which contains the necessary python scripts to run bAIes
* `predictions`, the directory where you can add your predictions. You will find an example subdirectory with a prediction of PaaA2.
* `force_field`, which contains the coil prior force field.
* `bAIes_simulation`, which contains an example output for the preparation of a bAIes simulation of PaaA2.

## Quickstart

### 0. Run a structure prediction

You can use AlphaFold2, AlphaFold3, Boltz-1, Boltz-2, ESMFold or ESMFold2.

Collect a predicted model (.pdb) and the distograms (.pkl) and deposit them in a subdirectory in the `predictions` directory.

If the predicted model is in mmCIF format (.cif), you can convert it to PDB by running the following on the terminal. Make sure to replace the elements in <> brackets by the desired file paths.

```bash
scripts/cif2pdb.py <path_to_cif_file>.cif <path_to_output_file>.pdb
```

Ensure that your structure predictor outputs the distograms [here](https://github.com/COSBlab/bAIes-IDP-OpenMM/blob/HEAD/tutorial/link)

### 1. Prepare the inputs

Run on the terminal:

```bash
./1-prepare_bAIes.bash -pdb <path_to_prediction_pdb>.pdb \
                       -pkl <path_to_distograms>.pkl \
                       -out <output_directory>
```

This command will read the pdb model, protonate it if necessary and deposit the protonated version in the <output_directory> under the same name. It will also read the distograms, analyze them and output the bAIes input file in the <output_directory>.

### 2. Run bAIes.

Run on the terminal:

```bash
./2-run_bAIes.bash -pdb <path_to_bAIes_pdb>.pdb \
                   -baies <path_to_baies_params>.dat \
                   -out <output_directory>
```

You can run a `coil` simulation by omitting the -baies argument.

## Example

To run a bAIes simulation of PaaA2, deposit the predicted PDB model (PaaA2.pdb) and associated distograms (PaaA2.pkl) in Predictions/PaaA2. Then run the following to prepare the inputs:

```bash
./1-prepare_bAIes.bash -pdb predictions/PaaA2/PaaA2.pdb \
                       -pkl predictions/PaaA2/PaaA2.pkl \
                       -out bAIes_simulations/PaaA2
```

Finally, after checking the simulation parameters, run the following:

```bash
./2-run_bAIes.bash -pdb bAIes_simulations/PaaA2/PaaA2.pdb \
                   -baies bAIes_simulations/PaaA2/bAIes_params.dat \
                   -out bAIes_simulations/PaaA2
```

## Tuning the simulation parameters

`2-run_bAIes.bash` calls the `bAIes_openMM.py` script in the `scripts` directory. You can directly tune the wanted number of steps, simulation time step and frame saving interval in `2-run_bAIes.bash`. You can also personalize bAIes_openMM.py with your desired features.
