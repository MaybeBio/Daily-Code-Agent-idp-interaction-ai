[![DOI](https://zenodo.org/badge/1219039183.svg)](https://doi.org/10.5281/zenodo.19709395)
# Command: enzywizard-batch

EnzyWizard-Batch is a command-line tool for running a complete EnzyWizard
analysis workflow. It performs multiple EnzyWizard modules, including residue property analysis,
hydrophobic cluster detection, energy evaluation, flexibility analysis, 
disordered region detection, optional conservation analysis, residue embedding generation,
binding pocket detection, optional substrate
feature generation, optional enzyme-substrate docking, molecular interaction calculation,
and final graph integration.
If substrate names or SMILES strings are provided, the full enzyme-substrate
workflow will be executed, including substrate preparation, docking,
enzyme-substrate molecular interaction detection, and strict graph integration when
all required side reports are available.
If no substrate input is provided, the program generates a protein-only integrated graph 
based on residue-level features and intra-enzyme molecular interactions.
If no MSA input is provided, conservation analysis is skipped and the integrated
residue nodes omit the conservation score field.
The final output is an integrated graph dataset that can be directly used for
graph-based analysis, machine learning, and enzyme function studies.


# Documentation index:

- example usage
- input parameters
- output files
- output report schema
- Process
- common errors and solutions
- dependencies
- references


# example usage:

The examples below use placeholder paths such as `path/to/input.cif`,
`path/to/alignment.sto`, and `path/to/output_dir/`; replace them with your own
cleaned protein structure file, optional matched MSA file, optional substrate
input, and output directory. The input structure must already be cleaned. MSA
input is optional: providing `-m` enables conservation analysis, cleaned MSA
output, and HMM profile output. Substrate input is optional: providing names or
SMILES strings enables substrate generation, docking, and protein-substrate
interaction detection. If substrate generation or docking fails, batch falls
back to the protein-only route.

Run the default protein-only workflow from a cleaned CIF structure.

```
enzywizard-batch -i path/to/input.cif -o path/to/output_default/
```

Run the default protein-only workflow from a cleaned PDB structure.

```
enzywizard-batch -i path/to/input.pdb -o path/to/output_pdb/
```

Run a protein-only workflow with a Stockholm MSA.

```
enzywizard-batch -i path/to/input.cif -m path/to/alignment.sto -o path/to/output_sto/
```

Run a protein-only workflow with an aligned FASTA MSA.

```
enzywizard-batch -i path/to/input.cif -m path/to/alignment.fasta -o path/to/output_fasta_msa/
```

Run a protein-only workflow with a gzip-compressed aligned FASTA MSA.

```
enzywizard-batch -i path/to/input.cif -m path/to/alignment.fasta.gz -o path/to/output_fasta_gz_msa/
```

Run a protein-only workflow with an A3M MSA.

```
enzywizard-batch -i path/to/input.cif -m path/to/alignment.a3m -o path/to/output_a3m/
```

Run the enzyme-substrate workflow for one named substrate. Substrate names are
resolved to SMILES through external chemical databases.

```
enzywizard-batch -i path/to/input.cif -s "glucose" -o path/to/output_glucose/
```

Run the enzyme-substrate workflow for two named substrates.

```
enzywizard-batch -i path/to/input.cif -m path/to/alignment.sto -s "glucose;fructose" -o path/to/output_named_substrates/
```

Run the full workflow with direct SMILES input. Direct SMILES skips external
name-to-SMILES lookup.

```
enzywizard-batch -i path/to/input.cif -m path/to/alignment.sto -s "CCO" -o path/to/output_smiles/
```

Run a mixed substrate input with one name and one SMILES string.

```
enzywizard-batch -i path/to/input.cif -m path/to/alignment.sto -s "glucose;CCO" -o path/to/output_mixed_substrates/
```

Use long option names for the same full workflow.

```
enzywizard-batch --cleaned_input_path path/to/input.cif --input_msa path/to/alignment.sto --substrate_names "glucose;fructose" --output_dir path/to/output_long_options/
```

Keep intermediate files such as cleaned MSA, HMM profile, substrate SDF files,
docked SDF files, and enzyme-substrate complex files. This is useful for
debugging or inspecting downstream inputs, but writes more files to disk.

```
enzywizard-batch -i path/to/input.cif -m path/to/alignment.sto -s "glucose;fructose" -o path/to/output_with_intermediates/ --save_extra_outputs
```

Use fewer energy minimization iterations. Smaller values reduce minimization
work and leave the structure closer to the starting conformation.

```
enzywizard-batch -i path/to/input.cif -m path/to/alignment.sto -o path/to/output_energy_20/ --energy_minimization_iteration 20
```

Use fewer substrate synonyms during name resolution. This can reduce API
requests and runtime, but may miss difficult or ambiguous substrate names that
need synonym-expanded matching.

```
enzywizard-batch -i path/to/input.cif -m path/to/alignment.sto -s "glucose;fructose" -o path/to/output_fast_lookup/ --substrate_max_synonyms 5
```

Use a catalytic residue as the docking box center. The residue index is the
cleaned protein residue index, and the residue CA atom coordinate is used as the
center. This skips PyVOL pocket detection and the global docking box fallback.
A smaller box focuses the search and can run faster, but may miss valid poses
outside the box. A larger box explores a broader region, but can increase
runtime and reduce search precision at the same exhaustiveness.

```
enzywizard-batch -i path/to/input.cif -m path/to/alignment.sto -s "glucose" -o path/to/output_catalytic_residue/ --dock_catalytic_residue 121 --dock_box_size 20,20,20
```

Use an explicit catalytic-site coordinate as the docking box center. This is
useful when the active-site coordinate is known from another analysis or a
reference structure.

```
enzywizard-batch -i path/to/input.cif -m path/to/alignment.sto -s "glucose" -o path/to/output_site_coord/ --dock_catalytic_site_coord 12.5,8.0,-3.2 --dock_box_size 18,18,18
```

Increase Vina exhaustiveness for a broader docking search. Larger values may
improve search coverage and docking robustness, but increase runtime. Smaller
values are faster but may miss better poses.

```
enzywizard-batch -i path/to/input.cif -m path/to/alignment.sto -s "glucose;fructose" -o path/to/output_high_exhaustiveness/ --dock_exhaustiveness 32
```

Disable docking early stop so batch continues after the first successful docking
result and searches other conformer and box combinations up to
`--dock_max_attempt_num`. This can improve the chance of finding a lower-energy
pose, but increases runtime.

```
enzywizard-batch -i path/to/input.cif -m path/to/alignment.sto -s "glucose;fructose" -o path/to/output_no_early_stop/ --dock_no_early_stop --dock_max_attempt_num 40
```

# input parameters:

-i, --cleaned_input_path
Required.
Path to the input cleaned protein structure file.
Supported file extensions: .cif, .pdb.

The file must:
- already be cleaned
- contain a valid single protein chain
- contain hydrogen atoms
- match the sequence used to generate the input MSA when `--input_msa` is provided

-m, --input_msa
Optional.
Path to the input MSA file.

Supported MSA formats include:
- Stockholm (.sto)
- aligned FASTA (.fa / .fasta / .afa / .fasta.gz)
- A3M

When provided, the MSA must be generated using the cleaned protein FASTA
sequence. When omitted, conservation analysis, cleaned MSA output, and HMM
profile output are skipped.

-s, --substrate_names
Optional.
Substrate names or SMILES strings.

Multiple substrates should be separated by ';'.

If provided, the following additional workflows will be executed:
- substrate feature generation
- substrate 3D structure generation
- enzyme-substrate docking
- enzyme-substrate molecular interaction calculation
- strict graph integration when an MSA is also provided

If substrate generation or docking fails, these additional workflows are skipped
for the final integration and the program continues with the protein-only route.

If not provided, substrate, docking, and enzyme-substrate molecular interaction steps
will be skipped.

-o, --output_dir
Required.
Directory to save batch outputs.

--save_extra_outputs
Optional.
Enable keeping intermediate and side output files.

By default, this option is disabled, and only the final integrated JSON outputs
and log.txt are kept.

When enabled, additional files such as cleaned MSA, HMM profile, substrate SDF
files, docked substrate files, and enzyme-substrate complex files may be saved.
Cleaned MSA and HMM files are generated only when `--input_msa` is provided.

--hydrocluster_cutoff
Optional.
Minimum contact area cutoff for hydrophobic cluster residue-residue connection.
Unit: square angstroms (A^2).
Default: 10.0.

--no_minimize_energy
Optional.
Disable energy minimization before energy evaluation.
By default, energy minimization is enabled.

--energy_minimization_iteration
Optional.
Maximum number of iterations for energy minimization.
Default: 100.

--flexibility_method
Optional.
Normal mode method for RMSF calculation.
Choices:
- ANM
- GNM
Default: ANM.

--flexibility_cutoff
Optional.
Distance cutoff used to determine residue-residue connections in ProDy.
Unit: angstroms (A).
Default: 15.0.

--flexibility_n_modes
Optional.
Number of low-frequency normal modes used for RMSF calculation.
Default: 20.

--disorder_window_size
Optional.
Sliding window size for FoldIndex-like disordered region score calculation.
Default: 11.

--disorder_min_region_length
Optional.
Minimum number of consecutive residues required to define a disordered region.
Default: 5.

--embedding_model_name
Optional.
ESM2 model used for residue embedding generation.

Choices:
- esm2_t6_8M_UR50D
- esm2_t12_35M_UR50D
- esm2_t30_150M_UR50D

Default: esm2_t6_8M_UR50D.

--pocket_min_rad
Optional.
Minimum probe radius used by PyVOL for binding pocket detection.
Unit: angstroms (A).
Default: 1.8.

--pocket_max_rad
Optional.
Maximum probe radius used by PyVOL for binding pocket detection.
Unit: angstroms (A).
Default: 6.2.

--pocket_min_volume
Optional.
Minimum binding pocket volume threshold.
Unit: cubic angstroms (A^3).
Default: 50.

--substrate_max_synonyms
Optional.
Maximum number of substrate synonyms retried when fetching SMILES from a
substrate name.
Default: 20.

--substrate_fp_radius
Optional.
Radius used for Morgan fingerprint generation.
Default: 2.

--substrate_n_bits
Optional.
Bit size of the Morgan fingerprint vector.
Default: 512.

--substrate_num_confs
Optional.
Maximum number of 3D conformers generated for each substrate.
Default: 5.

--substrate_prune_rms
Optional.
RMS threshold used to prune highly similar conformers during 3D conformer
generation.
Unit: angstroms (A).
Default: 0.5.

--dock_max_attempt_num
Optional.
Maximum number of docking attempts.
Default: 20.

--dock_no_early_stop
Optional.
Disable stopping immediately after the first successful docking result.

By default, early stopping is enabled.

--dock_exhaustiveness
Optional.
Exhaustiveness of AutoDock Vina search.
Default: 8.

--dock_cpu
Optional.
Number of CPUs used by AutoDock Vina.
Default: 0.

--dock_catalytic_residue
Optional.
Cleaned protein residue index used as the docking box center.

Example:
  121

--dock_catalytic_site_coord
Optional.
Catalytic site center coordinate separated by ','.
Unit: angstroms (A).

Example:
  12.5,8.0,-3.2

--dock_box_size
Optional.
Docking box size separated by ','.
Unit: angstroms (A).

Example:
  20,20,20

This parameter is required when --dock_catalytic_residue or
--dock_catalytic_site_coord is provided. All three values must be positive
numbers.

--hbond_bonded_h_min_distance
Optional.
Minimum bonded heavy atom-hydrogen distance used for hydrogen bond donor
detection.
Unit: angstroms (A).
Default: 0.8.

--hbond_bonded_h_max_distance
Optional.
Maximum bonded heavy atom-hydrogen distance used for hydrogen bond donor
detection.
Unit: angstroms (A).
Default: 1.3.

--hbond_da_max_distance
Optional.
Maximum donor-acceptor distance cutoff for hydrogen bond detection.
Unit: angstroms (A).
Default: 3.9.

--hbond_ha_max_distance
Optional.
Maximum hydrogen-acceptor distance cutoff for hydrogen bond detection.
Unit: angstroms (A).
Default: 2.5.

--hbond_angle
Optional.
Minimum donor-hydrogen-acceptor angle cutoff for hydrogen bond detection.
Unit: degrees.
Default: 90.0.

--ionic_distance_cutoff
Optional.
Maximum distance cutoff for ionic bond detection.
Unit: angstroms (A).
Default: 4.0.

--vdw_mu
Optional.
Mu parameter used in van der Waals interaction detection.
Unit: dimensionless.
Default: 0.01.

--ppstack_center_distance_cutoff
Optional.
Maximum ring-center distance cutoff for pi-pi stacking detection.
Unit: angstroms (A).
Default: 6.5.

--pication_distance_cutoff
Optional.
Maximum ring-cation distance cutoff for pi-cation interaction detection.
Unit: angstroms (A).
Default: 5.0.

--pication_angle_cutoff
Optional.
Maximum angle cutoff for pi-cation interaction detection.
Unit: degrees.
Default: 45.0.

--ssbond_max_distance
Optional.
Maximum sulfur-sulfur distance cutoff for disulfide bond detection.
Unit: angstroms (A).
Default: 2.5.


# output files:

The program always keeps the following files in the output directory:

`{protein_name}` is derived from the cleaned input structure file name.

1. An integrated JSON report
   - integrate_report_{protein_name}.json
     - Full integrated report containing overall statistics and integrated graph entries.

2. A node-only JSON file
   - integrate_nodes_{protein_name}.json
     - Array of integrated graph node records split from the integrated graph.

3. An edge-only JSON file
   - integrate_edges_{protein_name}.json
     - Array of integrated graph edge records split from the integrated graph.

4. A log file
   - log.txt
     - Processing log containing informational messages and errors.

When `--save_extra_outputs` is enabled, batch may also keep intermediate files
generated by the enabled workflow steps:

5. A cleaned Stockholm MSA file
   - cleaned_{msa_name}.sto
     - Cleaned MSA in Stockholm format. This file is generated only when
       `--input_msa` is provided.

6. A profile HMM file
   - hmm_profile_{msa_name}.hmm
     - HMM profile generated from the cleaned Stockholm MSA. This file is
       generated only when `--input_msa` is provided.

7. Substrate structure files in SDF format
   - {substrate_structure_name}.sdf
     - Generated 3D substrate conformer files. These files are generated only
       when `--substrate_names` is provided and substrate structure generation
       succeeds.
     - `{substrate_structure_name}` is cleaned for filesystem compatibility;
       characters such as commas and spaces are converted to underscores. The
       original resolved name is retained in the report's `substrate_name` field.

8. Docked substrate structure files in SDF format
   - docked_{substrate_name}.sdf
     - Docked SDF file for each substrate in the selected docking result. These
       files are generated only when `--substrate_names` is provided and docking
       succeeds.

9. Docked enzyme-substrate complex structure files
   - docked_{protein_name}_{substrate_names}.cif
     - Docked enzyme-substrate complex structure in CIF format.
   - docked_{protein_name}_{substrate_names}.pdb
     - Docked enzyme-substrate complex structure in PDB format.
     - These files are generated only when `--substrate_names` is provided and
       docking succeeds.


# output report schema:

The JSON report contains the following fields:

   - "report_type"
     - Data type: string
     - Expected value: "enzywizard_integrate"
     - Description: The field 'report_type' indicates the type of report ('report': http://purl.obolibrary.org/obo/IAO_0000088) generated by the EnzyWizard-Integrate software.

   - "overall_statistics"
     - Data type: object
     - Description: The field 'overall_statistics' indicates the overall summary statistics ('statistics': http://purl.obolibrary.org/obo/STATO_0000039) integrated from EnzyWizard reports ('report': http://purl.obolibrary.org/obo/IAO_0000088).

     The "overall_statistics" object may contain:

     - "sequence_length"
       - Data type: integer
       - Description: The field 'sequence_length' indicates the sequence length ('sequence length': http://edamontology.org/data_1249), measured as the number of amino acid residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) in the cleaned protein sequence ('protein sequence': http://edamontology.org/data_2976).

     - "total_molecular_weight"
       - Data type: number
       - Description: The field 'total_molecular_weight' indicates the total molecular weight, calculated as the sum of residue molecular weights ('molecular weight': https://goldbook.iupac.org/terms/view/R05271) across the protein sequence ('protein sequence': http://edamontology.org/data_2976). Unit: daltons (Da) ('dalton': http://qudt.org/vocab/unit/DA).

     - "total_net_charge"
       - Data type: number
       - Description: The field 'total_net_charge' indicates the total net charge, calculated as the sum of residue electric charges ('electric charge': https://goldbook.iupac.org/terms/view/E01923) across the protein sequence ('protein sequence': http://edamontology.org/data_2976). Unit: dimensionless ('dimensionless': http://qudt.org/vocab/unit/UNITLESS).

     - "total_residue_volume"
       - Data type: number
       - Description: The field 'total_residue_volume' indicates the total residue volume ('volume': http://purl.obolibrary.org/obo/PATO_0000918), calculated as the sum of residue volumes across the protein sequence ('protein sequence': http://edamontology.org/data_2976). Unit: cubic angstroms (Å^3) ('cubic angstrom': http://qudt.org/vocab/unit/ANGSTROM3).

     - "max_3d_diameter"
       - Data type: number
       - Description: The field 'max_3d_diameter' indicates the maximum three-dimensional diameter ('diameter': http://purl.obolibrary.org/obo/PATO_0001334), calculated as the maximum pairwise distance ('distance': http://purl.obolibrary.org/obo/PATO_0000040) between residue alpha-carbon coordinates. Unit: angstroms (Å) ('angstrom': http://qudt.org/vocab/unit/ANGSTROM).

     - "radius_of_gyration"
       - Data type: number
       - Description: The field 'radius_of_gyration' indicates the radius of gyration ('radius of gyration': https://goldbook.iupac.org/terms/view/R05121) calculated from residue alpha-carbon coordinates. Unit: angstroms (Å) ('angstrom': http://qudt.org/vocab/unit/ANGSTROM).

     - "asphericity"
       - Data type: number
       - Description: The field 'asphericity' indicates the asphericity ('asphericity': https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html) calculated from alpha-carbon coordinate covariance eigenvalues. Unit: dimensionless ('dimensionless': http://qudt.org/vocab/unit/UNITLESS).

     - "spherocity"
       - Data type: number
       - Description: The field 'spherocity' indicates the spherocity ('spherocity': https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html) calculated from alpha-carbon coordinate covariance eigenvalues. Unit: dimensionless ('dimensionless': http://qudt.org/vocab/unit/UNITLESS).

     - "principal_moment_ratio"
       - Data type: number
       - Description: The field 'principal_moment_ratio' indicates the ratio of the largest to the smallest principal moments ('moment of inertia': https://goldbook.iupac.org/terms/view/M04006) calculated from alpha-carbon coordinate covariance eigenvalues. Unit: dimensionless ('dimensionless': http://qudt.org/vocab/unit/UNITLESS).

     - "bounding_box_volume"
       - Data type: number
       - Description: The field 'bounding_box_volume' indicates the volume ('volume': http://purl.obolibrary.org/obo/PATO_0000918) of the axis-aligned bounding box ('bounding box': https://developer.mozilla.org/en-US/docs/Glossary/Bounding_box) enclosing all residue alpha-carbon coordinates. Unit: cubic angstroms (Å^3) ('cubic angstrom': http://qudt.org/vocab/unit/ANGSTROM3).

     - "mean_pairwise_ca_distance"
       - Data type: number
       - Description: The field 'mean_pairwise_ca_distance' indicates the mean pairwise alpha-carbon distance ('distance': http://purl.obolibrary.org/obo/PATO_0000040) between residue alpha-carbon coordinates. Unit: angstroms (Å) ('angstrom': http://qudt.org/vocab/unit/ANGSTROM).

     - "std_pairwise_ca_distance"
       - Data type: number
       - Description: The field 'std_pairwise_ca_distance' indicates the standard deviation ('standard deviation': http://purl.obolibrary.org/obo/STATO_0000237) of pairwise alpha-carbon distances ('distance': http://purl.obolibrary.org/obo/PATO_0000040) between residue alpha-carbon coordinates. Unit: angstroms (Å) ('angstrom': http://qudt.org/vocab/unit/ANGSTROM).

     - "residue_name_alanine_count"
       - Data type: integer
       - Description: The field 'residue_name_alanine_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is alanine ('alanine': http://purl.obolibrary.org/obo/CHEBI_16977).

     - "residue_name_cysteine_count"
       - Data type: integer
       - Description: The field 'residue_name_cysteine_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is cysteine ('cysteine': http://purl.obolibrary.org/obo/CHEBI_15356).

     - "residue_name_aspartic_acid_count"
       - Data type: integer
       - Description: The field 'residue_name_aspartic_acid_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is aspartic acid ('aspartic acid': http://purl.obolibrary.org/obo/CHEBI_22660).

     - "residue_name_glutamic_acid_count"
       - Data type: integer
       - Description: The field 'residue_name_glutamic_acid_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is glutamic acid ('glutamic acid': http://purl.obolibrary.org/obo/CHEBI_18237).

     - "residue_name_phenylalanine_count"
       - Data type: integer
       - Description: The field 'residue_name_phenylalanine_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is phenylalanine ('phenylalanine': http://purl.obolibrary.org/obo/CHEBI_28044).

     - "residue_name_glycine_count"
       - Data type: integer
       - Description: The field 'residue_name_glycine_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is glycine ('glycine': http://purl.obolibrary.org/obo/CHEBI_15428).

     - "residue_name_histidine_count"
       - Data type: integer
       - Description: The field 'residue_name_histidine_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is histidine ('histidine': http://purl.obolibrary.org/obo/CHEBI_27570).

     - "residue_name_isoleucine_count"
       - Data type: integer
       - Description: The field 'residue_name_isoleucine_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is isoleucine ('isoleucine': http://purl.obolibrary.org/obo/CHEBI_24898).

     - "residue_name_lysine_count"
       - Data type: integer
       - Description: The field 'residue_name_lysine_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is lysine ('lysine': http://purl.obolibrary.org/obo/CHEBI_25094).

     - "residue_name_leucine_count"
       - Data type: integer
       - Description: The field 'residue_name_leucine_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is leucine ('leucine': http://purl.obolibrary.org/obo/CHEBI_25017).

     - "residue_name_methionine_count"
       - Data type: integer
       - Description: The field 'residue_name_methionine_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is methionine ('methionine': http://purl.obolibrary.org/obo/CHEBI_16811).

     - "residue_name_asparagine_count"
       - Data type: integer
       - Description: The field 'residue_name_asparagine_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is asparagine ('asparagine': http://purl.obolibrary.org/obo/CHEBI_22653).

     - "residue_name_proline_count"
       - Data type: integer
       - Description: The field 'residue_name_proline_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is proline ('proline': http://purl.obolibrary.org/obo/CHEBI_17203).

     - "residue_name_glutamine_count"
       - Data type: integer
       - Description: The field 'residue_name_glutamine_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is glutamine ('glutamine': http://purl.obolibrary.org/obo/CHEBI_18050).

     - "residue_name_arginine_count"
       - Data type: integer
       - Description: The field 'residue_name_arginine_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is arginine ('arginine': http://purl.obolibrary.org/obo/CHEBI_29016).

     - "residue_name_serine_count"
       - Data type: integer
       - Description: The field 'residue_name_serine_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is serine ('serine': http://purl.obolibrary.org/obo/CHEBI_17822).

     - "residue_name_threonine_count"
       - Data type: integer
       - Description: The field 'residue_name_threonine_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is threonine ('threonine': http://purl.obolibrary.org/obo/CHEBI_16857).

     - "residue_name_valine_count"
       - Data type: integer
       - Description: The field 'residue_name_valine_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is valine ('valine': http://purl.obolibrary.org/obo/CHEBI_27266).

     - "residue_name_tryptophan_count"
       - Data type: integer
       - Description: The field 'residue_name_tryptophan_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is tryptophan ('tryptophan': http://purl.obolibrary.org/obo/CHEBI_27897).

     - "residue_name_tyrosine_count"
       - Data type: integer
       - Description: The field 'residue_name_tyrosine_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) whose residue name ('residue name': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) is tyrosine ('tyrosine': http://purl.obolibrary.org/obo/CHEBI_18186).

     - "residue_chemical_classification_uncharged_polar_count"
       - Data type: integer
       - Description: The field 'residue_chemical_classification_uncharged_polar_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) assigned to the uncharged polar residue chemical classification ('chemical classification': http://purl.obolibrary.org/obo/NCIT_C25161; 'uncharged polar': https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html).

     - "residue_chemical_classification_positively_charged_count"
       - Data type: integer
       - Description: The field 'residue_chemical_classification_positively_charged_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) assigned to the positively charged residue chemical classification ('chemical classification': http://purl.obolibrary.org/obo/NCIT_C25161; 'positively charged': https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html).

     - "residue_chemical_classification_negatively_charged_count"
       - Data type: integer
       - Description: The field 'residue_chemical_classification_negatively_charged_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) assigned to the negatively charged residue chemical classification ('chemical classification': http://purl.obolibrary.org/obo/NCIT_C25161; 'negatively charged': https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/IMGTclasses.html).

     - "residue_chemical_classification_hydrophobic_count"
       - Data type: integer
       - Description: The field 'residue_chemical_classification_hydrophobic_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) assigned to the hydrophobic residue chemical classification ('chemical classification': http://purl.obolibrary.org/obo/NCIT_C25161; 'hydrophobic': https://goldbook.iupac.org/terms/view/HT06964).

     - "residue_chemical_classification_aromatic_count"
       - Data type: integer
       - Description: The field 'residue_chemical_classification_aromatic_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) assigned to the aromatic residue chemical classification ('chemical classification': http://purl.obolibrary.org/obo/NCIT_C25161; 'aromatic': https://goldbook.iupac.org/terms/view/A00441).

     - "residue_chemical_classification_aliphatic_count"
       - Data type: integer
       - Description: The field 'residue_chemical_classification_aliphatic_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) assigned to the aliphatic residue chemical classification ('chemical classification': http://purl.obolibrary.org/obo/NCIT_C25161; 'aliphatic': https://goldbook.iupac.org/terms/view/A00217).

     - "residue_chemical_classification_heterocyclic_count"
       - Data type: integer
       - Description: The field 'residue_chemical_classification_heterocyclic_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) assigned to the heterocyclic residue chemical classification ('chemical classification': http://purl.obolibrary.org/obo/NCIT_C25161; 'heterocyclic': https://goldbook.iupac.org/terms/view/H02798).

     - "residue_chemical_classification_sulfur_containing_count"
       - Data type: integer
       - Description: The field 'residue_chemical_classification_sulfur_containing_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) assigned to the sulfur-containing residue chemical classification ('chemical classification': http://purl.obolibrary.org/obo/NCIT_C25161; 'sulfur': http://purl.obolibrary.org/obo/CHEBI_26833).

     - "secondary_structure_unassigned_count"
       - Data type: integer
       - Description: The field 'secondary_structure_unassigned_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) without assigned residue secondary structure ('secondary structure': http://edamontology.org/operation_1847).

     - "secondary_structure_alpha_helix_count"
       - Data type: integer
       - Description: The field 'secondary_structure_alpha_helix_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) assigned to alpha-helix residue secondary structures ('secondary structure': http://edamontology.org/operation_1847; 'alpha helix': https://manual.gromacs.org/current/onlinehelp/gmx-dssp.html).

     - "secondary_structure_beta_bridge_count"
       - Data type: integer
       - Description: The field 'secondary_structure_beta_bridge_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) assigned to beta-bridge residue secondary structures ('secondary structure': http://edamontology.org/operation_1847; 'beta bridge': https://manual.gromacs.org/current/onlinehelp/gmx-dssp.html).

     - "secondary_structure_extended_strand_count"
       - Data type: integer
       - Description: The field 'secondary_structure_extended_strand_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) assigned to extended-strand residue secondary structures ('secondary structure': http://edamontology.org/operation_1847; 'extended strand': https://manual.gromacs.org/current/onlinehelp/gmx-dssp.html).

     - "secondary_structure_three_ten_helix_count"
       - Data type: integer
       - Description: The field 'secondary_structure_three_ten_helix_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) assigned to 3-10 helix residue secondary structures ('secondary structure': http://edamontology.org/operation_1847; '3-10 helix': https://manual.gromacs.org/current/onlinehelp/gmx-dssp.html).

     - "secondary_structure_pi_helix_count"
       - Data type: integer
       - Description: The field 'secondary_structure_pi_helix_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) assigned to pi-helix residue secondary structures ('secondary structure': http://edamontology.org/operation_1847; 'pi helix': https://manual.gromacs.org/current/onlinehelp/gmx-dssp.html).

     - "secondary_structure_turn_count"
       - Data type: integer
       - Description: The field 'secondary_structure_turn_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) assigned to turn residue secondary structures ('secondary structure': http://edamontology.org/operation_1847; 'turn': https://manual.gromacs.org/current/onlinehelp/gmx-dssp.html).

     - "secondary_structure_bend_count"
       - Data type: integer
       - Description: The field 'secondary_structure_bend_count' indicates the count of residues ('residue': http://purl.obolibrary.org/obo/GENO_0000782) assigned to bend residue secondary structures ('secondary structure': http://edamontology.org/operation_1847; 'bend': https://manual.gromacs.org/current/onlinehelp/gmx-dssp.html).

     - "hydrophobic_cluster_count"
       - Data type: integer
       - Description: The field 'hydrophobic_cluster_count' indicates the count of hydrophobic clusters ('hydrophobic cluster': https://proteintools.uni-bayreuth.de/clusters/).

     - "max_hydrophobic_cluster_area"
       - Data type: number
       - Description: The field 'max_hydrophobic_cluster_area' indicates the maximum area ('area': http://purl.obolibrary.org/obo/PATO_0001323) of hydrophobic clusters ('hydrophobic cluster': https://proteintools.uni-bayreuth.de/clusters/). Unit: square angstroms (Å^2) ('angstrom': http://qudt.org/vocab/unit/ANGSTROM).

     - "total_hydrophobic_cluster_area"
       - Data type: number
       - Description: The field 'total_hydrophobic_cluster_area' indicates the total area ('area': http://purl.obolibrary.org/obo/PATO_0001323) of hydrophobic clusters ('hydrophobic cluster': https://proteintools.uni-bayreuth.de/clusters/). Unit: square angstroms (Å^2) ('angstrom': http://qudt.org/vocab/unit/ANGSTROM).

     - "disordered_region_count"
       - Data type: integer
       - Description: The field 'disordered_region_count' indicates the count of intrinsically disordered regions ('intrinsically disordered region': https://disprot.org/ontology).

     - "max_disordered_region_length"
       - Data type: integer
       - Description: The field 'max_disordered_region_length' indicates the maximum sequence length ('sequence length': http://edamontology.org/data_1249) of intrinsically disordered regions ('intrinsically disordered region': https://disprot.org/ontology).

     - "total_disordered_region_length"
       - Data type: integer
       - Description: The field 'total_disordered_region_length' indicates the total sequence length ('sequence length': http://edamontology.org/data_1249) of intrinsically disordered regions ('intrinsically disordered region': https://disprot.org/ontology).

     - "binding_pocket_count"
       - Data type: integer
       - Description: The field 'binding_pocket_count' indicates the count of binding pockets ('binding pocket': https://schlessinger-lab.github.io/pyvol/pocket_specification.html) calculated by PyVOL software ('PyVOL': https://bio.tools/PyVOL).

     - "max_binding_pocket_volume"
       - Data type: number
       - Description: The field 'max_binding_pocket_volume' indicates the maximum volume ('volume': http://purl.obolibrary.org/obo/PATO_0000918) of binding pockets ('binding pocket': https://schlessinger-lab.github.io/pyvol/index.html) calculated by PyVOL software ('PyVOL': https://bio.tools/PyVOL). Unit: cubic angstroms (Å^3) ('cubic angstrom': http://qudt.org/vocab/unit/ANGSTROM3).

     - "total_binding_pocket_volume"
       - Data type: number
       - Description: The field 'total_binding_pocket_volume' indicates the total volume ('volume': http://purl.obolibrary.org/obo/PATO_0000918) of binding pockets ('binding pocket': https://schlessinger-lab.github.io/pyvol/index.html) calculated by PyVOL software ('PyVOL': https://bio.tools/PyVOL). Unit: cubic angstroms (Å^3) ('cubic angstrom': http://qudt.org/vocab/unit/ANGSTROM3).

     - "total_potential_energy"
       - Data type: number
       - Description: The field 'total_potential_energy' indicates the total potential energy ('potential energy': https://goldbook.iupac.org/terms/view/P04778) calculated from the protein structure ('protein structure': http://edamontology.org/data_1537). Unit: kilojoules per mole (kJ/mol) ('kilojoule per mole': http://qudt.org/vocab/unit/KiloJ-PER-MOL).

     - "harmonic_bond_potential_energy"
       - Data type: number
       - Description: The field 'harmonic_bond_potential_energy' indicates the potential energy ('potential energy': https://goldbook.iupac.org/terms/view/P04778) contributed by the harmonic bond force term ('harmonic bond force term': https://docs.openmm.org/latest/userguide/theory/02_standard_forces.html#harmonicbondforce). Unit: kilojoules per mole (kJ/mol) ('kilojoule per mole': http://qudt.org/vocab/unit/KiloJ-PER-MOL).

     - "harmonic_angle_potential_energy"
       - Data type: number
       - Description: The field 'harmonic_angle_potential_energy' indicates the potential energy ('potential energy': https://goldbook.iupac.org/terms/view/P04778) contributed by the harmonic angle force term ('harmonic angle force term': https://docs.openmm.org/latest/userguide/theory/02_standard_forces.html#harmonicangleforce). Unit: kilojoules per mole (kJ/mol) ('kilojoule per mole': http://qudt.org/vocab/unit/KiloJ-PER-MOL).

     - "custom_bond_potential_energy"
       - Data type: number
       - Description: The field 'custom_bond_potential_energy' indicates the potential energy ('potential energy': https://goldbook.iupac.org/terms/view/P04778) contributed by the custom bond force term ('custom bond force term': https://docs.openmm.org/latest/userguide/theory/03_custom_forces.html#custombondforce). Unit: kilojoules per mole (kJ/mol) ('kilojoule per mole': http://qudt.org/vocab/unit/KiloJ-PER-MOL).

     - "custom_torsion_potential_energy"
       - Data type: number
       - Description: The field 'custom_torsion_potential_energy' indicates the potential energy ('potential energy': https://goldbook.iupac.org/terms/view/P04778) contributed by the custom torsion force term ('custom torsion force term': https://docs.openmm.org/latest/userguide/theory/03_custom_forces.html#customtorsionforce). Unit: kilojoules per mole (kJ/mol) ('kilojoule per mole': http://qudt.org/vocab/unit/KiloJ-PER-MOL).

     - "custom_nonbonded_potential_energy"
       - Data type: number
       - Description: The field 'custom_nonbonded_potential_energy' indicates the potential energy ('potential energy': https://goldbook.iupac.org/terms/view/P04778) contributed by the custom nonbonded force term ('custom nonbonded force term': https://docs.openmm.org/latest/userguide/theory/03_custom_forces.html#customnonbondedforce). Unit: kilojoules per mole (kJ/mol) ('kilojoule per mole': http://qudt.org/vocab/unit/KiloJ-PER-MOL).

     - "nonbonded_potential_energy"
       - Data type: number
       - Description: The field 'nonbonded_potential_energy' indicates the potential energy ('potential energy': https://goldbook.iupac.org/terms/view/P04778) contributed by the nonbonded force term ('nonbonded force term': https://docs.openmm.org/latest/userguide/theory/02_standard_forces.html#nonbondedforce). Unit: kilojoules per mole (kJ/mol) ('kilojoule per mole': http://qudt.org/vocab/unit/KiloJ-PER-MOL).

     - "periodic_torsion_potential_energy"
       - Data type: number
       - Description: The field 'periodic_torsion_potential_energy' indicates the potential energy ('potential energy': https://goldbook.iupac.org/terms/view/P04778) contributed by the periodic torsion force term ('periodic torsion force term': https://docs.openmm.org/latest/userguide/theory/02_standard_forces.html#periodictorsionforce). Unit: kilojoules per mole (kJ/mol) ('kilojoule per mole': http://qudt.org/vocab/unit/KiloJ-PER-MOL).

     - "cmap_torsion_potential_energy"
       - Data type: number
       - Description: The field 'cmap_torsion_potential_energy' indicates the potential energy ('potential energy': https://goldbook.iupac.org/terms/view/P04778) contributed by the CMAP torsion force term ('CMAP torsion force term': https://docs.openmm.org/latest/userguide/theory/02_standard_forces.html#cmaptorsionforce). Unit: kilojoules per mole (kJ/mol) ('kilojoule per mole': http://qudt.org/vocab/unit/KiloJ-PER-MOL).

     - "enzyme_substrate_binding_affinity"
       - Data type: number
       - Description: The field 'enzyme_substrate_binding_affinity' indicates the predicted binding affinity ('binding affinity': https://vina.scripps.edu/manual/#output) calculated by AutoDock Vina software ('AutoDock Vina': https://bio.tools/autodock_vina) from docking ('docking': https://goldbook.iupac.org/terms/view/11437) of the enzyme-substrate complex ('enzyme': https://purl.dsmz.de/schema/Enzyme; 'substrate': https://purl.dsmz.de/schema/Substrate; 'complex': https://goldbook.iupac.org/terms/view/C01203). Unit: kilocalories per mole (kcal/mol) ('kilocalorie': http://qudt.org/vocab/unit/KiloCAL; 'mole': http://qudt.org/vocab/unit/MOL).

     - "hydrogen_bond_count"
       - Data type: integer
       - Description: The field 'hydrogen_bond_count' indicates the count of hydrogen bonds ('hydrogen bond': https://goldbook.iupac.org/terms/view/H02899).

     - "ionic_bond_count"
       - Data type: integer
       - Description: The field 'ionic_bond_count' indicates the count of ionic bonds ('ionic bond': https://goldbook.iupac.org/terms/view/IT07058).

     - "van_der_waals_contact_count"
       - Data type: integer
       - Description: The field 'van_der_waals_contact_count' indicates the count of van der Waals contacts ('van der Waals forces': https://goldbook.iupac.org/terms/view/V06597).

     - "pi_pi_stacking_count"
       - Data type: integer
       - Description: The field 'pi_pi_stacking_count' indicates the count of pi-pi stacking interactions ('pi-pi stacking': https://goldbook.iupac.org/terms/view/13861).

     - "pi_cation_interaction_count"
       - Data type: integer
       - Description: The field 'pi_cation_interaction_count' indicates the count of pi-cation interactions ('cation-pi interaction': https://goldbook.iupac.org/terms/view/08154).

     - "disulfide_bond_count"
       - Data type: integer
       - Description: The field 'disulfide_bond_count' indicates the count of disulfide bonds ('disulfide bond': https://www.uniprot.org/help/disulfid).

   - "integrated_graph"
     - Data type: array
     - Description: The field 'integrated_graph' indicates the integrated graph ('graph': https://neo4j.com/docs/getting-started/appendix/graphdb-concepts/) containing molecular interactions ('molecular interaction': https://bioportal.bioontology.org/ontologies/MI) between source nodes and target nodes, and isolated nodes ('isolated node': https://mathworld.wolfram.com/IsolatedPoint.html) integrated from EnzyWizard reports ('report': http://purl.obolibrary.org/obo/IAO_0000088).

     Each item in "integrated_graph" must be one of the following forms:

     1. A molecular interaction entry object containing:

       - "molecular_interaction"
         - Data type: object
         - Description: The field 'molecular_interaction' indicates a molecular interaction ('molecular interaction': https://bioportal.bioontology.org/ontologies/MI) between the source node and the target node in the integrated graph ('graph': https://neo4j.com/docs/getting-started/appendix/graphdb-concepts/).

         The "molecular_interaction" object contains:

         - "molecular_interaction_type"
           - Data type: string
           - Allowed values: HBOND, IONIC, VDW, PIPISTACK, PICATION, SSBOND.
           - Description: The field 'molecular_interaction_type' indicates the type ('interaction type': http://purl.obolibrary.org/obo/MI_0190) of molecular interaction ('molecular interaction': https://bioportal.bioontology.org/ontologies/MI), using RING interaction codes ('RING interaction type': https://ring.biocomputingup.it/help/interactions): hydrogen bond ('hydrogen bond': https://goldbook.iupac.org/terms/view/H02899; value: HBOND), ionic bond ('ionic bond': https://goldbook.iupac.org/terms/view/IT07058; value: IONIC), van der Waals contact ('van der Waals forces': https://goldbook.iupac.org/terms/view/V06597; value: VDW), pi-pi stacking ('pi-pi stacking': https://goldbook.iupac.org/terms/view/13861; value: PIPISTACK), pi-cation interaction ('cation-pi interaction': https://goldbook.iupac.org/terms/view/08154; value: PICATION), and disulfide bond ('disulfide bond': https://www.uniprot.org/help/disulfid; value: SSBOND).

         - "molecular_interaction_one_hot_encoding"
           - Data type: array
           - Length: 6
           - Item data type: integer
           - Item enum: [0, 1]
           - Description: The field 'molecular_interaction_one_hot_encoding' indicates the one-hot encoding ('one-hot encoding': https://developers.google.com/machine-learning/glossary#one-hot_encoding) of the molecular interaction type ('interaction type': http://purl.obolibrary.org/obo/MI_0190).

         - "interaction_count"
           - Data type: integer
           - Description: The field 'interaction_count' indicates the count of molecular interactions ('molecular interaction': https://bioportal.bioontology.org/ontologies/MI) between the source node and the target node.

       - "source_node"
         - Data type: oneOf

       - "target_node"
         - Data type: oneOf

       The "source_node" and "target_node" objects are either residue nodes or substrate nodes.

       Residue node contains:

       - "node_index"
         - Data type: integer
         - Description: The field 'node_index' indicates the index ('index': http://purl.obolibrary.org/obo/NCIT_C25390) of the node ('node': https://neo4j.com/docs/getting-started/appendix/graphdb-concepts/#graphdb-node) in the integrated graph ('graph': https://neo4j.com/docs/getting-started/appendix/graphdb-concepts/).

       - "node_type"
         - Data type: string
         - Expected value: "residue"
         - Description: The field 'node_type' indicates the type of node ('node': https://neo4j.com/docs/getting-started/appendix/graphdb-concepts/#graphdb-node), with value 'residue' indicating a residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782).

       - "node_type_one_hot_encoding"
         - Data type: array
         - Length: 2
         - Item data type: integer
         - Item enum: [0, 1]
         - Description: The field 'node_type_one_hot_encoding' indicates the one-hot encoding ('one-hot encoding': https://developers.google.com/machine-learning/glossary#one-hot_encoding) of the node type.

       - "residue_index"
         - Data type: integer
         - Description: The field 'residue_index' indicates the index ('index': http://purl.obolibrary.org/obo/NCIT_C25390) of the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782).

       - "residue_name"
         - Data type: string
         - Pattern: ^[ACDEFGHIKLMNPQRSTVWY]$.
         - Description: The field 'residue_name' indicates the name of the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782), using one-letter code ('one-letter code': https://iupac.qmul.ac.uk/AminoAcid/A2021.html) to represent the amino acid residue.

       - "residue_name_one_hot_encoding"
         - Data type: array
         - Length: 20
         - Item data type: integer
         - Item enum: [0, 1]
         - Description: The field 'residue_name_one_hot_encoding' indicates the one-hot encoding ('one-hot encoding': https://developers.google.com/machine-learning/glossary#one-hot_encoding) of the residue name ('residue': http://purl.obolibrary.org/obo/GENO_0000782).

       - "residue_alpha_carbon_coordinate"
         - Data type: array
         - Length: 3
         - Item data type: number
         - Description: The field 'residue_alpha_carbon_coordinate' indicates the three-dimensional coordinate ('coordinate': http://purl.obolibrary.org/obo/NCIT_C44477) of the alpha carbon atom ('alpha carbon': https://www.rcsb.org/docs/general-help/glossary; 'atom': http://purl.obolibrary.org/obo/CHMO_0001075) in the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782). Unit: angstroms (Å) ('angstrom': http://qudt.org/vocab/unit/ANGSTROM).

       - "residue_chemical_classification"
         - Data type: string
         - Pattern: ^(uncharged_polar|positively_charged|negatively_charged|hydrophobic|aromatic|aliphatic|heterocyclic|sulfur_containing)(/(uncharged_polar|positively_charged|negatively_charged|hydrophobic|aromatic|aliphatic|heterocyclic|sulfur_containing))*$.
         - Description: The field 'residue_chemical_classification' indicates the chemical classification ('classification': http://purl.obolibrary.org/obo/NCIT_C25161) of the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782).

       - "residue_chemical_classification_multi_hot_encoding"
         - Data type: array
         - Length: 8
         - Item data type: integer
         - Item enum: [0, 1]
         - Description: The field 'residue_chemical_classification_multi_hot_encoding' indicates the multi-hot encoding ('multi-hot encoding': https://developers.google.com/machine-learning/crash-course/categorical-data/one-hot-encoding) of the chemical classification ('classification': http://purl.obolibrary.org/obo/NCIT_C25161) of the amino acid residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782).

       - "residue_secondary_structure"
         - Data type: string
         - Allowed values: -, H, B, E, G, I, T, S.
         - Description: The field 'residue_secondary_structure' indicates the secondary structure ('secondary structure': http://edamontology.org/operation_1847) assigned to the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782), using DSSP secondary-structure codes ('DSSP': https://manual.gromacs.org/current/onlinehelp/gmx-dssp.html).

       - "residue_secondary_structure_one_hot_encoding"
         - Data type: array
         - Length: 8
         - Item data type: integer
         - Item enum: [0, 1]
         - Description: The field 'residue_secondary_structure_one_hot_encoding' indicates the one-hot encoding ('one-hot encoding': https://developers.google.com/machine-learning/glossary#one-hot_encoding) of the residue secondary structure ('secondary structure': http://edamontology.org/operation_1847).

       - "residue_relative_solvent_accessibility"
         - Data type: number
         - Description: The field 'residue_relative_solvent_accessibility' indicates the relative solvent accessibility ('solvent accessibility': http://edamontology.org/data_1542) of the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782). Unit: dimensionless ('dimensionless': http://qudt.org/vocab/unit/UNITLESS).

       - "residue_backbone_phi_angle"
         - Data type: number
         - Description: The field 'residue_backbone_phi_angle' indicates the backbone phi torsion angle ('torsion angle': https://goldbook.iupac.org/terms/view/T06406) of the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782) in the protein backbone ('protein backbone': http://edamontology.org/operation_1825). Unit: degrees (°) ('degree': http://qudt.org/vocab/unit/DEG).

       - "residue_backbone_psi_angle"
         - Data type: number
         - Description: The field 'residue_backbone_psi_angle' indicates the backbone psi torsion angle ('torsion angle': https://goldbook.iupac.org/terms/view/T06406) of the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782) in the protein backbone ('protein backbone': http://edamontology.org/operation_1825). Unit: degrees (°) ('degree': http://qudt.org/vocab/unit/DEG).

       - "residue_net_charge"
         - Data type: number
         - Description: The field 'residue_net_charge' indicates the net electric charge ('net electric charge': https://goldbook.iupac.org/terms/view/N04111) of the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782). Unit: dimensionless ('dimensionless': http://qudt.org/vocab/unit/UNITLESS).

       - "residue_pka"
         - Data type: number
         - Description: The field 'residue_pka' indicates the pKa value ('pKa': https://goldbook.iupac.org/terms/view/15441) of the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782). Unit: dimensionless ('dimensionless': http://qudt.org/vocab/unit/UNITLESS).

       - "residue_volume"
         - Data type: number
         - Description: The field 'residue_volume' indicates the volume ('volume': http://purl.obolibrary.org/obo/PATO_0000918) of the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782). Unit: cubic angstroms (Å^3) ('cubic angstrom': http://qudt.org/vocab/unit/ANGSTROM3).

       - "residue_hydrophobicity"
         - Data type: number
         - Description: The field 'residue_hydrophobicity' indicates the hydrophobicity ('hydrophobicity': https://goldbook.iupac.org/terms/view/HT06964) of the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782). Unit: dimensionless ('dimensionless': http://qudt.org/vocab/unit/UNITLESS).

       - "residue_molecular_weight"
         - Data type: number
         - Description: The field 'residue_molecular_weight' indicates the molecular weight ('molecular weight': https://goldbook.iupac.org/terms/view/R05271) of the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782). Unit: daltons (Da) ('dalton': http://qudt.org/vocab/unit/DA).

       - "residue_isoelectric_point"
         - Data type: number
         - Description: The field 'residue_isoelectric_point' indicates the isoelectric point ('isoelectric point': https://goldbook.iupac.org/terms/view/I03275) of the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782). Unit: dimensionless ('dimensionless': http://qudt.org/vocab/unit/UNITLESS).

       - "residue_root_mean_square_fluctuation"
         - Data type: number
         - Description: The field 'residue_root_mean_square_fluctuation' indicates the root mean square fluctuation ('root mean square fluctuation': https://manual.gromacs.org/current/onlinehelp/gmx-rmsf.html) of the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782). Unit: angstroms (Å) ('angstrom': http://qudt.org/vocab/unit/ANGSTROM).

       - "residue_sequence_conservation_score"
         - Data type: number
         - Description: The field 'residue_sequence_conservation_score' indicates the sequence conservation score based on normalized Shannon information content ('Shannon entropy': https://mathworld.wolfram.com/Entropy.html; 'information content': https://www.ebsco.com/research-starters/library-and-information-science/information-content) of the residue position, calculated based on the 'normalized_emission_probability'. Unit: dimensionless ('dimensionless': http://qudt.org/vocab/unit/UNITLESS).

       - "residue_embedding"
         - Data type: array
         - Item data type: number
         - Description: The field 'residue_embedding' indicates the embedding ('embedding': https://developers.google.com/machine-learning/crash-course/embeddings) generated by the ESM-2 protein language model ('ESM-2': https://docs.nvidia.com/bionemo-framework/2.0/models/esm2/; 'protein language model': https://synbiointel.com/glossary/protein-language-model/) for the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782), represented as a numerical vector ('numerical vector': https://mathworld.wolfram.com/Vector.html).

       - "is_in_hydrophobic_cluster"
         - Data type: boolean
         - Description: The field 'is_in_hydrophobic_cluster' indicates whether the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782) is included in a hydrophobic cluster ('hydrophobic cluster': https://proteintools.uni-bayreuth.de/clusters/).

       - "is_in_disordered_region"
         - Data type: boolean
         - Description: The field 'is_in_disordered_region' indicates whether the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782) is included in an intrinsically disordered region ('intrinsically disordered region': https://disprot.org/ontology).

       - "is_in_binding_pocket"
         - Data type: boolean
         - Description: The field 'is_in_binding_pocket' indicates whether the residue ('residue': http://purl.obolibrary.org/obo/GENO_0000782) is included in a binding pocket ('binding pocket': https://schlessinger-lab.github.io/pyvol/index.html).

       Substrate node contains:

       - "node_index"
         - Data type: integer
         - Description: The field 'node_index' indicates the index ('index': http://purl.obolibrary.org/obo/NCIT_C25390) of the node ('node': https://neo4j.com/docs/getting-started/appendix/graphdb-concepts/#graphdb-node) in the integrated graph ('graph': https://neo4j.com/docs/getting-started/appendix/graphdb-concepts/).

       - "node_type"
         - Data type: string
         - Expected value: "substrate"
         - Description: The field 'node_type' indicates the type of node ('node': https://neo4j.com/docs/getting-started/appendix/graphdb-concepts/#graphdb-node), with value 'substrate' indicating a substrate ('substrate': https://purl.dsmz.de/schema/Substrate).

       - "node_type_one_hot_encoding"
         - Data type: array
         - Length: 2
         - Item data type: integer
         - Item enum: [0, 1]
         - Description: The field 'node_type_one_hot_encoding' indicates the one-hot encoding ('one-hot encoding': https://developers.google.com/machine-learning/glossary#one-hot_encoding) of the node type.

       - "substrate_index"
         - Data type: integer
         - Description: The field 'substrate_index' indicates the index ('index': http://purl.obolibrary.org/obo/NCIT_C25390) of the substrate ('substrate': https://purl.dsmz.de/schema/Substrate).

       - "substrate_name"
         - Data type: string
         - Description: The field 'substrate_name' indicates the name of the substrate ('substrate': https://purl.dsmz.de/schema/Substrate).

       - "substrate_smiles"
         - Data type: string
         - Description: The field 'substrate_smiles' indicates the SMILES representation ('SMILES': https://opensmiles.org/opensmiles.html) of the substrate ('substrate': https://purl.dsmz.de/schema/Substrate).

       - "substrate_atom_count"
         - Data type: integer
         - Description: The field 'substrate_atom_count' indicates the count of atoms ('atom': https://goldbook.iupac.org/terms/view/A00493) in the substrate ('substrate': https://purl.dsmz.de/schema/Substrate).

       - "substrate_molecular_weight"
         - Data type: number
         - Description: The field 'substrate_molecular_weight' indicates the molecular weight ('molecular weight': https://goldbook.iupac.org/terms/view/R05271) of the substrate ('substrate': https://purl.dsmz.de/schema/Substrate). Unit: daltons (Da) ('dalton': http://qudt.org/vocab/unit/DA).

       - "substrate_logp"
         - Data type: number
         - Description: The field 'substrate_logp' indicates the calculated logP value ('LogP': https://doktormike.gitlab.io/posts/navigating-logp-logd-pka-and-logs-a-physicists-guide/) of the substrate ('substrate': https://purl.dsmz.de/schema/Substrate). Unit: dimensionless ('dimensionless': http://qudt.org/vocab/unit/UNITLESS).

       - "substrate_tpsa"
         - Data type: number
         - Description: The field 'substrate_tpsa' indicates the topological polar surface area ('TPSA': https://www.rdkit.org/docs/GettingStartedInPython.html#descriptor-calculation) of the substrate ('substrate': https://purl.dsmz.de/schema/Substrate) calculated by RDKit software ('RDKit': https://www.rdkit.org/docs/index.html). Unit: square angstroms (Å^2) ('angstrom': http://qudt.org/vocab/unit/ANGSTROM).

       - "substrate_heavy_atom_count"
         - Data type: integer
         - Description: The field 'substrate_heavy_atom_count' indicates the count of heavy atoms ('atom': https://goldbook.iupac.org/terms/view/A00493) in the substrate ('substrate': https://purl.dsmz.de/schema/Substrate).

       - "substrate_hbond_donor_count"
         - Data type: integer
         - Description: The field 'substrate_hbond_donor_count' indicates the count of hydrogen bond donors ('hydrogen bond': https://goldbook.iupac.org/terms/view/H02899) in the substrate ('substrate': https://purl.dsmz.de/schema/Substrate) calculated by RDKit software ('RDKit': https://www.rdkit.org/docs/index.html).

       - "substrate_hbond_acceptor_count"
         - Data type: integer
         - Description: The field 'substrate_hbond_acceptor_count' indicates the count of hydrogen bond acceptors ('hydrogen bond': https://goldbook.iupac.org/terms/view/H02899) in the substrate ('substrate': https://purl.dsmz.de/schema/Substrate) calculated by RDKit software ('RDKit': https://www.rdkit.org/docs/index.html).

       - "substrate_rotatable_bond_count"
         - Data type: integer
         - Description: The field 'substrate_rotatable_bond_count' indicates the count of rotatable bonds ('bond': https://goldbook.iupac.org/terms/view/B00701) in the substrate ('substrate': https://purl.dsmz.de/schema/Substrate) calculated by RDKit software ('RDKit': https://www.rdkit.org/docs/index.html).

       - "substrate_molar_refractivity"
         - Data type: number
         - Description: The field 'substrate_molar_refractivity' indicates the molar refractivity ('molar refractivity': https://old.iupac.org/reports/1997/6905vandewaterbeemd/glossary.html) of the substrate ('substrate': https://purl.dsmz.de/schema/Substrate) calculated by RDKit software ('RDKit': https://www.rdkit.org/docs/index.html). Unit: cubic centimeters per mole (cm^3/mol) ('cubic centimeter': http://qudt.org/vocab/unit/CentiM3; 'mole': http://qudt.org/vocab/unit/MOL).

       - "substrate_structure_energy"
         - Data type: number
         - Description: The field 'substrate_structure_energy' indicates the energy ('energy': http://purl.obolibrary.org/obo/PATO_0001021) of a possible molecular structure ('molecular structure': http://edamontology.org/data_0883) generated for the substrate ('substrate': https://purl.dsmz.de/schema/Substrate). Unit: kilocalories per mole (kcal/mol) ('kilocalorie': http://qudt.org/vocab/unit/KiloCAL; 'mole': http://qudt.org/vocab/unit/MOL).

       - "substrate_structure_max_3d_diameter"
         - Data type: number
         - Description: The field 'substrate_structure_max_3d_diameter' indicates the maximum three-dimensional diameter ('diameter': http://purl.obolibrary.org/obo/PATO_0001334) of a possible molecular structure ('molecular structure': http://edamontology.org/data_0883) generated for the substrate ('substrate': https://purl.dsmz.de/schema/Substrate). Unit: angstroms (Å) ('angstrom': http://qudt.org/vocab/unit/ANGSTROM).

       - "substrate_structure_mean_pairwise_atom_distance"
         - Data type: number
         - Description: The field 'substrate_structure_mean_pairwise_atom_distance' indicates the mean pairwise atom distance ('distance': http://purl.obolibrary.org/obo/PATO_0000040) of a possible molecular structure ('molecular structure': http://edamontology.org/data_0883) generated for the substrate ('substrate': https://purl.dsmz.de/schema/Substrate). Unit: angstroms (Å) ('angstrom': http://qudt.org/vocab/unit/ANGSTROM).

       - "substrate_structure_std_pairwise_atom_distance"
         - Data type: number
         - Description: The field 'substrate_structure_std_pairwise_atom_distance' indicates the standard deviation ('standard deviation': http://purl.obolibrary.org/obo/STATO_0000237) of pairwise atom distances ('distance': http://purl.obolibrary.org/obo/PATO_0000040) of a possible molecular structure ('molecular structure': http://edamontology.org/data_0883) generated for the substrate ('substrate': https://purl.dsmz.de/schema/Substrate). Unit: angstroms (Å) ('angstrom': http://qudt.org/vocab/unit/ANGSTROM).

       - "substrate_structure_asphericity"
         - Data type: number
         - Description: The field 'substrate_structure_asphericity' indicates the asphericity ('asphericity': https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html) of a possible molecular structure ('molecular structure': http://edamontology.org/data_0883) generated for the substrate ('substrate': https://purl.dsmz.de/schema/Substrate) calculated by RDKit software ('RDKit': https://www.rdkit.org/docs/index.html). Unit: dimensionless ('dimensionless': http://qudt.org/vocab/unit/UNITLESS).

       - "substrate_structure_spherocity"
         - Data type: number
         - Description: The field 'substrate_structure_spherocity' indicates the spherocity index ('spherocity index': https://www.rdkit.org/docs/source/rdkit.Chem.rdMolDescriptors.html) of a possible molecular structure ('molecular structure': http://edamontology.org/data_0883) generated for the substrate ('substrate': https://purl.dsmz.de/schema/Substrate) calculated by RDKit software ('RDKit': https://www.rdkit.org/docs/index.html). Unit: dimensionless ('dimensionless': http://qudt.org/vocab/unit/UNITLESS).

       - "substrate_structure_principal_moment_ratio"
         - Data type: number
         - Description: The field 'substrate_structure_principal_moment_ratio' indicates the ratio of the largest to the smallest principal moments of inertia ('moment of inertia': https://goldbook.iupac.org/terms/view/M03954) of a possible molecular structure ('molecular structure': http://edamontology.org/data_0883) generated for the substrate ('substrate': https://purl.dsmz.de/schema/Substrate). Unit: dimensionless ('dimensionless': http://qudt.org/vocab/unit/UNITLESS).

       - "substrate_structure_radius_of_gyration"
         - Data type: number
         - Description: The field 'substrate_structure_radius_of_gyration' indicates the radius of gyration ('radius of gyration': https://goldbook.iupac.org/terms/view/R05121) of a possible molecular structure ('molecular structure': http://edamontology.org/data_0883) generated for the substrate ('substrate': https://purl.dsmz.de/schema/Substrate) calculated by RDKit software ('RDKit': https://www.rdkit.org/docs/index.html). Unit: angstroms (Å) ('angstrom': http://qudt.org/vocab/unit/ANGSTROM).

       - "docked_substrate_center_coordinate"
         - Data type: array
         - Length: 3
         - Item data type: number
         - Description: The field 'docked_substrate_center_coordinate' indicates the center coordinate ('coordinate': https://mathworld.wolfram.com/Coordinates.html) of the docked substrate ('substrate': https://purl.dsmz.de/schema/Substrate) in the enzyme-substrate complex ('enzyme': https://purl.dsmz.de/schema/Enzyme; 'substrate': https://purl.dsmz.de/schema/Substrate; 'complex': https://goldbook.iupac.org/terms/view/C01203). Unit: angstroms (Å) ('angstrom': http://qudt.org/vocab/unit/ANGSTROM).

       - "substrate_fingerprint_encoding"
         - Data type: array
         - Item data type: integer
         - Item enum: [0, 1]
         - Description: The field 'substrate_fingerprint_encoding' indicates the molecular fingerprint encoding ('molecular fingerprint': https://www.rdkit.org/docs/GettingStartedInPython.html#fingerprinting-and-molecular-similarity) of the substrate ('substrate': https://purl.dsmz.de/schema/Substrate) calculated by RDKit software ('RDKit': https://www.rdkit.org/docs/index.html).

     2. An isolated node entry object containing:

       - "isolated_node"
         - Data type: oneOf

       The "isolated_node" object is either a residue node or a substrate node, with the same node fields described above.

The node-only JSON file contains an array of residue node or substrate node objects with the same node fields described above.

The edge-only JSON file contains the following fields for each molecular interaction entry:

   - "molecular_interaction"
     - Data type: object
     - Description: The field 'molecular_interaction' indicates a molecular interaction ('molecular interaction': https://bioportal.bioontology.org/ontologies/MI) between the source node and the target node in the integrated graph ('graph': https://neo4j.com/docs/getting-started/appendix/graphdb-concepts/).

     The "molecular_interaction" object contains:

     - "molecular_interaction_type"
       - Data type: string
       - Allowed values: HBOND, IONIC, VDW, PIPISTACK, PICATION, SSBOND.
       - Description: The field 'molecular_interaction_type' indicates the type ('interaction type': http://purl.obolibrary.org/obo/MI_0190) of molecular interaction ('molecular interaction': https://bioportal.bioontology.org/ontologies/MI), using RING interaction codes ('RING interaction type': https://ring.biocomputingup.it/help/interactions): hydrogen bond ('hydrogen bond': https://goldbook.iupac.org/terms/view/H02899; value: HBOND), ionic bond ('ionic bond': https://goldbook.iupac.org/terms/view/IT07058; value: IONIC), van der Waals contact ('van der Waals forces': https://goldbook.iupac.org/terms/view/V06597; value: VDW), pi-pi stacking ('pi-pi stacking': https://goldbook.iupac.org/terms/view/13861; value: PIPISTACK), pi-cation interaction ('cation-pi interaction': https://goldbook.iupac.org/terms/view/08154; value: PICATION), and disulfide bond ('disulfide bond': https://www.uniprot.org/help/disulfid; value: SSBOND).

     - "molecular_interaction_one_hot_encoding"
       - Data type: array
       - Length: 6
       - Item data type: integer
       - Item enum: [0, 1]
       - Description: The field 'molecular_interaction_one_hot_encoding' indicates the one-hot encoding ('one-hot encoding': https://developers.google.com/machine-learning/glossary#one-hot_encoding) of the molecular interaction type ('interaction type': http://purl.obolibrary.org/obo/MI_0190).

     - "interaction_count"
       - Data type: integer
       - Description: The field 'interaction_count' indicates the count of molecular interactions ('molecular interaction': https://bioportal.bioontology.org/ontologies/MI) between the source node and the target node.

   - "source_node"
     - Data type: object
     - Description: The field 'source_node' indicates the source node ('node': https://neo4j.com/docs/getting-started/appendix/graphdb-concepts/#graphdb-node) in the integrated graph.

     The "source_node" object contains:

     - "node_index"
       - Data type: integer
       - Description: The field 'node_index' indicates the index ('index': http://purl.obolibrary.org/obo/NCIT_C25390) of the node ('node': https://neo4j.com/docs/getting-started/appendix/graphdb-concepts/#graphdb-node) in the integrated graph ('graph': https://neo4j.com/docs/getting-started/appendix/graphdb-concepts/).

   - "target_node"
     - Data type: object
     - Description: The field 'target_node' indicates the target node ('node': https://neo4j.com/docs/getting-started/appendix/graphdb-concepts/#graphdb-node) in the integrated graph.

     The "target_node" object contains:

     - "node_index"
       - Data type: integer
       - Description: The field 'node_index' indicates the index ('index': http://purl.obolibrary.org/obo/NCIT_C25390) of the node ('node': https://neo4j.com/docs/getting-started/appendix/graphdb-concepts/#graphdb-node) in the integrated graph ('graph': https://neo4j.com/docs/getting-started/appendix/graphdb-concepts/).

# Process:

This command processes the input cleaned structure and optional MSA as follows:

1. Validate input files
   - Check that cleaned_input_path exists.
   - If input_msa is provided, check that it exists.
   - Create output_dir if needed.

2. Resolve names
   - Extract protein_name from the cleaned structure filename.
   - If input_msa is provided, extract msa_name from the MSA filename.
   - Validate filename length.

3. Validate parameters
   - Check parameter ranges for hydrophobic cluster detection.
   - Check energy minimization parameters.
   - Check flexibility and disordered region parameters.
   - Check binding pocket detection parameters.
   - Check substrate generation parameters.
   - Check docking parameters.
   - Check interaction detection parameters.

4. Prepare output mode
   - If --save_extra_outputs is enabled, run directly in output_dir.
   - If disabled, run in a temporary directory and only copy final outputs.

5. Load cleaned structure
   - Read the cleaned CIF or PDB file.
   - Validate that it is a valid cleaned protein structure.
   - Check that hydrogen atoms are present.

6. Build identity clean report
   - Treat the input structure as already cleaned.
   - Build an identity residue mapping between old residues and new residues.
   - Generate an enzywizard_clean-style report for downstream integration.

7. Prepare OpenMM and sequence objects
   - Convert the cleaned structure into an OpenMM-compatible structure.
   - Build an OpenMM Modeller object.
   - Extract the cleaned protein sequence.

8. Run residue property analysis
   - Load DSSP information.
   - Calculate residue-level properties.
   - Generate the enzywizard_aaprops report.

9. Run hydrophobic cluster analysis
   - Detect hydrophobic clusters.
   - Generate the enzywizard_hydrocluster report.

10. Run energy analysis
   - Optionally minimize the protein structure.
   - Calculate energy terms using the selected force field.
   - Generate the enzywizard_energy report.

11. Run flexibility analysis
   - Calculate residue-level RMSF using ANM or GNM.
   - Generate the enzywizard_flexibility report.

12. Run disordered region analysis
   - Calculate FoldIndex-like disordered region scores.
   - Detect disordered regions.
   - Generate the enzywizard_disorder report.

13. Run conservation analysis when input_msa is provided
   - Load and validate the input MSA.
   - Automatically decompress the input MSA when the file is in .fasta.gz format.
   - Clean the MSA into Stockholm format.
   - Build an HMM profile.
   - Calculate residue-level conservation scores.
   - Generate the enzywizard_conservation report.
   - If input_msa is omitted, skip this step and omit conservation scores from
     the final integrated residue nodes.

14. Run embedding analysis
   - Generate residue-level ESM2 embeddings.
   - Generate the enzywizard_embedding report.

15. Run binding pocket analysis
   - Detect binding pockets using PyVOL.
   - Generate the enzywizard_pocket report.

16. Optionally run substrate analysis
   - Parse substrate names or SMILES strings.
   - Retrieve or complete SMILES information.
   - Generate substrate fingerprints and 3D conformers.
   - Save substrate structure files.
   - Generate the enzywizard_substrate report.
   - If substrate parsing, SMILES completion, conformer generation, structure saving,
     or report generation fails, log a warning and continue with the protein-only route.

17. Optionally run docking analysis
   - If --dock_catalytic_residue is provided, use the CA coordinate of that cleaned protein residue as the docking box center.
   - If --dock_catalytic_site_coord is provided, use that coordinate as the docking box center.
   - In either manual docking box mode, use --dock_box_size as the docking box size and skip docking-specific PyVOL pocket detection and the global docking box fallback.
   - If no manual docking box parameter is provided, dock generated substrate structures into automatically generated pocket and global fallback docking boxes.
   - Save docking results.
   - Generate the enzywizard_dock report.
   - Load docked substrate structures for molecular interaction analysis.
   - If docking, dock report generation, or docked substrate loading fails,
     log a warning and continue with the protein-only route.

18. Run molecular interaction analysis
   - If substrate input is provided and docking completed, filter valid docked substrates.
   - Calculate intra-enzyme molecular interactions.
   - Calculate enzyme-substrate molecular interactions when valid docked substrates exist.
   - Summarize molecular interaction counts.
   - Generate the enzywizard_interaction report.

19. Run graph integration
   - Collect all generated reports into report_dict.
   - Use strict integration when substrate input is provided, an MSA is provided,
     and substrate/docking workflows complete successfully.
   - Use non-strict integration when no substrate input is provided, no MSA is
     provided, or the workflow falls back to the protein-only route.
   - Merge residue-level, substrate-level, and molecular interaction-level information
     into a unified integrated graph.

20. Save integrated outputs
   - Write integrate_report_{protein_name}.json.
   - Split integrated_graph into node and edge lists.
   - Write integrate_nodes_{protein_name}.json.
   - Write integrate_edges_{protein_name}.json.

21. Finalize outputs
   - If --save_extra_outputs is disabled, copy only the final integrated JSON
     outputs and log.txt from the temporary directory to output_dir.
   - Finish the batch workflow.


# common errors and solutions:

- "the following arguments are required: -i/--cleaned_input_path, -o/--output_dir"
  - Cause: The required structure path or output directory was not provided.
  - Solution: Provide both `-i path/to/input.cif` and `-o path/to/output_dir/`. MSA and substrate input are optional.

- "Input cleaned protein file not found"
  - Cause: The structure file passed to `-i` or `--cleaned_input_path` does not exist, or the path points to the wrong location.
  - Solution: Check the structure path and use a cleaned `.cif` or `.pdb` file.

- "Filename too long"
  - Cause: The input structure or MSA file name without extension is longer than the supported filename limit.
  - Solution: Rename the input file to a shorter name and run the command again.

- "Unsupported format"
  - Cause: The input structure extension is not `.cif` or `.pdb`, or an optional MSA file uses an unsupported format.
  - Solution: Use a supported cleaned structure format and, when MSA is provided, use Stockholm, aligned FASTA, gzip-compressed aligned FASTA, or A3M.

- "Exception in loading structure for"
  - Cause: Biopython could not parse the structure file as a usable protein structure, or the file is empty, corrupted, or inconsistent with its extension.
  - Solution: Check that the file is valid and non-empty, then rerun with a cleaned CIF or PDB structure.

- "Input structure is not a valid cleaned structure."
  - Cause: The input is not a valid EnzyWizard-cleaned single-chain protein structure. Common causes include multiple chains, non-chain-A input, heterogens, insertion codes, non-standard residues, missing atoms, unexpected atoms, invalid occupancies, or non-continuous numbering.
  - Solution: Review the specific validation error above this summary in `log.txt`, run `enzywizard-clean`, and use its cleaned output.

- "Input cleaned structure does not contain hydrogen atoms."
  - Cause: Hydrogen atoms are missing from the cleaned protein structure.
  - Solution: Regenerate the cleaned structure with hydrogen addition enabled, then rerun batch.

- "Exception in loading dssp"
  - Cause: DSSP failed to run or failed to parse the cleaned structure.
  - Solution: Confirm that DSSP or `mkdssp` is installed and available, and check that the cleaned input structure is valid.

- "Failed to load OpenMM force field"
  - Cause: OpenMM could not load the force field used for energy calculation.
  - Solution: Check the OpenMM installation and force-field availability in the running environment.

- "Failed to create OpenMM system"
  - Cause: OpenMM could not build a molecular system from the cleaned topology and coordinates, often because the structure contains residues, atoms, or connectivity that the force field cannot parameterize.
  - Solution: Rerun structure cleaning, check the cleaned structure, and review the detailed OpenMM error in `log.txt`.

- "Failed to calculate RMSF by ProDy"
  - Cause: ProDy failed while building the elastic network or solving normal modes, often because of unsuitable coordinates, an extreme cutoff value, or an environment issue.
  - Solution: Review the detailed ProDy error in `log.txt`, check the cleaned structure, and try a standard cutoff such as `15.0`.

- "Input MSA file not found"
  - Cause: `-m` or `--input_msa` was provided, but the MSA file does not exist.
  - Solution: Check the MSA path, or omit `-m` to run without conservation analysis.

- "The first Stockholm MSA sequence does not match query_sequence after gap removal."
  - Cause: The first sequence in the Stockholm MSA is not the same as the cleaned input structure sequence after gaps are removed.
  - Solution: Put the cleaned protein sequence as the first MSA record and make sure it matches the structure passed to `-i`.

- "The first aligned FASTA MSA sequence does not match query_sequence after gap removal."
  - Cause: The first sequence in the aligned FASTA MSA is not the same as the cleaned input structure sequence after gaps are removed.
  - Solution: Put the cleaned protein sequence as the first MSA record and make sure all aligned sequences have consistent length.

- "The first A3M MSA sequence does not match query_sequence after removing lowercase insertions and gaps."
  - Cause: The first sequence in the A3M MSA is not the same as the cleaned input structure sequence after lowercase insertions and gaps are removed.
  - Solution: Put the cleaned protein sequence as the first A3M record and regenerate the alignment if needed.

- "hmmbuild failed"
  - Cause: HMMER `hmmbuild` failed while building the HMM profile from the cleaned Stockholm MSA.
  - Solution: Confirm that HMMER is installed and available, then check the cleaned MSA and earlier messages in `log.txt`.

- "HMM length"
  - Cause: The number of match emission rows parsed from the HMM profile does not match the cleaned protein sequence length.
  - Solution: Regenerate the MSA from the same cleaned protein sequence as the structure passed to `-i`, or omit `-m` if conservation scores are not needed.

- "Failed to load ESM2 model"
  - Cause: The selected ESM-2 model cannot be loaded, often because model files or dependencies are unavailable in the runtime environment.
  - Solution: Check the ESM installation and model cache, or use the default smaller model if a larger model is not available.

- "Failed to compute pockets."
  - Cause: PyVOL pocket detection failed on the cleaned structure or with the selected pocket parameters.
  - Solution: Check that PyVOL is installed and try standard pocket parameters before adjusting radius or volume thresholds.

- "Failed to obtain SMILES for substrate"
  - Cause: A substrate name could not be resolved to a SMILES string through the supported chemical lookup route.
  - Solution: Check the spelling, use a more specific substrate name, increase `--substrate_max_synonyms`, or provide the SMILES string directly with `-s`.

- "Invalid SMILES"
  - Cause: A direct SMILES input cannot be parsed by RDKit.
  - Solution: Check the SMILES syntax and use a valid canonical or isomeric SMILES string.

- "Failed to generate Mol(2D) for substrate"
  - Cause: RDKit could not convert the resolved SMILES into a valid 2D molecular object.
  - Solution: Check whether the SMILES string represents a supported small molecule and try a corrected substrate name or direct SMILES input.

- "Failed to save SDF file."
  - Cause: RDKit attempted to write a generated 3D conformation, but the SDF file was not created or was empty.
  - Solution: Check that the output directory is writable and that there is enough disk space.

- "mk_prepare_receptor.py failed with return code"
  - Cause: Meeko receptor preparation started but failed while converting the cleaned protein structure to PDBQT.
  - Solution: Review the output tail in `log.txt`, confirm Meeko is installed correctly, and check that the cleaned protein input is valid.

- "Vina docking failed for"
  - Cause: AutoDock Vina failed for a substrate combination and docking box, often because of receptor or ligand preparation issues, unsuitable docking box settings, or an unavailable Vina executable.
  - Solution: Confirm Vina is installed, review the docking error in `log.txt`, and try default docking settings or a manually defined box around the active site.

- "No valid docking results were found for any substrate combination and docking box."
  - Cause: Docking completed attempts but no valid pose could be parsed and accepted.
  - Solution: Check substrate SDF generation, box center and size, Vina installation, and consider increasing `--dock_max_attempt_num` or disabling early stop.

- "integrated_graph missing in integrate report."
  - Cause: The integration step returned an invalid report without the graph field expected by batch.
  - Solution: Check `log.txt` for the first earlier error, because this is usually caused by an upstream report-generation or integration failure.

- "Failed to save integrate JSON"
  - Cause: Batch could not write the integrated report, node-only JSON file, or edge-only JSON file because of a filesystem, permission, path, or disk-space problem.
  - Solution: Check that the output directory is writable and that there is enough disk space.

- Substrate output files are missing even though `--substrate_names` was provided.
  - Cause: Batch intentionally falls back to the protein-only workflow when substrate parsing, SMILES completion, substrate 3D generation, docking, dock report generation, or docked substrate validation fails.
  - Solution: Check `log.txt` for the first warning ending with `Falling back to protein-only workflow.`, verify substrate names or SMILES strings, and rerun with corrected substrate input or more suitable docking settings.

- Conservation scores are missing from residue nodes.
  - Cause: No MSA was provided, or MSA/HMM processing failed before conservation scores were generated.
  - Solution: Provide a matched MSA with `-m`, ensure it was generated from the cleaned protein sequence, and check `log.txt` for MSA/HMM errors.

- Output files are missing or fewer than expected.
  - Cause: The command failed before final integration files were written, or optional intermediate files were not requested with `--save_extra_outputs`.
  - Solution: Check `log.txt`, confirm the output directory passed to `-o`, and remember that cleaned MSA/HMM files require `-m`, while substrate, docked SDF, and complex files require successful substrate and docking workflows plus `--save_extra_outputs`.

- Output file names do not match the expected protein name.
  - Cause: Final output file names use the protein name derived from the cleaned input structure file name, after filename shortening if needed.
  - Solution: Check the cleaned input structure file name and look for `integrate_report_{protein_name}.json`, `integrate_nodes_{protein_name}.json`, `integrate_edges_{protein_name}.json`, and `log.txt` in the output directory.

# dependencies:

- Biopython
- NumPy
- OpenMM
- DSSP
- ProDy
- ESM
- HMMER
- PyVOL
- RDKit
- AutoDock Vina
- Meeko
- JSON


# references:

- Biopython:
  https://biopython.org/

- OpenMM:
  https://openmm.org/

- DSSP:
  https://github.com/PDB-REDO/dssp

- ProDy:
  http://prody.csb.pitt.edu/

- ESM:
  https://github.com/facebookresearch/esm

- HMMER:
  http://hmmer.org/

- PyVOL:
  https://github.com/schlessinger-lab/pyvol

- RDKit:
  https://www.rdkit.org/

- AutoDock Vina:
  https://vina.scripps.edu/

- Meeko:
  https://github.com/forlilab/Meeko

- JSON:
  https://www.json.org/
