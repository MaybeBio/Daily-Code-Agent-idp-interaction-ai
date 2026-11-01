# ssdata: Amino Acid Data and Mappings

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [soursop/data/phi_excluded_volume_tripeptides.pickle](soursop/data/phi_excluded_volume_tripeptides.pickle)
- [soursop/data/psi_excluded_volume_tripeptides.pickle](soursop/data/psi_excluded_volume_tripeptides.pickle)
- [soursop/ssdata.py](soursop/ssdata.py)

</details>



## Purpose and Scope

The `ssdata` module provides reference data and mappings for amino acid nomenclature, residue validation, and excluded volume calculations. This module contains well-defined biological constants and data structures used throughout SOURSOP for sequence analysis, residue validation, and conformational sampling assessment.

For information about how this data is used in the PENGUIN pipeline, see [Reference Models and EV Data](#5.4). For utilities that perform calculations using this data, see [sstools: Numerical and File Utilities](#7.2).

**Sources:** [soursop/ssdata.py:14-18]()

## Module Overview

The `ssdata` module serves as the central repository for amino acid-related constants and reference data. It contains no computational logic—only data structures that are imported and used by other modules throughout the codebase.

```mermaid
graph TB
    subgraph "ssdata.py Data Structures"
        THREE_TO_ONE["THREE_TO_ONE<br/>dict: 3-letter → 1-letter"]
        ONE_TO_THREE["ONE_TO_THREE<br/>dict: 1-letter → 3-letter"]
        VALID_RES["ALL_VALID_RESIDUE_NAMES<br/>list: recognized residues"]
        SIDECHAIN["DEFAULT_SIDECHAIN_VECTOR_ATOMS<br/>dict: residue → atom name"]
        EV_MAPPER["EV_RESIDUE_MAPPER<br/>dict: residue → EV group"]
        PSI_ANGLES["PSI_EV_ANGLES_DICT<br/>numpy array: EV reference"]
        PHI_ANGLES["PHI_EV_ANGLES_DICT<br/>numpy array: EV reference"]
    end
    
    subgraph "Consumer Modules"
        SSProtein["SSProtein<br/>sequence operations"]
        SSTrajectory["SSTrajectory<br/>residue validation"]
        SamplingQuality["SamplingQuality<br/>EV reference model"]
        SSTools["sstools<br/>helper functions"]
    end
    
    THREE_TO_ONE --> SSProtein
    ONE_TO_THREE --> SSProtein
    VALID_RES --> SSTrajectory
    SIDECHAIN --> SSProtein
    EV_MAPPER --> SamplingQuality
    PSI_ANGLES --> SamplingQuality
    PHI_ANGLES --> SamplingQuality
    
    VALID_RES --> SSTools
    THREE_TO_ONE --> SSTools
```

**Sources:** [soursop/ssdata.py:1-144]()

## Amino Acid Nomenclature Mappings

### Three-Letter to One-Letter Code

The `THREE_TO_ONE` dictionary maps standard three-letter amino acid codes to single-letter codes. This includes the 20 standard amino acids plus terminal capping groups.

| Three-Letter Code | One-Letter Code | Description |
|------------------|-----------------|-------------|
| ALA | A | Alanine |
| CYS | C | Cysteine |
| ASP | D | Aspartate |
| GLU | E | Glutamate |
| PHE | F | Phenylalanine |
| GLY | G | Glycine |
| HIS | H | Histidine |
| ILE | I | Isoleucine |
| LYS | K | Lysine |
| LEU | L | Leucine |
| MET | M | Methionine |
| ASN | N | Asparagine |
| PRO | P | Proline |
| GLN | Q | Glutamine |
| ARG | R | Arginine |
| SER | S | Serine |
| THR | T | Threonine |
| VAL | V | Valine |
| TRP | W | Tryptophan |
| TYR | Y | Tyrosine |

#### Terminal Groups and Caps

| Three-Letter Code | One-Letter Code | Description |
|------------------|-----------------|-------------|
| ACE | < | Acetyl N-terminal cap |
| NME | > | N-methylamide C-terminal cap |
| NAC | > | N-acetyl C-terminal cap |
| NH2 | ( | Amide C-terminal cap |
| FOR | ) | Formyl N-terminal cap |

**Sources:** [soursop/ssdata.py:22-46]()

### One-Letter to Three-Letter Code

The `ONE_TO_THREE` dictionary provides the inverse mapping from single-letter codes to three-letter codes. Note that this mapping is not perfectly bijective for terminal groups (e.g., both NME and NAC map to '>').

```mermaid
graph LR
    subgraph "Standard Amino Acids"
        A["A"] --> ALA["ALA"]
        K["K"] --> LYS["LYS"]
        R["R"] --> ARG["ARG"]
    end
    
    subgraph "Terminal Groups"
        LT["<"] --> ACE["ACE"]
        GT[">"] --> NME["NME"]
        LP["("] --> NH2["NH2"]
        RP[")"] --> FOR["FOR"]
    end
    
    subgraph "Usage"
        SSProtein["SSProtein.get_amino_acid_sequence()"]
        Conversion["sequence string conversion"]
    end
    
    ALA --> Conversion
    NME --> Conversion
    Conversion --> SSProtein
```

**Sources:** [soursop/ssdata.py:48-71]()

## Valid Residue Names

### ALL_VALID_RESIDUE_NAMES

The `ALL_VALID_RESIDUE_NAMES` list defines all residue names recognized by SOURSOP when validating whether a molecule is a protein. This list is used during trajectory loading to filter protein chains from other molecular species.

```python
ALL_VALID_RESIDUE_NAMES = [
    # Standard amino acids
    'ALA','CYS','ASP','ASH','GLU','GLH','PHE','GLY',
    'HIE','HIS','HID','HIP','ILE','LEU','LYS','LYD',
    'MET','ASN','PRO','GLN','ARG','SER','THR','VAL',
    'TRP','TYR',
    
    # Non-standard/modified amino acids
    'AIB',  # α-aminoisobutyric acid
    'ABA',  # α-aminobutyric acid
    'NVA',  # Norvaline
    'NLE',  # Norleucine
    'ORN',  # Ornithine
    'DAB',  # Diaminobutyric acid
    
    # Post-translational modifications
    'PTR',  # Phosphotyrosine
    'TPO',  # Phosphothreonine
    'SEP',  # Phosphoserine
    'KAC',  # Acetyllysine
    'KM1',  # Monomethyllysine
    'KM2',  # Dimethyllysine
    'KM3',  # Trimethyllysine
    
    # Terminal groups
    'ACE','NME','FOR','NH2'
]
```

### Protonation State Variants

Several residues have multiple protonation state variants that are all recognized:

| Base Residue | Variants | Description |
|-------------|----------|-------------|
| Histidine | HIS, HIE, HID, HIP | Different protonation states (epsilon, delta, both) |
| Aspartate | ASP, ASH | Deprotonated / protonated |
| Glutamate | GLU, GLH | Deprotonated / protonated |
| Lysine | LYS, LYD | Protonated / deprotonated |

**Sources:** [soursop/ssdata.py:111-112]()

## Sidechain Vector Atom Definitions

### DEFAULT_SIDECHAIN_VECTOR_ATOMS

The `DEFAULT_SIDECHAIN_VECTOR_ATOMS` dictionary maps each residue type to the default atom used for sidechain vector calculations. These atoms are typically at or near the end of the sidechain and are used for methods like `get_sidechain_alignment_angle()`.

```mermaid
graph TB
    subgraph "Atom Selection Strategy"
        Small["Small Residues<br/>ALA, CYS, SER, THR<br/>→ CB or first heavy atom"]
        Charged["Charged Residues<br/>LYS, ARG, ASP, GLU<br/>→ Terminal functional group"]
        Aromatic["Aromatic Residues<br/>PHE, TRP, TYR<br/>→ Ring center atom"]
        Branch["Branched Residues<br/>VAL, LEU, ILE<br/>→ Branch point or terminal"]
    end
    
    subgraph "Examples"
        LYS_Ex["LYS → NZ<br/>terminal amine"]
        ARG_Ex["ARG → CZ<br/>guanidinium carbon"]
        PHE_Ex["PHE → CZ<br/>ring center"]
        SER_Ex["SER → OG<br/>hydroxyl oxygen"]
    end
    
    Small --> SER_Ex
    Charged --> LYS_Ex
    Charged --> ARG_Ex
    Aromatic --> PHE_Ex
```

### Special Cases

| Residue | Atom | Notes |
|---------|------|-------|
| GLY | ERROR | Glycine has no sidechain |
| ACE, NME, FOR, NH2 | ERROR | Terminal groups have no standard sidechain |
| PRO | CG | Cyclic structure, CG chosen as representative |
| VAL | CB | Branch point is the sidechain vector |

### Post-Translational Modifications

Modified residues use the same atom as their unmodified counterparts:

| Modified Residue | Base | Atom | Modification |
|-----------------|------|------|--------------|
| PTR | TYR | CZ | Phosphorylation |
| TPO | THR | CB | Phosphorylation |
| SEP | SER | OG | Phosphorylation |
| KAC, KM1, KM2, KM3 | LYS | NZ | Acetylation/Methylation |

**Sources:** [soursop/ssdata.py:73-109]()

## Excluded Volume Reference Data

### EV_RESIDUE_MAPPER

The `EV_RESIDUE_MAPPER` dictionary maps amino acids into three groups for excluded volume (EV) calculations in the PENGUIN pipeline. This coarse-graining reduces the reference space while maintaining chemical realism.

```mermaid
graph TB
    subgraph "20 Amino Acids"
        PRO["PRO"]
        
        AlaGroup["GLY, ALA, CYS<br/>ASN, GLN, SER, THR"]
        
        LeuGroup["LEU, MET, ASP, GLU<br/>ARG, VAL, TRP, TYR<br/>PHE, HIS, ILE, LYS"]
    end
    
    subgraph "EV Groups"
        EV_PRO["PRO group<br/>unique cyclic structure"]
        EV_ALA["ALA group<br/>~alanine-like<br/>small sidechains"]
        EV_LEU["LEU group<br/>~leucine-like<br/>larger sidechains"]
    end
    
    PRO --> EV_PRO
    AlaGroup --> EV_ALA
    LeuGroup --> EV_LEU
    
    subgraph "Reference Data Usage"
        PSI["PSI_EV_ANGLES_DICT<br/>backbone psi angles"]
        PHI["PHI_EV_ANGLES_DICT<br/>backbone phi angles"]
    end
    
    EV_PRO --> PSI
    EV_ALA --> PSI
    EV_LEU --> PSI
    
    EV_PRO --> PHI
    EV_ALA --> PHI
    EV_LEU --> PHI
```

#### Group Assignments

**PRO Group** (Special):
- PRO: Unique due to cyclic backbone constraint

**ALA Group** (Approximately alanine-like):
- GLY: Lacks sidechain, treated as small
- ALA: Reference for small sidechains
- CYS: Small polar sidechain
- ASN, GLN: Amide sidechains
- SER, THR: Hydroxyl sidechains

**LEU Group** (Approximately leucine-like):
- LEU: Reference for larger sidechains
- All other residues with substantial sidechains
- Includes charged (ASP, GLU, ARG, LYS)
- Includes hydrophobic (VAL, ILE, MET)
- Includes aromatic (PHE, TYR, TRP, HIS)

**Sources:** [soursop/ssdata.py:115-140]()

### Angle Reference Distributions

#### PSI_EV_ANGLES_DICT and PHI_EV_ANGLES_DICT

These dictionaries contain precomputed backbone dihedral angle distributions from excluded volume tripeptide simulations. They are loaded from pickle files and used by the `SamplingQuality` class as reference distributions for assessing conformational sampling.

```mermaid
graph TB
    subgraph "Data Files"
        PsiFile["psi_excluded_volume_tripeptides.pickle<br/>soursop/data/"]
        PhiFile["phi_excluded_volume_tripeptides.pickle<br/>soursop/data/"]
    end
    
    subgraph "ssdata Loading"
        LoadPsi["np.load()<br/>PSI_EV_ANGLES_DICT"]
        LoadPhi["np.load()<br/>PHI_EV_ANGLES_DICT"]
    end
    
    subgraph "Dictionary Structure"
        Keys["Keys: tripeptide codes<br/>e.g., 'ALA-PRO-ALA'"]
        Values["Values: angle arrays<br/>from EV simulations"]
    end
    
    subgraph "SamplingQuality Usage"
        Mapper["Uses EV_RESIDUE_MAPPER<br/>to map sequence"]
        Reference["Loads reference angles<br/>from dictionaries"]
        Compare["Computes Hellinger distance<br/>between simulation and reference"]
    end
    
    PsiFile --> LoadPsi
    PhiFile --> LoadPhi
    LoadPsi --> Keys
    LoadPhi --> Keys
    Keys --> Values
    
    LoadPsi --> Mapper
    LoadPhi --> Mapper
    Mapper --> Reference
    Reference --> Compare
```

The reference data structure:

| Dictionary | Contains | Format |
|-----------|----------|--------|
| `PSI_EV_ANGLES_DICT` | Psi angle distributions for tripeptides | `dict[str, np.ndarray]` |
| `PHI_EV_ANGLES_DICT` | Phi angle distributions for tripeptides | `dict[str, np.ndarray]` |

Each key represents a tripeptide context (e.g., "ALA-PRO-ALA") using the EV group names, and each value is a NumPy array of angle values sampled from excluded volume simulations of that tripeptide.

**Sources:** [soursop/ssdata.py:142-144]()

## Integration with Other Modules

```mermaid
graph TB
    subgraph "ssdata.py Exports"
        T2O["THREE_TO_ONE"]
        O2T["ONE_TO_THREE"]
        VALID["ALL_VALID_RESIDUE_NAMES"]
        SC["DEFAULT_SIDECHAIN_VECTOR_ATOMS"]
        EVM["EV_RESIDUE_MAPPER"]
        PSI["PSI_EV_ANGLES_DICT"]
        PHI["PHI_EV_ANGLES_DICT"]
    end
    
    subgraph "SSTrajectory"
        ResCheck["_get_protein_residues()<br/>filters by ALL_VALID_RESIDUE_NAMES"]
        ChainID["identify protein chains<br/>during trajectory loading"]
    end
    
    subgraph "SSProtein"
        SeqConvert["get_amino_acid_sequence()<br/>uses THREE_TO_ONE"]
        SCVector["get_sidechain_alignment_angle()<br/>uses DEFAULT_SIDECHAIN_VECTOR_ATOMS"]
    end
    
    subgraph "SamplingQuality"
        EVMap["_map_sequence_to_ev()<br/>uses EV_RESIDUE_MAPPER"]
        RefLoad["PrecomputedDihedralInterface<br/>uses PSI/PHI_EV_ANGLES_DICT"]
        Hellinger["compute_hellinger_distance()<br/>compares to reference"]
    end
    
    VALID --> ResCheck
    ResCheck --> ChainID
    
    T2O --> SeqConvert
    SC --> SCVector
    
    EVM --> EVMap
    PSI --> RefLoad
    PHI --> RefLoad
    RefLoad --> Hellinger
```

### Usage Example Locations

| Module | Function/Method | ssdata Usage |
|--------|----------------|--------------|
| SSTrajectory | `__init__()` | Validates residues against `ALL_VALID_RESIDUE_NAMES` |
| SSProtein | `get_amino_acid_sequence()` | Converts residues using `THREE_TO_ONE` |
| SSProtein | `get_sidechain_alignment_angle()` | Looks up atoms in `DEFAULT_SIDECHAIN_VECTOR_ATOMS` |
| SamplingQuality | `_get_reference_distributions()` | Maps sequence using `EV_RESIDUE_MAPPER` |
| SamplingQuality | `PrecomputedDihedralInterface` | Loads data from `PSI_EV_ANGLES_DICT`, `PHI_EV_ANGLES_DICT` |

**Sources:** [soursop/ssdata.py:1-144]()

## Data File Locations

The excluded volume angle data is stored in the package data directory and accessed via the `soursop.get_data()` function:

| File | Location | Format | Size |
|------|----------|--------|------|
| psi_excluded_volume_tripeptides.pickle | soursop/data/ | NumPy pickle | Variable |
| phi_excluded_volume_tripeptides.pickle | soursop/data/ | NumPy pickle | Variable |

These files are included in the SOURSOP distribution and are loaded once when the `ssdata` module is imported. The data remains in memory for the lifetime of the Python session, avoiding repeated disk I/O.

**Sources:** [soursop/ssdata.py:142-144]()

---