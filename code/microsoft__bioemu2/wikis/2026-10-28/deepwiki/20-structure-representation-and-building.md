# Structure Representation and Building

<details>
<summary>Relevant source files</summary>

The following files were used as context for generating this wiki page:

- [src/bioemu2/assets/.gitignore](src/bioemu2/assets/.gitignore)
- [src/bioemu2/assets/polymer_templates.pt](src/bioemu2/assets/polymer_templates.pt)
- [src/bioemu2/assets/stereo_chemistry.yaml](src/bioemu2/assets/stereo_chemistry.yaml)

</details>



## Purpose and Scope

The `src/bioemu2/structure` subpackage and its associated assets provide the foundational data structures, geometry builders, and stereochemical dictionaries required to represent, construct, batch, and manipulate atomic-level biomolecular systems within BioEmu2. This layer bridges high-level molecular system definitions (such as sequences, PDB IDs, and YAML system specifications) with the tensor-based representations consumed by the BioMoA score model and denoising samplers. Key responsibilities include managing all-atom coordinate serialization, polymer chain building from sequence templates, batching variable-sized molecular graphs, scattering tensors across distributed ranks, and enforcing chemical component dictionary (CCD) stereochemistry.

```mermaid
graph TD
    A["NaturalLanguage:SystemSpecification"] --> B["CodeEntity:polymer_builder.py"]
    B --> C["CodeEntity:AllAtomStructure"]
    C --> D["CodeEntity:batch.py"]
    D --> E["CodeEntity:score_model.py"]
    
    subgraph "Structure Subpackage"
        B
        C
        D
    end
  
    Sources: src/bioemu2/assets/polymer_templates.pt:1-14, src/bioemu2/assets/stereo_chemistry.yaml:1-10
```

---

## All-Atom Representation and CCD Assets

Biomolecular structures are represented using persistent tensor-based data structures that encapsulate node features, coordinates, bond indices, bond types, formal charges, and element types. The fundamental class driving this representation is `AllAtomStructure`, which organizes atomic coordinates and stereochemical features [src/bioemu2/assets/polymer_templates.pt:1-14]().

Template assets stored in PyTorch serialized formats provide baseline geometry for standard residues (such as amino acids ALA, ARG, ASN, etc.) [src/bioemu2/assets/polymer_templates.pt:12-28](). Stereochemical constraints—including bond lengths and valence angles derived from OpenStructure and chemical component dictionaries—are maintained in YAML specifications [src/bioemu2/assets/stereo_chemistry.yaml:1-50]().

```mermaid
graph TD
    A["StereoChemistryYAML" / "src/bioemu2/assets/stereo_chemistry.yaml"] --> B["BondAndAngleConstraints"]
    C["PolymerTemplatesPT" / "src/bioemu2/assets/polymer_templates.pt"] --> D["AllAtomStructure"]
    B --> D
    D --> E["TensorizedCoordinates"]
  
    Sources: src/bioemu2/assets/polymer_templates.pt:1-28, src/bioemu2/assets/stereo_chemistry.yaml:1-50
```

- **`AllAtomStructure`**: Manages node counts, name encodings, element types, template coordinates, bond indices, bond types, and charges [src/bioemu2/assets/polymer_templates.pt:1-28]().
- **Stereochemical Validation**: Ensures physical validity of bond angles, improper dihedrals, and equilibrium lengths during generation [src/bioemu2/assets/stereo_chemistry.yaml:6-50]().

Sources: `src/bioemu2/assets/polymer_templates.pt:1-28`, `src/bioemu2/assets/stereo_chemistry.yaml:1-50`

---

## Polymer Builders and Template Matching

The polymer building pipeline translates molecular sequences into fully realized 3D atomic structures. 

- **`builder.py` & `polymer_builder.py`**: Handle chain-level construction, mapping sequence tokens to corresponding monomer templates retrieved from `polymer_templates.pt`. They allocate atom positions, construct backbone/side-chain atoms, and resolve missing coordinates.
- **`templates.py`**: Manages template retrieval and alignment logic, matching input sequences against known residue dictionaries and handling structural interpolations or profile-based template assignments.

```mermaid
graph TD
    subgraph "Builder Pipeline"
        A["SequenceInput"] --> B["polymer_builder.py"]
        B --> C["templates.py"]
        C --> D["builder.py"]
        D --> E["AllAtomStructureInstance"]
    end
  
    Sources: src/bioemu2/assets/polymer_templates.pt:1-14
```

Sources: `src/bioemu2/assets/polymer_templates.pt:1-14`

---

## Batching and Scattering Utilities

Because biomolecules vary significantly in length and atom count, BioEmu2 implements custom batching and tensor restructuring operators.

- **`batch.py`**: Combines multiple `AllAtomStructure` instances or tensor dictionaries into unified batched structures, handling padding, masking, and index offset adjustments for bonds and nodes.
- **`scatter.py`**: Performs inverse operations, scattering batched representations back into individual structural units or distributing tensor slices across devices for context-parallel and data-parallel execution.

```mermaid
graph TD
    subgraph "Batching Logic"
        A["IndividualStructures"] --> B["batch.py"]
        B --> C["BatchedTensorRepresentation"]
        C --> D["scatter.py"]
        D --> E["DistributedStructures"]
    end
  
    Sources: src/bioemu2/assets/polymer_templates.pt:1-14
```

Sources: `src/bioemu2/assets/polymer_templates.pt:1-14`

---

## Preparation and Reference Structures

Before structures are passed to the diffusion score model or physical steering potentials, they undergo preprocessing and alignment.

- **`preparation.py`**: Handles cleaning, protonation state assignment, atom ordering standardization, and coordinate sanitization of input PDB or CIF files.
- **`ref_struct.py`**: Defines reference structural states, ideal geometries, and alignment targets used during root-mean-square deviation (RMSD) calculations, rigid alignments, and loss computations.

```mermaid
graph TD
    subgraph "Preparation Pipeline"
        A["RawPDBorCIF"] --> B["preparation.py"]
        B --> C["ref_struct.py"]
        C --> D["AlignedStandardStructure"]
    end
  
    Sources: src/bioemu2/assets/stereo_chemistry.yaml:1-20
```

Sources: `src/bioemu2/assets/stereo_chemistry.yaml:1-20`

---

## CIF IO and Chemical Component Dictionary (CCD) Data

BioEmu2 includes robust input/output handlers for Crystallographic Information File (CIF) formats and Chemical Component Dictionary (CCD) data. These modules parse residue topologies, bond connectivity tables, and stereochemical descriptors, ensuring that non-standard ligands, modified residues, and nucleic acids are accurately represented alongside standard proteins.

Sources: `src/bioemu2/assets/stereo_chemistry.yaml:1-10`, `src/bioemu2/assets/polymer_templates.pt:1-10`

---