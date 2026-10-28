# google-deepmind/alphafold3 Overview
The AlphaFold 3 repository predicts the 3D structures of biological molecules and their complexes. It transforms diverse biological input data into structured features and employs a neural network to generate atomic coordinates and confidence metrics.

The system provides Docker images for configuring the development and runtime environment. These images ensure consistent dependency management, HMMER installation, and GPU optimization. Specialized scripts manage Solid State Drives for large biological databases, particularly in cloud environments. See [Setup and Environment Configuration](#setup-and-environment-configuration) for more information.

It parses and preprocesses various biological data formats, including CIF, FASTA, and Multiple Sequence Alignments. High-performance C++ implementations with Python bindings handle data extraction and serialization efficiently. This enables the conversion of raw biological data into a structured format. [Data Parsing and Preprocessing](#data-parsing-and-preprocessing) covers these processes.

The repository generates input features for the model from parsed biological data. This includes creating Multiple Sequence Alignments, identifying structural templates, and preparing molecular information using external tools like HMMER and RDKit. The featurization pipeline orchestrates these processes, often in parallel. More details are in [Input Feature Generation](#input-feature-generation).

Molecular structures are represented using an immutable, array-based data model. This robust system aggregates chains, residues, atoms, and bonds, facilitating structural analysis and manipulation. It specifically handles biological assemblies and chemical component data. [Molecular Structure Management](#molecular-structure-management) describes these internal representations.

The core model architecture processes input features through several neural network components:
*   An Evoformer/Pairformer refines single and pair representations.
*   Atom-level cross-attention mechanisms integrate atomic information.
*   A generative diffusion denoising process predicts atom positions.
These components are implemented using Haiku modules and JAX for efficient computation. Further details are in [Core Model Architecture and Pipeline](#core-model-architecture-and-pipeline).

JAX-compatible geometric primitives handle 3D vectors, rotation matrices, and rigid transformations. These primitives utilize a struct-of-arrays design, enabling efficient batch operations and seamless integration with JAX's capabilities. It supports calculations for distances, angles, and weighted means. See [Geometric Primitives and Operations](#geometric-primitives-and-operations) for more information.

The system includes functionalities for scoring and evaluating predicted molecular structures. This involves structural alignment, RMSD calculations, and chirality detection and comparison for ligands. It also generates pseudo-beta atom positions for quality assessment. [Model Scoring and Evaluation](#model-scoring-and-evaluation) provides more details.

A set of foundational utilities and constants support the codebase. This includes configuration management, secure data loading through restricted pickling, and access to atom types, chemical component definitions, and periodic table data. These ensure consistent data handling and reliable operations. [Common Utilities and Constants](#common-utilities-and-constants) describes these shared resources.

