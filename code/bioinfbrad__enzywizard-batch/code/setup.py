#!/usr/bin/env python
from setuptools import setup, find_packages
import os
import sys

# Read the version from version.py without importing the package
# This allows the package to be installed without its dependencies being available
version_file = os.path.join(os.path.dirname(__file__), 'src', 'enzywizard_batch', 'version.py')
with open(version_file) as f:
    exec(f.read())          # defines __version__

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="enzywizard-batch",
    version=__version__,                     # dynamically read from version.py
    author="bioinfbrad",
    description=(
        "Run a complete EnzyWizard analysis workflow from a cleaned protein structure "
        "and a matched MSA file. Integrates residue properties, hydrophobic clusters, "
        "energy, flexibility, disorder, conservation, embeddings, pockets, substrate "
        "docking, molecular interactions, and graph integration."
    ),
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/bioinfbrad/enzywizard-batch",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    python_requires=">=3.10",
    install_requires=[
        "biopython>=1.86",          # protein structure I/O, sequence handling
        "numpy>=1.23.5,<2",         # numerical backend
        "rdkit>=2026.3.1",          # cheminformatics for substrates
        "openmm>=8.5.0",            # molecular mechanics (energy, minimization)
        "prody>=2.6.1",             # elastic network models (flexibility)
        "fair-esm>=2.0.0",          # residue embeddings
        "bio-pyvol>=1.7.8",         # binding pocket detection
        "meeko>=0.7.1",             # ligand preparation for docking
        "pdbfixer>=1.12",           # structure cleaning
        "requests>=2.33.0",         # HTTP requests (API calls)
        "packaging>=26.1",          # version handling
        # External binaries are NOT listed here – they must be added in the Conda
        # recipe's run dependencies: hmmer, msms, dssp, vina.
    ],
    entry_points={
        "console_scripts": [
            "enzywizard-batch = enzywizard_batch.cli:main",
        ],
    },
    include_package_data=True,
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3.10",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
)
