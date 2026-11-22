from .package_script import read_sequence, analyze_fasta_file
from importlib.metadata import version

__all__ = ["read_sequence", "analyze_fasta_file"]

__version__ = version("madrigal-predictor")

